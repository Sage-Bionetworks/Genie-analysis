import os
import sys
import pandas as pd
from dotenv import load_dotenv
from tcga_gateway import TcgaGateway

# Analysis of TCGA data
class TcgaAnalysis:
    def __init__(self):
        load_dotenv(dotenv_path='/app/.env', verbose=True)

        self.gateway = TcgaGateway()
        self.pcawg_samples = self.parse_pcawg_samples()

    # Mutation frequency is the number of samples that have a mutation within a gene \
    # divided by the total number of samples for a given cancer type
    #
    # TERT is a special case here.
    # There are significantly more TERT mutations in GENIE data than TCGA data
    # Likely this is because regions covered during sequencing of TCGA data
    # As such, we limit samples for TERT mutations to WGS assays identified in the PCAWG dataset
    # 
    # Returns a df of mutation frequencies for a given cancer code
    def mutation_frequency_by_cancer_code(self, code):
        print(f"Beginning calculation of TCGA mutation frequency for {code}")
        tert_df = self.tert_mutation_frequency_by_cancer_code(code)
        try:
            adjusted_tert_mf = tert_df.at[0, 'tcga_mut_freq']
            tert_sample_count = tert_df.at[0, 'tcga_gene_sample_count']
            pcawg_sample_count = tert_df.at[0, 'tcga_total_sample_count']
        except:
            adjusted_tert_mf = False

        query = '''
            WITH genes AS (
                SELECT DISTINCT Hugo_Symbol FROM `project-genie-query-prod.consortium.mutation`
            ), tcga_data AS (
            SELECT COUNT(DISTINCT sample_barcode_tumor) AS unique_samples
                FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` tcga_mut
                WHERE tcga_mut.sample_barcode_tumor IN (SELECT samplebarcode FROM `isb-cgc.tcga_cohorts.*` WHERE _TABLE_SUFFIX = '{code}')
            )
            SELECT genes.Hugo_Symbol Hugo_Symbol,  
                    COUNT(DISTINCT sample_barcode_tumor)/(SELECT unique_samples FROM tcga_data) tcga_mut_fraq,
                    SAFE_MULTIPLY(COUNT(DISTINCT sample_barcode_tumor)/(SELECT unique_samples FROM tcga_data),100) tcga_mut_freq,
                    COUNT(DISTINCT sample_barcode_tumor) tcga_gene_sample_count,
                    (SELECT unique_samples FROM tcga_data) tcga_total_sample_count
                FROM genes, `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` tcga_mut
            WHERE genes.Hugo_Symbol = tcga_mut.Hugo_Symbol
                AND tcga_mut.Variant_Type = 'SNP'
                AND tcga_mut.sample_barcode_tumor IN (SELECT samplebarcode FROM `isb-cgc.tcga_cohorts.*` WHERE _TABLE_SUFFIX = '{code}')
            GROUP BY genes.Hugo_Symbol
            ORDER BY tcga_mut_fraq DESC
        '''.format(code=code)

        job = self.gateway.job(query)
        target_df = job.to_dataframe()
        if 'TERT' in target_df.values:
            try:
                original_tert_mf = target_df.loc[target_df['Hugo_Symbol'] == 'TERT']['tcga_mut_freq'].tolist()[0]
                original_tert_sample_count = target_df.loc[target_df['Hugo_Symbol'] == 'TERT']['tcga_gene_sample_count'].tolist()[0]
                total_sample_count = target_df.loc[target_df['Hugo_Symbol'] == 'TERT']['tcga_total_sample_count'].tolist()[0]
                print(f"Found {original_tert_sample_count} TERT samples out of {total_sample_count} samples")
            except:
                original_tert_mf = False
                print("Did not get an original tert mutation frequency...")

            if original_tert_mf and adjusted_tert_mf and adjusted_tert_mf != original_tert_mf:
                print("Found pcawg TERT samples. Limiting to pcawg subset for TERT...")
                print(f"TERT mutation frequency across all samples: {original_tert_mf}")
                print(f"TERT mutation frequency restricted to pcawg samples: {adjusted_tert_mf}")
                print(f"Number of pcawg samples in TCGA dataset for {code}: {pcawg_sample_count}")
                print(f"Number of TERT pcawg samples in TCGA dataset for {code}: {tert_sample_count}")
                target_df.loc[target_df['Hugo_Symbol'] == 'TERT', ['tcga_mut_fraq']] = adjusted_tert_mf
                target_df.loc[target_df['Hugo_Symbol'] == 'TERT', ['tcga_gene_sample_count']] = tert_sample_count
                target_df.loc[target_df['Hugo_Symbol'] == 'TERT', ['tcga_total_sample_count']] = pcawg_sample_count
        else:
            print("This sample did not have TERT mutations...")

        return target_df

    # Find mutation frequecy for samples containing TERT mutations within the PCAWG sample set
    def tert_mutation_frequency_by_cancer_code(self, code):
        query = '''
            WITH genes AS (
                SELECT DISTINCT Hugo_Symbol FROM `project-genie-query-prod.consortium.mutation`
                WHERE Hugo_Symbol = 'TERT'
            ), tcga_data AS (
            SELECT COUNT(DISTINCT sample_barcode_tumor) AS unique_samples
                FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` tcga_mut
                WHERE tcga_mut.sample_barcode_tumor IN (SELECT samplebarcode FROM `isb-cgc.tcga_cohorts.*` WHERE _TABLE_SUFFIX = '{code}')
                AND tcga_mut.sample_barcode_tumor IN UNNEST ({pcawg_samples})
            )
            SELECT genes.Hugo_Symbol Hugo_Symbol,  
                    COUNT(DISTINCT sample_barcode_tumor)/(SELECT unique_samples FROM tcga_data) tcga_mut_fraq,
                    SAFE_MULTIPLY(COUNT(DISTINCT sample_barcode_tumor)/(SELECT unique_samples FROM tcga_data),100) tcga_mut_freq,
                    COUNT(DISTINCT sample_barcode_tumor) tcga_gene_sample_count,
                    (SELECT unique_samples FROM tcga_data) tcga_total_sample_count
                FROM genes, `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` tcga_mut
            WHERE genes.Hugo_Symbol = tcga_mut.Hugo_Symbol
                AND tcga_mut.Variant_Type = 'SNP'
                AND tcga_mut.sample_barcode_tumor IN (SELECT samplebarcode FROM `isb-cgc.tcga_cohorts.*` WHERE _TABLE_SUFFIX = '{code}')
            GROUP BY genes.Hugo_Symbol
            ORDER BY tcga_mut_fraq DESC
        '''.format(code=code,pcawg_samples=self.pcawg_samples)

        try:
            job = self.gateway.job(query).to_dataframe()
        except:
            job = pd.DataFrame() 
        return job

    # Create a whitelist of samples to analyze for TERT mutations
    def parse_pcawg_samples(self):
        self.pcawg_df = pd.read_csv(str(os.getenv("PCAWG_SUPPLEMENTARY_TABLE")), comment="#")
        spec_id = set(self.pcawg_df["submitted_specimen_id"].to_list())
        donor_id = set(self.pcawg_df["submitted_donor_id"].to_list())
        return list(spec_id.union(donor_id))