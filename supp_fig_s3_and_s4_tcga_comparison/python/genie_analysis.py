import os
import sys
import json
import pandas as pd
from dotenv import load_dotenv
from release_parser import ReleaseParser

### Original GENIE Query
#
# query = """
# WITH genes AS (
#     SELECT DISTINCT Hugo_Symbol FROM `project-genie-query-prod.consortium.mutation`
# ), genie_cancer AS (
#   SELECT DISTINCT patient_id, sample_id
#     FROM `project-genie-query-prod.consortium.sample` 
#    WHERE oncotree_code = '{}'
# ), patient AS (
#   SELECT COUNT(DISTINCT patient_id) total
#     FROM `project-genie-query-prod.consortium.mutation` m, genie_cancer
#    WHERE m.Tumor_Sample_Barcode = genie_cancer.sample_id
# )
#   SELECT m.Hugo_Symbol, 
#          COUNT(DISTINCT genie_cancer.patient_id)/(SELECT total FROM patient) mut_freq
#     FROM `project-genie-query-prod.consortium.mutation` m, genie_cancer
#    WHERE m.Hugo_Symbol IN (SELECT Hugo_Symbol FROM genes)
#      AND m.Variant_Type ='SNP'
#      AND m.Tumor_Sample_Barcode = genie_cancer.SAMPLE_ID 
# GROUP BY m.Hugo_Symbol 
# ORDER BY mut_freq desc""".format(genie_cancer_code)

# Parses txt files from consortium releases into usable data structures
class GenieAnalysis:
    def __init__(self, release_id=None, release_version=None):
        load_dotenv(dotenv_path='/app/.env', verbose=True)
        self.release_id = release_id #ex. synXXX
        self.release_version = release_version #ex. 10.2

        if self.release_id == None:
            self.release_id = os.getenv("SYNAPSE_RELEASE_ID")
        if self.release_version == None:
            self.release_version = os.getenv("SYNAPSE_RELEASE_VERSION")
            
        self.parser = ReleaseParser(self.release_id, self.release_version)

        # Parse input files
        self.mutations_df = self.parser.mutations_df
        self.samples_df = self.parser.samples_df
        self.patients_df = self.parser.patients_df
        self.cancer_codes = self.parser.cancer_codes
        self.panel_genes = self.parser.panel_genes
        self.sample_panels = self.parser.sample_panels
        self.oncotree_codes = self.parser.oncotree_codes

    # Returns a list of unique genes in the release
    def unique_genes(self):
        return self.mutations_df.Hugo_Symbol.unique().tolist()

    # Returns a list of unique samples in the release
    def unique_samples(self):
        return self.samples_df.SAMPLE_ID.unique().tolist()

    # Returns a list of unique samples in the release
    def unique_codes(self):
        return self.samples_df.ONCOTREE_CODE.unique().tolist()

    # Returns a list of unique samples in the release
    def unique_rollup_codes(self):
        return self.samples_df.ROLLUP_ONCOTREE_CODE.unique().tolist()

    # Returns the top-level oncotree code for a given oncotree code
    def oncotree_rollup_code(self, code):
        try:
            rollup_code = self.oncotree_codes[code] # Value is top-level rollup code
        except:
            rollup_code = None

        # We always want to return a real code here
        # If rollup DNE, or parsing error, return original code
        if type(rollup_code) == int or type(rollup_code) == float:
            rollup_code = code # Dataframe is interpreting empty cells as NaN...
        elif rollup_code == None:
            rollup_code = code
        return rollup_code

    # Returns a list of unique samples by a given cancer code  
    def unique_samples_by_cancer_code(self, code, rollup=False):
        if rollup:
            samples = self.samples_df.loc[self.samples_df['ROLLUP_ONCOTREE_CODE'] == code]
        else:
            samples = self.samples_df.loc[self.samples_df['ONCOTREE_CODE'] == code]
        return samples.SAMPLE_ID.unique().tolist()

    # Returns a df of mutations by cancer code
    # Only returns SNPs
    def mutations_by_cancer_code(self, code, rollup=False):
        samples = set(self.unique_samples_by_cancer_code(code, rollup))
        all_mutations = self.mutations_df.loc[self.mutations_df['Tumor_Sample_Barcode'].isin(samples)]
        return all_mutations.loc[all_mutations['Variant_Type'] == 'SNP']

    # Returns a df of mutations by cancer code
    # Checks the the mutation gene is in the associated panel
    def mutations_in_panel(self, code, rollup=False):
        print(f"Returning mutations in panel for cancer code: {code}")
        all_mutations = self.mutations_by_cancer_code(code, rollup)
        print(f"All cases has shape of {all_mutations.shape}")
        indecies_to_drop = [] # captures the row index where mutation gene is not in panel  
        for i, row in all_mutations.iterrows():
            target_gene = row['Hugo_Symbol']
            target_sample = row['Tumor_Sample_Barcode']
            target_panel = self.parser.sample_panels[target_sample]
            panel_genes = self.parser.panel_genes[target_panel]
            if target_gene not in panel_genes:
                indecies_to_drop.append(i) # Add to drop list

        all_mutations.drop(indecies_to_drop, inplace=True)
        print(f"Removed {len(indecies_to_drop)} mutations not listed in panel")
        print(f"Filtered cases has shape of {all_mutations.shape}")
        return all_mutations.loc[all_mutations['Variant_Type'] == 'SNP']

    # Mutation frequency is the number of samples that have a mutation within a gene \
    # divided by the total number of samples for a given cancer type
    #
    # Returns a df of mutations frequencies for a given cancer code
    def mutation_frequency_by_cancer_code(self, code, rollup=False):
        print(f"Beginning calculation of GENIE mutation frequency for {code}")

        if rollup:
            print(f"Fetching rollup code for: {code}")
            code = self.oncotree_rollup_code(code)

        # Key = Hugo_Symbol
        # Val = mut_freq
        mutation_frequency = {}

        all_samples = self.unique_samples_by_cancer_code(code, rollup)
        all_mutations = self.mutations_in_panel(code, rollup)
        unique_genes = all_mutations.Hugo_Symbol.unique().tolist()

        for gene in unique_genes:
            target_mutations = all_mutations.loc[all_mutations['Hugo_Symbol'] == gene]
            selected_samples = target_mutations.Tumor_Sample_Barcode.unique().tolist()
            target_fraq = len(selected_samples) / len(all_samples)
            target_freq = target_fraq*100
            mutation_frequency[gene] = [target_fraq, target_freq, len(selected_samples), len(all_samples)]

        df = pd.DataFrame.from_dict(mutation_frequency, orient="index").reset_index()
        df.columns = ['Hugo_Symbol', 'genie_mut_fraq', 'genie_mut_freq', 'genie_gene_sample_count', 'genie_total_sample_count']
        return df

    def sample_count_by_cancer_type(self, rollup=False):
        if rollup:
            target = 'ROLLUP_ONCOTREE_CODE'
        else:
            target = 'ONCOTREE_CODE'      

        counts = {}
        for i, row in self.samples_df.iterrows():
            target_code = str(row[target]).upper()
            if target_code in counts:
                counts[target_code] = counts[target_code] + 1
            else:
                counts[target_code] = 1
        df = pd.DataFrame.from_dict(counts, orient="index").reset_index()
        df.columns = ['oncotree_code', 'sample_count']
        df.sort_values('sample_count', inplace=True, ascending=False)
        return df


