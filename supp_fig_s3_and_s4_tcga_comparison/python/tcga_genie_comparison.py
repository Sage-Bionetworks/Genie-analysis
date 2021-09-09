import os
import glob
import time
import math
from pathlib import Path
from numpy.core.numeric import roll
import pandas as pd
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from adjustText import adjust_text
from tcga_gateway import TcgaGateway
from genie_analysis import GenieAnalysis
from tcga_analysis import TcgaAnalysis
from plot import Plot

# Compare TCGA and GENIE data
class TcgaGenieComparison:
    def __init__(self):
        load_dotenv(dotenv_path='/app/.env', verbose=True)

        self.gateway = TcgaGateway()
        self.genie_analysis = GenieAnalysis(str(os.getenv('SYNAPSE_RELEASE_ID')), str(os.getenv('SYNAPSE_RELEASE_VERSION')))
        self.tcga_analysis = TcgaAnalysis()

    # Driver code to generate plots
    def execute(self):
        start = time.time()
        self.create_infrastructure()
        self.process_data(rollup=True)
        self.process_data(rollup=False)
        self.handle_deviations()
        end = time.time()
        print(f"Process took {end - start} seconds")

    # Create output directory structure
    def create_infrastructure(self):
        targets = ["mutation_frequencies", 
                    "raw_data", 
                    "sample_counts_by_cancer_type", 
                    "sample_counts_by_gene",
                    "rmsd_plots",
                    "deviations"]
        for target in targets:
            Path(f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/direct_comparison/{target}").mkdir(parents=True, exist_ok=True)
            Path(f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/direct_comparison/no_tert/{target}").mkdir(parents=True, exist_ok=True)
            Path(f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/rollup/{target}").mkdir(parents=True, exist_ok=True)
            Path(f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/rollup/no_tert/{target}").mkdir(parents=True, exist_ok=True)

    # Calculate counts, MF, and draw plots
    def process_data(self, rollup=False):
        outpath = self.find_outpath(rollup)
        print(f"Outpath is {outpath}")
        sample_counts_comparison = []
        rmsd_comparison = [] #tcga_cancer_code, genie_cancer_code, rmsd, wrmsd
        cancer_codes = self.genie_analysis.cancer_codes["cancer_codes"]
        for i, arr in enumerate(cancer_codes):
            tcga_cancer_code = arr["tcga_cancer_code"]
            genie_cancer_code = arr["genie_cancer_code"]
            cancer_type_label = arr["cancer_type_label"]

            print(f"Generating results for TCGA code: {tcga_cancer_code}, GENIE code: {genie_cancer_code}")
            result = self.frequency_count_data(tcga_cancer_code, genie_cancer_code, rollup)
            rmsd, wrmsd = self.rmsd_by_cancer_type(result)
            rmsd_comparison.append([tcga_cancer_code, genie_cancer_code, rmsd, wrmsd])
            result.to_csv(f"{outpath}/raw_data/{genie_cancer_code}_results.tsv", sep='\t', index=True)

            try:
                sample_counts_comparison.append([genie_cancer_code, 
                        result.tcga_total_sample_count.unique().tolist()[0], 
                        result.genie_total_sample_count.unique().tolist()[0]])
            except:
                print(f"Unable to parse out sample counts for {tcga_cancer_code}")

            # Generate plots
            self.mutation_frequencies(result, tcga_cancer_code, genie_cancer_code, cancer_type_label, rmsd, wrmsd, rollup)
            self.sample_counts_by_gene(result, tcga_cancer_code, genie_cancer_code, cancer_type_label, rollup)
            print(" ")

        # rmsd by gene calcs and plots
        rmsd_by_gene = self.rmsd_by_gene(outpath)
        rmsd_gene_df = pd.DataFrame(rmsd_by_gene, columns = ['gene', 'rmsd', 'wrmsd', 'error_sum', 'cancer_type_count'])
        rmsd_gene_df.sort_values('rmsd', inplace=True, ascending=False)
        rmsd_gene_df.to_csv(f"{outpath}/raw_data/rmsd_by_gene_raw.tsv", sep='\t', index=False)
        # Filter to genes in at least 3 cancer types
        print(rmsd_gene_df)
        drop_targets = rmsd_gene_df[rmsd_gene_df['cancer_type_count'] < 3].index
        rmsd_gene_df.drop(drop_targets, inplace = True)
        print(rmsd_gene_df)
 
        # Sort by RMSD
        rmsd_gene_df.sort_values('rmsd', inplace=True, ascending=False)
        Plot.rmsd_by_gene(outpath, rmsd_gene_df, 'rmsd_by_gene')
        # Sort by wRMSD
        rmsd_gene_df.sort_values('wrmsd', inplace=True, ascending=False)
        Plot.rmsd_by_gene(outpath, rmsd_gene_df, 'wrmsd_by_gene')
        # Drop TERT
        rmsd_gene_df.drop(rmsd_gene_df[rmsd_gene_df['gene'] == 'TERT'].index, inplace = True)
        # Sort by RMSD
        rmsd_gene_df.sort_values('rmsd', inplace=True, ascending=False)
        Plot.rmsd_by_gene(outpath, rmsd_gene_df, 'rmsd_by_gene_no_tert')
        # Sort by wRMSD
        rmsd_gene_df.sort_values('wrmsd', inplace=True, ascending=False)
        Plot.rmsd_by_gene(outpath, rmsd_gene_df, 'wrmsd_by_gene_no_tert')
        # Sort by Error Sum
        rmsd_gene_df.sort_values('error_sum', inplace=True, ascending=False)
        Plot.error_sum_by_gene(outpath, rmsd_gene_df)


        # Calculate rmsd by cancer type
        rmsd_df = pd.DataFrame(rmsd_comparison, columns = ['tcga_cancer_code', 'genie_cancer_code', 'rmsd', 'wrmsd']) 
        rmsd_df.sort_values('rmsd', inplace=True, ascending=False)
        rmsd_df.to_csv(f"{outpath}/raw_data/rmsd_by_cancer_type.tsv", sep='\t', index=False)

        # Aggregate sample plots
        print("Generateing aggregate plots...")
        self.sample_counts_by_cancer_type(sample_counts_comparison, rollup)
        self.genie_oncotree_distribution(rollup)
        print(f"Done processing data for all cancer codes. Rollup was {str(rollup)}")


    def rmsd_by_gene(self, outpath):
        print(f"Calculating RMSD and wRMSD by gene.")
        gene_dict = {}
        for cancer_results_path in glob.glob(f"{outpath}/raw_data/*"):
            if 'rmsd' in cancer_results_path or 'RMSD' in cancer_results_path:
                continue
            with open(cancer_results_path, "r") as f:
                for line in f:
                    spl = line.split("\t")
                    gene = spl[1]
                    tcga_mut_freq = spl[3]
                    genie_mut_freq = spl[7]

                    if gene in gene_dict:
                        gene_dict[gene].append([tcga_mut_freq, genie_mut_freq])
                    else:
                        gene_dict[gene] = [[tcga_mut_freq, genie_mut_freq]]
        rmsd_by_gene = []
        for gene, mut_freqs in gene_dict.items():
            if gene == "Hugo_Symbol":
                continue
            absolute_sum = []
            err = []
            weighted_err = []
            for arr in mut_freqs:
                mut_freq_x = float(arr[0]) #tcga
                mut_freq_y = float(arr[1]) #genie
                
                # calculate absolute distance from point to y=x
                if mut_freq_y == mut_freq_x:
                    err.append(0)
                    weighted_err.append(0)
                    absolute_sum.append(0)
                elif mut_freq_x > mut_freq_y:
                    hyp = mut_freq_x - mut_freq_y
                    opp = hyp * math.sin(math.radians(45))
                    err.append(opp)
                    absolute_sum.append(opp)
                    weighted_err.append(opp*mut_freq_x/100)
                elif mut_freq_y > mut_freq_x:
                    hyp = mut_freq_y - mut_freq_x
                    adj = hyp * math.cos(math.radians(45))
                    err.append(adj)
                    absolute_sum.append(adj)
                    weighted_err.append(adj*mut_freq_y/100)
            
            error_sum = sum(absolute_sum)
            cancer_type_count = len(absolute_sum)
            mse = sum([x**2 for x in err]) / len(err)
            rmsd = round(mse**0.5, 2)
            weighted_mse = sum([x**2 for x in weighted_err]) / len(weighted_err)
            wrmsd = round(weighted_mse**0.5, 2)
            rmsd_by_gene.append([gene, rmsd, wrmsd, error_sum, cancer_type_count])
        return rmsd_by_gene

    def rmsd_by_cancer_type(self, result):
        err = []
        weighted_err = []
        for i, row in result.iterrows():
            mut_freq_x = row['tcga_mut_freq']
            mut_freq_y = row['genie_mut_freq']

            # calculate absolute distance from point to y=x
            if mut_freq_y == mut_freq_x:
                err.append(0)
                weighted_err.append(0)
            elif mut_freq_x > mut_freq_y:
                hyp = mut_freq_x - mut_freq_y
                opp = hyp * math.sin(math.radians(45))
                err.append(opp)
                weighted_err.append(opp*mut_freq_x/100)
            elif mut_freq_y > mut_freq_x:
                hyp = mut_freq_y - mut_freq_x
                adj = hyp * math.cos(math.radians(45))
                err.append(adj)
                weighted_err.append(adj*mut_freq_y/100)
        # mse = (1/n) * sum[ (acual - pred)**2 ]
        # rmsd = sqrt(mse)
        mse = sum([x**2 for x in err]) / len(err)
        rmsd = round(mse**0.5, 2)
        print(f"Calculated RMSD of {rmsd}")
        weighted_mse = sum([x**2 for x in weighted_err]) / len(weighted_err)
        wrmsd = round(weighted_mse**0.5, 2)
        print(f"Calculated wRMSD of {wrmsd}")

        return rmsd, wrmsd

    # Get GENIE && TCGA mutation frequency and count data
    def frequency_count_data(self, tcga_cancer_code, genie_cancer_code, rollup=False):
        print(f"Fetching MF/count results for TCGA cancer code: {tcga_cancer_code}")
        genie_mutation_frequencies = self.genie_analysis.mutation_frequency_by_cancer_code(genie_cancer_code, rollup)
        tcga_mutation_frequencies = self.tcga_analysis.mutation_frequency_by_cancer_code(tcga_cancer_code)
        return pd.merge(tcga_mutation_frequencies, genie_mutation_frequencies, on='Hugo_Symbol')

    # Calculate mutation frequences and generate comparison plots for all cancer codes
    def mutation_frequencies(self, result, tcga_cancer_code, genie_cancer_code, cancer_type_label, rmsd, wrmsd, rollup=False):
        outpath = self.find_outpath(rollup)

        # Generate plot...
        print(f"Plotting mutation frequency results for TCGA cancer code: {tcga_cancer_code}")
        Plot.reset()
        Plot.mutation_frequencies(outpath, str(os.getenv('SYNAPSE_RELEASE_VERSION')), result, genie_cancer_code, tcga_cancer_code, cancer_type_label, rmsd, wrmsd)
        # Plots without TERT...
        try:
            result.drop(result[result['Hugo_Symbol'] == 'TERT'].index, inplace = True)
            Plot.reset()
            Plot.mutation_frequencies(f"{outpath}/no_tert", str(os.getenv('SYNAPSE_RELEASE_VERSION')), result, genie_cancer_code, tcga_cancer_code, cancer_type_label, rmsd, wrmsd)
        except Exception as e:
            print(e)
            print(f"Error in creating plot without TERT...")

    # Calculate mutation sample counts and generate comparison plots for all cancer codes
    # Limits analysis to > 5 mutations in either sample
    # Limits analysis to no more than 40 genes
    def sample_counts_by_gene(self, result, tcga_cancer_code, genie_cancer_code, cancer_type_label, rollup=False):
        outpath = self.find_outpath(rollup)
        indecies_to_drop = [] # captures the row index where mutation gene is not in panel  
        result.sort_values('tcga_gene_sample_count', inplace=True, ascending=False)
        for i, row in result.iterrows():
            tcga_count = row['tcga_gene_sample_count']
            genie_count = row['genie_gene_sample_count']

            if tcga_count <= 5 and genie_count <= 5:
                indecies_to_drop.append(i) # Add to drop list
                continue

            if i >= 40:
                indecies_to_drop.append(i) # Add to drop list
                continue

        result.drop(indecies_to_drop, inplace=True)
        # Generate plot...
        print(f"Plotting sample counts by gene results for TCGA cancer code: {tcga_cancer_code}")
        Plot.reset()
        Plot.sample_counts_by_gene(outpath, 
            str(os.getenv('SYNAPSE_RELEASE_VERSION')), 
            result, 
            genie_cancer_code, 
            cancer_type_label)

    # Plot total sample counts by cancer type for GENIE and TCGA cohorts
    def sample_counts_by_cancer_type(self, result, rollup=False):
        outpath = self.find_outpath(rollup)
        print(f"Plotting total sample count results...")
        Plot.reset()
        Plot.sample_counts_by_cancer_type(outpath, 
            str(os.getenv('SYNAPSE_RELEASE_VERSION')), 
            result)

    # Plot total sample counts by cancer type for GENIE and TCGA cohorts
    # Limit to top 40 cancer types
    def genie_oncotree_distribution(self, rollup=False):
        outpath = self.find_outpath(rollup)
        result = self.genie_analysis.sample_count_by_cancer_type(rollup)
        indecies_to_drop = [] # captures the row index where mutation gene is not in panel  
        result.sort_values('sample_count', inplace=True, ascending=False)
        for i, row in result.iterrows():
            if i >= 40:
                indecies_to_drop.append(i) # Add to drop list
                continue
        result.drop(indecies_to_drop, inplace=True)
        print(f"Plotting GENIE sample count bargraph...")
        Plot.reset()
        Plot.genie_oncotree_distribution(outpath, 
            str(os.getenv('SYNAPSE_RELEASE_VERSION')), 
            result)

    # Returns the output path depending on if a rollup oncotree code is used
    def find_outpath(self, rollup=True):
        if rollup:
            return f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/rollup"
        else:
            return f"/app/outputs/{os.getenv('SYNAPSE_RELEASE_VERSION')}/direct_comparison"

    # One-off cases that deviate from regular analysis pipeline
    def handle_deviations(self):
        # Create chart for UCEC rollup without UCS
        # Create charts splitting ductal and lobular breast
        # Create charts splitting uterine serous vs endometrioid uterine
        # uterine_serous_oncotree_code = "USC"
        # uterine_endometrioid_oncotree_code = "UEC"
        # lobular_breast_oncotree_code = "ILC"
        # ductal_breast_oncotree_code = "IDC"
        targets = [
            ["UCEC", "Uterine Corpus Endometrial Carcinoma", "UEC", False],
            ["UCEC", "Uterine Corpus Endometrial Carcinoma", "USC", False],
            ["BRCA", "Breast Invasive Carcinoma", "IDC", False],
            ["BRCA", "Breast Invasive Carcinoma", "ILC", False]
        ]

        for arr in targets:
            tcga_cancer_code = arr[0]
            cancer_type_label = arr[1]
            genie_cancer_code = arr[2]
            rollup = arr[3]
            outpath = f"{self.find_outpath(rollup)}/deviations"
            print(f"Outpath is {outpath}")

            genie_mutation_frequencies = self.genie_analysis.mutation_frequency_by_cancer_code(genie_cancer_code, rollup)
            tcga_mutation_frequencies = self.tcga_analysis.mutation_frequency_by_cancer_code(tcga_cancer_code)
            results = pd.merge(tcga_mutation_frequencies, genie_mutation_frequencies, on='Hugo_Symbol')
            results.to_csv(f"{outpath}/{genie_cancer_code}_results.tsv", sep='\t', index=True)
            results.drop(results[results['Hugo_Symbol'] == 'TERT'].index, inplace = True)
            Plot.reset()

            # Get sample counts for axes title
            try:
                tcga_total_sample_count = results.tcga_total_sample_count.unique().tolist()[0]
                genie_total_sample_count = results.genie_total_sample_count.unique().tolist()[0]
            except:
                tcga_total_sample_count = 'UNK'
                genie_total_sample_count = 'UNK'

            # Fetch mutation frequency results from df
            tcga_mut_freq = results.tcga_mut_freq.tolist()
            genie_mut_freq = results.genie_mut_freq.tolist()

            # General plot details
            sup_title_text = f"{cancer_type_label} Mutation Frequency"
            title_text = f"Oncotree: {genie_cancer_code} TCGA: {tcga_cancer_code}"
            plt.style.use('ggplot')
            plt.style.use('seaborn-whitegrid')
            plt.figure(figsize=(10,10))
            plt.xlim(0, 100)
            plt.ylim(0, 100)
            plt.title(title_text, fontsize='x-large')
            plt.suptitle(sup_title_text, fontsize='xx-large') 
            plt.xlabel(f"TCGA Mutation Frequency (%)\nN={tcga_total_sample_count}", fontsize='x-large')
            plt.ylabel(f"GENIE Mutation Frequency (%)\nN={genie_total_sample_count}", fontsize='x-large')

            # Dividing line
            div_x = [0,100]
            div_y = [0,100]
            plt.plot(div_x, div_y, color='lightgray')

            # Scatter
            plt.scatter(tcga_mut_freq, genie_mut_freq, marker='o', color='k', alpha=0.9)

            # Handle point labels
            spent_index = []
            annotations = []
            genes = []
            delta_dict = {}
            leftovers = {}
            for i, txt in enumerate(results.Hugo_Symbol):
                genes.append(txt)
                # Add a label if:
                #   a). TCGA or GENIE mut_freq meets min threshold
                #   a). Delta signifies cohort discrep
                tcga_mut_freq = results.tcga_mut_freq.iat[i]
                genie_mut_freq = results.genie_mut_freq.iat[i]
                delta = abs(genie_mut_freq - tcga_mut_freq)
                # Parse out secondary targets, and then leftovers
                if (genie_mut_freq >= 5 and tcga_mut_freq >= 5) or (genie_mut_freq >= 1 and delta > 10) or (tcga_mut_freq >= 1 and delta > 10):
                    delta_dict[i] = delta
                else:
                    leftovers[i] = delta

                # Handle outliers AKA primary targets
                if (tcga_mut_freq > 30 or 
                    genie_mut_freq > 30 or 
                    (genie_mut_freq > 20 and tcga_mut_freq > 10) or
                    (tcga_mut_freq > 20 and genie_mut_freq > 10) or
                    delta > 25):
                    # Old annotation method...
                    #   plt.annotate(txt, (results.mut_fraq_x.iat[i] + 1, results.mut_fraq_y.iat[i]))
                    # 
                    # New annotation method with adjust_text spacing...
                    spent_index.append(i)
                    annotations.append(plt.text(results.tcga_mut_freq.iat[i], results.genie_mut_freq.iat[i], txt))

            # Get a min number of labels based upon delta threshold
            label_min = 8
            axis_min = 3
            if (len(annotations) <= label_min and len(delta_dict.items()) > 0):
                index_list = []
                delta_list = []
                for key, value in delta_dict.items():
                    index_list.append(key)
                    delta_list.append(value)

                index_list, delta_list = zip(*sorted(zip(index_list, delta_list)))
                index_list, delta_list = (list(x) for x in zip(*sorted(zip(index_list, delta_list))))
                while (len(annotations) <= label_min):
                    target_i = index_list.pop(0)
                    if (target_i not in spent_index):
                        spent_index.append(i)
                        annotations.append(plt.text(results.tcga_mut_freq.iat[target_i], results.genie_mut_freq.iat[target_i], genes[target_i]))
                    if len(index_list) == 0:
                        break
                        
            # If you got to this point, all the points are hugging the axis.
            if (len(annotations) <= axis_min and len(leftovers.items()) > 0):
                index_list = []
                delta_list = []
                for key, value in leftovers.items():
                    index_list.append(key)
                    delta_list.append(value)

                index_list, delta_list = zip(*sorted(zip(index_list, delta_list)))
                index_list, delta_list = (list(x) for x in zip(*sorted(zip(index_list, delta_list))))
                while (len(annotations) <= axis_min):
                    target_i = index_list.pop(0)
                    if (target_i not in spent_index):
                        spent_index.append(i)
                        annotations.append(plt.text(results.tcga_mut_freq.iat[target_i], results.genie_mut_freq.iat[target_i], genes[target_i]))
                    if len(index_list) == 0:
                        break

            # Use adjustText to help point label overlap...
            # Arrows can be removed to just have text labels
            # Sample arrow code:
            #   arrowprops=dict(arrowstyle="-", color='dimgray', lw=0.5)
            adjust_text(annotations, arrowprops=dict(arrowstyle="-", color='dimgray', lw=0.5))

            # Generate plot
            plt.tight_layout()
            cancer_underscore = cancer_type_label.replace(" ", "_")
            plt.savefig(f"{outpath}/{genie_cancer_code}_{cancer_underscore}.png")