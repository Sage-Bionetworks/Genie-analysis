import os
import sys
import json
import pandas as pd
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import gc
from adjustText import adjust_text

# Generates plots
class Plot:

    @classmethod
    def reset(self):
        plt.cla() 
        plt.clf()
        plt.close('all')
        gc.collect()

    # Creatings a plot comparing mutation frequency between TCGA and GENIE samples
    @classmethod
    def mutation_frequencies(self, outpath, synapse_release_version, results, genie_cancer_code, tcga_cancer_code, cancer_type_label, rmsd, wrmsd):
        # Get sample counts for axes title
        try:
            tcga_total_sample_count = results.tcga_total_sample_count.unique().tolist()[0]
            genie_total_sample_count = results.genie_total_sample_count.unique().tolist()[0]
        except:
            tcga_total_sample_count = 'UNK'
            genie_total_sample_count = 'UNK'

        # Fetch mutation frequency results from df
        mut_freq_x = results.tcga_mut_freq.tolist()
        mut_freq_y = results.genie_mut_freq.tolist()

        import sys
        print( f"code: {tcga_cancer_code}", file=sys.stderr)
#        if genie_cancer_code == "LGG":
        if tcga_cancer_code == "LGG":
            title_text = f"Oncotree: LGGNOS,DIFG TCGA: {tcga_cancer_code}"
        else:
            title_text = f"Oncotree: {genie_cancer_code} TCGA: {tcga_cancer_code}"

        # General plot details
        sup_title_text = f"{cancer_type_label} Mutation Frequency"
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
        plt.scatter(mut_freq_x, mut_freq_y, marker='o', color='k', alpha=0.9)

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
            tcga_mut_fraq = results.tcga_mut_fraq.iat[i]
            genie_mut_fraq = results.genie_mut_fraq.iat[i]
            delta = abs(genie_mut_fraq - tcga_mut_fraq)
            # Parse out secondary targets, and then leftovers
            if (genie_mut_fraq >= 5 and tcga_mut_fraq >= 5) or (genie_mut_fraq >= 1 and delta > 10) or (tcga_mut_fraq >= 1 and delta > 10):
                delta_dict[i] = delta
            else:
                leftovers[i] = delta

            # Handle outliers AKA primary targets
            if (tcga_mut_fraq > 30 or 
                genie_mut_fraq > 30 or 
                (genie_mut_fraq > 20 and tcga_mut_fraq > 10) or
                (tcga_mut_fraq > 20 and genie_mut_fraq > 10) or
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
        # Just add in until axis_min
        # TODO: DRY up; pull this code into a helper function
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

        # Add in the rmsd label
        props = dict(boxstyle='Square', facecolor='None', edgecolor="dimgray", alpha=0.5)
        plt.text(5,95,f"RMSD: {rmsd}\nwRMSD: {wrmsd}",bbox=props)

        # Generate plot
        plt.tight_layout()
        cancer_underscore = cancer_type_label.replace(" ", "_")
        plt.savefig(f"{outpath}/mutation_frequencies/{cancer_underscore}.png")

    # Creates a plot comparing sample counts by gene between TCGA and GENIE samples
    @classmethod
    def sample_counts_by_gene(self, outpath, synapse_release_version, data, genie_cancer_code, cancer_type_label):
        plt.style.use('ggplot')
        plt.style.use('seaborn-whitegrid')
        plt.figure(figsize=(10,10))

        # Plot Cell Frequency Data
        x = data.Hugo_Symbol.tolist()
        y_tcga = data.tcga_gene_sample_count.tolist()
        y_genie = data.genie_gene_sample_count.tolist()

        # Plotting...
        tcga = plt.scatter(x, y_tcga, marker='o', color="tab:red", alpha=0.9, label="TCGA");  
        genie = plt.scatter(x, y_genie, marker='o', color="tab:blue", alpha=0.9, label="GENIE");  
        plt.title(F"{genie_cancer_code} Sample Counts")
        plt.xlabel("Gene")
        plt.ylabel("Sample Count")

        # ax.plot([0, 1], [0, 1], transform=ax.transAxes)
        plt.xticks(rotation=90)
        plt.legend(handles=[tcga, genie], loc='upper right', fontsize='xx-small', frameon=True)
        plt.savefig(f"{outpath}/sample_counts_by_gene/{genie_cancer_code}_sample_counts.png")

    # Creates a plot comparing sample counts by gene between TCGA and GENIE samples
    @classmethod
    def sample_counts_by_cancer_type(self, outpath, synapse_release_version, data):
        plt.style.use('ggplot')
        plt.style.use('seaborn-whitegrid')
        plt.tight_layout()
        plt.figure(figsize=(10, 10))

        # Plot Cell Frequency Data
        x = [] # Cancer type
        y_genie = []
        y_tcga = []
        for arr in data:
            x.append(arr[0])
            y_tcga.append(arr[1])
            y_genie.append(arr[2])

        # Plotting...
        tcga = plt.scatter(x, y_tcga, marker='o', color="tab:red", alpha=0.9, label="TCGA");  
        genie = plt.scatter(x, y_genie, marker='o', color="tab:blue", alpha=0.9, label="GENIE");  

        # Point Labels...
        for i, txt in enumerate(x):
            if abs(float(y_genie[i]) - float(y_tcga[i])) > 500: # Only show labels where discrepency is greater than 50 counts
                plt.annotate(y_genie[i], (x[i], y_genie[i]), xytext=(i+0.1, y_genie[i]), rotation=30)
                plt.annotate(y_tcga[i], (x[i], y_tcga[i]), xytext=(i+0.1, y_tcga[i]), rotation=30)

        plt.title(F"Sample Counts by Oncotree Code")
        plt.xlabel("Oncotree Code")
        plt.ylabel("Sample Count")
        plt.xticks(rotation=90)
        plt.legend(handles=[tcga, genie], loc='upper right', fontsize='xx-small', frameon=True)
        plt.gcf().subplots_adjust(bottom=0.3)
        plt.savefig(f"{outpath}/sample_counts_by_cancer_type/sample_counts_comparison.png")
        
    # Creates a barchart show distribution of cancer types in GENIE data
    @classmethod
    def genie_oncotree_distribution(self, outpath, synapse_release_version, data):
        plt.style.use('ggplot')
        plt.style.use('seaborn-whitegrid')

        x = data.oncotree_code.tolist()
        y_genie = data.sample_count.tolist()

        # Plotting...
        genie = plt.bar(x, y_genie, color="tab:blue");  
        
        plt.title(F"GENIE Sample Counts by Oncotree Code")
        plt.xlabel("Oncotree Code")
        plt.ylabel("Sample Count")
        plt.xticks(rotation=90)
        plt.savefig(f"{outpath}/sample_counts_by_cancer_type/genie_sample_counts.png")

    @classmethod
    def rmsd_by_gene(self, outpath, data, filename):
        plt.style.use('ggplot')
        plt.style.use('seaborn-whitegrid')
        plt.figure(figsize=(10,10))

        x = data.gene.tolist()[0:50]
        y_rmsd = data.rmsd.tolist()[0:50]
        y_wrmsd = data.wrmsd.tolist()[0:50]

        # Plotting...
        rmsd = plt.scatter(x, y_rmsd, marker='o', color="tab:red", alpha=0.9, label="RMSD");  
        wrmsd = plt.scatter(x, y_wrmsd, marker='o', color="tab:blue", alpha=0.9, label="wRMSD");  
        plt.title(F"Root Mean Square Deviation by Gene")
        plt.xlabel("Gene")
        plt.ylabel("RMSD")
        plt.xticks(rotation=90)
        plt.legend(handles=[rmsd, wrmsd], loc='upper right', fontsize='xx-small', frameon=True)
        plt.tight_layout()
        plt.savefig(f"{outpath}/rmsd_plots/{filename}.png")

    @classmethod
    def error_sum_by_gene(self, outpath, data):
        plt.style.use('ggplot')
        plt.style.use('seaborn-whitegrid')
        plt.figure(figsize=(10,10))

        x = data.gene.tolist()[0:50]
        y = data.error_sum.tolist()[0:50]

        # Plotting...
        error = plt.scatter(x, y, marker='o', color="tab:red", alpha=0.9, label="Error Sum");  
        plt.title(F"Aggregate Error Sum by Gene")
        plt.xlabel("Gene")
        plt.ylabel("Error Sum")
        plt.xticks(rotation=90)
        plt.legend(handles=[error], loc='upper right', fontsize='xx-small', frameon=True)
        plt.tight_layout()
        plt.savefig(f"{outpath}/rmsd_plots/error_sum_by_gene.png")
