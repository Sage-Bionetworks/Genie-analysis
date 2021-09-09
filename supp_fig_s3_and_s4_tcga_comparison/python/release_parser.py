import os
import sys
import json
import glob
import pandas as pd
from dotenv import load_dotenv

# Parses txt files from consortium releases into usable data structures
class ReleaseParser:
    def __init__(self, release_id=None, release_version=None):
        load_dotenv(dotenv_path='/app/.env', verbose=True)
        self.release_id = release_id #ex. synXXX
        self.release_version = release_version #ex. 10.2
        self.release_root = f"/app/releases/{release_id}"

        if self.release_id == None:
            self.release_id = os.getenv("SYNAPSE_RELEASE_ID")
        if self.release_version == None:
            self.release_version = os.getenv("SYNAPSE_RELEASE_VERSION")

        # filepaths
        self.release_dir = f"/app/releases/{self.release_id}"
        self.mutations_path = f"/app/releases/{self.release_id}/data_mutations_extended_{self.release_version}.txt"
        self.samples_path = f"/app/releases/{self.release_id}/data_clinical_sample_{self.release_version}.txt"
        self.patients_path = f"/app/releases/{self.release_id}/data_clinical_patient_{self.release_version}.txt"

        self.parse_all()

    # Help function to parse all files into dfs
    def parse_all(self):
        self.parse_oncotree_codes()
        self.parse_cancer_codes()
        self.parse_mutations()
        self.parse_patients()
        self.parse_samples()
        self.parse_panel_genes()
        self.create_sample_panel_dict() # dependent on samples df

    # Returns a json object of cancer codes in TCGA and GENIE 
    def parse_cancer_codes(self):
        print("Parsing cancer codes...")
        with open(os.getenv("CANCER_CODES_PATH")) as cc:
            self.cancer_codes = json.load(cc)
            return self.cancer_codes

    # Returns a dict of cancer codes
    # Key code
    # Value top level code
    # level_1	level_2	level_3	level_4	level_5	level_6	level_7	metamaintype	metacolor	metanci	metaumls	history
    def parse_oncotree_codes(self):
        print("Parsing oncotree codes...")
        self.oncotree_codes = {}
        with open(os.getenv("CANCER_CODES_PATH")) as cc:
            targets = json.load(cc)['oncotree_rollup']
            for target in targets:
                genie_code = target['genie_cancer_code']
                rollup_codes = target['rollup_codes']
                for code in rollup_codes:
                    self.oncotree_codes[code] = genie_code
        return self.oncotree_codes

    # Returns a df of mutations in the release
    def parse_mutations(self):
        print("Parsing mutations...")
        self.mutations_df = pd.read_csv(self.mutations_path, sep='\t', comment="#")
        return self.mutations_df

    # Returns a df of samples in the release
    def parse_samples(self):
        print("Parsing samples...")
        self.samples_df = pd.read_csv(self.samples_path, sep='\t', comment="#")
        rollup_codes = []
        for i, row in self.samples_df.iterrows():
            try:
                oncotree_code = row['ONCOTREE_CODE']
                rollup_code = self.oncotree_codes[oncotree_code]
            except:
                rollup_code = None

            rollup_codes.append(rollup_code)
        self.samples_df["ROLLUP_ONCOTREE_CODE"] = rollup_codes
        return self.samples_df

    # Returns a df of patients in the release
    def parse_patients(self):
        print("Parsing patients...")
        self.patients_df = pd.read_csv(self.patients_path, sep='\t', comment="#")
        return self.patients_df
    
    # Returns a dictionary of gene panels
    # Key = stable_id
    # Value = gene_list
    def parse_panel_genes(self):
        print("Parsing panel genes...")
        self.panel_genes = {}
        panel_files = glob.glob(f"{self.release_dir}/data_gene_panel*.txt")
        for panel_file in panel_files:
            with open(panel_file) as f:
                panel = ""
                genes = []
                for line in f:
                    # keys are 'stable_id', 'description', and 'gene_list'
                    spl = line.split(":")
                    if spl[0] == 'stable_id':
                        panel = spl[1].strip().upper()

                    if spl[0] == 'gene_list':
                        genes = spl[1].strip().upper().split("\t")

                self.panel_genes[panel] = set(genes)
        return self.panel_genes

    # Returns a dict of samples and their associated panels
    # Key = SAMPLE_ID
    # Value = SEQ_ASSAY_ID
    def create_sample_panel_dict(self):
        print("Parsing sample panels...")
        self.sample_panels = {}
        for i, row in self.samples_df.iterrows():
            self.sample_panels[row['SAMPLE_ID']] = row['SEQ_ASSAY_ID']
        return self.sample_panels

    def save_parsed_dataframes(self):
        print("Saving parsed dataframes and dicts to file...")
        df_list = [
            ["mutations_df", self.mutations_df],
            ["samples_df", self.samples_df],
            ["patients_df", self.patients_df]
        ]
    
        os.makedirs(f"/app/outputs/{self.release_version}", exist_ok=True)
        
        for arr in df_list:
            name = arr[0]
            df = arr[1]
            df.to_csv(f"/app/outputs/{self.release_version}/{name}.tsv", sep='\t', index=False)
