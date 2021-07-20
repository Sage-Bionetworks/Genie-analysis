"""Example analysis code
Must have synapseclient installed
`pip install synapseclient`
Highly recommend pandas when working with dataframes,
but it is not necessary `pip install pandas`
"""
import pandas
import synapseclient


def main():
    # Need to be logged into synapse. Set up ~/.synapseConfig
    syn = synapseclient.login()

    # Navigate to GENIE site and get release 8.0 clinical sample
    # Synapse Id (synapse.org/genie)
    sample_ent = syn.get("syn22228695")
    print(sample_ent.path)

    # If you have pandas installed:
    sampledf = pandas.read_csv(sample_ent.path,
                               sep="\t", comment="#")
    # obtain oncotree code distribution
    oncotree_counts = sampledf.ONCOTREE_CODE.value_counts()

    # Write files out
    oncotree_counts.to_csv("oncotree_code_dist.csv",
                           header=['count'])

if __name__ == "__main__":
    main()
