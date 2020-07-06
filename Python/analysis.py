"""Example analysis code"""
import synapseclient

import utils

RELEASE = "8.0-public"

def main():
    # Need to be logged into synapse. Set up ~/.synapseConfig
    syn = synapseclient.login()

    # Do not modify these lines
    fileview_synid = utils.get_genie_fileview_synid(syn)
    releases = utils.get_available_releases(syn, fileview_synid)
    if RELEASE in releases:
        genie_file_map = utils.get_all_genie_files(syn, releases[RELEASE])
    else:
        raise ValueError(f"Release does not exist: {RELEASE}")

    # Begin your analysis

    # Write files out

if __name__ == "__main__":
    main()
