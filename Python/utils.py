import synapseclient
from synapseclient import Synapse
import synapseutils


def get_available_releases(syn: Synapse, fileview_synid: str) -> dict:
    """Get available releases and their synapse ids

    Returns:
        Mapping between release and Synapse id of its folder
    """
    releases = syn.tableQuery(f"select name, id from {fileview_synid}")
    releasesdf = releases.asDataFrame()
    release_map = {}
    for _, row in releasesdf.iterrows():
        print(row['name'])
        release_name = row['name'].split(" ")
        release = release_name[0]
        if row['name'] != "case_lists":
            # This logic is because the consortium release folders are
            # in 'Release X' folder names.
            if len(release_name) > 1:
                release = release_name[1]
            if "." in release or "-" in release:
                release_map[release] = row['id']
    print("Choose from these releases: {}".format(
        ", ".join(release_map.keys())
    ))
    return release_map


def get_genie_fileview_synid(syn: Synapse) -> str:
    """Depending on the Synapse user, get the genie file view.
    If part of the GENIE consortium release provide GENIE internal file view
    """
    current_user = syn.getUserProfile()
    # Determine if part of GENIE consortium
    try:
        syn.restGET(f"/team/3326313/member/{current_user.ownerId}")
        fileview_synid = "syn17019650"
    except synapseclient.core.exceptions.SynapseHTTPError:
        fileview_synid = "syn22233011"
    return fileview_synid


def get_genie_files(syn, folderid, download_location=None) -> dict:
    """Download all GENIE files for a specific release

    Returns:
        {filename: filepath}
    """
    files = synapseutils.syncFromSynapse(
        syn, folderid,
        downloadLocation=download_location, followLink=True
    )
    genie_file_map = {file_ent.name: file_ent.path for file_ent in files}
    return genie_file_map
