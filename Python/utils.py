"""Using GENIE data utility functions"""
import synapseclient
from synapseclient import Synapse
import synapseutils


def get_available_releases(syn: Synapse, fileview_synid: str) -> dict:
    """Get available releases and their synapse ids

    Args:
        syn: Synapse connection
        fileview_synid: Synapse id of GENIE Release File View

    Returns:
        Mapping between release and Synapse id of its folder

    """
    releases = syn.tableQuery(f"select name, id from {fileview_synid} "
                              "where name <> 'case_lists'")
    releasesdf = releases.asDataFrame()
    release_map = {}
    for _, row in releasesdf.iterrows():
        print(row['name'])
        release_name = row['name'].split(" ")
        release = release_name[0]
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

    Args:
        syn: Synapse connection

    Returns:
        Synapse id of GENIE Release File View

    """
    current_user = syn.getUserProfile()
    # Determine if part of GENIE consortium
    try:
        syn.restGET(f"/team/3326313/member/{current_user.ownerId}")
        fileview_synid = "syn17019650"
    except synapseclient.core.exceptions.SynapseHTTPError:
        fileview_synid = "syn22233011"
    return fileview_synid


def get_all_genie_files(syn: Synapse, folderid: str,
                        download_location: str = None) -> dict:
    """Download all GENIE files for a specific release

    Args:
        syn: Synapse connection
        folderid: Synapse id of GENIE release folder
        download_location: Download

    Returns:
        {entity name: filepath}
    """
    files = synapseutils.syncFromSynapse(
        syn, folderid,
        downloadLocation=download_location, followLink=True
    )
    genie_file_map = {file_ent.name: file_ent.path for file_ent in files}
    return genie_file_map


def get_genie_file(syn: Synapse, folderid: str, filetype: str,
                   download_location: str = None) -> dict:
    """Download all GENIE files for a specific release

    Args:
        syn: Synapse connection
        folderid: Synapse id of GENIE release folder
        filetype: Download specific filetypes.
                  ['mutation', 'clinical', 'cna', 'fusion']
        download_location: Download

    Returns:
        {entity name: filepath}
    """
    allowed_filetypes = ['mutation', 'clinical', 'cna', 'fusion']
    if filetype not in allowed_filetypes:
        raise ValueError("filetype not one of {}".format(
            ", ".join(allowed_filetypes)
        ))
    files = syn.getChildren(folderid)
    genie_file_map = {}
    for file_info in files:
        filename = file_info['name']
        if filetype in filename and not filename.startswith("meta"):
            file_ent = syn.get(file_info['id'])
            genie_file_map[filename] = file_ent.path

    return genie_file_map


# TODO: can also download them as csv files for people to just read
def get_genie_bpc_tables(syn) -> dict:
    """Get BPC tables

    Returns:
        {
            'table name': {
                'id': 'syn1234',
                'form': 'form',
                'primary_key': 'key1, key2',
                'double_curated': False,
                'createdOn': 123,
                'modifiedOn': 234
            }
        }
    """
    table_view = syn.tableQuery("select * from syn21446696")
    table_viewdf = table_view.asDataFrame()
    table_viewdf.index = table_viewdf['name']
    del table_viewdf['name']
    return table_viewdf.to_dict('index')
