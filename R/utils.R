library(glue)
library(synapser)
library(synapserutils)
synLogin()

#' Get available releases and their synapse ids
#' 
#' @param fileview_synid: Synapse id of GENIE Release File View
#' @return Mapping between release and Synapse id of its folder
#' @examples
#' get_available_releases("syn22233011")
get_available_releases <- function(fileview_synid) {
  releases = synTableQuery(glue::glue("select name, id from {fileview_synid} where name <> 'case_lists'"))
  releasesdf = releases$asDataFrame()
  release_mapping = apply(releasesdf, 1, function(each_row) {
    release_name = unlist(strsplit(each_row['name'], " "))
    release_map = list()
    release = release_name[1]
    if (length(release_name) > 1) {
      release = release_name[2]
    }
    if (grepl("[.]", release) || grepl("-", release)) {
      release_map[release] = each_row['id']
    }
    release_map
  })
  unlist(release_mapping)
}


#' Depending on the Synapse user, get the genie file view.
#' If part of the GENIE consortium release provide GENIE internal file view
#' 
#' @return Synapse id of GENIE Release File View
#' @examples
#' get_genie_fileview_synid()
get_genie_fileview_synid <- function() {
  current_user = synGetUserProfile()$ownerId
  current_user = 5
  fileview_synid = tryCatch({
      # this rest call will throw an error if member is not part of the team
      synRestGET(glue::glue("/team/3326313/member/{current_user}"))
      fileview_synid = "syn17019650"
    }, error = function(e) {
      fileview_synid = "syn22233011"
    }
  )
  fileview_synid
}


#' Download all GENIE files for a specific release
#' 
#' @param folderid Synapse id of GENIE release folder
#' @param download_location Download path

#' @return c(entity name = filepath)
#' @examples
#' get_all_genie_files("syn7844527")
#' get_all_genie_files("syn7844527", "./")
get_all_genie_files <- function(folderid) {
  files = syncFromSynapse(folderid, followLink = T)
  genie_file_map = c()
  for (file_ent in files) {
    genie_file_map[file_ent$properties$name] = file_ent$path
  }
  genie_file_map
}


#' Download GENIE files of specific filetypes from a specific release
#'
#' @param folderid Synapse id of GENIE release folder
#' @param filetype Download specific filetypes. c('mutation', 'clinical', 'cna', 'fusion')
#' @param download_location Download path
#' @return c(entity name = filepath)
#' @examples
#' get_genie_file("syn7844527", "clinical")
#' get_all_genie_files("syn7844527", "clinical", download_location=./")
get_genie_file <- function(folderid, filetype,
                           download_location = NULL) {
  allowed_filetypes = c('mutation', 'clinical', 'cna', 'fusion')
  if (!filetype %in% allowed_filetypes) {
    filetype_str = paste(allowed_filetypes, collapse = ", ")
    stop(glue::glue("filetype not one of {filetype_str}"))
  }
  files = as.list(synGetChildren(folderid))
  genie_file_map = c()
  for (file_info in files) {
    filename = file_info$name
    if (grepl(filetype, filename) & (!startsWith(filename, "meta"))) {
      file_ent = synGet(file_info$id,
                        downloadLocation = download_location)
      genie_file_map[filename] = file_ent$path
    }
  }
  genie_file_map
}

