# Description: Create plots by center and oncotree code for tumor-mutation burden
#               distributions for 100k GENIE samples.  
# Author: Haley Hunter-Zinck
# Date: 2022-01-14

# pre-setup  ---------------------------

library(optparse)

waitifnot <- function(cond, msg) {
  if (!cond) {
    
    for (str in msg) {
      message(str)
    }
    message("Press control-C to exit and try again.")
    
    while(T) {}
  }
}

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of input file"),
  make_option(c("-o", "--synid_folder_output"), type = "character",
              help="Synapse ID of output folder"),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Output script messages to the user.")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input) && !is.null(opt$synid_folder_output),
          msg = "Rscript template.R -h")

synid_file_input <- opt$synid_file_input
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(ggplot2)
library(synapser)
synLogin()

# files
synid_table_sample <- "syn7517674"
outplot <- "tmb.pdf"

# constants
BINS <- c("Low (<2)", "Mid (2-16)", "High (>16)")

# functions ----------------------------

#' Download and load data stored in csv or other delimited format on Synapse
#' into an R data frame.
#' 
#' @param synapse_id Synapse ID
#' @version Version of the Synapse entity to download.  NA will load current
#' version
#' @param set Delimiter for file
#' @param na.strings Vector of strings to be read in as NA values
#' @param header TRUE if the file contains a header row; FALSE otherwise.
#' @param check_names TRUE if column names should be modified for compatibility 
#' with R upon reading; FALSE otherwise.
#' @return data frame
get_synapse_entity_data_in_csv <- function(synapse_id, 
                                           version = NA,
                                           sep = ",", 
                                           na.strings = c("NA"), 
                                           header = T,
                                           check_names = F) {
  
  if (is.na(version)) {
    entity <- synGet(synapse_id)
  } else {
    entity <- synGet(synapse_id, version = version)
  }
  
  data <- read.csv(entity$path, stringsAsFactors = F, 
                   na.strings = na.strings, sep = sep, check.names = check_names,
                   header = header)
  return(data)
}

#' Store a file on Synapse with options to define provenance.
#' 
#' @param path Path to the file on the local machine.
#' @param parent_id Synapse ID of the folder or project to which to load the file.
#' @param file_name Name of the Synapse entity once loaded
#' @param prov_name Provenance short description title
#' @param prov_desc Provenance long description
#' @param prov_used Vector of Synapse IDs of data used to create the current
#' file to be loaded.
#' @param prov_exec String representing URL to script used to create the file.
#' @return Synapse ID of entity representing file
save_to_synapse <- function(path, 
                            parent_id, 
                            file_name = NA, 
                            prov_name = NA, 
                            prov_desc = NA, 
                            prov_used = NA, 
                            prov_exec = NA) {
  
  if (is.na(file_name)) {
    file_name = path
  } 
  file <- File(path = path, parentId = parent_id, name = file_name)
  
  if (!is.na(prov_name) || !is.na(prov_desc) || !is.na(prov_used) || !is.na(prov_exec)) {
    act <- Activity(name = prov_name,
                    description = prov_desc,
                    used = prov_used,
                    executed = prov_exec)
    file <- synStore(file, activity = act)
  } else {
    file <- synStore(file)
  }
  
  return(file$properties$id)
}

# read ----------------------------

data <- get_synapse_entity_data_in_csv(synid_file_input, sep = "\t")

query <- glue("SELECT SAMPLE_ID, ONCOTREE_CODE FROM {synid_table_sample}")
sam <- as.data.frame(synTableQuery(query, includeRowIdAndRowVersion = F))

# main ----------------------------

# join TMB and sample info
mut <- data %>% 
  mutate(center = unlist(lapply(strsplit(SAMPLE_ID, split = "-"), function(x) {return(x[[2]])}))) %>%
  left_join(sam, by = "SAMPLE_ID")

# TMB bin by center
df_center <- mut %>%
  group_by(tmb_bin, center) %>%
  count()

# TMB bin by code
common_code <- as.character(unlist(mut %>% 
  group_by(ONCOTREE_CODE) %>%
  count() %>% 
  filter(n > 1000 & !is.na(ONCOTREE_CODE)) %>%
  select(ONCOTREE_CODE)))
df_cancer <- mut %>%
  filter(is.element(ONCOTREE_CODE, common_code)) %>%
  group_by(tmb_bin, ONCOTREE_CODE) %>%
  count()

# plots -----------------------------

pdf(outplot)

# plot TMB bin by center
ggplot(df_center, aes(fill = factor(tmb_bin, levels = BINS), y = n, x = factor(center, levels = rev(sort(unique(center)))))) +
  geom_bar(position = "fill", stat = "identity") +
  labs(fill = "TMB bin") +
  ylab("Fraction of samples") + 
  xlab("Center") + 
  coord_flip()

# plot TMB bin by cancer
ggplot(df_cancer, aes(fill = factor(tmb_bin, levels = BINS), y = n, x = factor(ONCOTREE_CODE, levels = rev(sort(unique(ONCOTREE_CODE)))))) +
  geom_bar(position = "fill", stat = "identity") +
  labs(fill = "TMB bin") +
  ylab("Fraction of samples") + 
  xlab("OncoTree code") + 
  coord_flip()

graphics.off()

# write -------------------------------

synid_file_output <- save_to_synapse(path = outplot, 
                parent_id = synid_folder_output,
                prov_name = "tmb plots", 
                prov_desc = "plot tmb for genie release 9.1-public by center and oncotree code", 
                prov_used = synid_file_input, 
                prov_exec = "https://github.com/hhunterzinck/genie_requests/blob/main/2022-01-14_100k_tmb_by_site.R")

file.remove(outplot)

# close out ----------------------------

print(glue("Plots loaded to {synid_file_output} as '{outplot}'."))

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
