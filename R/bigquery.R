# This script is an example of how to interact with a GENIE BigQuery
# Table.  There are other more advanced SQL commands that you can
# leverage.
library(bigrquery)

bq_auth()

# I think you can use this billing project, but am unsure
billing = "project-genie-query-prod"

# 4 tables available
# project-genie-query-prod.consortium.genomic_information
# project-genie-query-prod.consortium.mutation
# project-genie-query-prod.consortium.patient
# project-genie-query-prod.consortium.sample

sql <- "SELECT distinct(SEQ_ASSAY_ID) FROM `project-genie-query-prod.consortium.sample`"

tb <- bq_project_query(billing, sql)

# bq_table_download(tb, max_results = 10)
seq_assay_ids = bq_table_download(tb)


