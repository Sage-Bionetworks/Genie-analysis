# AACR/GENIE

The AACR GENIE project is a publically accessible curation of cancer data. This codebase serves as a means for PSJH to query and analyze published GENIE data, and submit it's genomics own data to the [project](https://www.aacr.org/professionals/research/aacr-project-genie/).

### Docker Build

Execute, "./docker-build.sh VERSION" to build the image locally.  

### ENV

```
## General
CANCER_CODES_PATH=/app/references/cancer_codes.json
ONCOTREE_CODES_PATH=/app/references/oncotree_cancer_codes.tsv
PCAWG_SUPPLEMENTARY_TABLE=/app/releases/pcawg/supplementary_table.csv

## GCP
GBQ_KEY_PATH=/path/to/key.json
SERVICE_GCP_PROJECT=

## Synapse
SYNAPSE_USERNAME=
SYNAPSE_PASSWORD=
SYNAPSE_RELEASE_ID=syn25451657
SYNAPSE_RELEASE_VERSION=9.1-public
```

### Pull a release from Synapse
```
./docker-entrypoint VERSION
cd /app/python
python
>>> from synapse_gateway import SynapseGateway
>>> x = SynapseGateway()
>>> x.login()
Welcome, galvinma!

>>> x.fetch_release("syn25451657")
Downloading  [####################]100.00%   1.2MB/1.2MB (1.9MB/s) 9.1-public.html Done...
Downloading  [####################]100.00%   17.8kB/17.8kB (178.8kB/s) assay_information_9.1-public.txt Done...
Downloading  [####################]100.00%   5.7MB/5.7MB (6.0MB/s) data_clinical_patient_9.1-public.txt Done...
```

### TCGA vs. GENIE Mutation Frequency Analysis

To generate mutation frequency TSVs and plots execute the following:

```
./docker-entrypoint VERSION
cd /app/python
python
from tcga_genie_comparison import TcgaGenieComparison
x = TcgaGenieComparison()
x.exeute()
```

### References
- Wiki: https://github.com/EACRI/biocoor/wiki/GENIE
- Synapse: https://www.synapse.org/
- Synapse GENIE Data Portal: https://www.synapse.org/#!Synapse:syn7222066/files/
- PHS GENIE repository: https://www.synapse.org/#!Synapse:syn11703429
- GENIE Data Issue Tracking: https://github.com/Sage-Bionetworks/GENIE-Data-QC-Tracking
