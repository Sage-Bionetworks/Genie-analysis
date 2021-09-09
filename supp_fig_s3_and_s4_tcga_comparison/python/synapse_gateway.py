import os
import sys
from pathlib import Path
from shutil import copyfile
import synapseclient 
import synapseutils
from dotenv import load_dotenv

# Push and pull from synapse
class SynapseGateway:
	def __init__(self):
	        load_dotenv(dotenv_path='/app/.env', verbose=True)

	        self.synapse_username = str(os.getenv('SYNAPSE_USERNAME'))
	        self.synapse_password = str(os.getenv('SYNAPSE_PASSWORD'))

	def login(self):
	    try:
	    	self.syn = synapseclient.Synapse()
	    	self.syn.login(self.synapse_username, self.synapse_password)
	    except Exception as err:
	        print("Unable to login to Synapse...")
	        print(err)

	def fetch_release(self, release):
		path = f"/app/releases/{release}"
		synapseutils.syncFromSynapse(self.syn, release, followLink=True, path=path) 

		# Check for 'PHS-TRISEQ-V2'
		target = Path(f"{path}/data_gene_panel_PHS-TRISEQ-V2.txt")
		if not target.is_file():
			copyfile("/app/references/data_gene_panel_PHS-TRISEQ-V2.txt", f"{path}/data_gene_panel_PHS-TRISEQ-V2.txt")
			

