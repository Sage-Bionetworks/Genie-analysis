import os
import sys
from google.cloud import bigquery
from google.oauth2 import service_account
from dotenv import load_dotenv

# Push and pull from TCGA GBQ
class TcgaGateway:
    def __init__(self):
        load_dotenv(dotenv_path='/app/.env', verbose=True)
        self.login()

    def login(self):
        try:
            credentials = service_account.Credentials.from_service_account_file(
                str(os.getenv('GBQ_KEY_PATH')), scopes=["https://www.googleapis.com/auth/cloud-platform"],
            )
            self.client = bigquery.Client(credentials=credentials, project=os.getenv('SERVICE_GCP_PROJECT'))
        except Exception as err:
            print("GBQ auth failed...")
            print(err)

    def job(self, query):
        return self.client.query(query)