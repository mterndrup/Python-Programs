import logging
import sys
from datetime import datetime
from ete3 import NCBITaxa
import pandas as pd
import os
import requests

# ======================
# Logging setup
# ======================
log_name = f"resistance_ingest_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    filename=log_name,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
console = logging.StreamHandler(sys.stdout)
console.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

logging.info("==== Starting Resistance Ingestion Pipeline ====")

# ======================
# Helper Functions
# ======================
ncbi = NCBITaxa()  # Initialize ete3 NCBI interface

def classify_organism(name):
    """Classify organism type and dynamically fetch associated diseases using NCBI taxonomy + MeSH"""
    try:
        # Get NCBI taxid
        taxid = ncbi.get_name_translator([name]).get(name)
        if not taxid:
            logging.warning(f"Organism not found in NCBI: {name}")
            return "Unknown", []

        taxid = taxid[0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        kingdom_tid = next((tid for tid, rank in ranks.items() if rank == "superkingdom"), None)
        kingdom = ncbi.get_taxid_translator([kingdom_tid])[kingdom_tid] if kingdom_tid else "Unknown"

        # Fetch associated diseases dynamically from MeSH via NCBI E-utilities
        disease_list = fetch_diseases_from_ncbi(name)

        logging.info(f"Organism: {name} → Type: {kingdom}, Diseases: {', '.join(disease_list) if disease_list else 'Unknown'}")
        return kingdom, disease_list
    except Exception as e:
        logging.error(f"Error classifying organism {name}: {e}")
        return "Unknown", []

def fetch_diseases_from_ncbi(name):
    """Query NCBI PubMed/MeSH for diseases associated with the organism"""
    try:
        # Search MeSH via E-utilities
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "mesh",
            "term": f"{name}[Mesh]",
            "retmode": "json",
            "retmax": 5
        }
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        # Fetch MeSH terms for IDs
        url_summary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {"db": "mesh", "id": ",".join(ids), "retmode": "json"}
        resp = requests.get(url_summary, params=params, timeout=10)
        resp.raise_for_status()
        summary_data = resp.json()
        diseases = [summary_data["result"][i]["ds_name"] for i in ids if i in summary_data["result"]]
        return diseases
    except Exception as e:
        logging.warning(f"Could not fetch diseases for {name}: {e}")
        return []

def read_local_file():
    """Prompt user for local resistance/mutation data file and read it"""
    file_path = input("Enter path to local resistance/mutation data file (CSV/Excel): ").strip()
    if not os.path.exists(file_path):
        logging.error(f"File not found: {file_path}")
        return None
    try:
        if file_path.endswith(".csv"):
            df = pd.read_csv(file_path)
        elif file_path.endswith((".xls", ".xlsx")):
            df = pd.read_excel(file_path)
        else:
            logging.error("Unsupported file format. Please provide CSV or Excel.")
            return None
        logging.info(f"Successfully read file: {file_path}, {len(df)} records loaded")
        return df
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return None

def resistance_pipeline():
    """Main resistance ingestion workflow"""
    organism = input("Enter organism name: ").strip()
    org_type, diseases = classify_organism(organism)

    print(f"Organism type: {org_type}")
    print(f"Associated disease(s): {', '.join(diseases) if diseases else 'Unknown'}")

    # Prompt user for local resistance/mutation data
    df = read_local_file()
    if df is not None:
        print(f"Loaded {len(df)} records from local file for {organism}")
        logging.info(f"Loaded data for {organism}: {df.head(5).to_dict(orient='records')}")
    else:
        print("No valid data loaded.")

# ======================
# Run pipeline
# ======================
if __name__ == "__main__":
    resistance_pipeline()
