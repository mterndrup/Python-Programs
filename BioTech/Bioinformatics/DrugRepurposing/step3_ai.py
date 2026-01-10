import logging
import sys
from datetime import datetime
from ete3 import NCBITaxa
import pandas as pd
import os
import requests
import re
from collections import Counter
import xml.etree.ElementTree as ET

# ======================
# Logging setup
# ======================
log_name = f"resistance_ingest_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    filename=log_name,
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
console = logging.StreamHandler(sys.stdout)
console.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

logging.info("==== Starting Resistance Ingestion Pipeline ====")

# ======================
# Initialize NCBITaxa
# ======================
ncbi = NCBITaxa()


# ======================
# Helper Functions
# ======================
def classify_organism(name):
    """Classify organism type and fetch associated diseases and biological processes."""
    try:
        taxid_dict = ncbi.get_name_translator([name])
        if name not in taxid_dict:
            logging.warning(f"Organism not found in NCBI: {name}")
            return "Unknown", [], []

        taxid = taxid_dict[name][0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)

        # Try superkingdom → phylum → class as fallback
        for rank in ["superkingdom", "phylum", "class"]:
            rank_tid = next((tid for tid, r in ranks.items() if r == rank), None)
            if rank_tid:
                kingdom = ncbi.get_taxid_translator([rank_tid])[rank_tid]
                break
        else:
            kingdom = "Unknown"

        # Fetch diseases and biological processes dynamically
        diseases, biological_processes = get_diseases_and_processes_for_organism(name)

        # Log separately
        logging.info(f"Organism: {name} → Type: {kingdom}")
        logging.info(f"Diseases: {', '.join(diseases) if diseases else 'Unknown'}")
        logging.info(f"Biological Processes: {', '.join(biological_processes) if biological_processes else 'Unknown'}")

        return kingdom, diseases, biological_processes

    except Exception as e:
        logging.error(f"Error classifying organism {name}: {e}")
        return "Unknown", [], []


def get_diseases_and_processes_for_organism(name):
    """Get diseases and biological processes associated with organism using multiple dynamic approaches."""
    all_findings = []

    # Method 1: Search PubMed for disease associations
    pubmed_findings = search_pubmed_for_diseases(name)
    if pubmed_findings:
        all_findings.extend(pubmed_findings)
        logging.info(f"Found {len(pubmed_findings)} findings from PubMed for {name}")

    # Method 2: Search NCBI Taxonomy for pathogenicity information
    taxonomy_findings = get_diseases_from_taxonomy(name)
    if taxonomy_findings:
        all_findings.extend(taxonomy_findings)
        logging.info(f"Found {len(taxonomy_findings)} findings from NCBI Taxonomy for {name}")

    # Method 3: Search MeSH terms more effectively
    mesh_findings = search_mesh_for_diseases(name)
    if mesh_findings:
        all_findings.extend(mesh_findings)
        logging.info(f"Found {len(mesh_findings)} findings from MeSH for {name}")

    # Categorize into diseases vs biological processes
    diseases = []
    biological_processes = []

    if all_findings:
        # Define biological processes
        biological_process_terms = {
            'necrosis', 'apoptosis', 'phagocytosis', 'mitosis', 'meiosis',
            'necroptosis', 'pyroptosis', 'ferroptosis', 'autophagy',
            'transcription', 'translation', 'replication', 'metabolism',
            'pathogenesis', 'morphogenesis', 'angiogenesis', 'carcinogenesis',
            'fibrosis', 'sclerosis', 'hyperplasia', 'dysplasia', 'metaplasia',
            'proliferation', 'differentiation', 'activation', 'inhibition',
            'signaling', 'pathway', 'mechanism', 'process', 'response',
            'expression', 'synthesis', 'degradation', 'transport', 'secretion',
            'adhesion', 'migration', 'invasion', 'cytotoxicity', 'lysis',
            'phosphorylation', 'ubiquitination', 'glycosylation', 'acetylation'
        }

        # Define disease indicators
        disease_indicators = {
            'infection', 'disease', 'syndrome', 'disorder', 'condition',
            'tuberculosis', 'pneumonia', 'sepsis', 'meningitis', 'endocarditis',
            'gastroenteritis', 'diarrhea', 'dysentery', 'malaria', 'cholera',
            'aids', 'hiv', 'uti', 'mrsa', 'clostridium', 'salmonella'
        }

        # Valid disease suffixes
        disease_suffixes = ['itis', 'osis', 'emia', 'pathy']

        # Count frequencies for ranking
        finding_counts = Counter([f.lower().strip() for f in all_findings])

        # Categorize each finding
        for finding in all_findings:
            finding_clean = finding.lower().strip()

            # Skip very short or generic terms
            if len(finding_clean) <= 3:
                continue

            # Check if it's a biological process
            is_biological_process = any(process in finding_clean for process in biological_process_terms)

            # Check if it's a disease
            is_disease = (
                    any(indicator in finding_clean for indicator in disease_indicators) or
                    (any(finding_clean.endswith(suffix) for suffix in disease_suffixes) and
                     not any(exclude in finding_clean for exclude in ['ptosis', 'lysis', 'genesis'])) or
                    finding_clean in ['tb', 'aids', 'hiv', 'uti', 'mrsa', 'copd']
            )

            # Categorize
            if is_biological_process:
                if finding_clean not in [p.lower() for p in biological_processes]:
                    biological_processes.append(finding.title())
            elif is_disease:
                if finding_clean not in [d.lower() for d in diseases]:
                    diseases.append(finding.title())

        # Return top 5 of each category
        return diseases[:5], biological_processes[:5]

    return [], []


def search_pubmed_for_diseases(organism_name):
    """Search PubMed abstracts for disease associations."""
    diseases = []
    try:
        # Multiple search strategies
        search_queries = [
            f'("{organism_name}"[Title/Abstract]) AND (pathogen OR disease OR infection OR clinical)',
            f'("{organism_name}"[Title/Abstract]) AND (syndrome OR illness OR disorder)',
            f'("{organism_name}"[Title/Abstract]) AND (treatment OR therapy OR antimicrobial)',
        ]

        for query in search_queries:
            pmids = search_pubmed(query, max_results=15)
            if pmids:
                abstracts = fetch_pubmed_abstracts(pmids)
                extracted_diseases = extract_diseases_from_abstracts(abstracts, organism_name)
                diseases.extend(extracted_diseases)

                # Also try to get MeSH terms from these articles
                mesh_terms = get_mesh_terms_from_pmids(pmids)
                diseases.extend(mesh_terms)

        return diseases

    except Exception as e:
        logging.warning(f"Error searching PubMed for {organism_name}: {e}")
        return []


def search_pubmed(query, max_results=15):
    """Search PubMed and return PMIDs."""
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": max_results,
            "sort": "relevance",
            "datetype": "pdat",
            "reldate": "1825"  # Last 5 years for more recent associations
        }

        response = requests.get(url, params=params, timeout=20)
        response.raise_for_status()

        data = response.json()
        pmids = data.get("esearchresult", {}).get("idlist", [])
        return pmids

    except Exception as e:
        logging.error(f"Error searching PubMed with query '{query}': {e}")
        return []


def fetch_pubmed_abstracts(pmids):
    """Fetch abstracts for given PMIDs."""
    try:
        if not pmids:
            return ""

        # Limit to prevent timeout
        pmids = pmids[:8]

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "rettype": "abstract",
            "retmode": "text"
        }

        response = requests.get(url, params=params, timeout=25)
        response.raise_for_status()

        return response.text

    except Exception as e:
        logging.error(f"Error fetching abstracts: {e}")
        return ""


def get_mesh_terms_from_pmids(pmids):
    """Extract MeSH terms from PubMed articles."""
    mesh_diseases = []
    try:
        if not pmids:
            return []

        # Limit to prevent timeout
        pmids = pmids[:5]

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml"
        }

        response = requests.get(url, params=params, timeout=20)
        response.raise_for_status()

        root = ET.fromstring(response.text)

        # Extract MeSH headings
        for mesh_heading in root.findall(".//MeshHeading"):
            descriptor = mesh_heading.find("DescriptorName")
            if descriptor is not None:
                mesh_term = descriptor.text
                # Filter for disease-related terms
                if any(keyword in mesh_term.lower() for keyword in
                       ['infection', 'disease', 'syndrome', 'disorder', 'osis', 'itis', 'emia']):
                    mesh_diseases.append(mesh_term)

        return mesh_diseases

    except Exception as e:
        logging.warning(f"Error extracting MeSH terms: {e}")
        return []


def get_diseases_from_taxonomy(organism_name):
    """Try to get disease information from NCBI Taxonomy comments."""
    try:
        taxid_dict = ncbi.get_name_translator([organism_name])
        if organism_name not in taxid_dict:
            return []

        taxid = taxid_dict[organism_name][0]

        # Search NCBI Taxonomy database for pathogenicity info
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "taxonomy",
            "id": str(taxid),
            "retmode": "xml"
        }

        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()

        root = ET.fromstring(response.text)
        diseases = []

        # Look for disease mentions in various fields
        for prop in root.findall(".//Property"):
            prop_name = prop.find("PropName")
            prop_value = prop.find("PropValueString")

            if prop_name is not None and prop_value is not None:
                if "disease" in prop_name.text.lower() or "pathogen" in prop_name.text.lower():
                    diseases.append(prop_value.text)

        # Also check comments
        for comment in root.findall(".//Comment"):
            if comment.text:
                extracted = extract_diseases_from_text(comment.text)
                diseases.extend(extracted)

        return diseases

    except Exception as e:
        logging.warning(f"Error getting taxonomy diseases for {organism_name}: {e}")
        return []


def search_mesh_for_diseases(organism_name):
    """Search MeSH database more effectively for disease associations."""
    try:
        # Search for the organism in MeSH
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "mesh",
            "term": f'"{organism_name}"',
            "retmode": "json",
            "retmax": 10
        }

        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()

        data = response.json()
        mesh_ids = data.get("esearchresult", {}).get("idlist", [])

        if not mesh_ids:
            return []

        # Get detailed info for MeSH terms
        url_fetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "mesh",
            "id": ",".join(mesh_ids),
            "retmode": "xml"
        }

        response = requests.get(url_fetch, params=params, timeout=15)
        response.raise_for_status()

        # Safely parse XML with error handling
        try:
            root = ET.fromstring(response.text)
        except ET.ParseError as parse_error:
            logging.warning(f"XML parsing error for MeSH data: {parse_error}")
            return []

        diseases = []

        # Extract related disease terms
        for descriptor in root.findall(".//DescriptorRecord"):
            # Look at concept relations and see-also references
            for see_related in descriptor.findall(".//SeeRelatedList/SeeRelatedDescriptor"):
                related_term = see_related.find("DescriptorName/String")
                if related_term is not None and any(keyword in related_term.text.lower()
                                                    for keyword in
                                                    ['infection', 'disease', 'syndrome', 'tuberculosis', 'pneumonia']):
                    # Filter out biological processes
                    if not any(exclude in related_term.text.lower() for exclude in
                               ['diagnosis', 'therapy', 'genetics', 'metabolism', 'physiology']):
                        diseases.append(related_term.text)

        return diseases

    except Exception as e:
        logging.warning(f"Error searching MeSH for {organism_name}: {e}")
        return []


def extract_diseases_from_abstracts(text, organism_name):
    """Extract disease names from abstracts using comprehensive pattern matching."""
    diseases = []

    if not text:
        return diseases

    # Comprehensive disease patterns
    disease_patterns = [
        # Suffix-based diseases (more specific)
        r'\b(\w*tuberculosis)\b',
        r'\b(\w+osis)(?!\s*(?:diagnosis|prognosis|osmosis))\b',  # tuberculosis, candidiasis, etc. but exclude diagnosis
        r'\b(\w+itis)(?!\s*(?:arthritis|dermatitis))\b',  # meningitis, pneumonitis, etc.
        r'\b(\w+emia)(?!\s*(?:anemia|leukemia))\b',  # septicemia, bacteremia, etc.
        r'\b(\w+pathy)\b',  # nephropathy, etc.

        # Specific diseases
        r'\b(pneumonia)\b',
        r'\b(sepsis|septicemia)\b',
        r'\b(meningitis)\b',
        r'\b(endocarditis)\b',
        r'\b(gastroenteritis)\b',
        r'\b(diarrhea|diarrhoea)\b',
        r'\b(dysentery)\b',
        r'\b(malaria)\b',
        r'\b(cholera)\b',
        r'\b(tuberculosis|TB)\b',
        r'\b(AIDS|HIV)\b',

        # Infection types
        r'\b(urinary tract infection|UTI)\b',
        r'\b(hospital[- ]acquired infection)\b',
        r'\b(bloodstream infection)\b',
        r'\b(respiratory tract infection)\b',
        r'\b(skin and soft tissue infection)\b',
        r'\b(wound infection)\b',
        r'\b(surgical site infection)\b',
        r'\b(catheter[- ]related infection)\b',

        # Clinical syndromes
        r'\b(invasive disease)\b',
        r'\b(opportunistic infection)\b',
        r'\b(nosocomial infection)\b',
        r'\b(community[- ]acquired pneumonia)\b',

        # Multi-word diseases
        r'\b(pelvic inflammatory disease)\b',
        r'\b(inflammatory bowel disease)\b',
        r'\b(necrotizing fasciitis)\b',
        r'\b(toxic shock syndrome)\b',
    ]

    # Terms to exclude (biological processes, not diseases)
    exclude_terms = {
        ''
    }

    text_lower = text.lower()

    for pattern in disease_patterns:
        matches = re.findall(pattern, text_lower, re.IGNORECASE)
        for match in matches:
            if isinstance(match, tuple):
                match = [m for m in match if m][0]  # Get first non-empty match

            # Clean up the match
            disease = match.strip().title()

            # Filter out excluded terms and very generic or short terms
            if (len(disease) > 2 and
                    disease.lower() not in exclude_terms and
                    disease.lower() not in ['the', 'and', 'with', 'from'] and
                    not disease.lower().endswith('sis') or disease.lower() in ['sepsis', 'tuberculosis']):
                diseases.append(disease)

    return diseases


def extract_diseases_from_text(text):
    """Extract disease names from general text."""
    if not text:
        return []

    disease_keywords = [
        'tuberculosis', 'pneumonia', 'sepsis', 'meningitis', 'endocarditis',
        'gastroenteritis', 'infection', 'disease', 'syndrome', 'disorder',
        'malaria', 'cholera', 'dysentery', 'diarrhea'
    ]

    diseases = []
    text_lower = text.lower()

    for keyword in disease_keywords:
        if keyword in text_lower:
            diseases.append(keyword.title())

    return diseases


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
    org_type, diseases, biological_processes = classify_organism(organism)

    print(f"Organism type: {org_type}")
    print(f"Associated disease(s): {', '.join(diseases) if diseases else 'Unknown'}")
    print(
        f"Associated biological process(es): {', '.join(biological_processes) if biological_processes else 'Unknown'}")

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