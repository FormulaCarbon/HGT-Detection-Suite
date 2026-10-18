from pathlib import Path
import os
import urllib.request
from email.utils import parsedate_to_datetime
import shutil

from dotenv import load_dotenv
import yaml
import pandas as pd

from iohelpers import header
from taxahelpers import gen_taxa_db
from dlhelpers import download_dataset, split_by_filetype

load_dotenv()

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
NCBI_API_KEY = os.getenv('NCBI_API_KEY')

def create_db(db_path: Path, max_workers: int = 10, update_ncbi_db = False):
    db_path.mkdir(parents=True, exist_ok=True)
    config_path = db_path / "config.yaml"
    
    with urllib.request.urlopen(REFSEQ_URL) as req:
        refseq_timestamp = parsedate_to_datetime(req.headers["Last-Modified"]).date()
    
    config = {
        "last_updated": str(refseq_timestamp),
        "proteomes_path": "PROT/",
        "nucleotides_path": "NUCL/",
        "gbff_path": "GBFF/",
        "lineage_path": "BACTERIAL_LINEAGES.txt",
        "min_lineage_path": "BACTERIAL_LINEAGES_MIN.txt",
        "refseq_path": "assembly_summary_refseq.txt",
        "log_path": "log.txt"
    }
    with open(config_path, 'w') as conf:
        yaml.dump(config, conf, default_flow_style = False)
        
    header("Downloading RefSeq...")
    refseq_path = db_path / config['refseq_path']
    urllib.request.urlretrieve(REFSEQ_URL, refseq_path)
    
    header("Creating Lineage File...")
    header(f"Importing RefSeq from {refseq_path}...")
    refseq = pd.read_csv(refseq_path, sep = "\t", header=0, skiprows=1, low_memory=False)
    refseq.columns = refseq.columns.str.lstrip("#").str.strip()

    header("Filtering Bacteria Only...")
    refseq = refseq[(refseq["group"] == "bacteria") & (refseq["version_status"] == "latest") & (refseq["assembly_level"] == "Complete Genome")]
    refseq = refseq.loc[refseq.groupby(['assembly_accession', 'infraspecific_name'])['seq_rel_date'].idxmax()]

    header("Adding Taxonomy Information...")
    taxa_db = gen_taxa_db(refseq, update_ncbi_db)
    
    lineage_path = db_path / config['lineage_path']

    header(f"Saving to {lineage_path}...")
    taxa_db.to_csv(lineage_path, sep = "\t", index = False, mode = "w", float_format='%.0f')

    header("Lineage File Complete, Starting Download...")
    with open(db_path / "accessions.txt", 'w') as acc_file:
        for accession in taxa_db["accession"]:
            acc_file.write(f"{accession}\n")

    
    (db_path / "temp").mkdir(parents=True, exist_ok=True)
    download_dataset(db_path / "accessions.txt", db_path / "ncbi_dataset.zip", db_path / "temp", max_workers, NCBI_API_KEY)
    header("Sorting by Filetype...")
    split_by_filetype(db_path / "temp" / "ncbi_dataset" / "data", db_path / config['proteomes_path'], db_path / config['nucleotides_path'], db_path / config['gbff_path'])

    shutil.rmtree(db_path / "temp")
    header("Done!")
def update_db():
    pass

def add_to_db():
    pass
