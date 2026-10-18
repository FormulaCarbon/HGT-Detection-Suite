import argparse
import json
import os
import shutil

from pathlib import Path
import pandas as pd
from ete3 import NCBITaxa
from tqdm import tqdm
import numpy as np
from Bio import Entrez
from dotenv import load_dotenv

from iohelpers import header

load_dotenv()

Entrez.email = os.getenv('EMAIL')

def gen_taxa_db(target: pd.DataFrame, update_db = False) -> pd.DataFrame:
    Path("temp").mkdir(parents=True, exist_ok=True)
    ncbi = NCBITaxa()
    if update_db:
        ncbi.update_taxonomy_database()
    
    header("Cleaning Data")
    clean = target.dropna(subset = ["assembly_accession", "ftp_path", "species_taxid"])
    
    taxIDs = clean["species_taxid"].astype(int).tolist()
    
    lineages = []
    bad_taxids = []
    clean_taxids = []
    for taxID in tqdm(taxIDs, desc="Getting Lineage Information"):
        try:
            lineages.append(ncbi.get_lineage(taxID))
            clean_taxids.append(taxID)
        except ValueError as e:
            print(f"{e}, skipping")
            bad_taxids.append(taxID)
    taxIDs = clean_taxids
    clean = clean[clean["species_taxid"].isin(taxIDs)]
            
    
    accs = clean["assembly_accession"].astype(str).tolist()
    ftps = clean['ftp_path'].astype(str).tolist()
    
    # convert to fasta link
    urls = [f"{genomeUrl.rstrip("/")}/{genomeUrl.rstrip("/").split('/')[-1]}_protein.faa.gz" for genomeUrl in ftps]
    
    strain = [items.split(', /')[0].split('=')[1] if items != 'na' else '' for items in clean['infraspecific_name'].astype(str).tolist()]
    
    data = {
        "accession": accs,
        "tax_id": taxIDs,
        "kingdom" : [],
        "phylum" : [],
        "class" : [],
        "order" : [],
        "family" : [],
        "genus" : [],
        "species": [],
        "strain": strain,
        "kingdom_id" : [],
        "phylum_id" : [],
        "class_id" : [],
        "order_id" : [],
        "family_id" : [],
        "genus_id" : [],
        "species_id": [],
        "ftp_path": urls
    }
    
    for lineage in tqdm(lineages, desc="Creating Table"):
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(list(ranks.keys()))
        ranks = {v: k for k, v in ranks.items()}
        
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]: 
            if rank in ranks:
                data[rank].append(names.get(ranks[rank], np.nan))
                data[rank + "_id"].append(int(ranks[rank]))
            else:
                data[rank].append(np.nan)
                data[rank + "_id"].append(np.nan)
        
    # may not need this anymore        
    for k in ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]:
        cleaned = []
        
        for item in data[k]:
            if item is not np.nan:
                for c in ["<", ">", ":", '"', "/", "\\", "|", "?", "*"]:
                    item = item.replace(c, " ") 
            cleaned.append(item)
        data[k] = cleaned
          
    out = pd.DataFrame.from_dict(data)
    # Keep rows needed for APP bucketing; do not drop purely for missing higher-rank labels
    out = out.dropna(subset=["accession", "tax_id", "species_id", "genus_id", "family_id", "ftp_path"])
    shutil.rmtree("temp")
    return out

def graph(lineages: pd.DataFrame, outfile: Path):
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(lineages['tax_id'].astype(str).tolist())
    
    for node in tree.traverse():
        if node.name.isdigit():
            names = ncbi.get_taxid_translator([int(node.name)])
            if names:
                node.name = names[int(node.name)]
    
    tree.write(outfile=str(outfile))
    
def accession2taxid(acc: str, db="nucleotide") -> str:
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0] # type: ignore
    handle = Entrez.esummary(db=db, id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    return str(taxid)

def check_accession_updates(old_accs: list[str], new_accs: list[str]) -> dict[str, list[str]]:
    out = {
        "update": [],
        "new": [],
    }
    old_accs_nv = {acc.split('.')[0] for acc in old_accs}
    out["new"] = [acc for acc in new_accs if acc.split('.')[0] not in old_accs_nv]
    
    out["update"] = [acc for acc in set(new_accs) - set(old_accs) if acc not in out["new"]]
    
    return out
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("lineagePath")
    parser.add_argument("outFile")
    args = parser.parse_args()
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)

    
    lineages = pd.read_csv(args.lineagePath, sep = "\t", header=0, low_memory=False)
    graph(lineages, Path(args.outFile))
        
