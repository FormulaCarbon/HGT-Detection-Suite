from pathlib import Path

import pandas as pd
import yaml
from itertools import product

from apphelpers import *
from dlhelpers import download_url, unzip_fastas

def run(
    accession: str, 
    db_path: Path, 
    working_dir: Path, 
    conf: Path, 
    program: str, 
    threads: int,
    species_id_override: int | None = None, 
    use_min: bool = False, 
    mge: bool = False, 
    faa: Path | None = None,
    add_new_to_db: bool = False
    ):
    with open(db_path / "config.yaml", 'r') as conff:
        db_config = yaml.safe_load(conff)
        
    with open(conf, 'r') as conff:
        config = yaml.safe_load(conff)
        
    lineage_path = db_path / (db_config['min_lineage_path'] if use_min else db_config['lineage_path'])
    proteomes_dir = db_path / db_config['proteomes_path']
    
    lineages = pd.read_csv(lineage_path, sep = "\t")
    
    query = lineages[lineages["accession"] == accession].iloc[0].to_dict()
    
    temp_dir = working_dir / "temp"
    
    if faa is None:
        url = '_'.join(query["ftp_path"].split('_')[:-1]) + "_translated_cds.faa.gz"
        download_url(url, working_dir / "query.faa.gz")
        unzip_fastas(working_dir)
        print("Processing FASTA")
        query["file_path"] = process_fasta(working_dir / "query.faa")
    else:
    
        print("Processing FASTA")
        query["file_path"] = process_fasta(faa, out = str(working_dir / "query.faa.modified"))
        
    # Header files for type 1 inclusion
    acc_s1 = working_dir / "species_1_Accession_no_species.txt"
    acc_g1 = working_dir / "genus_1_Accession_no_species.txt"
    acc_f1 = working_dir / "family_1_Accession_no_species.txt"

    # Header files for type 2 inclusion
    acc_s2 = working_dir / "species_2_Accession_no_species.txt"
    acc_g2 = working_dir / "genus_2_Accession_no_species.txt"
    acc_f2 = working_dir / "family_2_Accession_no_species.txt"

    species_id = query["species_id"]
    tax_id = query["tax_id"]

    inclusion_flag = blast_all_levels(
        query['accession'],
        query['file_path'],
        lineages,
        program,
        threads,
        proteomes_dir,
        working_dir,
        species_id_override if species_id_override else 0,
        add_new_to_db
    )
    
    # Process hits
    
    for level, t in product(['species', 'genus', 'fasta'], [1, 2]):
        if inclusion_flag != 1 and (level, t) == ("species", 1): continue
        if inclusion_flag != 0 and (level, t) == ("species", 2): continue

        header(f"Processing hits for {level} {t}")
        process_hits(
            query['file_path'],
            working_dir / f"{level}_{t}_output_blasthits",
            working_dir / f"{level}_{t}_Accession_no_species.txt",
            level,
            t,
            config['identity_thresholds'][level],
            working_dir,
            config
        )
        
    thresh_pident_species = config["identity_thresholds"]["species"]
    thresh_pident_genus = config["identity_thresholds"]["genus"]
    thresh_pident_family = config["identity_thresholds"]["family"]

    calc_s1 = working_dir / f"{thresh_pident_species}_calculations_species_1.txt"
    calc_g1 = working_dir / f"{thresh_pident_genus}_calculations_genus_1.txt"
    calc_f1 = working_dir / f"{thresh_pident_family}_calculations_family_1.txt"

    calc_s2 = working_dir / f"{thresh_pident_species}_calculations_species_2.txt"
    calc_g2 = working_dir / f"{thresh_pident_genus}_calculations_genus_2.txt"
    calc_f2 = working_dir / f"{thresh_pident_family}_calculations_family_2.txt"

    classify_genes(
        query['file_path'],
        inclusion_flag,
        calc_s1,
        calc_g1,
        calc_f1,
        calc_s2,
        calc_g2,
        calc_f2,
        working_dir
    )