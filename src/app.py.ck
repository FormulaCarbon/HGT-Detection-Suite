from pathlib import Path

import pandas as pd
import yaml
from itertools import product

from apphelperscopied import *
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
    print(query)
    temp_dir = working_dir / "temp"
    
    working_dir.mkdir(parents=True, exist_ok=True)
    if faa is None:
        url = '_'.join(query["ftp_path"].split('_')[:-1]) + "_translated_cds.faa.gz"
        download_url(url, working_dir / "query.faa.gz")
        unzip_fastas(working_dir)
        print("Processing FASTA")
        query["file_path"] = process_fasta(working_dir / "query.faa")
    else:
    
        print("Processing FASTA")
        query["file_path"] = process_fasta(faa, out = str(working_dir / "query.faa.modified"))
    
    #query["file_path"] = working_dir / "query.faa.modified"
    # Header files for type 1 inclusion
    acc_s1 = working_dir / "s_1_Accession_no_species.tsv"
    acc_g1 = working_dir / "g_1_Accession_no_species.tsv"
    acc_f1 = working_dir / "f_1_Accession_no_species.tsv"

    # Header files for type 2 inclusion
    acc_s2 = working_dir / "s_2_Accession_no_species.tsv"
    acc_g2 = working_dir / "g_2_Accession_no_species.tsv"
    acc_f2 = working_dir / "f_2_Accession_no_species.tsv"

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
    
    #inclusion_flag=0
    
    # Process hits
    
    for level, t in product(['s', 'g', 'f'], [1, 2]):
        if inclusion_flag != 1 and (level, t) == ("s", 1): continue
        if inclusion_flag != 0 and (level, t) == ("s", 2): continue

        header(f"Processing hits for {level} {t}")
        process_hits(
            query['file_path'],
            working_dir / f"{level}_{t}_output_blasthits",
            working_dir / f"{level}_{t}_Accession_no_species.tsv",
            level,
            t,
            config['identity_thresholds']['species' if level == 's' else 'genus' if level=='g' else 'family'],
            working_dir,
            config
        )
        
    thresh_pident_species = config["identity_thresholds"]["species"]
    thresh_pident_genus = config["identity_thresholds"]["genus"]
    thresh_pident_family = config["identity_thresholds"]["family"]

    calc_s1 = working_dir / f"{thresh_pident_species}_calculations_s_1.tsv"
    calc_g1 = working_dir / f"{thresh_pident_genus}_calculations_g_1.tsv"
    calc_f1 = working_dir / f"{thresh_pident_family}_calculations_f_1.tsv"

    calc_s2 = working_dir / f"{thresh_pident_species}_calculations_s_2.tsv"
    calc_g2 = working_dir / f"{thresh_pident_genus}_calculations_g_2.tsv"
    calc_f2 = working_dir / f"{thresh_pident_family}_calculations_f_2.tsv"

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
