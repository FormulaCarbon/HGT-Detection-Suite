from collections import defaultdict
from pathlib import Path
import re
import subprocess
from typing import Counter

from tqdm import tqdm
import pandas as pd

from dlhelpers import download_url, unzip_fastas
from db import add_to_db
from iohelpers import header, error

def process_fasta(faa: Path | str, out = "") -> Path:
    if out == "":
        out = f"{faa}.modified"
    with open(faa, 'r') as infile, open(out, 'w') as outfile:
        counter = 0
        for line in infile:
            if line.startswith(">"):
                if counter > 0:
                    outfile.write("\n")
                counter += 1
                line = line.replace("-->", '').replace("->",'')
                header_info = line[1:].split("] [")
                found = False

                for key in header_info:
                    if "protein_id" in key:
                        acc_no = key.split("=")[1].strip()
                        outfile.write(f">{acc_no}:g{counter}\n")
                
                        found = True
                        break
                if not found:
                    outfile.write(f"{line}:g{counter}")
            else:
                outfile.write(line.strip())
    return Path(out)

def compile_multifasta(files: list[Path], outfile: Path, map_outfile: Path) -> Path:
    with open(outfile, 'w', encoding='utf-8') as out, open(map_outfile, "w", encoding="utf-8") as map_out:
        for i, file in tqdm(enumerate(files), total=len(files)):
            filename = file.stem.replace('.', '_')
            with open(file, 'r', encoding='utf-8') as infile:
                for line in infile:
                    if line.startswith('>'):
                        protein_id = line[1:].split()[0]
                        out.write(f">{protein_id}:g{i}\n")
                        map_out.write(f"{protein_id}:g{i}\t{filename}\n")
                    else:
                        out.write(line)
                        
    return outfile

def run_blastp(query_file: Path, multifasta: Path, working_dir: Path, threads: int, level: str, t: int, size: int):
    header(f"Creating BLAST DB from {multifasta}")
    subprocess.run([
        "makeblastdb",
        "-in", str(multifasta),
        "-parse_seqids",
        "-dbtype", "prot"
    ], check = True )
    
    header(f"Running BLASTp with {threads} threads...")
    subprocess.run([
        "blastp",
        "-num_threads", str(threads),
        "-query", str(query_file),
        "-db", str(multifasta),
        "-out", str(working_dir / f"{level}_{t}_output_blasthits"),
        "-outfmt", '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore',
        "-max_target_seqs", str(size)
    ], check=True)
    header(f"Done")
    
def run_diamond(query_file: Path, multifasta: Path, outfile: Path, working_dir: Path, threads: int, level: str, t: int, size: int):
    
    header(f"Creating DIAMOND DB at {outfile}.dmnd from {multifasta}")
    subprocess.run([
        "diamond",
        "makedb",
        "--in",
        str(multifasta),
        "--db",
        str(outfile),
        "--verbose"
    ], check = True )
    
    header(f"Running DIAMOND with {threads} threads...")
    subprocess.run([
        "diamond", "blastp",
        "-q",  str(query_file),
        "-d", str(outfile),
        "-o", str(working_dir / f"{level}_{t}_output_blasthits"),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "qcovhsp", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "--threads", str(threads),
        "--max-target-seqs", str(size),
        "--ultra-sensitive",
        "--comp-based-stats", "4"
    ], check=True)
    header(f"Done")

def get_fastas(accs: list[str], genomes_dir: Path, working_dir: Path, lineages: pd.DataFrame, add_new_to_db: bool) -> list[Path]:
    out = []
    for acc in accs:
        file = genomes_dir / f"{acc}.faa"
        if file.is_file():
            out.append(file)
        else:
            if add_new_to_db:
                add_to_db() # TODO: Implement this when add_to_db is done
                out.append(file)
            else:
                if (not (working_dir / "temp" / f"{acc}.faa").is_file()) and download_url(lineages.loc[lineages['accession']==acc, 'ftp_path'].values[0], working_dir / "temp" / f"{acc}.faa.gz") == 0:
                    out.append(working_dir / "temp" / f"{acc}.faa")
                else:
                    print('error downloading ', acc)
                    continue
    unzip_fastas(working_dir / "temp")
    return out
                
def blast_single_level(query_file:Path, id_at_level: int, id_to_exclude: int| None, lineages: pd.DataFrame, level: str, t: int, program: str, threads: int, genomes_dir: Path, working_dir: Path, add_new_to_db: bool) -> int:
    exclude_col, target_col = \
        ("taxid", "species_id") if level == "s" else ("species_id", "genus_id") if level == "g" else \
        ("genus_id", "family_id") if level == "f" else ("family_id", "order_id") if level == 'o' else (None, None)
    
    if target_col is None: raise ValueError("Valid level not provided. Level should be 's', 'g', 'f', or 'o' ")
    
    matches = lineages[lineages[target_col].astype(str) == str(id_at_level)]
    matches = matches[lineages[exclude_col].astype(str) != str(id_to_exclude)] if t == 2 else matches
    
    accs = matches['accession'].dropna().drop_duplicates().tolist()
    
    files = get_fastas(accs, genomes_dir, working_dir, lineages, add_new_to_db)
    
    if len(files) < 1:
        error(f"No {level} t{t} files")
        return 1
    else:
        header(f"{level} t{t}: {len(files)} matches")
        
    mfa = compile_multifasta(
        files, 
        working_dir / f"{level}_{t}_raw_database.fasta",
        working_dir / f"{level}_{t}_Accession_no_species.tsv"
    )
    
    if program == "blastp":
        run_blastp(
            query_file,
            mfa, 
            working_dir,
            threads, 
            level, 
            t,
            len(files)
        )
    elif program == "diamond":
        run_diamond(
            query_file,
            mfa, 
            working_dir / f"{level}_{t}_db",
            working_dir,
            threads,
            level,
            t,
            len(files)
        )
    return 0
        
def blast_all_levels(acc: str, query_file: Path, lineages: pd.DataFrame, program: str, threads: int, genomes_dir: Path, working_dir: Path, sid: int = 0, add_new_to_db: bool = False) -> int:
    kingdom, phylum, clas, order, family, genus, species = \
        lineages.loc[lineages["accession"] == acc, ["kingdom_id", "phylum_id", "class_id", "order_id", "family_id", "genus_id", "species_id"]].iloc[0] 
    
    if sid != 0:
        species = sid

    inclusion_flag = 0
    
    header("Species (Inclusive)")
    blast_single_level(
        query_file,
        species,
        None,
        lineages,
        's',
        1,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    header("Genus (Inclusive)")
    blast_single_level(
        query_file,
        species,
        None,
        lineages,
        'g',
        1,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    header("Family (Inclusive)")
    family_exists = blast_single_level(
        query_file,
        species,
        None,
        lineages,
        'f',
        1,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    if family_exists == 1:
        error("No Family (Inclusive) Relatives Found")
        header("Order (Inclusive)")
        blast_single_level(
            query_file,
            species,
            None,
            lineages,
            'o',
            1,
            program,
            threads,
            genomes_dir,
            working_dir,
            add_new_to_db
        )
        
    header("Species (Exclusive)")
    inclusion_flag = blast_single_level(
        query_file,
        species,
        None,
        lineages,
        's',
        2,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    header("Genus (Exclusive)")
    blast_single_level(
        query_file,
        species,
        None,
        lineages,
        'g',
        2,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    header("Family (Exclusive)")
    family_exists = blast_single_level(
        query_file,
        species,
        None,
        lineages,
        'f',
        2,
        program,
        threads,
        genomes_dir,
        working_dir,
        add_new_to_db
    )
    
    if family_exists == 1:
        error("No Family (Exclusive) Relatives Found")
        header("Order (Exclusive)")
        blast_single_level(
            query_file,
            species,
            None,
            lineages,
            'o',
            2,
            program,
            threads,
            genomes_dir,
            working_dir,
            add_new_to_db
        )
        
    return inclusion_flag

def process_hits(query_file: Path, hits: Path, map_file: Path, level: str, t: int, pident_thresh: int, working_dir: Path, config: dict):
    names = "qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore".split()
    blasthits = pd.read_csv(hits, sep="\t", header=None, names=names)
    
    mappings = (
        pd.read_csv(map_file, sep="\t", header=None, names=["prot_acc", "gen_acc"])
          .dropna()
          .set_index("prot_acc")["gen_acc"]
          .to_dict()
    )
    
    with open(query_file, "r", encoding="utf-8", errors="replace") as f:
        sequences = f.readlines()
        
    main_out_path = working_dir / f"Main_blast_output_counts_pident_{pident_thresh}_{level}_{t}.tsv"
    calc_out_path = working_dir / f"{pident_thresh}_calculations_{level}_{t}.tsv"
    temp_out_path = working_dir / f"FOR_INPUT_{level}_{t}.tsv"
    unique_genomes = sorted(set(mappings.values()))
    
    accs = []
    for s in sequences:
        s = s.strip()
        if s.startswith(">"):
            accs.append(s[1:])
            
    blasthits_counter = defaultdict(Counter)
    for gene in accs:
        for genome in unique_genomes:
            blasthits_counter[gene][genome] = 0
            
    for hit in blasthits.itertuples(index=False):
        query = str(hit.qseqid)
        subject = str(hit.sseqid)
        pident = float(hit.pident)
        qcov = float(hit.qcovs)
        
        target = mappings.get(subject)
        if target is None:
            continue

        if pident >= pident_thresh and qcov >= config['coverage_threshold']:
            blasthits_counter[query][target] = 1
    
    header_cols = ["genes"] + unique_genomes
    
    with open(main_out_path, "w", encoding="utf-8") as main, \
         open(temp_out_path, "w", encoding="utf-8") as temp:

        

        main.write("\t".join(header_cols) + "\n")
        temp.write("\t".join(header_cols) + "\n")

        n_genomes = len(unique_genomes)

        for gene in sorted(blasthits_counter.keys()):
            row_vals = [str(blasthits_counter[gene][g]) for g in unique_genomes]
            line = "\t".join([gene] + row_vals) + "\n"

            main.write(line)
            temp.write(line)
            
    with open(temp_out_path, "r", encoding="utf-8") as infile, \
         open(calc_out_path, "w", encoding="utf-8") as outfile:

        outfile.write("Gene\ttotal_count\tpresent_count\tfrac\n")
        
        for line in infile:
            line = line.strip()
            
            if re.match(r"^genes", line):
                continue
            
            arr = line.split('\t')
            gene = arr[0].strip()
            
            gene_counter = 0
            
            for i in range(0, len(arr) - 1):
                try:
                    if int(arr[i + 1]) == 1:
                        gene_counter += 1
                except ValueError:
                    pass
                    
            size = len(arr)-1
            frac = (float(gene_counter) / float(size)) * 100
            
            outfile.write(f"{gene}\t{size}\t{gene_counter}\t{frac}\n")

def calculate_pident(calc_files: Path):
    out = {}
    with open(calc_files, 'r', encoding="utf-8") as inputfile:
        for line in inputfile:
            line = line.strip()
            if re.match(r"^Gene", line):
                continue
            arr = line.split("\t")
            acc = arr[0].strip()
            pdt = arr[3].strip()
            out[acc] = {"pdt": float(pdt)}
    
    return out

def classify_genes(
    query_file: Path,
    flag: int,
    calc_s1: Path,
    calc_g1: Path,
    calc_f1: Path,
    calc_s2: Path,
    calc_g2: Path,
    calc_f2: Path,
    working_dir: Path
    ):
    full_out_path = working_dir / "phyletic_pattern_analysis_native_alien_genes_full_record.tsv"
    out_path = working_dir / "APP_Alien_genes.tsv"
    
    pdts_s = calculate_pident(calc_s1) if flag == 1 else calculate_pident(calc_s2)
    pdts_g1 = calculate_pident(calc_g1)
    pdts_f1 = calculate_pident(calc_f1)
    pdts_g2 = calculate_pident(calc_g2)
    pdts_f2 = calculate_pident(calc_f2)
    
    accession_nos = {}
    counter = 1
    with open(query_file, 'r') as q:
        for line in q:
            if line.startswith(">"):
                g = line.split(">")[1].strip()
                accession_nos[counter] = g
                counter += 1
                
    out = {}
    
    for acc in pdts_s:
        pdt_species = float(pdts_s.get(acc, {}).get("pdt", 0.0))
        pdt_genus_t1 = float(pdts_g1.get(acc, {}).get("pdt", 0.0))
        pdt_family_t1 = float(pdts_f1.get(acc, {}).get("pdt", 0.0))
        pdt_genus_t2 = float(pdts_g2.get(acc, {}).get("pdt", 0.0))
        pdt_family_t2 = float(pdts_f2.get(acc, {}).get("pdt", 0.0))
        
        if pdt_species <= 30:
            out[acc] = {"type": 0, "mode": 0}
        elif pdt_species >= 80:
            if pdt_genus_t2 >= 70:
                if pdt_family_t2 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t2 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t2 > 30 and pdt_family_t2 < 40:
                    out[acc] = {"type": "Ambiguous1", "mode": 2}
                else:
                    out[acc] = {"type": "Error1"}
            elif pdt_genus_t2 <= 70 and pdt_family_t2 >= 30:
                if pdt_family_t1 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t1 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t1 > 30 and pdt_family_t1 < 40:
                    out[acc] = {"type": "Ambiguous2", "mode": 2}
                else:
                    out[acc] = {"type": "Error2"}
            elif pdt_genus_t2 < 30:
                out[acc] = {"type": 0, "mode": 1}
            else:
                out[acc] = {"type": "Error3"}
        elif pdt_species < 80 and pdt_species > 30:
            if pdt_genus_t1 <= 30:
                out[acc] = {"type": 0, "mode": 1}
            elif pdt_genus_t1 >= 70:
                out[acc] = {"type": 1, "mode": 2}
            elif pdt_genus_t1 > 30 and pdt_genus_t1 < 70:
                if pdt_family_t1 <= 30:
                    out[acc] = {"type": 0, "mode": 1}
                elif pdt_family_t1 >= 40:
                    out[acc] = {"type": 1, "mode": 2}
                elif pdt_family_t1 > 30 and pdt_family_t1 < 40:
                    out[acc] = {"type": "Ambiguous3", "mode": 2}
                else:
                    out[acc] = {"type": "Error4"}
        else:
            out[acc] = {"type": "Error"}
            
        with open(full_out_path, 'w') as full, open(out_path, 'w') as outfile:
            for a in sorted(out.keys()):
                b = pdts_s.get(a, {}).get("pdt", "")
                c = pdts_g1.get(a, {}).get("pdt", "")
                d = pdts_f1.get(a, {}).get("pdt", "")
                f = pdts_g2.get(a, {}).get("pdt", "")
                g = pdts_f2.get(a, {}).get("pdt", "")
                h = out[a].get("type", "")
                t = out[a].get("mode", "")
                full.write(f"{a}\t{b}\t{c}\t{d}\t{f}\t{g}\t{h}\t{t}\n")
                if h == 0:
                    if t == 1:
                        outfile.write(f"{a}\t{t}\tAncient\n")
                    elif t == 0:
                        outfile.write(f"{a}\t{t}\tRecent\n")

def marker_gene_enrichment():
    pass
    
if __name__ == "__main__":
    process_fasta("query.faa")