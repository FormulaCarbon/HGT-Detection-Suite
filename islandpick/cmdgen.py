#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import shutil
import urllib.request
from pathlib import Path
import json


DB_PATH = Path("/work/11587/siddharthk86/ls6/APP/GENOME_DB")
OUTPUT_BASE = Path("/scratch/11587/siddharthk86/APP/118genometest/outputs")
DOWNLOAD_DIR = Path("/scratch/11587/siddharthk86/APP/118genometest/fastas")


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "accessions",
        help="File containing one GCF accession per line"
    )

    parser.add_argument(
        "output",
        help="Output file containing generated python3 commands"
    )

    parser.add_argument(
        "--db",
        default=str(DB_PATH),
        help="APP GENOME_DB directory"
    )

    parser.add_argument(
        "--download-dir",
        default=str(DOWNLOAD_DIR),
        help="Directory for FASTAs not already in the database"
    )

    return parser.parse_args()


def load_summary(db_path):
    summary = db_path / "BACTERIAL_LINEAGES.txt"

    database = {}

    with open(summary, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            accession = row["accession"].strip()
            database[accession] = row

    return database


def extract_url(value):
    """
    Handles either:

        https://ftp.ncbi.nlm.nih.gov/...

    or Markdown:

        [https://ftp.ncbi.nlm.nih.gov/...](https://ftp.ncbi.nlm.nih.gov/...)
    """

    value = value.strip()

    match = re.search(r"\]\((https?://[^)]+)\)", value)

    if match:
        return match.group(1)

    match = re.search(r"https?://\S+", value)

    if match:
        return match.group(0).rstrip(")")

    return value


def download_from_url(url, accession, download_dir):
    download_dir.mkdir(parents=True, exist_ok=True)

    # Get filename from URL.
    filename = Path(url.split("?")[0]).name

    if not filename:
        filename = f"{accession}_protein.faa.gz"

    compressed = download_dir / filename
    fasta = download_dir / (str(accession) + ".faa")

    # Don't download again if it already exists.
    if fasta.exists():
        return fasta

    if not compressed.exists():
        print(f"Downloading {accession}...")
        
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": "APP-genome-downloader/1.0"
            }
        )

        with urllib.request.urlopen(request) as response:
            with open(compressed, "wb") as out:
                shutil.copyfileobj(response, out)

    # Convert .gz -> .faa
    print(f"Extracting {accession}...")

    with gzip.open(compressed, "rb") as src:
        with open(fasta, "wb") as dst:
            shutil.copyfileobj(src, dst)

    return fasta


def download_missing_from_ncbi(accession, download_dir):
    """
    Download a protein FASTA for an accession that is completely absent
    from BACTERIAL_LINEAGES.tsv.

    Uses the standard NCBI RefSeq directory structure.
    """

    # GCF_000005825.2
    print(accession)
    match = re.match(r"^(GCF)_(\d{3})(\d{3})(\d{3})\.(\d+)$", accession)
    if not match:
        raise ValueError(
            f"Cannot determine NCBI directory structure for {accession}"
        )

    _, a, b, c, version = match.groups()

    base = f"{accession}"

    # NCBI assembly directories use the assembly name after the accession.
    # Since we don't know the assembly name for an accession absent from
    # BACTERIAL_LINEAGES.tsv, retrieve it from NCBI's assembly_summary.
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/"
        "assembly_summary.txt"
    )

    # Download the summary only once.
    download_dir.mkdir(parents=True, exist_ok=True)
    summary_file = download_dir / "assembly_summary_refseq_bacteria.txt"

    if not summary_file.exists():
        print("Downloading NCBI bacterial assembly summary...")

        request = urllib.request.Request(
            summary_url,
            headers={"User-Agent": "APP-genome-downloader/1.0"}
        )

        with urllib.request.urlopen(request) as response:
            with open(summary_file, "wb") as out:
                shutil.copyfileobj(response, out)

    ftp_path = None
    species_id = None

    with open(summary_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if not fields:
                continue

            if fields[0] != accession:
                continue

            # RefSeq assembly_summary columns:
            # 0 accession
            # 5 organism_name
            # 6 infraspecific_name
            # ...
            # 20 ftp_path
            if len(fields) > 20:
                ftp_path = fields[19]

            # Taxid is column 6 in current assembly_summary format.
            if len(fields) > 5:
                species_id = fields[5]

            break

    if not ftp_path:
        raise RuntimeError(
            f"Could not find {accession} in NCBI assembly_summary"
        )

    # The FTP path looks like:
    #
    # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/...
    #
    # Find the actual protein filename.
    ftp_path = ftp_path.rstrip("/")

    assembly_name = ftp_path.split("/")[-1]

    protein_url = (
        f"{ftp_path}/{assembly_name}_protein.faa.gz"
    )

    fasta = download_from_url(
        protein_url,
        accession,
        download_dir
    )

    return fasta, species_id


def main():
    args = parse_args()

    db_path = Path(args.db)
    download_dir = Path(args.download_dir)

    database = load_summary(db_path)

    with open(args.accessions, "r", encoding="utf-8") as f:
        accessions = list(json.load(f).values())

    commands = []

    for accession in accessions:
        if not accession.startswith("GCF"):
            continue

        if (OUTPUT_BASE / accession / "APP_Alien_genes.tsv").is_file():
            continue
        print(f"Processing {accession}...")

        # ---------------------------------------------------------
        # CASE 1:
        # Genome exists in the local GENOMES database.
        # ---------------------------------------------------------

        database_fasta = (
            db_path /
            "GENOMES" /
            f"{accession}.faa"
        )

        if database_fasta.exists():

            row = database.get(accession)

            if row is None:
                print(
                    f"{database_fasta} exists but {accession} "
                    "is missing from BACTERIAL_LINEAGES.tsv"
                )
                continue

            species_id = row["species_id"]

            fasta = database_fasta

            print(f"  Using database FASTA: {fasta}")

        # ---------------------------------------------------------
        # CASE 2:
        # Genome is in the summary but FASTA isn't in GENOMES.
        # ---------------------------------------------------------

        elif accession in database:

            row = database[accession]

            species_id = row["species_id"]

            url = extract_url(row["ftp_path"])

            print(f"  Downloading from database ftp_path")

            fasta = download_from_url(
                url,
                accession,
                download_dir
            )

        # ---------------------------------------------------------
        # CASE 3:
        # Genome isn't in the database summary at all.
        # ---------------------------------------------------------

        else:

            print(
                f"  {accession} not in BACTERIAL_LINEAGES.tsv; "
                "looking up on NCBI"
            )
            try:

                fasta, species_id = download_missing_from_ncbi(
                accession,
                download_dir
                )

                if not species_id:
                    raise RuntimeError(
                    f"Could not determine species ID for {accession}"
                    )
            except:
                print('err')
                continue
        # ---------------------------------------------------------
        # Generate requested command.
        # ---------------------------------------------------------

        output_dir = OUTPUT_BASE / accession

        command = (
            f"python3 src/main.py app "
            f"{accession} "
            f"{db_path}/ "
            f"-d {output_dir} "
            f"-t 64 "
            f"-M "
            f"--sid {species_id} "
            f"-f {fasta}"
        )

        commands.append(command)

    # -------------------------------------------------------------
    # Write ONLY the python3 lines to the output file.
    # -------------------------------------------------------------

    with open(args.output, "w", encoding="utf-8") as f:
        for command in commands:
            f.write(command + "\n")

    print(
        f"\nWrote {len(commands)} commands to {args.output}"
    )


if __name__ == "__main__":
    main()
