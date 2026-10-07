#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mapping", help="TSV mapping: NUCCORE<TAB>REFSEQ")
parser.add_argument("folder", help="Folder containing .faa files")
args = parser.parse_args()

folder = Path(args.folder)

# Read mapping
refseq_to_nuccores = defaultdict(list)

with open(args.mapping) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        nuccore, refseq = line.split("\t")[:2]
        refseq_to_nuccores[refseq].append(nuccore)

# Process each RefSeq accession
for refseq, nuccores in refseq_to_nuccores.items():
    out_file = folder / f"{refseq}.faa"

    first = True
    for nc in nuccores:
        in_file = folder / f"{nc}.faa"

        if not in_file.exists():
            print(f"Missing: {in_file.name}")
            continue

        if first:
            shutil.move(in_file, out_file)
            first = False
        else:
            with open(out_file, "ab") as out, open(in_file, "rb") as inp:
                shutil.copyfileobj(inp, out)
            in_file.unlink()

    if not first:
        print(f"Created {out_file.name} ({len(nuccores)} source files)")
