#!/usr/bin/env python3
"""
generate_app_commands.py

For every FASTA file in a directory, writes one line to an output file
containing the APP command to run on that genome:

  python3 src/main.py app {acc} {genome_db} -d {outdir_base}/{acc} -t {threads} -M --sid {sid} -f {fasta_path}

  - {acc}   = RefSeq accession, taken from the FASTA filename (the filename
              minus its extension, e.g. GCF_000009665.1.faa -> GCF_000009665.1)
  - {sid}   = taxonomic ID, looked up from a TSV table (RefSeq accession in
              column 1, taxid in column 2; any further columns are ignored)
  - {fasta_path} = path to the FASTA file itself (absolute by default)

Any FASTA file whose accession isn't found in the TSV is skipped and
reported at the end, rather than silently producing a broken command line.

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
python3 generate_app_commands.py \
    --fasta-dir genomes/ \
    --taxid-tsv taxid_table.tsv \
    --output app_commands.txt \
    --genome-db /work/11587/siddharthk86/ls6/APP/GENOME_DB/ \
    --outdir-base /scratch/11587/siddharthk86/APP/18genometest4/outputs \
    --threads 64

All of --genome-db / --outdir-base / --threads / --main-script have the
defaults shown in the example command, so you can omit any you don't need
to change.
"""

import argparse
import csv
import sys
from pathlib import Path

FASTA_EXTENSIONS = ('.fasta', '.fa', '.faa', '.fna')


def load_taxid_table(tsv_path: str) -> dict:
    """RefSeq accession (col 1) -> taxid (col 2). Extra columns ignored."""
    mapping = {}
    with open(tsv_path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row or not row[0].strip():
                continue
            acc = row[0].strip()
            if len(row) < 2 or not row[1].strip():
                continue
            taxid = row[1].strip()
            mapping[acc] = taxid
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta-dir', required=True, help='Directory containing FASTA files (one per genome)')
    ap.add_argument('--taxid-tsv', required=True, help='TSV: refseq_accession<TAB>taxid<TAB>...(ignored)')
    ap.add_argument('--output', required=True, help='Output file to write, one command per line')
    ap.add_argument('--pattern', default=None,
                     help='Glob pattern within --fasta-dir (default: *.fasta, *.fa, *.faa, *.fna)')
    ap.add_argument('--genome-db', default='/work/11587/siddharthk86/ls6/APP/GENOME_DB/',
                     help='GENOME_DB path (positional arg after acc in the APP command)')
    ap.add_argument('--outdir-base', default='/scratch/11587/siddharthk86/APP/18genometest6/outputs',
                     help='Base output directory; each command writes to <outdir-base>/<acc>')
    ap.add_argument('--threads', type=int, default=64, help='-t value')
    ap.add_argument('--main-script', default='src/main.py', help='Path to APP\'s main.py')
    ap.add_argument('--absolute-fasta-paths', action='store_true', default=True,
                     help='Use absolute paths for -f (default: on)')
    ap.add_argument('--relative-fasta-paths', dest='absolute_fasta_paths', action='store_false',
                     help='Use the fasta paths exactly as found (relative to --fasta-dir) instead of resolving to absolute')
    args = ap.parse_args()

    fasta_dir = Path(args.fasta_dir)
    if args.pattern:
        files = sorted(fasta_dir.glob(args.pattern))
    else:
        files = sorted(set(p for ext in FASTA_EXTENSIONS for p in fasta_dir.glob(f'*{ext}')))
    if not files:
        sys.exit(f"No FASTA files found in {fasta_dir}.")

    taxid_map = load_taxid_table(args.taxid_tsv)
    print(f"Loaded {len(taxid_map)} accession -> taxid entries from {args.taxid_tsv}")
    print(f"Found {len(files)} FASTA file(s) in {fasta_dir}")

    lines = []
    skipped = []
    completed = []
    for f in files:
        acc = f.stem  # filename minus extension, e.g. GCF_000009665.1.faa -> GCF_000009665.1
        if (Path(args.outdir_base) / acc / "APP_Alien_genes.tsv").is_file():
            completed.append(acc)
            continue
        sid = taxid_map.get(acc)
        if sid is None:
            skipped.append(acc)
            continue
        fasta_path = f.resolve() if args.absolute_fasta_paths else f
        cmd = (
            f"python3 {args.main_script} app {acc} {args.genome_db} "
            f"-d {args.outdir_base}/{acc} -t {args.threads} -M --sid {sid} -f {fasta_path}"
        )
        lines.append(cmd)

    with open(args.output, 'w') as out_fh:
        for line in lines:
            out_fh.write(line + '\n')

    print(f"Wrote {len(lines)} command(s) to {args.output}")
    if skipped:
        print(f"WARNING: {len(skipped)} FASTA file(s) had no matching taxid in the TSV and were skipped:")
        for acc in skipped:
            print(f"  {acc}")
    if completed:
        print(f"WARNING: {len(completed)} FASTA file(s) were already completed:")
        for acc in completed:
            print(f"  {acc}")

if __name__ == '__main__':
    main()
