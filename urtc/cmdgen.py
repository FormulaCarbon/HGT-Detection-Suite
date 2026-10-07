#!/usr/bin/env python3
"""
generate_full_dataset_cmds.py

Builds cmds.txt for the 700+ genome dataset: one APP command line per
genome, in the same style as before --

  python3 src/main.py app {acc} /work/11587/siddharthk86/ls6/APP/GENOME_DB/ \
      -d /scratch/11587/siddharthk86/APP/urtc/outputs/{acc} -t 64 -M \
      --sid {species_id} -f {fasta_path}

  - {acc}         = refseq_accession (GCF) from the mapping file produced
                    by prepare_full_dataset.py
  - {species_id}  = species_id from that same mapping file
  - {fasta_path}  = resolved from the mapping file's "source" column:
                      existing_db     -> <db-genomes-dir>/{acc}.faa (or
                                          <db-genomes-dir>/{acc}* if that
                                          exact name isn't found)
                      downloaded_new  -> <new-genomes-dir>/{acc}.faa
                      anything else (unresolved / download_failed)
                                      -> skipped, reported at the end

Before writing each line, checks whether
    <outdir-base>/{acc}/APP_Alien_genes.tsv
already exists -- if so, that genome is treated as already completed and
is NOT written to cmds.txt.

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
python3 generate_full_dataset_cmds.py \
    --mapping-tsv genome_mapping.tsv \
    --db-genomes-dir /path/to/dbpath/GENOMES/ \
    --new-genomes-dir new_genomes/ \
    --output cmds.txt

All of --genome-db / --outdir-base / --threads / --main-script /
--alien-filename have the defaults shown above and can be overridden.
"""

import argparse
import csv
import sys
from pathlib import Path


def resolve_fasta_path(acc: str, source: str, db_genomes_dir: Path, new_genomes_dir: Path):
    if source == 'existing_db':
        exact = db_genomes_dir / f"{acc}.faa"
        if exact.exists():
            return exact
        matches = sorted(db_genomes_dir.glob(f"{acc}*"))
        return matches[0] if matches else None
    elif source == 'downloaded_new':
        candidate = new_genomes_dir / f"{acc}.faa"
        return candidate if candidate.exists() else None
    else:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mapping-tsv', required=True,
                     help='genome_mapping.tsv from prepare_full_dataset.py '
                          '(nuccore_accession, refseq_accession, species_id, genome_id, organism, source)')
    ap.add_argument('--db-genomes-dir', required=True, help='Existing dbpath/GENOMES/ directory')
    ap.add_argument('--new-genomes-dir', required=True, help='Folder of newly downloaded genomes (from prepare_full_dataset.py)')
    ap.add_argument('--output', default='cmds.txt', help='Output commands file (default: cmds.txt)')
    ap.add_argument('--genome-db', default='/work/11587/siddharthk86/ls6/APP/GENOME_DB/')
    ap.add_argument('--outdir-base', default='/scratch/11587/siddharthk86/APP/urtc/outputs')
    ap.add_argument('--threads', type=int, default=64)
    ap.add_argument('--main-script', default='src/main.py')
    ap.add_argument('--alien-filename', default='APP_Alien_genes.tsv',
                     help='Filename checked under <outdir-base>/<acc>/ to decide whether a genome is already done')
    ap.add_argument('--absolute-fasta-paths', action='store_true', default=True)
    ap.add_argument('--relative-fasta-paths', dest='absolute_fasta_paths', action='store_false',
                     help='Use fasta paths as given instead of resolving to absolute')
    args = ap.parse_args()

    db_genomes_dir = Path(args.db_genomes_dir)
    new_genomes_dir = Path(args.new_genomes_dir)
    outdir_base = Path(args.outdir_base)

    with open(args.mapping_tsv, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'refseq_accession', 'species_id', 'source'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{args.mapping_tsv} is missing column(s): {missing}. Found: {reader.fieldnames}")
        rows = list(reader)

    print(f"Loaded {len(rows)} genome(s) from {args.mapping_tsv}")

    lines = []
    n_already_done = n_no_fasta = n_no_species_id = n_bad_source = 0

    for row in rows:
        acc = row['refseq_accession'].strip()
        species_id = row['species_id'].strip()
        source = row['source'].strip()

        if not acc or source in ('unresolved', 'download_failed'):
            n_bad_source += 1
            continue

        alien_file = outdir_base / acc / args.alien_filename
        if alien_file.exists():
            n_already_done += 1
            continue

        if not species_id:
            print(f"  WARNING: {acc} has no species_id, skipping")
            n_no_species_id += 1
            continue

        fasta_path = resolve_fasta_path(acc, source, db_genomes_dir, new_genomes_dir)
        if fasta_path is None:
            print(f"  WARNING: could not find a FASTA file for {acc} (source={source}), skipping")
            n_no_fasta += 1
            continue
        if args.absolute_fasta_paths:
            fasta_path = fasta_path.resolve()

        cmd = (
            f"python3 {args.main_script} app {acc} {args.genome_db} "
            f"-d {args.outdir_base}/{acc} -t {args.threads} -M --sid {species_id} -f {fasta_path}"
        )
        lines.append(cmd)

    with open(args.output, 'w') as out_fh:
        for line in lines:
            out_fh.write(line + '\n')

    print(f"\nWrote {len(lines)} command(s) to {args.output}")
    print(f"  Already completed (alien-genes file exists) : {n_already_done}")
    print(f"  Skipped -- unresolved/download_failed source : {n_bad_source}")
    print(f"  Skipped -- no FASTA file found                : {n_no_fasta}")
    print(f"  Skipped -- no species_id                      : {n_no_species_id}")


if __name__ == '__main__':
    main()
