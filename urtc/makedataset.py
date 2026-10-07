#!/usr/bin/env python3
"""
prepare_full_dataset.py

Takes the paper's genome table (Table S1: GenomeID, Division, Family,
Organism, Accession, Replicon, Length -- "Accession" is a nuccore
accession, e.g. NC_008009) and reconciles it against an existing local
APP reference database:

    <dbpath>/BACTERIAL_LINEAGES.tsv   (accession[GCF] <tab> tax_id <tab> ... <tab> species_id <tab> ... <tab> ftp_path)
    <dbpath>/GENOMES/                 (protein FASTAs, one per GCF accession already in the db)

For every genome in Table S1:
  1. Resolve its nuccore accession to a GCF assembly accession AND a
     species-level taxid (NCBI's "SpeciesTaxid" for that assembly) in one
     combined Entrez Direct lookup (esearch -> esummary -> xtract).
  2. Check whether that GCF accession already exists in
     BACTERIAL_LINEAGES.tsv AND has a FASTA file under <dbpath>/GENOMES/.
       - If yes: use the species_id already recorded in the db (kept
         authoritative over a fresh lookup), don't touch the db or
         download anything.
       - If no: download that genome's protein FASTA from NCBI, save it
         as <new-genomes-dir>/GCF_ACCESSION.faa (a separate folder -- the
         existing db is never written to by this script), and use the
         freshly-resolved species-level taxid.
  3. Write one mapping row per genome either way:
         nuccore_accession, refseq_accession (GCF), species_id, genome_id,
         organism, source (existing_db / downloaded_new / unresolved)

This mapping file is what feeds generate_app_commands.py's --taxid-tsv
(GCF in column 1, species_id in column 2), and --fasta-dir would then be
either <dbpath>/GENOMES/ or <new-genomes-dir> depending on each genome's
"source" -- see the printed summary for how many fall into each bucket.

------------------------------------------------------------------------
REQUIREMENTS
------------------------------------------------------------------------
NCBI Entrez Direct (esearch/esummary/xtract) on PATH.

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
python3 prepare_full_dataset.py \
    --table-s1 TableS1.tsv \
    --db-lineage-tsv /path/to/dbpath/BACTERIAL_LINEAGES.tsv \
    --db-genomes-dir /path/to/dbpath/GENOMES/ \
    --new-genomes-dir new_genomes/ \
    --mapping-output genome_mapping.tsv \
    --email you@example.com

Table S1 can be .tsv/.csv/.xlsx -- detected by extension.
"""

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


# --------------------------------------------------------------------
# Entrez Direct: resolve a nuccore accession -> (GCF accession, species taxid)
# in a single esearch -> esummary -> xtract round trip.
# --------------------------------------------------------------------

def get_assembly_info(nc_accession):
    search = subprocess.Popen(
        ['esearch', '-db', 'assembly', '-query', nc_accession],
        stdout=subprocess.PIPE
    )
    summary = subprocess.Popen(
        ['esummary'],
        stdin=search.stdout,
        stdout=subprocess.PIPE
    )
    xtract = subprocess.Popen(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'RefSeq', 'SpeciesTaxid'],
        stdin=summary.stdout,
        stdout=subprocess.PIPE,
        text=True
    )
    output, _ = xtract.communicate()
    return output


class AssemblyResolver:
    def __init__(self, cache_path, delay=0.35):
        self.cache_path = Path(cache_path)
        self.cache = json.loads(self.cache_path.read_text()) if self.cache_path.exists() else {}
        self.delay = delay

    def save_cache(self):
        self.cache_path.write_text(json.dumps(self.cache, indent=2, sort_keys=True))

    def resolve(self, nc_accession):
        """Returns (gcf_accession_or_None, species_taxid_or_None, note)."""
        if nc_accession in self.cache:
            e = self.cache[nc_accession]
            return e.get('gcf'), e.get('species_taxid'), e.get('note', 'cached')

        gcf, species_taxid, note = None, None, ''
        try:
            time.sleep(self.delay)
            raw = get_assembly_info(nc_accession)
            lines = [l for l in raw.splitlines() if l.strip()]
            if not lines:
                note = 'no assembly found (empty xtract output)'
            else:
                parts = lines[0].split('\t')
                gcf = parts[0].strip() if len(parts) > 0 and parts[0].strip() else None
                species_taxid = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
                note = 'resolved via esearch/esummary/xtract'
                if len(lines) > 1:
                    note += f' ({len(lines)} assemblies matched, used first)'
                if not gcf:
                    note = 'esearch/esummary matched but RefSeq (GCF) field was empty'
        except FileNotFoundError:
            note = 'esearch/esummary/xtract not found on PATH -- is Entrez Direct installed?'
        except Exception as e:
            note = f'lookup error: {e}'

        self.cache[nc_accession] = {'gcf': gcf, 'species_taxid': species_taxid, 'note': note}
        self.save_cache()
        return gcf, species_taxid, note


def get_gcf_accession(nc_accession):
    """Kept for standalone reuse elsewhere (e.g. rename_merge_by_gcf.py's
    interface); prefer AssemblyResolver.resolve() within this script since
    it fetches species_taxid in the same round trip."""
    search = subprocess.Popen(['esearch', '-db', 'assembly', '-query', nc_accession], stdout=subprocess.PIPE)
    summary = subprocess.Popen(['esummary'], stdin=search.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'RefSeq'],
                               stdin=summary.stdout, stdout=subprocess.PIPE, text=True)
    output, _ = xtract.communicate()
    return output


# --------------------------------------------------------------------
# Downloading missing genomes
# --------------------------------------------------------------------

def fetch_protein_fasta(nc_accession, email, api_key, tool='full-dataset-prep', retries=3):
    from urllib.parse import urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError

    params = {'db': 'nuccore', 'id': nc_accession, 'rettype': 'fasta_cds_aa', 'retmode': 'text', 'tool': tool}
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urlencode(params)}"

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={'User-Agent': f'{tool} ({email or "no-email-provided"})'})
            with urlopen(req, timeout=60) as resp:
                data = resp.read().decode('utf-8', errors='replace')
            if not data.startswith('>'):
                raise ValueError(f"Unexpected response (not FASTA) for {nc_accession}: {data[:200]}")
            return data
        except (HTTPError, URLError, ValueError) as e:
            last_err = e
            time.sleep(1.5 * attempt)
    raise RuntimeError(f"Failed to fetch {nc_accession} after {retries} attempts: {last_err}")


# --------------------------------------------------------------------
# Local database lookups
# --------------------------------------------------------------------

def load_lineage_db(tsv_path):
    """GCF accession -> full row dict (must include a 'species_id' column)."""
    db = {}
    with open(tsv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'species_id' not in reader.fieldnames or 'accession' not in reader.fieldnames:
            sys.exit(f"{tsv_path} is missing an 'accession' or 'species_id' column. "
                      f"Found columns: {reader.fieldnames}")
        for row in reader:
            db[row['accession'].strip()] = row
    return db


def genome_exists_in_db(gcf_accession, genomes_dir: Path):
    """Existing db files may be named exactly GCF_x.faa or with an extra
    assembly-name suffix (as in the ftp_path convention); match either."""
    if (genomes_dir / f"{gcf_accession}.faa").exists():
        return True
    return any(genomes_dir.glob(f"{gcf_accession}*"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table-s1', required=True, help='Table S1 file (.tsv/.csv/.xlsx) with the GenomeID/.../Accession columns')
    ap.add_argument('--db-lineage-tsv', required=True, help='Existing dbpath/BACTERIAL_LINEAGES.tsv')
    ap.add_argument('--db-genomes-dir', required=True, help='Existing dbpath/GENOMES/ directory')
    ap.add_argument('--new-genomes-dir', required=True, help='Separate output folder for genomes not already in the db')
    ap.add_argument('--mapping-output', required=True, help='Output TSV: nuccore_accession, refseq_accession, species_id, genome_id, organism, source')
    ap.add_argument('--email', default=None, help='Your email, for NCBI efetch (recommended)')
    ap.add_argument('--api-key', default=None, help='NCBI API key (optional, raises rate limits)')
    ap.add_argument('--cache', default='assembly_resolve_cache.json')
    ap.add_argument('--dry-run', action='store_true', help='Resolve and report, but do not download anything')
    args = ap.parse_args()

    table_path = Path(args.table_s1)
    if table_path.suffix.lower() in ('.xlsx', '.xls'):
        table = pd.read_excel(table_path, dtype=str)
    else:
        table = pd.read_csv(table_path, sep=None, engine='python', dtype=str)

    required_cols = {'GenomeID', 'Organism', 'Accession'}
    missing_cols = required_cols - set(table.columns)
    if missing_cols:
        sys.exit(f"Table S1 is missing expected column(s): {missing_cols}. Found: {list(table.columns)}")

    print(f"Loaded {len(table)} genome(s) from {table_path}")

    lineage_db = load_lineage_db(args.db_lineage_tsv)
    print(f"Loaded {len(lineage_db)} existing genome(s) from {args.db_lineage_tsv}")

    db_genomes_dir = Path(args.db_genomes_dir)
    new_genomes_dir = Path(args.new_genomes_dir)
    new_genomes_dir.mkdir(parents=True, exist_ok=True)

    resolver = AssemblyResolver(args.cache)

    rows_out = []
    n_existing = n_downloaded = n_unresolved = n_download_failed = 0

    for i, row in table.iterrows():
        nc_acc = str(row['Accession']).strip()
        genome_id = str(row['GenomeID']).strip()
        organism = str(row['Organism']).strip()

        gcf, species_taxid_fresh, note = resolver.resolve(nc_acc)
        print(f"[{i+1}/{len(table)}] {nc_acc} ({organism}) -> GCF={gcf}, species_taxid={species_taxid_fresh}  ({note})")

        if gcf is None:
            n_unresolved += 1
            rows_out.append({'nuccore_accession': nc_acc, 'refseq_accession': '', 'species_id': '',
                              'genome_id': genome_id, 'organism': organism, 'source': 'unresolved'})
            continue

        in_db = gcf in lineage_db and genome_exists_in_db(gcf, db_genomes_dir)
        if in_db:
            species_id = lineage_db[gcf].get('species_id', '').strip()
            n_existing += 1
            rows_out.append({'nuccore_accession': nc_acc, 'refseq_accession': gcf, 'species_id': species_id,
                              'genome_id': genome_id, 'organism': organism, 'source': 'existing_db'})
            continue

        # Not in the existing db -- needs downloading.
        species_id = species_taxid_fresh or ''
        out_fasta = new_genomes_dir / f"{gcf}.faa"
        if not args.dry_run:
            if out_fasta.exists() and out_fasta.stat().st_size > 0:
                print(f"    {out_fasta.name} already downloaded previously, skipping fetch")
            else:
                try:
                    data = fetch_protein_fasta(nc_acc, args.email, args.api_key)
                    out_fasta.write_text(data)
                    print(f"    downloaded -> {out_fasta.name} ({len(data):,} bytes)")
                except Exception as e:
                    print(f"    DOWNLOAD FAILED: {e}")
                    n_download_failed += 1
                    rows_out.append({'nuccore_accession': nc_acc, 'refseq_accession': gcf, 'species_id': species_id,
                                      'genome_id': genome_id, 'organism': organism, 'source': 'download_failed'})
                    continue
        n_downloaded += 1
        rows_out.append({'nuccore_accession': nc_acc, 'refseq_accession': gcf, 'species_id': species_id,
                          'genome_id': genome_id, 'organism': organism, 'source': 'downloaded_new'})

    with open(args.mapping_output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['nuccore_accession', 'refseq_accession', 'species_id',
                                                 'genome_id', 'organism', 'source'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"\n{'=' * 60}")
    print(f"Total genomes           : {len(table)}")
    print(f"Already in existing db  : {n_existing}")
    print(f"Downloaded as new       : {n_downloaded}" + (" (dry-run: not actually fetched)" if args.dry_run else ""))
    print(f"Download failed         : {n_download_failed}")
    print(f"Could not resolve GCF   : {n_unresolved}")
    print(f"Mapping written to      : {args.mapping_output}")
    if n_unresolved or n_download_failed:
        print(f"\nReview 'unresolved'/'download_failed' rows in the mapping file -- "
              f"these genomes are not yet usable with APP.")


if __name__ == '__main__':
    main()
