import pandas as pd
from pathlib import Path

accs = pd.read_csv("accs.tsv", sep="\t")

with open('cmds.txt', 'w') as f:
	for row in accs.itertuples():
		acc = row.accession
		if Path(f"/scratch/11587/siddharthk86/APP/18genometest3x/outputs/{acc}/APP_Alien_genes.tsv").is_file():
			continue
		sid = row.sid
		f.write(f"python3 src/main.py app {acc} /work/11587/siddharthk86/ls6/APP/GENOME_DB/ -d /scratch/11587/siddharthk86/APP/18genometest3x/outputs/{acc} -t 64 -M --sid {sid} -f fastas/{acc}.faa\n")
