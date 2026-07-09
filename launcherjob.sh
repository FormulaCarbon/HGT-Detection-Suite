#!/bin/bash
#SBATCH   -p development
#SBATCH   -J APP18genome
#SBATCH   -o APP18genome.o%j
#SBATCH   -e APP18genome.o%j
#SBATCH   -N 3
#SBATCH   -t 2:00:00


module load python/3.12.11
module load pylauncher
source .venv/bin/activate
pip install paramiko
python3 launcher.py
