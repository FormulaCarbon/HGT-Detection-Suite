#!/bin/bash
#SBATCH   -p normal
#SBATCH   -J APPIP
#SBATCH   -o APPIP.o%j
#SBATCH   -e APPIP.o%j
#SBATCH   -N 8
#SBATCH   -t 4:00:00


module load python/3.12.11
module load pylauncher
source .venv/bin/activate
pip install paramiko
python3 islandpick/launcher.py
