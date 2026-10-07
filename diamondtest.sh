#!/bin/bash
#SBATCH   -p development
#SBATCH   -J APP_dmnd_test
#SBATCH   -o APP_dmnd_test.o%j
#SBATCH   -e APP_dmnd_test.o%j
#SBATCH   -N 1
#SBATCH   -t 2:00:00


module load python/3.12.11
module load pylauncher
source .venv/bin/activate
pip install paramiko


python3 launcher.py
