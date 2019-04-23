#!/bin/bash

#SBATCH --job-name=all
#SBATCH --mem-per-cpu=32G
#SBATCH --account=def-jpviroid
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /project/6003961/vana2406/bin/virtenv/bin/activate

python /project/6003961/vana2406/Scripts/processingG4/Parser_Fasta $1
