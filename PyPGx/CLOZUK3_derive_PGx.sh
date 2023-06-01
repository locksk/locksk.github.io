#!/bin/bash

#SBATCH --account=scw2063
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=02:00:00

echo "Setting up"
# set up
module load python/3.9.2
echo "..."
source ve_pypgx/bin/activate
echo "..."
pip install -r ve_pypgx/requirements.txt
echo "complete!"
echo
echo
echo "Running PyPGx..."
echo
echo "CYP1A2 starting..."
echo
# run cyp1a2
pypgx run-chip-pipeline CYP1A2 CYP1A2-pipeline data/CLOZUK3.chr15.dose.pgx4.vcf.gz

mv CYP1A2-pipeline output/
cd output/CYP1A2-pipeline 
bash ../../unzip_script.sh
cp results/*/data.tsv ../results/cyp1a2.tsv
cd ..
cd ..
echo "complete!"
echo
echo



