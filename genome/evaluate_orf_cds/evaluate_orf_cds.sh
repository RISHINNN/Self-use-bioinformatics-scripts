input=helixer.bgi.gff.cds

## CPC2
# https://github.com/gao-lab/CPC2_standalone
/01_software/CPC/CPC2_standalone-1.0.1/bin/CPC2.py -i ${input} -o ${input}.CPC.txt

## PSAURON
# https://github.com/salzberg-lab/PSAURON
# Note: internal stop codons are ignored by PSAURON. A high PSAURON score does not guarantee a sequence contains a valid ORF. This is intended behavior, as alternate frame scores are used by default to boost the power of the model.
export MAMBA_EXE='/01_soft/mamba/bin/micromamba'
export MAMBA_ROOT_PREFIX='/home_micromamba'
micromamba activate
source activate /home_micromamba/envs/psauron

psauron -i ${input} -o ${input}.PSAURON.csv
