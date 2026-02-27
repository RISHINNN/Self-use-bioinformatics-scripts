#!/bin/bash

genome=Stichopus_variegatus.fa
trans=transcripts.fa
pep=proteins.faa
cpu=15

/01_software/EviAnn-2.0.2/bin/eviann.sh -t $cpu -g $genome -e $PWD/$trans -p $PWD/$pep --partial --debug -l

[ -d out_final ] || mkdir out_final
mv ${genome}.pseudo_label.gff out_final
mv ${genome}.transcripts.fasta out_final
mv ${genome}.proteins.fasta out_final
cp *.sh.o* out_final/log.txt

#rm ${genome}.* tissue* *stringtie*.sh *success broken* makeblastdb.out blastp2.out combine.out makeblastdb.sex2mex.out blastp5.out proteins.faa.uniq miniprot.err check_cds.out
