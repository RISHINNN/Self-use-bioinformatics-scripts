#!/bin/bash
set +o posix

genome=Stichopus_variegatus.fa
bgi_gff=Stichopus_variegatus.bgi.gff
chr_id=chr.id
repeatmask_out=HiTE.update.out

## get chromosome information
seqtk subseq ${genome} ${chr_id} | seqkit seq -w 100 > tmp_chr.fa
fishInWinter.pl -bf table -ff table ${chr_id} ${bgi_gff} > tmp_chr.gff

## stat
iTools Fatools stat -InPut ${genome} -OutPut ${genome}.chrlist
grep -v "#" ${genome}.chrlist | grep chr | awk '{print "chr - "$1" "$1" 0 "$2" "$1}' > karyotype.txt
awk '$3=="mRNA"' tmp_chr.gff | awk '{print $1"\t"$4"\t"$5}' | sort -k 1V,1 -k 2n,2 > gene.bed
cut -f 3,6 -d " " karyotype.txt | awk '{print $1"\t"$2}' > chr.length

## make windows (no more 2500 windows)
bedtools makewindows -g chr.length -n 70 > chr.window

## stat
bedtools coverage -a chr.window -b gene.bed | cut -f 1-4 | sort -k 1V,1 -k 2n,2 > gene_density.txt
bedtools nuc -fi tmp_chr.fa -bed chr.window | cut -f 1,2,3,5 | sed '1d' > GC_content.txt

## deal TE
grep "DNA/" ${repeatmask_out} | awk '{print $5"\t"$6"\t"$7}' | fishInWinter.pl -bf table -ff table ${chr_id} - |sort -k 1V,1 -k 2n,2 | bedtools merge -i - > DNA.TE.bed
grep "LINE/" ${repeatmask_out} | awk '{print $5"\t"$6"\t"$7}' | fishInWinter.pl -bf table -ff table ${chr_id} - |sort -k 1V,1 -k 2n,2 | bedtools merge -i - > LINE.TE.bed
grep "SINE/" ${repeatmask_out} | awk '{print $5"\t"$6"\t"$7}' | fishInWinter.pl -bf table -ff table ${chr_id} - |sort -k 1V,1 -k 2n,2 | bedtools merge -i - > SINE.TE.bed
grep "LTR/" ${repeatmask_out} | awk '{print $5"\t"$6"\t"$7}' | fishInWinter.pl -bf table -ff table ${chr_id} - |sort -k 1V,1 -k 2n,2 | bedtools merge -i - > LTR.TE.bed
bedtools coverage -a chr.window -b DNA.TE.bed | cut -f 1-4 | sort -k 1V,1 -k 2n,2 > DNA_TE_density.txt
bedtools coverage -a chr.window -b LINE.TE.bed | cut -f 1-4 | sort -k 1V,1 -k 2n,2 > LINE_TE_density.txt
bedtools coverage -a chr.window -b SINE.TE.bed | cut -f 1-4 | sort -k 1V,1 -k 2n,2 > SINE_TE_density.txt
bedtools coverage -a chr.window -b LTR.TE.bed | cut -f 1-4 | sort -k 1V,1 -k 2n,2 > LTR_TE_density.txt

rm tmp_chr.fa tmp_chr.gff DNA.TE.bed LINE.TE.bed SINE.TE.bed LTR.TE.bed *.fai

## circos
cp /02_pilple/circos/circos.conf .
cp /02_pilple/circos/ticks.conf .

/usr/bin/singularity exec --bind $PWD:$PWD /01_soft/singularity_all/circos.sif circos -conf circos.conf
