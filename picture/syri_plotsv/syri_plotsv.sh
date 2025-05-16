## https://github.com/nschan/nf-plotsv 

ref_genome=Chr.male.fasta
ref_genome_chr_num=19
query_genome=manli_telo.chr_genome.fa
query_genome_chr_num=19
cpu_for_minimap=4
mem_for_minimap=20

export PATH="/dellfsqd2/ST_OCEAN/USER/lishuo1/00_tools/:$PATH"

## rename chromesome
iTools Fatools stat -InPut ${ref_genome} -OutPut ${ref_genome}.chrlen
iTools Fatools stat -InPut ${query_genome} -OutPut ${query_genome}.chrlen
grep -v "#" ${ref_genome}.chrlen | sort -k 2nr,2 | head -n ${ref_genome_chr_num} | awk '{print $1}' | seqtk subseq ${ref_genome} - | seqtk rename - Chr | seqkit seq -w 100 - > ref.rename.fa
grep -v "#" ${query_genome}.chrlen | sort -k 2nr,2 | head -n ${query_genome_chr_num} | awk '{print $1}' | seqtk subseq ${query_genome} - | seqtk rename - Chr | seqkit seq -w 100 - > query.rename.fa

## make samplesheet
echo "name,fasta" >> samplesheet.csv
echo "genome1,$PWD/ref.rename.fa" >> samplesheet.csv
echo "genome2,$PWD/query.rename.fa" >> samplesheet.csv

## run
cp base.demo.config base.config
sed -i "s/ABCD/${cpu_for_minimap}/g;s/EFGH/${mem_for_minimap}/g" base.config
nextflow -config $PWD/base.config run /dellfsqd2/ST_OCEAN/USER/lishuo11/01_soft/nf-plotsv --samplesheet samplesheet.csv -profile local --reference ref_genome --ref_genome $PWD/ref.rename.fa --subset_pattern Chr[1-9] --reorient true

## rm tmp
rm -rf .nextflow* work ref.rename.fa query.rename.fa *chrlen *.chrlist  samplesheet.csv base.config
