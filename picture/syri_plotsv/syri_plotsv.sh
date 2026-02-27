## https://github.com/nschan/nf-plotsv 

## Use relative paths
ref_genome=final-father-chr-mela.fa
ref_genome_chr_num=24
query_genome=final-mother-chr-mela.fa
query_genome_chr_num=24
cpu_for_minimap=20
mem_for_minimap=100
plot_all_chr_num=24
plotsr_Space_for_homologous_chromosome=0.7
plotsr_height=10
plotsr_width=7
plotsr_font_size=8
plotsr_minimum_size_of_SR_to_be_plotted=5000

export PATH="/00_tools/:$PATH"

## rename chromesome
iTools Fatools stat -InPut ${ref_genome} -OutPut ${ref_genome}.chrlen
iTools Fatools stat -InPut ${query_genome} -OutPut ${query_genome}.chrlen
grep -v "#" ${ref_genome}.chrlen | sort -k 2nr,2 | head -n ${ref_genome_chr_num} | awk '{print $1}' | seqtk subseq ${ref_genome} - | seqtk rename - Chr | seqkit seq -w 100 - > ref.rename.fa
grep -v "#" ${query_genome}.chrlen | sort -k 2nr,2 | head -n ${query_genome_chr_num} | awk '{print $1}' | seqtk subseq ${query_genome} - | seqtk rename - Chr | seqkit seq -w 100 - > query.rename.fa

## make samplesheet
echo "name,fasta" >> samplesheet.csv
## the names also be used in the final picture.
echo "ref,$PWD/ref.rename.fa" >> samplesheet.csv
echo "query,$PWD/query.rename.fa" >> samplesheet.csv

## run
cp /01_soft/nf-plotsv/configs/base.demo.config base.config
sed -i "s/ABCD/${cpu_for_minimap}/g;s/EFGH/${mem_for_minimap}/g" base.config
nextflow -config $PWD/base.config run /01_soft/nf-plotsv --samplesheet samplesheet.csv -profile local --reference ref_genome --ref_genome $PWD/ref.rename.fa --subset_pattern Chr[1-9] --reorient true

## re-draw
cp work/*/*/plotsr_infile.tsv plotsv/syri_pairwise/plotsr_infile.tsv
cp work/*/*/files.txt plotsv/syri_pairwise/files.txt
cp /nf-plotsv/assets/plotsr_config.conf plotsv/syri_pairwise/plotsr_config.conf
for i in `seq 1 ${plot_all_chr_num}` ; do echo Chr${i} >> plotsv/syri_pairwise/chr.order ; done

cd plotsv/syri_pairwise/
cp ../align_pairwise/*fa ./
sr=$(cat files.txt)
singularity run /01_soft/singularity_all/fixchr-syri-plotsr.sif plotsr --genomes plotsr_infile.tsv ${sr} --cfg plotsr_config.conf -o replot.pdf -S ${plotsr_Space_for_homologous_chromosome} -W ${plotsr_width} -H ${plotsr_height} -f ${plotsr_font_size} -s ${plotsr_minimum_size_of_SR_to_be_plotted} --chrord chr.order
pigz --best -p ${cpu_for_minimap} *syri.out
pigz --best -p ${cpu_for_minimap} *syri.vcf
pigz --best -p ${cpu_for_minimap} *.fa
cd ../../

## rm tmp
rm -rf .nextflow* work ref.rename.fa query.rename.fa *chrlen *.chrlist  samplesheet.csv base.config plotsv/prepare_genomes plotsv/align_pairwise
