#!/bin/bash
set +o posix
set -eo pipefail

export PATH="/01_software/lastal/bin/:$PATH"
export PATH="/01_software/latex/bin/x86_64-linux/:$PATH"
python=/01_soft/mamba/env/jcvi/bin/python

## use ${sp1}.bed ${sp2}.bed ${sp1}.pep ${sp2}.pep ${sp1}.genome.fa ${sp2}.genome.fa
sp1=O_curv
sp2=O_mela

### dotplot
${python} -m jcvi.compara.catalog ortholog --dbtype prot --cpus=1 --no_strip_names ${sp1} ${sp2}

### make .simple
${python} -m jcvi.compara.synteny screen --simple ${sp1}.${sp2}.anchors ${sp1}.${sp2}.anchors.new
simpletolink.py ${sp1}.${sp2}.anchors.simple

### use ${sp1}.genome.fa ${sp2}.genome.fa
iTools Fatools stat -InPut ${sp1}.genome.fa -OutPut ${sp1}.genome.fa.chrlist
iTools Fatools stat -InPut ${sp2}.genome.fa -OutPut ${sp2}.genome.fa.chrlist

grep -v "#" ${sp1}.genome.fa.chrlist | cut -f 1,2 | fishInWinter.pl -bf table -ff table <( cut -f 1 ${sp1}.bed | awk '!a[$0]++' ) - | awk '{print "chr - " "'"$sp1"'" "_" $1 " " "'"$sp1"'" "_" $1 " 0 " $2 " " "chr"NR}' > karyotype_sp1.txt
grep -v "#" ${sp2}.genome.fa.chrlist | cut -f 1,2 | fishInWinter.pl -bf table -ff table <( cut -f 1 ${sp2}.bed | awk '!a[$0]++' ) - | awk '{print "chr - " "'"$sp2"'" "_" $1 " " "'"$sp2"'" "_" $1 " 0 " $2 " " "chr"NR}' > karyotype_sp2.txt

awk '{print "'"$sp1"'" "_" $1 "\t" $2 "\t" $3 "\t" "'"$sp2"'" "_" $4 "\t" $5 "\t" $6}' ${sp1}.${sp2}.anchors.simple_link.txt > anchors.simple_link.rename.txt

generate_circos_configs.py anchors.simple_link.rename.txt
cp circos_config/* .
/usr/bin/singularity exec --bind $PWD:$PWD /01_soft/singularity_all/circos.sif circos -conf circos_config_output/circos.conf

rm *.ssp *.tis *.sds *.des *.prj *.suf *.bck *.chrlist
