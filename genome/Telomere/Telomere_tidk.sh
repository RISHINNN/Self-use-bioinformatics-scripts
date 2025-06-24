genome=Stichopus_variegatus.fa
out_dir=test1

## download tidk_database.csv to ~/.local/share/tidk
#/software/miniconda3/bin/tidk build

[ -d ${out_dir} ] || mkdir ${out_dir}

## explore
/software/miniconda3/bin/tidk explore --distance 0.05 --minimum 5 --maximum 7 -t 30 ${genome} > ${out_dir}/candicate_TR_unit.tsv

## find
for i in $(cat ${out_dir}/candicate_TR_unit.tsv | sed '1d' | awk '{print $1}'); do /software/miniconda3/bin/tidk search -s ${i} -o find_${i} -d ${out_dir} -w 50000 ${genome} ; done

## check
for i in $(ls ${out_dir}/*_windows.tsv) ; do awk '$3>50 || $4>50' ${i} > ${i}.check.tsv ; done
