
indel=H_ota.indel.vcf.gz
genome=Hexagrammos_otakii_genomic.fna
male_id=male.list
female_id=female.list
out1=sex_specific_indels.txt
out2=primers.txt
tolerance=0.2                       # 容错率 (0.1 表示允许 10% 的样本不符合规律)
maxgap=400                          # indel聚类距离
min_diff=30                         # 引物最小净长度差异

python indel_find.py --vcf ${indel} --males ${male_id} --females ${female_id} --out ${out1} -t ${tolerance}
python design_primers.py -i ${out1} -g ${genome} -o ${out2} --max_gap ${maxgap} --min_net_diff ${min_diff}
