## 主要功能
- 支持三种主流多序列比对工具：PRANK、MUSCLE和MAFFT
- 基于蛋白质序列引导的DNA密码子感知比对
- 自动提取4D兼并位点（四重简并位点）
- 多种重复物种处理策略：基于比对质量、最长序列、首个序列等
- 多种缺失数据处理策略：缺口填充、物种排除、基因排除
- 自动翻译CDS序列到蛋白质序列
- 支持TrimAl修剪比对结果
- 并行处理以提高效率
- 详细的日志和统计信息输出

## 完整参数列表
- 必要参数:  
  --input_dir DIR          包含CDS文件的目录 (默认: ".")
   
  --output_dir DIR         输出文件目录 (默认: "./output")
    
  --aligner {prank,muscle,mafft}    
                           使用的比对工具 (默认: "prank")
  
  --prank_path PATH        PRANK可执行文件的绝对路径 (选择prank时必需)
  
  --muscle_path PATH       MUSCLE可执行文件的绝对路径 (选择muscle时必需)
  
  --mafft_path PATH        MAFFT可执行文件的绝对路径 (选择mafft时必需)   
  
- 常用选项:  
  --supergene_output FILE  超基因输出文件名 (默认: "supergene_4d.fasta")
  
  --threads N              并行处理的线程数 (默认: 4)
  
  --no_codon_aware         禁用密码子感知比对 (默认启用)
   
  --duplicate_strategy {longest,first,rename,alignment_quality}
    
                           处理重复物种的策略 (默认: alignment_quality)
   
  --skip_existing          如果比对文件已存在则跳过处理
  
  --min_coverage_pct N     物种必须存在的最低基因百分比 (默认: 50.0%)
  
  --log_level {DEBUG,INFO,WARNING,ERROR}
                           设置日志级别 (默认: INFO)   

- TrimAl相关选项:
  --use_trimal            使用TrimAl修剪蛋白质比对
  
  --trimal_path PATH      TrimAl可执行文件的绝对路径
  
  --trimal_automated      使用TrimAl自动化修剪方法 (默认: True)
  
  --gap_threshold N       TrimAl最小缺口阈值
  
  --consistency_threshold N  TrimAl一致性阈值  
  
  --conservation_threshold N   TrimAl保守性阈值
   
  --trim_supergene        对最终蛋白质超基因应用TrimAl   
  
- 高级选项:  
  --f N                   PRANK插入开放概率 (默认: 0.2)
  
  --gaprate N             PRANK缺口开放率
  
  --gapext N              PRANK缺口扩展概率
  
  --use_logs              在PRANK中使用对数计算(大数据集)
  
  --penalize_terminal_gaps   
                          在PRANK中正常惩罚末端缺口   

  --clean_temp            处理后清理临时文件 (默认: True)
  
  --create_protein_msa    创建蛋白质序列的多序列比对   

## 主要结果文件说明
4d_sites/supergene_4d_*.fasta: 4D兼并位点超基因，用于系统发育分析  
  
full_cds/supergene_full_*.fasta: 完整CDS序列的超基因  
  
proteins/supergene_protein_*.fasta: 翻译的蛋白质序列超基因   
   
stats/species_coverage_matrix_*.tsv: 物种/基因覆盖矩阵，显示每个物种在各基因中的存在情况   

## 输出
output_dir/  
├── 4d_sites/               # 4D位点序列  
│   ├── supergene_4d_gaps.fasta             # 用缺口填充策略的超基因  
│   ├── supergene_4d_exclude_species.fasta   # 排除物种策略的超基因  
│   └── supergene_4d_exclude_genes.fasta     # 排除基因策略的超基因  
├── alignments/             # 各基因的比对结果  
│   ├── gene1.best.fas  
│   ├── gene2.best.fas  
│   └── ...  
├── full_cds/               # 完整CDS序列的超基因   
│   ├── supergene_full_gaps.fasta   
│   ├── supergene_full_exclude_species.fasta  
│   └── supergene_full_exclude_genes.fasta   
├── proteins/               # 翻译的蛋白质序列   
│   ├── supergene_protein_gaps.fasta   
│   ├── supergene_protein_exclude_species.fasta   
│   └── supergene_protein_exclude_genes.fasta   
├── protein_msa/            # 蛋白质多序列比对结果（如果使用--create_protein_msa,这个参数废弃，结果已经在proteins/ 中）  
│   ├── gene1_protein_msa.fasta  
│   ├── gene2_protein_msa.fasta  
│   ├── supergene_protein_msa_gaps.fasta  
│   └── ...  
├── stats/                  # 统计信息  
│   ├── species_coverage_matrix_gaps.tsv  
│   ├── species_coverage_matrix_exclude_species.tsv  
│   └── species_coverage_matrix_exclude_genes.tsv  
└── temp/                   # 临时文件（处理完成后可能被删除）   
