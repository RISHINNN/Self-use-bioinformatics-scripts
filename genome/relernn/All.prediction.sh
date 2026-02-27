source activate /micromamba/envs/ReLERNN
SIMULATE="ReLERNN_SIMULATE"
SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="6e-9"
GENTIME="1"
URTR="1"
DIR="./ABCD_output/"
VCF="./ABCD.vcf"
GENOME="./S_maxi.fa.fai.bed"
CPU="10"
Maxwinsize="500"
Minsites="10"
batchsize="10"
#MASK="./accessibility_mask.bed"

## prepare
#awk '$2>5000000' S_maxi.fa.fai | awk '{print $1"\t0\t"$2}' > S_maxi.fa.fai.bed
bcftools view -S ABCD -m 2 -M 2 S_maxi.SNP.filter.vcf.gz | vcftools --vcf - --maf 0.05 --max-maf 0.95 --max-missing 0.1 --stdout --recode | bcftools annotate --remove QUAL,FILTER,INFO,^FORMAT/GT  | grep -v contig > ABCD.vcf

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --nTrain 13000 \
    --nVali 2000 \
    --nTest 100 \
    --forceDiploid  \
    --maxSites ${Maxwinsize} \
    -t ${CPU} \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --nEpochs 2 \
    --nValSteps 2 \
    -t ${CPU} \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED} \
    --minSites ${Minsites} \
    --batchSizeOverride ${batchsize} \
    --phased

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nSlice 2 \
    --nReps 2 \
    --seed ${SEED} \
    -t ${CPU}

## remove tmp
rm ${VCF}
rsync --delete-before -a /rsync_tmp/ ${DIR}/train/
rsync --delete-before -a rsync_tmp/ ${DIR}/splitVCFs/
rsync --delete-before -a /rsync_tmp/ ${DIR}/vali/
rsync --delete-before -a /rsync_tmp/ ${DIR}/test/
rsync --delete-before -a /rsync_tmp/ ${DIR}/networks/
rm -rf ${DIR}/train/ ${DIR}/splitVCFs/ ${DIR}/vali/ ${DIR}/test/ ${DIR}/networks/
