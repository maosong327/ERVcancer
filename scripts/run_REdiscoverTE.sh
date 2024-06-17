#!/bin/bash

# Extract fastq files from TCGA BAM file
bam=/home/public/Project/TCGA/2.READ/db82542f-db4b-4e59-8555-1765713e32cf/6a29d5a3-a6a7-4bf2-b0f8-973fa4bc1ce0_gdc_realn_rehead.bam
base=$(basename $bam)
sample=${base/_gdc_realn_rehead.bam/}

if [ ! -d ${sample} ];then
  echo "mkdir output dir"
  mkdir -p ${sample}
fi

# test library type
paired=$(samtools view -H $bam | grep "_2.fastq" | wc -w)
FASTQ_ARGS=""
if [[ "$paired" != "0" ]];then
  echo "${bam} is paired end"
  FASTQ_ARGS+="-1 ${sample}/${sample}_R1.fq.gz -2 ${sample}/${sample}_R2.fq.gz -0 /dev/null -s /dev/null"
else
  echo "${bam} is single end"
  FASTQ_ARGS+="-o ${sample}/${sample}_R1.fq.gz"
fi
echo $FASTQ_ARGS

samtools sort -@ 16 -n $bam | samtools fastq -@ 16 $FASTQ_ARGS -n -

# run RE_discoverTE
RE_discoverTE_PATH='/home/maosong/biosoft/REdiscoverTE'
mkdir -p ${PWD}/${sample}/result
R_1="$PWD"/"${sample}"/"${sample}"_R1.fq.gz
R_2="$PWD"/"${sample}"/"${sample}"_R2.fq.gz
echo $R_1
echo $R_2
cp "${RE_discoverTE_PATH}"/Makefile_backup "${RE_discoverTE_PATH}"/Makefile
# editting makefile
sed -i "s|FASTQ_READS_1=SIMULATED_FASTQS/input_R1.fq.gz|FASTQ_READS_1=$R_1|g" "${RE_discoverTE_PATH}"/Makefile
sed -i "s|FASTQ_READS_2=SIMULATED_FASTQS/input_R2.fq.gz|FASTQ_READS_2=$R_2|g" "${RE_discoverTE_PATH}"/Makefile
sed -i "s|SALMON_COUNTS_DIR=Step_2_salmon_counts|SALMON_COUNTS_DIR=${PWD}/${sample}/result|g" "${RE_discoverTE_PATH}"/Makefile
sed -i "s|ROLLUP_RESULTS_DIR=Step_4_rollup|ROLLUP_RESULTS_DIR=${PWD}/${sample}/result|g" "${RE_discoverTE_PATH}"/Makefile
if [[ "$paired" != "0" ]];then
  echo "fastq files are paired end"
else
  echo "fastq file is single end"
  sed -i "s|FASTQ_PAIRED_TYPE=PAIRED|FASTQ_PAIRED_TYPE=SINGLE|g" "${RE_discoverTE_PATH}"/Makefile  
fi
ORIG_DIR=$PWD

cd $RE_discoverTE_PATH
make all
mv "${RE_discoverTE_PATH}"/Makefile "$ORIG_DIR"/"${sample}"/result

# featureCounts 
cd $ORIG_DIR
COUNT_ARGS=""
if [[ "$paired" != "0" ]];then
  COUNT_ARGS+="-p"
fi
echo $COUNT_ARGS
# featureCounts quantify with unique mapping reads
featureCounts $COUNT_ARGS -a /home/reference/Hsapiens/UCSC/hg38/Annotation/repeat.hg38.gtf -T 8 -o ${sample}.unique.repeats.counts $bam
awk 'BEGIN{FS=OFS="\t"}!/#/{print $1,$6,$7}' ${sample}.unique.repeats.counts >${sample}.unique.repeats.counts.simple

# featureCounts quantify multiple mapping reads
featureCounts $COUNT_ARGS -M --fraction -a /home/reference/Hsapiens/UCSC/hg38/Annotation/repeat.hg38.gtf -T 8 -o ${sample}.multi.repeats.counts $bam
awk 'BEGIN{FS=OFS="\t"}!/#/{print $1,$6,$7}' ${sample}.multi.repeats.counts >${sample}.multi.repeats.counts.simple