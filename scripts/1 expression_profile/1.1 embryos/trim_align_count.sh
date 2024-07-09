#!/bin/bash

######################################### 0.init
function printHelp(){
    echo "Version: 1.0.0"
    echo
    echo
    echo "Required Arguments:"
    echo
    echo "--sample=File                      Sample file, separated by \",\",  (Default: sample.txt)"
    echo "--paired=Boolean                   Whether paired end, Possible values: {true or false} (Default: true)"
    echo
    echo
    echo "Optional Arguments:"
    echo
    echo "--species=character                'hs' for human. (Default: hs)"
    echo "--suffix=character                 Suffix of fastq files, 'fastq.gz' or 'fq.gz'. (Default: fastq.gz)"
    echo "--split=character                  Rule of how to name paired reads, a.1.fq.gz or a_1.fq.gz, '.' or '_'. (Default: _)"
}

for p in "$@"
do
    if [[ "$p" = *"help"* || "$p" = *"-h" ]]; then
        printHelp
        exit 88
    fi

    if [[ "$p" = "--sample"* ]]; then
        SAMPLE=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--paired"* ]]; then
        PAIRED=`echo $p | cut -d "=" -f 2`
        #echo $PAIRED
    fi

    if [[ "$p" = "--suffix"* ]]; then
       suffix=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--split"* ]]; then
      split=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--species"* ]]; then
        species=`echo $p | cut -d "=" -f 2`
    fi

    if [[ "$p" = "--thread"* ]]; then
        nThread=`echo $p | cut -d "=" -f 2`
    fi
done


######################################### 1.default variant
if [ ! -n "$SAMPLE" ]; then
    SAMPLE="sample.txt"
fi
if [ ! -n "$nThread" ]; then
    nThread=16
fi
if [ ! -n "$suffix" ]; then
    #suffix=".fastq.gz"
    suffix=".fq.gz"
fi
if [ ! -n "$sp" ]; then
    sp="_"
fi
if [ ! -n "$PAIRED" ]; then
   PAIRED="true"
fi
if [ ! -n "$species" ]; then
    species="hs"
fi

if [ "$species" = "hs" ]; then
    stargenome_R='/home/reference/Hsapiens/UCSC/hg38/Sequence/STARIndex/repeats'
    myanno_R='/home/reference/Hsapiens/UCSC/hg38/Annotation/repeat.hg38.gtf'
fi

######################################### 2.Test sample infomation
if [ ! -e $SAMPLE ];then
    printHelp
    echo "Erro: Sample information file is required, default file is sample.txt, default format is"
    echo "      \"sampleName\", if file name is sampleName.fq.gz or sampleName.fastq.gz"
    exit 99
else
    sampleLine=`cat $SAMPLE | wc -l`
    if [[ "$sampleLine" = "0" ]]; then
        echo "$SAMPLE is empty!"
        exit 100
    fi
fi

#echo sample infomation
echo
echo Start time: `date`
echo Working directory: `pwd`
#dir=`pwd`
#echo Fastq files directory: `pwd`/fastq

for sample in $(cat $SAMPLE)
do
    if [[ "$PAIRED" = "true" ]]; then
        if ([ ! -e ${sample}${sp}1${suffix} ] && [ ! -e ${sample}${sp}2${suffix} ]);then
            echo "${sample}${sp}1${suffix} or ${sample}${sp}2${suffix} do not exist!"
            exit 110
        fi
    else
        if [ ! -e ${sample}${suffix} ];then
            echo "${sample}${suffix} do not exist!"
            exit 110
        fi            
    fi
done

echo
echo "Data Processing..."

################################3.mkdir
#if [ ! -d fastqc ];then
#        mkdir -p fastqc
#fi
if [ ! -d trim ];then
    echo "mkdir trim"
    mkdir -p trim
fi
if [ ! -d starOut ];then
    echo "mkdir starOut"
    mkdir -p starOut
fi

if [ ! -d repeats.count ];then
    echo "mkdir repeats.count"
    mkdir -p repeats.count
fi

#################################

echo "Input species: $species "
echo "Input repeats STAR index: $stargenome_R "
echo "Input repeats annotation gtf: $myanno_R "

for line in $(cat $SAMPLE)
do
    sample=$(basename $line)
    echo "Start analysising for ${sample}..."

    TRIM_ARGS=""
    MAP_G_ARGS=""
    MAP_R_ARGS=""
    COUNT_ARGS=""

    if [[ "$PAIRED" = "true" ]]; then
        TRIM_ARGS+="--paired ${line}${sp}1${suffix} ${line}${sp}2${suffix}"
        MAP_R_ARGS+="--genomeDir $stargenome_R --readFilesIn trim/${sample}_1_val_1.fq.gz trim/${sample}_2_val_2.fq.gz --outFileNamePrefix starOut/${sample}.repeats."
        COUNT_ARGS+="-p"
    else
        TRIM_ARGS+="${line}${suffix}"
        MAP_R_ARGS+="--genomeDir $stargenome_R --readFilesIn trim/${sample}_trimmed.fq.gz --outFileNamePrefix starOut/${sample}.repeats."
    fi

    if ([[ "$TRIM_ARGS" = "--paired"* ]] && [ ! -e trim/${sample}_1_val_1.fq.gz ] && [ ! -e trim/${sample}_2_val_2.fq.gz ]) || ([[ "$TRIM_ARGS" != "--paired"* ]] && [ ! -e trim/${sample}_trimmed.fq.gz ]); then      
		echo "Trimming adaptor and low-quality reads..."
		trim_galore $TRIM_ARGS --gzip -o trim
    fi

    if [  ! -e starOut/${sample}.repeats.Aligned.sortedByCoord.out.bam ];then
        echo "Mapping to reference repeats genome..."
        STAR $MAP_R_ARGS \
--outFilterMismatchNoverLmax 0.04 \
--runThreadN $nThread \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 50000000000 \
--outFilterMultimapNmax 500 \
--outMultimapperOrder Random \
--outSAMmultNmax 1 \
--readFilesCommand zcat
    fi
    
    if [ ! -e repeats.count/${sample}.repeats.count ];then
        featureCounts $COUNT_ARGS -M -a $myanno_R -T $nThread -o repeats.count/${sample}.repeats.count starOut/${sample}.repeats.Aligned.sortedByCoord.out.bam
        awk 'BEGIN{FS=OFS="\t"}!/#/{print $1,$6,$7}' repeats.count/${sample}.repeats.count > repeats.count/${sample}.repeats.count.simple
    fi    

done
