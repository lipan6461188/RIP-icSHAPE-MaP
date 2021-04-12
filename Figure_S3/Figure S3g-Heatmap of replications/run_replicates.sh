#!/bin/env bash

if [ $# != 4 ]; then
    echo -e "Usage: $0 inputdir dmso_file nai_file outputdir"
    exit;
fi

IN_DIR=$1
DMSO_FILE=$2
NAI_FILE=$3
PROCESS_DIR=$4

if [ ! -f $IN_DIR/${DMSO_FILE} ] || [ ! -f $IN_DIR/${NAI_FILE} ]; then
    echo -e "$IN_DIR/${DMSO_FILE} and $IN_DIR/${NAI_FILE} should be in folder"
    exit;
fi

mkdir -p $PROCESS_DIR
mkdir -p $PROCESS_DIR/mapping
mkdir -p $PROCESS_DIR/parseMut
mkdir -p $PROCESS_DIR/splitMut
mkdir -p $PROCESS_DIR/final_shape
mkdir -p $PROCESS_DIR/shape_files

dmso_samples=(${DMSO_FILE%%.*})
nai_samples=(${NAI_FILE%%.*})

errorFn=$PROCESS_DIR/error
logFn=$PROCESS_DIR/log

function job_id() { echo $2 | gawk 'match($0, /<([0-9]+)>/, id){print id[1]}'; }

###  1. Trim adaptor

jobIDs=()

for dmso_sample in ${dmso_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "cutadapt -a AGATCGGAAGAGCACACGTCTG -m 20 $IN_DIR/${dmso_sample}.fastq.gz -o $PROCESS_DIR/${dmso_sample}.trimmed.fastq.gz")
    jobIDs+=("$(job_id $info)")
done

for nai_sample in ${nai_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "cutadapt -a AGATCGGAAGAGCACACGTCTG -m 20 $IN_DIR/${nai_sample}.fastq.gz -o $PROCESS_DIR/${nai_sample}.trimmed.fastq.gz")
    jobIDs+=("$(job_id $info)")
done

for jobid in ${jobIDs[@]}; 
do
    echo -e "cutadapt bwait $jobid"
    bwait $jobid
done

### 2. Quanlity Trim

jobIDs=()

trimmomaticBIN=/Share2/home/zhangqf/shaodi/app/Trimmomatic-0.33/trimmomatic-0.33.jar
for dmso_sample in ${dmso_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "java -Xmx30g -jar ${trimmomaticBIN} SE -threads 20 -phred33 $PROCESS_DIR/${dmso_sample}.trimmed.fastq.gz $PROCESS_DIR/${dmso_sample}.trimmed.qual.fastq TRAILING:20 MINLEN:20")
    jobIDs+=("$(job_id $info)")
done

for nai_sample in ${nai_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "java -Xmx30g -jar ${trimmomaticBIN} SE -threads 20 -phred33 $PROCESS_DIR/${nai_sample}.trimmed.fastq.gz $PROCESS_DIR/${nai_sample}.trimmed.qual.fastq TRAILING:20 MINLEN:20")
    jobIDs+=("$(job_id $info)")
done

for jobid in ${jobIDs[@]}; 
do
    echo -e "trimmomatic bwait $jobid"
    bwait $jobid
done

### 3. Unique

jobIDs=()

readCollapseBin=/Share2/home/zhangqf7/lipan/usr/icSHAPE/scripts/readCollapse.pl
for dmso_sample in ${dmso_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "perl ${readCollapseBin} -U $PROCESS_DIR/${dmso_sample}.trimmed.qual.fastq -o $PROCESS_DIR/${dmso_sample}.rmdup.trimmed.qual.fastq")
    jobIDs+=("$(job_id $info)")
done

for nai_sample in ${nai_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "perl ${readCollapseBin} -U $PROCESS_DIR/${nai_sample}.trimmed.qual.fastq -o $PROCESS_DIR/${nai_sample}.rmdup.trimmed.qual.fastq")
    jobIDs+=("$(job_id $info)")
done

for jobid in ${jobIDs[@]}; 
do
    echo -e "readCollapse bwait $jobid"
    bwait $jobid
done

### 4. Head crop

jobIDs=()

for dmso_sample in ${dmso_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "cutadapt -u 8 -m 20 $PROCESS_DIR/${dmso_sample}.rmdup.trimmed.qual.fastq -o $PROCESS_DIR/${dmso_sample}.clean.fastq.gz")
    jobIDs+=("$(job_id $info)")
done

for nai_sample in ${nai_samples[@]};
do
    info=$(bsub -q Z-ZQF -e ${errorFn} -o ${logFn} "cutadapt -u 8 -m 20 $PROCESS_DIR/${nai_sample}.rmdup.trimmed.qual.fastq -o $PROCESS_DIR/${nai_sample}.clean.fastq.gz")
    jobIDs+=("$(job_id $info)")
done

for jobid in ${jobIDs[@]}; 
do
    echo -e "cutadapt-2 bwait $jobid"
    bwait $jobid
done


### 5. Mapping

function mapping_NAI_gz()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    STAR --readFilesIn $INPUT \
        --outFileNamePrefix $OUTPUT \
        --genomeDir $REF \
        --runThreadN 20 \
        --genomeLoad NoSharedMemory \
        --runMode alignReads \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 \
        --scoreGap -1000 \
        --alignEndsType Local \
        --outSAMattributes All \
        --outFilterMultimapNmax 10 \
        --outMultimapperOrder Random \
        --outFilterMismatchNmax 3 \
        --scoreDelBase -1 \
        --scoreInsBase -1 \
        --limitBAMsortRAM 999999999999 \
        --outFilterMismatchNoverLmax 0.3 \
        --readFilesCommand "gzip -dc"
}

function mapping_DMSO_gz()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    STAR --readFilesIn $INPUT \
        --outFileNamePrefix $OUTPUT \
        --genomeDir $REF \
        --runThreadN 20 \
        --genomeLoad NoSharedMemory \
        --runMode alignReads \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 \
        --scoreGap -1000 \
        --alignEndsType Local \
        --outSAMattributes All \
        --limitBAMsortRAM 999999999999 \
        --outFilterMultimapNmax 10 \
        --outMultimapperOrder Random \
        --outFilterMismatchNmax 2 \
        --outFilterMismatchNoverLmax 0.3 \
        --readFilesCommand "gzip -dc"
}

export -f mapping_NAI_gz
export -f mapping_DMSO_gz
REF=/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/

jobIDs=()

for dmso_sample in ${dmso_samples[@]};
do
    info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "mapping_DMSO_gz $REF $PROCESS_DIR/${dmso_sample}.clean.fastq.gz $PROCESS_DIR/mapping/${dmso_sample}.")
    jobIDs+=("$(job_id $info)")
done

for nai_sample in ${nai_samples[@]};
do
    info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "mapping_NAI_gz $REF $PROCESS_DIR/${nai_sample}.clean.fastq.gz $PROCESS_DIR/mapping/${nai_sample}.")
    jobIDs+=("$(job_id $info)")
done

for jobid in ${jobIDs[@]}; 
do
    echo -e "STAR bwait $jobid"
    bwait $jobid
done


### 6. Merge BAM files

mv $PROCESS_DIR/mapping/${DMSO_FILE%%.*}.Aligned.sortedByCoord.out.bam $PROCESS_DIR/RIP_DMSO_repX.bam
mv $PROCESS_DIR/mapping/${NAI_FILE%%.*}.Aligned.sortedByCoord.out.bam $PROCESS_DIR/RIP_NAI_repX.bam

### 7. Parse mutations

samtools view -h $PROCESS_DIR/RIP_DMSO_repX.bam > $PROCESS_DIR/parseMut/RIP_DMSO_repX.sam 
samtools view -h $PROCESS_DIR/RIP_NAI_repX.bam > $PROCESS_DIR/parseMut/RIP_NAI_repX.sam

jobIDs=()

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "shapemapper_mutation_parser --min_mapq 0 --max_internal_match 7 --min_qual 30 -i $PROCESS_DIR/parseMut/RIP_DMSO_repX.sam -o $PROCESS_DIR/parseMut/RIP_DMSO_repX.step1 --input_is_unpaired")
jobIDs+=("$(job_id $info)")

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "shapemapper_mutation_parser --min_mapq 0 --max_internal_match 7 --min_qual 30 -i $PROCESS_DIR/parseMut/RIP_NAI_repX.sam -o $PROCESS_DIR/parseMut/RIP_NAI_repX.step1 --input_is_unpaired")
jobIDs+=("$(job_id $info)")

for jobid in ${jobIDs[@]}; 
do
    echo -e "shapemapper_mutation_parser bwait $jobid"
    bwait $jobid
done

jobIDs=()

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "samtools view $PROCESS_DIR/RIP_DMSO_repX.bam | paste -d \"\t\" - $PROCESS_DIR/parseMut/RIP_DMSO_repX.step1 | cut -f 3,20- > $PROCESS_DIR/parseMut/RIP_DMSO_repX.step2")
jobIDs+=("$(job_id $info)")

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "samtools view $PROCESS_DIR/RIP_NAI_repX.bam | paste -d \"\t\" - $PROCESS_DIR/parseMut/RIP_NAI_repX.step1 | cut -f 3,20- > $PROCESS_DIR/parseMut/RIP_NAI_repX.step2")
jobIDs+=("$(job_id $info)")

for jobid in ${jobIDs[@]}; 
do
    echo -e "paste bwait $jobid"
    bwait $jobid
done

rm $PROCESS_DIR/parseMut/RIP_DMSO_repX.sam
rm $PROCESS_DIR/parseMut/RIP_NAI_repX.sam
rm $PROCESS_DIR/parseMut/RIP_DMSO_repX.step1
rm $PROCESS_DIR/parseMut/RIP_NAI_repX.step1

### 8. 拆分mutation parsed文件

jobIDs=()

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "splitMutateParser.py --max_mut 2 $PROCESS_DIR/parseMut/RIP_DMSO_repX.step2 $PROCESS_DIR/splitMut RIP_DMSO_repX")
jobIDs+=("$(job_id $info)")

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "splitMutateParser.py --max_mut 2 $PROCESS_DIR/parseMut/RIP_NAI_repX.step2 $PROCESS_DIR/splitMut RIP_NAI_repX")
jobIDs+=("$(job_id $info)")

for jobid in ${jobIDs[@]}; 
do
    echo -e "splitMutateParser bwait $jobid"
    bwait $jobid
done

### 9. run final shapemapper

info=$(bsub -q Z-ZQF -n 20 -e ${errorFn} -o ${logFn} "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Sample/shapemapper_final.py $PROCESS_DIR/splitMut $PROCESS_DIR/final_shape /Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa")
echo -e "Wait final"
bwait $(job_id $info)


### 10. Clean

finalshapes=$(ls --color=never $PROCESS_DIR/final_shape/*/RIP_NAI_repX.shape)
for file in ${finalshapes[@]};
do
    path=${file%/*}
    tid=${path##*/}
    cp ${file} $PROCESS_DIR/shape_files/${tid}.shape
done

#rm $PROCESS_DIR/*fastq.gz
#rm $PROCESS_DIR/*.fastq
#rm $PROCESS_DIR/*.fa
#rm $PROCESS_DIR/*.bam
#rm -r $PROCESS_DIR/mapping
#rm -r $PROCESS_DIR/parseMut
#rm -r $PROCESS_DIR/splitMut
#rm -r $PROCESS_DIR/final_shape


