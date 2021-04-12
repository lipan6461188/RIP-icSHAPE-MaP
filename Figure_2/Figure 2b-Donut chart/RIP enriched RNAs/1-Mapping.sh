
############################
##### 1. Mapping
############################

function mapping_DMSO_gz()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 20 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
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
            --outFilterMultimapNmax 20 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 170084483246 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.4 \
            --readFilesCommand \"gzip -dc\""
}

function mapping_DMSO()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 20 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
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
            --outFilterMultimapNmax 20 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 170084483246 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.4"
}

#IP.W        RIP_DMSO_repX_20190514
#IP.C        RIP_repX_20190606
#input.W     INPUT_repX_20190514
#input.C     INPUT_repX_20190606

cd /Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/1.mapping

mapping_DMSO_gz ${SMALLRNAINDEX} ${RIP_DMSO_rep1_20190514} RIP_DMSO_rep1_20190514
mapping_DMSO_gz ${SMALLRNAINDEX} ${RIP_DMSO_rep2_20190514} RIP_DMSO_rep2_20190514

mapping_DMSO_gz ${SMALLRNAINDEX} ${INPUT_rep1_20190514} INPUT_rep1_20190514
mapping_DMSO_gz ${SMALLRNAINDEX} ${INPUT_rep2_20190514} INPUT_rep2_20190514
mapping_DMSO_gz ${SMALLRNAINDEX} ${INPUT_rep3_20190514} INPUT_rep3_20190514

mapping_DMSO ${SMALLRNAINDEX} ${FSS_rep1} FSS_rep1
mapping_DMSO ${SMALLRNAINDEX} ${FSS_rep2} FSS_rep2
mapping_DMSO ${SMALLRNAINDEX} ${FSS_rep3} FSS_rep3
mapping_DMSO ${SMALLRNAINDEX} ${FSS_rep4} FSS_rep4

mapping_DMSO ${SMALLRNAINDEX} ${MnCl2_rep1} MnCl2_rep1
mapping_DMSO ${SMALLRNAINDEX} ${MnCl2_rep2} MnCl2_rep2
mapping_DMSO ${SMALLRNAINDEX} ${MnCl2_rep3} MnCl2_rep3
mapping_DMSO ${SMALLRNAINDEX} ${MnCl2_rep4} MnCl2_rep4

mapping_DMSO ${SMALLRNAINDEX} ${OEMutDicer_rep1} OEMutDicer_rep1
mapping_DMSO ${SMALLRNAINDEX} ${OEMutDicer_rep2} OEMutDicer_rep2
mapping_DMSO ${SMALLRNAINDEX} ${OEMutDicer_rep3} OEMutDicer_rep3
mapping_DMSO ${SMALLRNAINDEX} ${OEMutDicer_rep4} OEMutDicer_rep4

mapping_DMSO ${SMALLRNAINDEX} ${OEWTDicer_rep1} OEWTDicer_rep1
mapping_DMSO ${SMALLRNAINDEX} ${OEWTDicer_rep2} OEWTDicer_rep2
mapping_DMSO ${SMALLRNAINDEX} ${OEWTDicer_rep3} OEWTDicer_rep3
mapping_DMSO ${SMALLRNAINDEX} ${OEWTDicer_rep4} OEWTDicer_rep4

############################
##### 2. Count reads
############################

bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/RIP_DMSO_rep1_20190514Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/RIP_DMSO_rep1_20190514"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/RIP_DMSO_rep2_20190514Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/RIP_DMSO_rep2_20190514"

bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/INPUT_rep1_20190514Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/INPUT_rep1_20190514"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/INPUT_rep2_20190514Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/INPUT_rep2_20190514"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/INPUT_rep3_20190514Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/INPUT_rep3_20190514"


bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/FSS_rep1Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/FSS_rep1"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/FSS_rep2Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/FSS_rep2"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/FSS_rep3Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/FSS_rep3"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/FSS_rep4Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/FSS_rep4"

bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/MnCl2_rep1Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/MnCl2_rep1"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/MnCl2_rep2Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/MnCl2_rep2"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/MnCl2_rep3Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/MnCl2_rep3"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/MnCl2_rep4Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/MnCl2_rep4"

bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEMutDicer_rep1Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEMutDicer_rep1"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEMutDicer_rep2Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEMutDicer_rep2"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEMutDicer_rep3Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEMutDicer_rep3"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEMutDicer_rep4Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEMutDicer_rep4"

bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEWTDicer_rep1Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEWTDicer_rep1"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEWTDicer_rep2Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEWTDicer_rep2"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEWTDicer_rep3Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEWTDicer_rep3"
bsub -q Z-ZQF -e error "samtools view -F 20 1.mapping/OEWTDicer_rep4Aligned.sortedByCoord.out.bam | cut -f 3 | uniq -c > 2.countreads/OEWTDicer_rep4"






