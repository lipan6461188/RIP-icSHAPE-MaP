
DMSO_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/dmso_rep1.clean.fastq.gz
NAI_100mm_10min_rep1_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_100mm_10min_rep1.clean.fastq.gz
NAI_100mm_10min_rep2_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_100mm_10min_rep2.clean.fastq.gz
NAI_100mm_5min_rep1_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_100mm_5min_rep1.clean.fastq.gz
NAI_100mm_5min_rep2_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_100mm_5min_rep2.clean.fastq.gz
NAI_50mm_10min_rep1_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_50mm_10min_rep1.clean.fastq.gz
NAI_50mm_10min_rep2_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_50mm_10min_rep2.clean.fastq.gz
NAI_50mm_5min_rep1_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_50mm_5min_rep1.clean.fastq.gz
NAI_50mm_5min_rep2_20190412=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190404/preprocessing/preprocessing/in_vitro_50mm_5min_rep2.clean.fastq.gz

ROOT=/Share/home/zhangqf7/jinsong_zhang/shapemap/20190422/preprocessing/preprocessing/
DMSO_20190422=$ROOT/dmso_rep1.clean.fastq.gz
NAI_100mm_10min_rep1_20190422=$ROOT/in_vitro_100mm_10min_rep1.clean.fastq.gz
NAI_100mm_10min_rep2_20190422=$ROOT/in_vitro_100mm_10min_rep2.clean.fastq.gz
NAI_100mm_5min_rep1_20190422=$ROOT/in_vitro_100mm_5min_rep1.clean.fastq.gz
NAI_100mm_5min_rep2_20190422=$ROOT/in_vitro_100mm_5min_rep2.clean.fastq.gz
NAI_50mm_10min_rep1_20190422=$ROOT/in_vitro_50mm_10min_rep1.clean.fastq.gz
NAI_50mm_10min_rep2_20190422=$ROOT/in_vitro_50mm_10min_rep2.clean.fastq.gz
NAI_50mm_5min_rep1_20190422=$ROOT/in_vitro_50mm_5min_rep1.clean.fastq.gz
NAI_50mm_5min_rep2_20190422=$ROOT/in_vitro_50mm_5min_rep2.clean.fastq.gz


################
## 创建目录
################

cd /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/1.STAR_mapping

mkdir -p DMSO_dicer_1_CIRL_BENR_SSII_rep1
mkdir -p DMSO_dicer_1_CIRL_BENR_SSII_rep2
mkdir -p DMSO_dicer_2_CIRL_BENR_SSII_rep1
mkdir -p DMSO_dicer_2_CIRL_BENR_SSII_rep2
mkdir -p NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1
mkdir -p NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2
mkdir -p NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1
mkdir -p NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2

mkdir DMSO_SMR_SSII_rep1
mkdir DMSO_SMR_SSII_rep2
mkdir NAI_100mm_vivo_SMR_SSII_rep1
mkdir NAI_100mm_vivo_SMR_SSII_rep2

mkdir DMSO_SMR_BENR_SSII
mkdir NAI_100mm_vitro_SMR_BENR_SSII
mkdir NAI_100mm_vitro_SMR_CENR_SSII
mkdir NAI_100mm_vitro_SMR_SSII_rep1
mkdir NAI_100mm_vitro_SMR_SSII_rep2
mkdir NAI_200mm_vitro_SMR_BENR_SSII

mkdir NAI_100mm_vivo_CIRL_CENR_SSII_rep1
mkdir NAI_100mm_vivo_CIRL_CENR_SSII_rep2
mkdir NAI_100mm_vivo_CIRL_SSII
mkdir NAI_100mm_vivo_CIRL_TGIII

mkdir DMSO_20190412
mkdir NAI_100mm_10min_rep1_20190412
mkdir NAI_100mm_10min_rep2_20190412
mkdir NAI_100mm_5min_rep1_20190412
mkdir NAI_100mm_5min_rep2_20190412
mkdir NAI_50mm_10min_rep1_20190412
mkdir NAI_50mm_10min_rep2_20190412
mkdir NAI_50mm_5min_rep1_20190412
mkdir NAI_50mm_5min_rep2_20190412

mkdir DMSO_20190422
mkdir NAI_100mm_10min_rep1_20190422
mkdir NAI_100mm_10min_rep2_20190422
mkdir NAI_100mm_5min_rep1_20190422
mkdir NAI_100mm_5min_rep2_20190422
mkdir NAI_50mm_10min_rep1_20190422
mkdir NAI_50mm_10min_rep2_20190422
mkdir NAI_50mm_5min_rep1_20190422
mkdir NAI_50mm_5min_rep2_20190422

cd /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome

mkdir -p DMSO_dicer_1_CIRL_BENR_SSII_rep1
mkdir -p DMSO_dicer_1_CIRL_BENR_SSII_rep2
mkdir -p DMSO_dicer_2_CIRL_BENR_SSII_rep1
mkdir -p DMSO_dicer_2_CIRL_BENR_SSII_rep2
mkdir -p NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1
mkdir -p NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2
mkdir -p NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1
mkdir -p NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2

mkdir DMSO_SMR_SSII_rep1
mkdir DMSO_SMR_SSII_rep2
mkdir NAI_100mm_vivo_SMR_SSII_rep1
mkdir NAI_100mm_vivo_SMR_SSII_rep2

mkdir DMSO_SMR_BENR_SSII
mkdir NAI_100mm_vitro_SMR_BENR_SSII
mkdir NAI_100mm_vitro_SMR_CENR_SSII
mkdir NAI_100mm_vitro_SMR_SSII_rep1
mkdir NAI_100mm_vitro_SMR_SSII_rep2
mkdir NAI_200mm_vitro_SMR_BENR_SSII

mkdir NAI_100mm_vivo_CIRL_CENR_SSII_rep1
mkdir NAI_100mm_vivo_CIRL_CENR_SSII_rep2
mkdir NAI_100mm_vivo_CIRL_SSII
mkdir NAI_100mm_vivo_CIRL_TGIII

mkdir DMSO_20190412
mkdir NAI_100mm_10min_rep1_20190412
mkdir NAI_100mm_10min_rep2_20190412
mkdir NAI_100mm_5min_rep1_20190412
mkdir NAI_100mm_5min_rep2_20190412
mkdir NAI_50mm_10min_rep1_20190412
mkdir NAI_50mm_10min_rep2_20190412
mkdir NAI_50mm_5min_rep1_20190412
mkdir NAI_50mm_5min_rep2_20190412

mkdir DMSO_20190422
mkdir NAI_100mm_10min_rep1_20190422
mkdir NAI_100mm_10min_rep2_20190422
mkdir NAI_100mm_5min_rep1_20190422
mkdir NAI_100mm_5min_rep2_20190422
mkdir NAI_50mm_10min_rep1_20190422
mkdir NAI_50mm_10min_rep2_20190422
mkdir NAI_50mm_5min_rep1_20190422
mkdir NAI_50mm_5min_rep2_20190422


## Data 20190514
mkdir -p INPUT_rep1_20190514
mkdir -p INPUT_rep2_20190514
mkdir -p INPUT_rep3_20190514
mkdir -p RIP_DMSO_rep1_20190514
mkdir -p RIP_DMSO_rep2_20190514
mkdir -p RIP_DMSO_rep3_20190514
mkdir -p RIP_NAIN3_rep1_20190514
mkdir -p RIP_NAIN3_rep2_20190514
mkdir -p RIP_NAIN3_rep3_20190514


## Data 20190606
mkdir -p INPUT_rep1_20190606
mkdir -p INPUT_rep2_20190606
mkdir -p RIP_rep1_20190606
mkdir -p RIP_rep2_20190606

################
## 低要求Mapping
################

function mapping_NAI()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 10 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --scoreGap -1000 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 50 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 5 \
            --limitBAMsortRAM 5049772506 \
            --scoreDelBase -1 \
            --scoreInsBase -1 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3"
}

function mapping_DMSO()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 10 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --scoreGap -1000 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 50 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 5049772506 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3"
}

function mapping_NAI_gz()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 5 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --scoreGap -1000 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 50 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 5 \
            --limitBAMsortRAM 5049772506 \
            --scoreDelBase -1 \
            --scoreInsBase -1 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3 \
            --readFilesCommand \"gzip -dc\""
}

function mapping_DMSO_gz()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-ZQF -n 5 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --scoreGap -1000 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 50 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 5049772506 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3 \
            --readFilesCommand \"gzip -dc\""
}


lpd

REF=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/ref/
OUT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/1.STAR_mapping

mapping_DMSO $REF $DMSO_dicer_1_CIRL_BENR_SSII_rep1 $OUT/DMSO_dicer_1_CIRL_BENR_SSII_rep1/
mapping_DMSO $REF $DMSO_dicer_1_CIRL_BENR_SSII_rep2 $OUT/DMSO_dicer_1_CIRL_BENR_SSII_rep2/

mapping_DMSO $REF $DMSO_dicer_2_CIRL_BENR_SSII_rep1 $OUT/DMSO_dicer_2_CIRL_BENR_SSII_rep1/
mapping_DMSO $REF $DMSO_dicer_2_CIRL_BENR_SSII_rep2 $OUT/DMSO_dicer_2_CIRL_BENR_SSII_rep2/

mapping_NAI $REF $NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1 $OUT/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2 $OUT/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2/

mapping_NAI $REF $NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1 $OUT/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2 $OUT/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2/


mapping_DMSO $REF $DMSO_SMR_SSII_rep1 $OUT/DMSO_SMR_SSII_rep1/
mapping_DMSO $REF $DMSO_SMR_SSII_rep2 $OUT/DMSO_SMR_SSII_rep2/

mapping_DMSO $REF $NAI_100mm_vivo_SMR_SSII_rep1 $OUT/NAI_100mm_vivo_SMR_SSII_rep1/
mapping_DMSO $REF $NAI_100mm_vivo_SMR_SSII_rep2 $OUT/NAI_100mm_vivo_SMR_SSII_rep2/

mapping_DMSO $REF $DMSO_SMR_BENR_SSII $OUT/DMSO_SMR_BENR_SSII/
mapping_NAI $REF $NAI_100mm_vitro_SMR_BENR_SSII $OUT/NAI_100mm_vitro_SMR_BENR_SSII/
mapping_NAI $REF $NAI_100mm_vitro_SMR_CENR_SSII $OUT/NAI_100mm_vitro_SMR_CENR_SSII/
mapping_NAI $REF $NAI_100mm_vitro_SMR_SSII_rep1 $OUT/NAI_100mm_vitro_SMR_SSII_rep1/
mapping_NAI $REF $NAI_100mm_vitro_SMR_SSII_rep2 $OUT/NAI_100mm_vitro_SMR_SSII_rep2/
mapping_NAI $REF $NAI_200mm_vitro_SMR_BENR_SSII $OUT/NAI_200mm_vitro_SMR_BENR_SSII/

mapping_NAI $REF $NAI_100mm_vivo_CIRL_CENR_SSII_rep1 $OUT/NAI_100mm_vivo_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $NAI_100mm_vivo_CIRL_CENR_SSII_rep2 $OUT/NAI_100mm_vivo_CIRL_CENR_SSII_rep2/
mapping_NAI $REF $NAI_100mm_vivo_CIRL_SSII $OUT/NAI_100mm_vivo_CIRL_SSII/
mapping_NAI $REF $NAI_100mm_vivo_CIRL_TGIII $OUT/NAI_100mm_vivo_CIRL_TGIII/


mapping_DMSO_gz $REF $DMSO_20190412 $OUT/DMSO_20190412/
mapping_NAI_gz $REF $NAI_100mm_10min_rep1_20190412 $OUT/NAI_100mm_10min_rep1_20190412/
mapping_NAI_gz $REF $NAI_100mm_10min_rep2_20190412 $OUT/NAI_100mm_10min_rep2_20190412/
mapping_NAI_gz $REF $NAI_100mm_5min_rep1_20190412 $OUT/NAI_100mm_5min_rep1_20190412/
mapping_NAI_gz $REF $NAI_100mm_5min_rep2_20190412 $OUT/NAI_100mm_5min_rep2_20190412/
mapping_NAI_gz $REF $NAI_50mm_10min_rep1_20190412 $OUT/NAI_50mm_10min_rep1_20190412/
mapping_NAI_gz $REF $NAI_50mm_10min_rep2_20190412 $OUT/NAI_50mm_10min_rep2_20190412/
mapping_NAI_gz $REF $NAI_50mm_5min_rep1_20190412 $OUT/NAI_50mm_5min_rep1_20190412/
mapping_NAI_gz $REF $NAI_50mm_5min_rep2_20190412 $OUT/NAI_50mm_5min_rep2_20190412/


mapping_DMSO_gz $REF $DMSO_20190422 $OUT/DMSO_20190422/
mapping_NAI_gz $REF $NAI_100mm_10min_rep1_20190422 $OUT/NAI_100mm_10min_rep1_20190422/
mapping_NAI_gz $REF $NAI_100mm_10min_rep2_20190422 $OUT/NAI_100mm_10min_rep2_20190422/
mapping_NAI_gz $REF $NAI_100mm_5min_rep1_20190422 $OUT/NAI_100mm_5min_rep1_20190422/
mapping_NAI_gz $REF $NAI_100mm_5min_rep2_20190422 $OUT/NAI_100mm_5min_rep2_20190422/
mapping_NAI_gz $REF $NAI_50mm_10min_rep1_20190422 $OUT/NAI_50mm_10min_rep1_20190422/
mapping_NAI_gz $REF $NAI_50mm_10min_rep2_20190422 $OUT/NAI_50mm_10min_rep2_20190422/
mapping_NAI_gz $REF $NAI_50mm_5min_rep1_20190422 $OUT/NAI_50mm_5min_rep1_20190422/
mapping_NAI_gz $REF $NAI_50mm_5min_rep2_20190422 $OUT/NAI_50mm_5min_rep2_20190422/


## Data 20190514
mapping_DMSO_gz $REF $INPUT_rep1_20190514 $OUT/INPUT_rep1_20190514/
mapping_DMSO_gz $REF $INPUT_rep2_20190514 $OUT/INPUT_rep2_20190514/
mapping_DMSO_gz $REF $INPUT_rep3_20190514 $OUT/INPUT_rep3_20190514/
mapping_DMSO_gz $REF $RIP_DMSO_rep1_20190514 $OUT/RIP_DMSO_rep1_20190514/
mapping_DMSO_gz $REF $RIP_DMSO_rep2_20190514 $OUT/RIP_DMSO_rep2_20190514/
mapping_DMSO_gz $REF $RIP_DMSO_rep3_20190514 $OUT/RIP_DMSO_rep3_20190514/
mapping_NAI_gz $REF $RIP_NAIN3_rep1_20190514 $OUT/RIP_NAIN3_rep1_20190514/
mapping_NAI_gz $REF $RIP_NAIN3_rep2_20190514 $OUT/RIP_NAIN3_rep2_20190514/
mapping_NAI_gz $REF $RIP_NAIN3_rep3_20190514 $OUT/RIP_NAIN3_rep3_20190514/

## Data 20190606

mapping_DMSO_gz $REF $INPUT_rep1_20190606 $OUT/INPUT_rep1_20190606/
mapping_DMSO_gz $REF $INPUT_rep2_20190606 $OUT/INPUT_rep2_20190606/
mapping_DMSO_gz $REF $RIP_rep1_20190606 $OUT/RIP_rep1_20190606/
mapping_DMSO_gz $REF $RIP_rep2_20190606 $OUT/RIP_rep2_20190606/


################
## 剩余的Reads map到基因组上
################

IN=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/1.STAR_mapping
REF=/150T/zhangqf/GenomeAnnotation/INDEX/STAR/hg38_NCBI/
OUT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome

function mapping_NAI()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-HNODE -n 5 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "/150T/zhangqf/GenomeAnnotation/INDEX/bin/STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 2 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 5049772506 \
            --scoreDelBase -1 \
            --scoreInsBase -1 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outSAMstrandField intronMotif \
            --outSJfilterOverhangMin 30 12 12 12"
}
function mapping_DMSO()
{
    REF=$1
    INPUT=$2
    OUTPUT=$3
    
    bsub -q Z-HNODE -n 5 -e ${OUTPUT}_error -o ${OUTPUT}_log \
        "/150T/zhangqf/GenomeAnnotation/INDEX/bin/STAR --readFilesIn $INPUT \
            --outFileNamePrefix $OUTPUT \
            --genomeDir $REF \
            --runThreadN 20 \
            --genomeLoad NoSharedMemory \
            --runMode alignReads \
            --outSAMtype BAM Unsorted \
            --outSAMmultNmax 1 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 2 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 2 \
            --limitBAMsortRAM 5049772506 \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNoverLmax 0.3 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outSAMstrandField intronMotif \
            --outSJfilterOverhangMin 30 12 12 12"
}

mapping_DMSO $REF $IN/DMSO_dicer_1_CIRL_BENR_SSII_rep1/Unmapped.out.mate1 $OUT/DMSO_dicer_1_CIRL_BENR_SSII_rep1/
mapping_DMSO $REF $IN/DMSO_dicer_1_CIRL_BENR_SSII_rep2/Unmapped.out.mate1 $OUT/DMSO_dicer_1_CIRL_BENR_SSII_rep2/

mapping_DMSO $REF $IN/DMSO_dicer_2_CIRL_BENR_SSII_rep1/Unmapped.out.mate1 $OUT/DMSO_dicer_2_CIRL_BENR_SSII_rep1/
mapping_DMSO $REF $IN/DMSO_dicer_2_CIRL_BENR_SSII_rep2/Unmapped.out.mate1 $OUT/DMSO_dicer_2_CIRL_BENR_SSII_rep2/

mapping_NAI $REF $IN/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1/Unmapped.out.mate1 $OUT/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $IN/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2/Unmapped.out.mate1 $OUT/NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_rep2/

mapping_NAI $REF $IN/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1/Unmapped.out.mate1 $OUT/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $IN/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2/Unmapped.out.mate1 $OUT/NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_rep2/


mapping_DMSO $REF $IN/DMSO_SMR_SSII_rep1/Unmapped.out.mate1 $OUT/DMSO_SMR_SSII_rep1/
mapping_DMSO $REF $IN/DMSO_SMR_SSII_rep2/Unmapped.out.mate1 $OUT/DMSO_SMR_SSII_rep2/

mapping_NAI $REF $IN/NAI_100mm_vivo_SMR_SSII_rep1/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_SMR_SSII_rep1/
mapping_NAI $REF $IN/NAI_100mm_vivo_SMR_SSII_rep2/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_SMR_SSII_rep2/


mapping_DMSO $REF $IN/DMSO_SMR_BENR_SSII/Unmapped.out.mate1 $OUT/DMSO_SMR_BENR_SSII/
mapping_NAI $REF $IN/NAI_100mm_vitro_SMR_BENR_SSII/Unmapped.out.mate1 $OUT/NAI_100mm_vitro_SMR_BENR_SSII/
mapping_NAI $REF $IN/NAI_100mm_vitro_SMR_CENR_SSII/Unmapped.out.mate1 $OUT/NAI_100mm_vitro_SMR_CENR_SSII/
mapping_NAI $REF $IN/NAI_100mm_vitro_SMR_SSII_rep1/Unmapped.out.mate1 $OUT/NAI_100mm_vitro_SMR_SSII_rep1/
mapping_NAI $REF $IN/NAI_100mm_vitro_SMR_SSII_rep2/Unmapped.out.mate1 $OUT/NAI_100mm_vitro_SMR_SSII_rep2/
mapping_NAI $REF $IN/NAI_200mm_vitro_SMR_BENR_SSII/Unmapped.out.mate1 $OUT/NAI_200mm_vitro_SMR_BENR_SSII/


mapping_NAI $REF $IN/NAI_100mm_vivo_CIRL_CENR_SSII_rep1/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_CIRL_CENR_SSII_rep1/
mapping_NAI $REF $IN/NAI_100mm_vivo_CIRL_CENR_SSII_rep2/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_CIRL_CENR_SSII_rep2/
mapping_NAI $REF $IN/NAI_100mm_vivo_CIRL_SSII/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_CIRL_SSII/
mapping_NAI $REF $IN/NAI_100mm_vivo_CIRL_TGIII/Unmapped.out.mate1 $OUT/NAI_100mm_vivo_CIRL_TGIII/


mapping_DMSO $REF $IN/DMSO_20190412/Unmapped.out.mate1 $OUT/DMSO_20190412/
mapping_NAI $REF $IN/NAI_100mm_10min_rep1_20190412/Unmapped.out.mate1 $OUT/NAI_100mm_10min_rep1_20190412/
mapping_NAI $REF $IN/NAI_100mm_10min_rep2_20190412/Unmapped.out.mate1 $OUT/NAI_100mm_10min_rep2_20190412/
mapping_NAI $REF $IN/NAI_100mm_5min_rep1_20190412/Unmapped.out.mate1 $OUT/NAI_100mm_5min_rep1_20190412/
mapping_NAI $REF $IN/NAI_100mm_5min_rep2_20190412/Unmapped.out.mate1 $OUT/NAI_100mm_5min_rep2_20190412/
mapping_NAI $REF $IN/NAI_50mm_10min_rep1_20190412/Unmapped.out.mate1 $OUT/NAI_50mm_10min_rep1_20190412/
mapping_NAI $REF $IN/NAI_50mm_10min_rep2_20190412/Unmapped.out.mate1 $OUT/NAI_50mm_10min_rep2_20190412/
mapping_NAI $REF $IN/NAI_50mm_5min_rep1_20190412/Unmapped.out.mate1 $OUT/NAI_50mm_5min_rep1_20190412/
mapping_NAI $REF $IN/NAI_50mm_5min_rep2_20190412/Unmapped.out.mate1 $OUT/NAI_50mm_5min_rep2_20190412/


mapping_DMSO $REF $IN/DMSO_20190422/Unmapped.out.mate1 $OUT/DMSO_20190422/
mapping_NAI $REF $IN/NAI_100mm_10min_rep1_20190422/Unmapped.out.mate1 $OUT/NAI_100mm_10min_rep1_20190422/
mapping_NAI $REF $IN/NAI_100mm_10min_rep2_20190422/Unmapped.out.mate1 $OUT/NAI_100mm_10min_rep2_20190422/
mapping_NAI $REF $IN/NAI_100mm_5min_rep1_20190422/Unmapped.out.mate1 $OUT/NAI_100mm_5min_rep1_20190422/
mapping_NAI $REF $IN/NAI_100mm_5min_rep2_20190422/Unmapped.out.mate1 $OUT/NAI_100mm_5min_rep2_20190422/
mapping_NAI $REF $IN/NAI_50mm_10min_rep1_20190422/Unmapped.out.mate1 $OUT/NAI_50mm_10min_rep1_20190422/
mapping_NAI $REF $IN/NAI_50mm_10min_rep2_20190422/Unmapped.out.mate1 $OUT/NAI_50mm_10min_rep2_20190422/
mapping_NAI $REF $IN/NAI_50mm_5min_rep1_20190422/Unmapped.out.mate1 $OUT/NAI_50mm_5min_rep1_20190422/
mapping_NAI $REF $IN/NAI_50mm_5min_rep2_20190422/Unmapped.out.mate1 $OUT/NAI_50mm_5min_rep2_20190422/


## Data 20190514
mapping_DMSO $REF $IN/INPUT_rep1_20190514/Unmapped.out.mate1 $OUT/INPUT_rep1_20190514/
mapping_DMSO $REF $IN/INPUT_rep2_20190514/Unmapped.out.mate1 $OUT/INPUT_rep2_20190514/
mapping_DMSO $REF $IN/INPUT_rep3_20190514/Unmapped.out.mate1 $OUT/INPUT_rep3_20190514/
mapping_DMSO $REF $IN/RIP_DMSO_rep1_20190514/Unmapped.out.mate1 $OUT/RIP_DMSO_rep1_20190514/
mapping_DMSO $REF $IN/RIP_DMSO_rep2_20190514/Unmapped.out.mate1 $OUT/RIP_DMSO_rep2_20190514/
mapping_DMSO $REF $IN/RIP_DMSO_rep3_20190514/Unmapped.out.mate1 $OUT/RIP_DMSO_rep3_20190514/
mapping_NAI $REF $IN/RIP_NAIN3_rep1_20190514/Unmapped.out.mate1 $OUT/RIP_NAIN3_rep1_20190514/
mapping_NAI $REF $IN/RIP_NAIN3_rep2_20190514/Unmapped.out.mate1 $OUT/RIP_NAIN3_rep2_20190514/
mapping_NAI $REF $IN/RIP_NAIN3_rep3_20190514/Unmapped.out.mate1 $OUT/RIP_NAIN3_rep3_20190514/

## Data 20190606
mapping_DMSO $REF $IN/INPUT_rep1_20190606/Unmapped.out.mate1 $OUT/INPUT_rep1_20190606/
mapping_DMSO $REF $IN/INPUT_rep2_20190606/Unmapped.out.mate1 $OUT/INPUT_rep2_20190606/
mapping_DMSO $REF $IN/RIP_rep1_20190606/Unmapped.out.mate1 $OUT/RIP_rep1_20190606/
mapping_DMSO $REF $IN/RIP_rep2_20190606/Unmapped.out.mate1 $OUT/RIP_rep2_20190606/


################
## 统计Reads Map到哪些基因上
################

######## Python Script ########
######## Python Script ########
######## Python Script ########


