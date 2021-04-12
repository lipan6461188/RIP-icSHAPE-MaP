
cd /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref
cp /150T/zhangqf/GenomeAnnotation/Gencode/yeast_transcriptome.fa .

cd-hit-est -i yeast_transcriptome.fa -o yeast.uniq.fa -c 0.95 -n 10 -d 0 -M 16000 -T 8


REF_rRNA=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref/rRNA/rRNA
REF=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref

######################
#### 1. Clean rRNA
######################

ROOT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/raw_data/DMS-MapSeq
IN_REP1=${ROOT}/GSM2241644-Scer_TGIRT_DMSvivo_Rep1_genomeWide/SRR3929621.fastq.gz
IN_REP2=${ROOT}/GSM2241645-Scer_TGIRT_DMSvivo_Rep2_genomeWide/SRR3929622.fastq.gz
IN_REP2_add=${ROOT}/GSM2241646-Scer_TGIRT_DMSvivo_Rep2_genomeWide_additionalSequencing/SRR3929623.fastq.gz

cd /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/1.Map_rRNA

bsub -q Z-ZQF -n 20 -e error -o log \
    "icSHAPE-pipe cleanFq -i $IN_REP1 -o REP1/rep1.fq -x $REF_rRNA -p 20 --mode Local --sam REP1/rep1.sam"

bsub -q Z-ZQF -n 20 -e error -o log \
    "icSHAPE-pipe cleanFq -i $IN_REP2 -o REP2/rep2.fq -x $REF_rRNA -p 20 --mode Local --sam REP2/rep2.sam"

bsub -q Z-ZQF -n 20 -e error -o log \
    "icSHAPE-pipe cleanFq -i $IN_REP2_add -o REP2_add/rep2_add.fq -x $REF_rRNA -p 20 --mode Local --sam REP2_add/rep2_add.sam"


######################
#### 2. Map to transcriptome
######################

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
            --outSAMtype BAM SortedByCoordinate \
            --outSAMmultNmax 1 \
            --scoreGap -1000 \
            --alignEndsType Local \
            --outSAMattributes All \
            --outFilterMultimapNmax 5 \
            --outMultimapperOrder Random \
            --outFilterMismatchNmax 3 \
            --limitBAMsortRAM 5049772506 \
            --scoreDelBase -1 \
            --scoreInsBase -1 \
            --outFilterMismatchNoverLmax 0.3"
}

IN=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/1.Map_rRNA
OUT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/2.Map_Transcript

mapping_NAI $REF $IN/REP1/rep1.fq $OUT/REP1/
mapping_NAI $REF $IN/REP2/rep2.fq $OUT/REP2/


##########################
### 3. Marge Bam
##########################

ln -s ../2.Map_Transcript/REP1/Aligned.sortedByCoord.out.bam REP1.bam
ln -s ../2.Map_Transcript/REP2/Aligned.sortedByCoord.out.bam REP2.bam  

##########################
###  4. Parse mutations
##########################

IN=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/3.merge_bam
OUT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/4.parse_mutation

declare samSets=()
declare samSets=$(ls --color=never $IN/*bam)

function job_id() { echo $2 | gawk 'match($0, /<([0-9]+)>/, id){print id[1]}'; }

for bamFullFile in ${samSets[@]};
do
    bamFile=${bamFullFile##*/}
    name=${bamFile%.*}

    mkfifo $OUT/${name}.sam
    BJOB1=$(bsub -q Z-ZQF -e error "samtools view -h $IN/${name}.bam > $OUT/${name}.sam | shapemapper_mutation_parser --min_mapq 0 --max_internal_match 7 --min_qual 30 -i $OUT/${name}.sam -o $OUT/${name}.step1 --input_is_unpaired")
    BJOB2=$(bsub -q Z-ZQF -e error -w "done($(job_id $BJOB1))" "samtools view $IN/${name}.bam > $OUT/${name}.sam | paste -d \"\t\" $OUT/${name}.sam $OUT/${name}.step1 | cut -f 3,20- > $OUT/${name}.step2")
    BJOB3=$(bsub -q Z-ZQF -e error -w "done($(job_id $BJOB2))" "rm $OUT/${name}.sam;rm $OUT/${name}.step1")
done

##########################
###  5. Split
##########################

OUT=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/5.split_mutation
IN=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/4.parse_mutation

bsub -q Z-ZQF -e error "splitMutateParser.py --max_mut 2 ${IN}/REP1.step2 ${OUT} REP1"
bsub -q Z-ZQF -e error "splitMutateParser.py --max_mut 2 ${IN}/REP2.step2 ${OUT} REP2"














