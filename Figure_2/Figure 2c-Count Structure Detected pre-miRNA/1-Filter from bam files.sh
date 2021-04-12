

############################
######  1.获得所有pre-miRNA的coverage
############################

INPATH=/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/2.combine_BAM
OUTPARH=/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/how_many_premiRNA

samtools view -h ${INPATH}/DMSO_SMR_SSII_repX.bam | grep miRNA | samtools view -bh > ${OUTPARH}/DMSO_SMR_SSII_repX.miRNA.bam
samtools view -h ${INPATH}/RIP_DMSO_repX_20190514.bam | grep miRNA | samtools view -bh > ${OUTPARH}/RIP_DMSO_repX_20190514.miRNA.bam

samtools view -F 20 -bh DMSO_SMR_SSII_repX.miRNA.bam > DMSO_SMR_SSII_repX.pos.bam 
samtools view -F 20 -bh RIP_DMSO_repX_20190514.miRNA.bam > RIP_DMSO_repX_20190514.pos.bam 

samtools depth -a DMSO_SMR_SSII_repX.pos.bam > DMSO_SMR_SSII_repX.bedGraph
samtools depth -a RIP_DMSO_repX_20190514.pos.bam > RIP_DMSO_repX_20190514.bedGraph


samtools view DMSO_SMR_SSII_repX.pos.bam | cut -f 3 | uniq -c > DMSO_SMR_SSII_repX.count
samtools view RIP_DMSO_repX_20190514.pos.bam | cut -f 3 | uniq -c > RIP_DMSO_repX_20190514.count








