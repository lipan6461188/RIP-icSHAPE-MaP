

root=/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/2.countreads
#Dicer_RIP_icSHAPE-MaP_DMSO_rep1.R1.fastq.gz 
awk '{print $2"\t"$1}' $root/RIP_DMSO_rep1_20190514 > Dicer_RIP_icSHAPE-MaP_DMSO_rep1_readcount.txt
#Dicer_RIP_icSHAPE-MaP_DMSO_rep2.R1.fastq.gz
awk '{print $2"\t"$1}' $root/RIP_DMSO_rep2_20190514 > Dicer_RIP_icSHAPE-MaP_DMSO_rep2_readcount.txt
#Dicer_RIP_input_icSHAPE-MaP_rep1.R1.fastq.gz
awk '{print $2"\t"$1}' $root/INPUT_rep1_20190514 > Dicer_RIP_input_icSHAPE-MaP_rep1_readcount.txt
#Dicer_RIP_input_icSHAPE-MaP_rep2.R1.fastq.gz
awk '{print $2"\t"$1}' $root/INPUT_rep2_20190514 > Dicer_RIP_input_icSHAPE-MaP_rep2_readcount.txt
#Dicer_RIP_input_icSHAPE-MaP_rep3.R1.fastq.gz
awk '{print $2"\t"$1}' $root/INPUT_rep3_20190514 > Dicer_RIP_input_icSHAPE-MaP_rep3_readcount.txt
#sRNA_seq_OE_WT_Dicer_rep1.fastq.gz
awk '{print $2"\t"$1}' $root/OEWTDicer_rep1 > sRNA_seq_OE_WT_Dicer_rep1_readcount.txt
#sRNA_seq_OE_WT_Dicer_rep2.fastq.gz
awk '{print $2"\t"$1}' $root/OEWTDicer_rep2 > sRNA_seq_OE_WT_Dicer_rep2_readcount.txt
#sRNA_seq_OE_WT_Dicer_rep3.fastq.gz
awk '{print $2"\t"$1}' $root/OEWTDicer_rep3 > sRNA_seq_OE_WT_Dicer_rep3_readcount.txt
#sRNA_seq_OE_WT_Dicer_rep4.fastq.gz
awk '{print $2"\t"$1}' $root/OEWTDicer_rep4 > sRNA_seq_OE_WT_Dicer_rep4_readcount.txt
#sRNA_seq_OE_Mut_Dicer_rep1.fastq.gz
awk '{print $2"\t"$1}' $root/OEMutDicer_rep1 > sRNA_seq_OE_Mut_Dicer_rep1_readcount.txt
#sRNA_seq_OE_Mut_Dicer_rep2.fastq.gz
awk '{print $2"\t"$1}' $root/OEMutDicer_rep2 > sRNA_seq_OE_Mut_Dicer_rep2_readcount.txt
#sRNA_seq_OE_Mut_Dicer_rep3.fastq.gz
awk '{print $2"\t"$1}' $root/OEMutDicer_rep3 > sRNA_seq_OE_Mut_Dicer_rep3_readcount.txt
#sRNA_seq_OE_Mut_Dicer_rep4.fastq.gz
awk '{print $2"\t"$1}' $root/OEMutDicer_rep4 > sRNA_seq_OE_Mut_Dicer_rep4_readcount.txt


