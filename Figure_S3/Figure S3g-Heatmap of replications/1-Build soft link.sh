
ROOT=/Share/home/zhangqf7/jinsong_zhang/mirna_biogenesis/RIP_data/2019-05-14
RIP_DMSO_rep1=${ROOT}/RIP.DMSO.rep1.fastq.gz
RIP_DMSO_rep2=${ROOT}/RIP.DMSO.rep2.fastq.gz
RIP_DMSO_rep3=${ROOT}/RIP.DMSO.rep3.fastq.gz
RIP_NAI_rep1=${ROOT}/RIP.NAIN3.rep1.fastq.gz
RIP_NAI_rep2=${ROOT}/RIP.NAIN3.rep2.fastq.gz
RIP_NAI_rep3=${ROOT}/RIP.NAIN3.rep3.fastq.gz

cd /Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Replicates
mkdir -p rep1 rep2 rep3

ln -s $RIP_DMSO_rep1 rep1/RIP_DMSO_rep1.fastq.gz
ln -s $RIP_DMSO_rep2 rep2/RIP_DMSO_rep2.fastq.gz
ln -s $RIP_DMSO_rep3 rep3/RIP_DMSO_rep3.fastq.gz
ln -s $RIP_NAI_rep1 rep1/RIP_NAI_rep1.fastq.gz
ln -s $RIP_NAI_rep2 rep2/RIP_NAI_rep2.fastq.gz
ln -s $RIP_NAI_rep3 rep3/RIP_NAI_rep3.fastq.gz


./run_replicates.sh rep1 RIP_DMSO_rep1.fastq.gz RIP_NAI_rep1.fastq.gz Out_rep1
./run_replicates.sh rep2 RIP_DMSO_rep2.fastq.gz RIP_NAI_rep2.fastq.gz Out_rep2
./run_replicates.sh rep3 RIP_DMSO_rep3.fastq.gz RIP_NAI_rep3.fastq.gz Out_rep3




