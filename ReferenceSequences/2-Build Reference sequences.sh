cp /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa .
cat losted_ref_seq.fa high_exp_RNAs.fa > small_ref_prepare.fa

inFa=/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/ref_genome.fa
faformat -in ${inFa} -out small_ref_prepare.fa -fp_chrid mRNA -append

#cd-hit-est -i small_ref_prepare.fa -o small_ref.fa -c 1.0

gaper = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed")
def format_fasta(inFa,outFa):
    OUT = open(outFa, 'w')
    for line in open(inFa):
        if line[0]=='>':
            line = line[1:-1]
            if 'mRNA' in line:
                OUT.writelines('>mRNA_'+line+"\n")
            else:
                OUT.writelines('>'+line+'\n')
        else:
            OUT.writelines(line)
    OUT.close()

format_fasta("small_ref_prepare.fa", "small_ref_format.fa")