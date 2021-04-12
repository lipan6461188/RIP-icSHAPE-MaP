
#########################
#### 合并replicates比较哪些RNAs是Enrich的
#########################

inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/genecount.txt"
genecount = pd.read_csv(inFn, sep="\t", index_col=0)

genecount = genecount.loc[genecount.max(axis=1)>100,:]

rRNA_Filter = np.array([ not (it.startswith('rRNA') or it.startswith('chrM')) for it in genecount.index.tolist()])
genecount = genecount.loc[rRNA_Filter,:]

INPUT = genecount['INPUT_rep1']+genecount['INPUT_rep2']+genecount['INPUT_rep3']
RIP = genecount['RIP_1']+genecount['RIP_2']

INPUT_NORMFAC = INPUT.sum()
RIP_NORMFAC = RIP.sum()

NORM_INPUT = INPUT / INPUT_NORMFAC
NORM_RIP = RIP / RIP_NORMFAC

log2FC = np.log2(NORM_RIP/(NORM_INPUT+1e-15))

df = pd.concat([INPUT, RIP, NORM_INPUT, NORM_RIP, log2FC], axis=1)
df.columns = ['INPUT', 'RIP', 'NORM_INPUT','NORM_RIP', 'log2FC']
df = df.sort_values(by=['log2FC'],ascending=False)

df.to_csv("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/5.foldChange/foldChange.txt", sep="\t")


"""
总结所有的基因类型：
awk 'NR==1{print $0}$6>2{print $0}' foldChange.txt | awk '{split($1,arr,"_"); print arr[1]}' | sort | uniq -c
"""
