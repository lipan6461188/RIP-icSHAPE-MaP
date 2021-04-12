
importCommon()

def readReadCount(rcFn):
    readcounts =  {}
    for count,tid in [ line.strip().split() for line in open(rcFn).readlines() ]:
        readcounts[tid] = int(count)
    return readcounts

root = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/2.countreads/"

INPUT_rep1_20190514 = readReadCount(root+"INPUT_rep1_20190514")
INPUT_rep2_20190514 = readReadCount(root+"INPUT_rep2_20190514")
INPUT_rep3_20190514 = readReadCount(root+"INPUT_rep3_20190514")
RIP_DMSO_rep1_20190514 = readReadCount(root+"RIP_DMSO_rep1_20190514")
RIP_DMSO_rep2_20190514 = readReadCount(root+"RIP_DMSO_rep2_20190514")

count_table = []
gene_names = []
for tid in set(INPUT_rep1_20190514)|set(INPUT_rep2_20190514)|set(INPUT_rep3_20190514)|set(RIP_DMSO_rep1_20190514)|set(RIP_DMSO_rep2_20190514):
    r1 = INPUT_rep1_20190514.get(tid,0)
    r2 = INPUT_rep2_20190514.get(tid,0)
    r3 = INPUT_rep3_20190514.get(tid,0)
    r4 = RIP_DMSO_rep1_20190514.get(tid,0)
    r5 = RIP_DMSO_rep2_20190514.get(tid,0)
    count_table.append([r1,r2,r3,r4,r5])
    gene_names.append(tid)

col_names = ['INPUT_rep1','INPUT_rep2','INPUT_rep3','RIP_1','RIP_2']
count_table = pd.DataFrame(count_table, columns=col_names, index=gene_names)

rRNA_Filter = np.array([ not (it.startswith('rRNA') or it.startswith('chrM')) for it in count_table.index.tolist()])
count_table = count_table.loc[rRNA_Filter,:]

outFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/genecount.txt"
count_table.to_csv(outFn, sep="\t")


########################
#### For FSS_repX, MnCl2_repX, OEMutDicer_repX, OEWTDicer_repX
########################

FSS_rep1,FSS_rep2,FSS_rep3,FSS_rep4 = readReadCount(root+"FSS_rep1"),readReadCount(root+"FSS_rep2"),readReadCount(root+"FSS_rep3"),readReadCount(root+"FSS_rep4")
MnCl2_rep1,MnCl2_rep2,MnCl2_rep3,MnCl2_rep4 = readReadCount(root+"MnCl2_rep1"),readReadCount(root+"MnCl2_rep2"),readReadCount(root+"MnCl2_rep3"),readReadCount(root+"MnCl2_rep4")
OEMutDicer_rep1,OEMutDicer_rep2,OEMutDicer_rep3,OEMutDicer_rep4 = readReadCount(root+"OEMutDicer_rep1"),readReadCount(root+"OEMutDicer_rep2"),readReadCount(root+"OEMutDicer_rep3"),readReadCount(root+"OEMutDicer_rep4")
OEWTDicer_rep1,OEWTDicer_rep2,OEWTDicer_rep3,OEWTDicer_rep4 = readReadCount(root+"OEWTDicer_rep1"),readReadCount(root+"OEWTDicer_rep2"),readReadCount(root+"OEWTDicer_rep3"),readReadCount(root+"OEWTDicer_rep4")

INPUT_rep1_20190514 = readReadCount(root+"INPUT_rep1_20190514")
INPUT_rep2_20190514 = readReadCount(root+"INPUT_rep2_20190514")
INPUT_rep3_20190514 = readReadCount(root+"INPUT_rep3_20190514")
RIP_DMSO_rep1_20190514 = readReadCount(root+"RIP_DMSO_rep1_20190514")
RIP_DMSO_rep2_20190514 = readReadCount(root+"RIP_DMSO_rep2_20190514")

FSS_repX = set(FSS_rep1)|set(FSS_rep2)|set(FSS_rep3)|set(FSS_rep4)
MnCl2_repX = set(MnCl2_rep1)|set(MnCl2_rep2)|set(MnCl2_rep3)|set(MnCl2_rep4)
OEMutDicer_repX = set(OEMutDicer_rep1)|set(OEMutDicer_rep2)|set(OEMutDicer_rep3)|set(OEMutDicer_rep4)
OEWTDicer_repX = set(OEWTDicer_rep1)|set(OEWTDicer_rep2)|set(OEWTDicer_rep3)|set(OEWTDicer_rep4)
INPUT_repX_20190514 = set(INPUT_rep1_20190514)|set(INPUT_rep2_20190514)|set(INPUT_rep3_20190514)
RIP_DMSO_repX_20190514 = set(RIP_DMSO_rep1_20190514)|set(RIP_DMSO_rep2_20190514)

count_table = []
gene_names = []
for tid in FSS_repX|MnCl2_repX|OEMutDicer_repX|OEWTDicer_repX|INPUT_repX_20190514|RIP_DMSO_repX_20190514:
    tmp = []
    for Dict in [FSS_rep1,FSS_rep2,FSS_rep3,FSS_rep4]+[MnCl2_rep1,MnCl2_rep2,MnCl2_rep3,MnCl2_rep4]+[OEMutDicer_rep1,OEMutDicer_rep2,OEMutDicer_rep3,OEMutDicer_rep4]+[OEWTDicer_rep1,OEWTDicer_rep2,OEWTDicer_rep3,OEWTDicer_rep4]+[INPUT_rep1_20190514,INPUT_rep2_20190514,INPUT_rep3_20190514,RIP_DMSO_rep1_20190514,RIP_DMSO_rep2_20190514]:
        tmp.append( Dict.get(tid,0) )
    count_table.append(tmp)
    gene_names.append(tid)

col_names = ['FSS_rep1','FSS_rep2','FSS_rep3','FSS_rep4'] + \
            ['MnCl2_rep1','MnCl2_rep2','MnCl2_rep3','MnCl2_rep4'] + \
            ['OEMutDicer_rep1','OEMutDicer_rep2','OEMutDicer_rep3','OEMutDicer_rep4'] + \
            ['OEWTDicer_rep1','OEWTDicer_rep2','OEWTDicer_rep3','OEWTDicer_rep4'] + \
            ['INPUT_rep1_20190514','INPUT_rep2_20190514','INPUT_rep3_20190514'] + \
            ['RIP_DMSO_rep1_20190514','RIP_DMSO_rep2_20190514']

count_table = pd.DataFrame(count_table, columns=col_names, index=gene_names)

rRNA_Filter = np.array([ not (it.startswith('rRNA') or it.startswith('chrM')) for it in count_table.index.tolist()])
count_table = count_table.loc[rRNA_Filter,:]

outFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/RNA-Seq-genecount.txt"
count_table.to_csv(outFn, sep="\t")



