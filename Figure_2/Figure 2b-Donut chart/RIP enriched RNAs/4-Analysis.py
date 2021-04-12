
importCommon()

def read_shape(inFolder):
    SHAPE = {}; Sequence = {}
    inFolder = inFolder.rstrip('/')+'/'
    files = os.listdir(inFolder)
    for file in files:
        if file.endswith(".shape"):
            tid = file.rstrip('.shape')
            data = General.load_SHAPEMap(inFolder+file)
            shape = data['shape_pro_list']
            seq = data['seq']
            if shape.count('NULL')/len(shape)>0.6: continue
            SHAPE[tid] = shape
            Sequence[tid] = seq
    return SHAPE, Sequence

DE_table = pd.read_csv("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/4.limmaDE/limma.txt",sep="\t")
DE_table.index = DE_table['Gene']
DE_table = DE_table.drop(columns=['Gene'])

GC_table = pd.read_csv("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/genecount.txt",sep="\t",index_col=0)

DE_sig = DE_table.loc[(DE_table['adj.P.Val']<0.05)&(DE_table['logFC']>1)]

######################
#### 查看那些显著表达基因他们的readcount是怎样的
######################

miRNAs = [ it for it in DE_sig.index if it.startswith('miRNA') ]
mRNAIntron = [ it for it in DE_sig.index if it.startswith('mRNAIntron') ]
mRNAExon = [ it for it in DE_sig.index if it.startswith('mRNAExon') ]
snoRNA = [ it for it in DE_sig.index if it.startswith('snoRNA') ]

GC_table.loc[miRNAs,:]
GC_table.loc[mRNAIntron,:]
GC_table.loc[mRNAExon,:]
GC_table.loc[snoRNA,:]

######################
#### 显著的RNA和结构之间的overlap
######################

SHAPE, Sequence = read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/Out/shape_files")

set(SHAPE)&set(DE_sig.index)



miRNAs = os.popen("awk '$2>0&&$6<0.05{print $0}' /Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/4.limmaDE/limma.txt | grep miRNA | cut -f 1").readlines()
miRNAs = [it.strip() for it in miRNAs]

miRNAs_dct = os.listdir("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/Out/shape_files")
miRNAs_dct = [ it.rstrip('.shape') for it in miRNAs_dct if it.startswith('miRNA') ]

len(set(miRNAs)&set(miRNAs_dct))



