
"""

主要数：
1. 有结构分数的各类RNA的数量：Figure S2F
2. Dicer富集的各类RNA的数量：Figure S2C
3. 有结构分数且Dicer富集的各类RNA的数量：Figure 2B
4. Dicer富集的各类RNA和Cell2014的overlap：Supplementary Table 1

"""


importCommon()
import DicerMiRNA

def readEnrichedRNAs():
    inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/4.limmaDE/limma.txt"
    df = pd.read_csv(inFn, sep="\t")
    enriched_RNAs_1 = df.loc[(df['logFC']>1)&(df['adj.P.Val']<0.05),'Gene'].tolist()
    
    inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/5.foldChange/foldChange.txt"
    df = pd.read_csv(inFn, sep="\t", index_col=0)
    #### 2019/8/28这里使用0作为enrich的条件
    enriched_RNAs_2 = df.loc[(df['log2FC']>0),:].index.tolist()
    
    return enriched_RNAs_1, enriched_RNAs_2

def readCell2014DicerCLIPTable():
    inFn1 = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/Supplementary/table1_hg38.bed"
    inFn2 = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/Supplementary/table2_hg38.bed"
    Cell_2014_1 = [ line.strip().split() for line in open(inFn1) ]
    Cell_2014_2 = [ line.strip().split() for line in open(inFn2) ]
    for data in Cell_2014_1:
        data[1],data[2],data[3] = int(data[1]),int(data[2]),data[5]
        del data[4]
        del data[4]
    for data in Cell_2014_2:
        data[1],data[2],data[3] = int(data[1]),int(data[2]),data[5]
        del data[4]
        del data[4]
    return Cell_2014_1,Cell_2014_2

def compare_Seqs_And_Locus(Seqs, Locus, tolerance=100):
    """
    Seqs        -- {tid1:seq1, tid2:seq2, tid3:seq3...}
    Locus       -- [ ['chr1', pos1, pos2, strand],... ]
    """
    ## Overlap EnrichedRNAs and Cell2014DicerCLIPTable
    ### 1. Find all sequences for each locus
    i = 0
    locus_with_Seq = []
    for chrID,start,end,strand in Locus:
        start -= 500
        end += 500
        try:
            LongSeq = Seqer.fetch(chrID,start,end,strand)
        except:
            continue
        locus_with_Seq.append( [chrID,start,end,strand,LongSeq] )
    
    overlap = []
    for tid in Seqs:
        seq = smallRNASeqs[tid]
        for chrID,start,end,strand,LongSeq in locus_with_Seq:
            findpos = LongSeq.find(seq)
            if findpos == -1:
                continue
            else:
                if findpos<len(LongSeq)-500+tolerance and findpos+len(seq)>500-tolerance:
                    has_find =  True
                    overlap.append([tid, (chrID,start,end,strand)])
                    break
    return overlap

def filter_wellcovtrans(shapeDict, min_cov_num=30, min_cov_ratio=0.5):
    valid_tid_list = []
    for tid in shapeDict:
        shapemap =shapeDict[tid]
        validnum = len(shapemap)-shapemap.count("NULL")
        if validnum >= min_cov_num and 1.0*validnum/len(shapemap) >= min_cov_ratio:
            valid_tid_list.append(tid)
    return valid_tid_list


smallRNASeqs = General.load_fasta("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/0.ref_seqs/index/smallRNAWithGenome.fa")
Seqer = Seq.seqClass("/150T/zhangqf/GenomeAnnotation/genome/hg38.fa")
Gaper = GAP.init("/150T/zhangqf/GenomeAnnotation/NCBI/hg38.genomeCoor.bed")

#### We use enriched_RNAs_2 finally, 这是一个更松的条件
enriched_RNAs_1,enriched_RNAs_2 = readEnrichedRNAs()
Cell2014_peaks_1,Cell2014_peaks_2 = readCell2014DicerCLIPTable()

enrich_withCellTable = compare_Seqs_And_Locus(enriched_RNAs_2, Cell2014_peaks_1, tolerance=0)
print("Enriched RNAs:{}; Cell2014_peaks:{}; Overlap:{}".format(len(enriched_RNAs_2), len(Cell2014_peaks_1), len(enrich_withCellTable)))

shape_1000, seq_1000 = DicerMiRNA.read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000")
wc_tid_1000 = filter_wellcovtrans(shape_1000); print(len(wc_tid_1000))
enrich_withStructure_1000 = set(wc_tid_1000)&set(enriched_RNAs_2); print( len(enrich_withStructure_1000) )

def count_RNA_type(RNA_ids):
    RNAtypes = [ it.split('_')[0] for it in RNA_ids ]
    print("Gene type", "Count", sep="\t")
    Count = []
    for Type in set(RNAtypes):
        if Type.startswith('rRNA') or Type=='chrM': continue
        Count.append([Type, RNAtypes.count(Type)])
    Count.sort(key=lambda x:x[1],reverse=True)
    for rnatype,count in Count:
        print(rnatype, count, sep="\t")



count_RNA_type(wc_tid_1000)   ####   Figure S2F
count_RNA_type(enriched_RNAs_2)   ####   Figure S2C
count_RNA_type(enrich_withStructure_1000)  ####   Figure 2B
count_RNA_type( set([it[0] for it in enrich_withCellTable]) )

##########  Save enriched genes with structures

OUT = open("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/6.enrich_with_structure/enrich_with_structure.txt", 'w')
for tid in enrich_withStructure_1000:
    print(tid, file=OUT)

OUT.close()

################################
####  Table S1
####  Some functions to process Cell Peaks
################################

def WhatIsIt(chrID,chrStart,chrEnd,strand,GAPer):
    transList = GAPer.genomeCoor2transCoor(chrID,chrStart,chrEnd,strand)
    transIDList = [ it[3] for it in transList ]
    ID = ""
    if transList:
        for tid in transIDList:
            ft = GAPer.getTransFeature(tid)
            if ft['gene_type'] in ('guide_RNA','snRNA','snoRNA','tRNA','precursor_RNA','rRNA','scRNA','misc_RNA'):
                RNAType=ft['gene_type'].replace('guide_RNA','guideRNA')
                RNAType=RNAType.replace('precursor_RNA','miRNA')
                RNAType=RNAType.replace('misc_RNA','miscRNA')
                ID = RNAType+"_"+ft['gene_name']
                break
        if ID == "":
            for tid in transIDList:
                ft = GAPer.getTransFeature(tid)
                if ft['gene_type'] == 'mRNA':
                    ID = "mRNAExon_"+ft['gene_name']
                    break
        if ID == "":
            for tid in transIDList:
                ft = GAPer.getTransFeature(tid)
                if ft['gene_type'] == 'lnc_RNA':
                    ID = "lncRNAExon_"+ft['gene_name']
                    break
        if ID == "":
            tid = transIDList[0]
            ft = GAPer.getTransFeature(tid)
            RNAType=ft['gene_type'].replace('guide_RNA','guideRNA')
            RNAType=RNAType.replace('precursor_RNA','miRNA')
            RNAType=RNAType.replace('misc_RNA','miscRNA')
            RNAType=RNAType.replace('lnc_RNA','lncRNA')
            ID = RNAType+"Exon_"+ft['gene_name']
    else:
        geneList = GAPer.genomeCoor2geneCoor(chrID,chrStart,chrEnd,strand)
        geneIDList = [ it[3] for it in geneList ]
        if geneList:
            for gid in geneIDList:
                ft = GAPer.getGeneParser()[gid]
                if 'mRNA' in ft['gene_type']:
                    ID = "mRNAIntron_"+ft['gene_name']
            if ID == "":
                for gid in geneIDList:
                    ft = GAPer.getGeneParser()[gid]
                    if 'lnc_RNA' in ft['gene_type']:
                        ID = "lncRNAIntron_"+ft['gene_name']
            if ID == "":
                gid = geneIDList[0]
                ft = GAPer.getGeneParser()[gid]
                RNAType=ft['gene_type'][0].replace('guide_RNA','guideRNA')
                RNAType=RNAType.replace('precursor_RNA','miRNA')
                RNAType=RNAType.replace('lnc_RNA','lncRNA')
                RNAType=RNAType.replace('misc_RNA','miscRNA')
                ID = RNAType+"Intron_"+ft['gene_name']
        else:
            ID = "intergenic"
    
    return ID

def annotate_peaks(peaks,gaper,seqer):
    annotPeaks = []
    for chrID,start,end,strand in peaks:
        try:
            gene_type = WhatIsIt(chrID,start,end,strand,gaper)
            seq = seqer.fetch(chrID,start,end,strand)
        except KeyError:
            continue
        genome_locus = [chrID,start,end,strand]
        annotPeaks.append([gene_type,seq,genome_locus])
    return annotPeaks

def inFastaDB(querySeq, fastaDB, overlap=100):
    for tid in fastaDB:
        start = querySeq.find(fastaDB[tid])
        if start == -1: continue
        end = start + len(fastaDB[tid])
        if end>overlap and start<len(querySeq)-overlap:
            return tid
    return False

def annotate_with_smallRNA(annotPeak, small_RNA, Seqer):
    for item in annotPeak:
        chrID,start,end,strand = item[2]
        longer_seq = Seqer.fetch(chrID,max(start-500,1),end+500,strand)
        tid = inFastaDB(longer_seq, small_RNA, 500)
        if tid:
            print(item[0],"=>",tid,sep="  ")
            item[0] = tid

Annot_cell_peaks = annotate_peaks(Cell2014_peaks_1, Gaper, Seqer)
annotate_with_smallRNA(Annot_cell_peaks, smallRNASeqs, Seqer)

count_RNA_type([it[0] for it in Annot_cell_peaks])



