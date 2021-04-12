
importCommon()
import DicerMiRNA

#########################
## Functions: About SHAPE and Structure
#########################

class TransInfo:
    def __init__(self, name):
        self.full_shape = []
        self.full_seq = ""
        self.stem_loop = [0,0,0,0]
        self.full_dot = ""
        self.name = name
        self.enrich_score = 0
        self.cleavage_score = 0
        
        self.stemloop_seq = ""
        self.stemloop_shape = []
        self.stemloop_dot = []
        self.ago2 = 0
        self.ago3 = 0
        
        self.ago2_normed = 0
        self.ago3_normed = 0
        
        self.rnatype = ""
        
        self.cluster_id = -1
    def stemLength(self):
        return ((self.stem_loop[1]-self.stem_loop[0])+(self.stem_loop[3]-self.stem_loop[2]))/2

def PredictStructure(dataframe, Sequence, SHAPE):
    TransList = []
    
    for tid in dataframe.index:
        if tid not in SHAPE:
            continue
        print("Process ", tid, "...")
        
        trans = TransInfo(tid)
        seq = Sequence[tid]
        shape = SHAPE[tid]
        if len(seq)>250:
            continue
        if tid not in dataframe.index:
            continue
        enrich_score = round(dataframe.loc[tid, 'RIP/INPUT'],3)
        cleavage_score = round(dataframe.loc[tid, 'OEMut/OEWT'],3)
        ago2 = round(dataframe.loc[tid, 'AGO2/INPUT'],3)
        ago3 = round(dataframe.loc[tid, 'AGO3/INPUT'],3)
        dot = Structure.predict_structure(seq, shape)
        #stem_loops = Structure.find_stem_loop(dot, max_loop_len=15, max_stem_gap=5, min_stem_len=5)
        if enrich_score>1:
            stem_loops = Structure.find_stem_loop(dot, max_loop_len=15, max_stem_gap=5, min_stem_len=10)
        else:
            stem_loops = Structure.find_stem_loop(dot, max_loop_len=10, max_stem_gap=3, min_stem_len=4)
        if len(stem_loops) == 0:
            continue
        
        max_sl = []
        if enrich_score>1:
            aucs = []
            for stemloop in stem_loops:
                subdot = dot[stemloop[0]-1:stemloop[3]]
                subshape = shape[stemloop[0]-1:stemloop[3]]
                if len(subshape)-subshape.count('NULL')<10:
                    continue
                auc = General.calc_AUC_v2(subdot, subshape)
                aucs.append([stemloop,auc])
            aucs.sort(key=lambda x: x[1],reverse=True)
            if len(aucs)==0:
                continue
            max_sl,sl_auc = aucs[0]
        else:
            mid = int( len(stem_loops)/2 )
            max_sl = stem_loops[mid]
        
        trans.full_shape = shape
        trans.full_seq = seq
        trans.stem_loop = max_sl
        trans.full_dot = dot
        trans.enrich_score = enrich_score
        trans.cleavage_score = cleavage_score
        trans.ago2_normed = ago2
        trans.ago3_normed = ago3
        
        TransList.append(trans)
    
    print ("Raw transcripts:", len(SHAPE))
    print ("Transcripts with stemloop:", len(TransList))
    
    return TransList

def RNAType_To_Colors(RNAType_list):
    RNAColor = {'mRNAExon': Colors.RGB['amber'],\
        'mRNAIntron': Colors.RGB['amber'],\
        'intergenic': Colors.RGB['amber'],\
        'mRNA-intergenic-lncRNA': Colors.RGB['amber'],\
        'SNORA': Colors.RGB['light_green'],\
        'SNORD': Colors.RGB['blue'],\
        'miRNA': Colors.RGB['green'],\
        'tRNA': Colors.RGB['deep_orange'],\
        'snRNA': Colors.RGB['brown']}
    return [ RNAColor.get(rt,'#f1f1f1') for rt in RNAType_list ]

def is_valid_RNA(RNA_id):
    if RNA_id.startswith('mRNA') or \
        RNA_id.startswith('intergenic') or \
        RNA_id.startswith('lncRNA') or \
        RNA_id.startswith('snoRNA') or \
        RNA_id.startswith('miRNA') or \
        RNA_id.startswith('tRNA') or \
        RNA_id.startswith('miscRNA') or \
        RNA_id.startswith('snRNA'):
        return True
    return False

def read_depth(inFn):
    Depth = {}
    for line in open(inFn):
        tid,pos,cov = line.strip().split()
        if tid not in Depth:
            Depth[tid] = []
        Depth[tid].append(int(cov))
    return Depth

def stemLoopGetSHAPE(TransList, SHAPE):
    TransWithStemLoop = []
    for trans in TransList:
        tid = trans.name
        sl = trans.stem_loop
        dot = trans.full_dot
        seq = trans.full_seq
        shape = trans.full_shape
        
        ls,le,rs,re = sl
        sl_mid = int((le+rs)/2)
        pseudo_shape = ['NULL']*60+[min(max(float(it),0),1) if it!='NULL' else 'NULL' for it in shape]+['NULL']*60
        pseudo_dot = '.'*60+dot+'.'*60
        pseudo_seq = 'N'*60+seq+'N'*60
        start = sl_mid+60-30
        end = sl_mid+60+30
        
        seg_seq = pseudo_seq[start:end]
        seg_shape = pseudo_shape[start:end]
        seg_dot = pseudo_dot[start:end]
        
        if seg_shape.count('NULL')>5:
            continue
        if not is_valid_RNA(trans.name):
            print("Warning: filter ", trans.name, " Not a valid RNA")
            continue
        
        trans.stemloop_seq = seg_seq
        trans.stemloop_shape = seg_shape
        trans.stemloop_dot = seg_dot
        TransWithStemLoop.append(trans)
    
    print("Raw number of transcript:",len(TransList))
    print("Number of transcripts with structures:",len(TransWithStemLoop))
    return TransWithStemLoop

def stemLoopGetAgo(TransList, AGO2, AGO3):
    for trans in TransList:
        tid = trans.name
        sl = trans.stem_loop
        left,right = sl[0],sl[3]
        if tid in AGO2:
            trans.ago2 = max(AGO2[tid][left:right])
        if tid in AGO3:
            trans.ago3 = max(AGO3[tid][left:right])

def stemLoopGetEnrich(TransList, normed_df_detected):
    for trans in TransList:
        tid = trans.name
        if tid in normed_df_detected.index:
            trans.enrich_score = round(normed_df_detected.loc[tid,'RIP/INPUT'],3)
            trans.cleavage_score = round(normed_df_detected.loc[tid,'OEMut/OEWT'],3)
        else:
            trans.enrich_score = trans.cleavage_score = np.nan

def renameTrans(trans_list, Gaper):
    for trans in trans_list:
        name_items = trans.name.split('_')
        if name_items[1].startswith('ENST'):
            gene_name = Gaper.getTransFeature(name_items[1])['gene_name']
            trans.name = name_items[0]+"_"+name_items[1]+"_"+gene_name
        trans.rnatype = trans.name.split('_')[0]
        if 'SNORA' in trans.name:
            trans.rnatype = 'SNORA'
        if 'SNORD' in trans.name:
            trans.rnatype = 'SNORD'
        if trans.rnatype in ('mRNAExon','mRNAIntron','intergenic','lncRNA','miscRNA'):
            trans.rnatype = 'mRNA-intergenic-lncRNA'

#########################
## Read shape and sequence
#########################

shape_1000, seq_1000 = DicerMiRNA.read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000")
Gaper = GAP.init("/150T/zhangqf/GenomeAnnotation/Gencode/hg38.genomeCoor.bed")

AGO2 = read_depth("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/3.mapBWA/AGO2.depth")
AGO3 = read_depth("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/3.mapBWA/AGO3.depth")

#########################
## Read RNA-Seq data & Get SHAPE & Add Ago
#########################

def get_normed_rnaseq_table():
    #########################
    ## Read data and combine
    #########################
    genecountFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/3.organized_table/RNA-Seq-genecount.txt"
    df = pd.read_csv(genecountFn,sep="\t",index_col=0)
    
    AGO2 = read_depth("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/3.mapBWA/AGO2.depth")
    AGO3 = read_depth("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Cell_2014/3.mapBWA/AGO3.depth")
    df['AGO2'] = [0] * df.shape[0]
    df['AGO3'] = [0] * df.shape[0]
    for tid in df.index:
        if tid in AGO2:
            df.loc[tid, 'AGO2'] = sum(AGO2[tid])
    for tid in df.index:
        if tid in AGO3:
            df.loc[tid, 'AGO3'] = sum(AGO3[tid])
    
    FSS = df['FSS_rep1']+df['FSS_rep2']+df['FSS_rep3']+df['FSS_rep4']
    OEWT = df['OEWTDicer_rep1']+df['OEWTDicer_rep2']+df['OEWTDicer_rep3']+df['OEWTDicer_rep4']
    OEMut = df['OEMutDicer_rep1']+df['OEMutDicer_rep2']+df['OEMutDicer_rep3']+df['OEMutDicer_rep4']
    INPUT = df['INPUT_rep1_20190514']+df['INPUT_rep2_20190514']+df['INPUT_rep3_20190514']
    RIP = df['RIP_DMSO_rep1_20190514']+df['RIP_DMSO_rep2_20190514']
    AGO2 = df['AGO2']
    AGO3 = df['AGO3']
    concat_df = pd.concat([INPUT, RIP, FSS, OEWT, OEMut, AGO2, AGO3],axis=1)
    concat_df.columns = ['INPUT', 'RIP', 'FSS', 'OEWT', 'OEMut', 'AGO2', 'AGO3']
    RNASeq_high = concat_df.loc[:,('FSS','OEWT','OEMut')].sum(axis=1)>30
    
    #########################
    ## Normalize Data
    #########################
    libsize = concat_df.sum()
    
    normed_df = concat_df/libsize
    normed_df = normed_df.loc[(concat_df['RIP']>10)&(RNASeq_high),:]
    normed_df['RIP/INPUT'] = np.log2(normed_df.RIP/(normed_df.INPUT+1e-15))
    normed_df['FSS/OEWT'] = np.log2(normed_df.FSS/(normed_df.OEWT+1e-15))
    normed_df['OEMut/OEWT'] = np.log2(normed_df.OEMut/(normed_df.OEWT+1e-15))
    normed_df['AGO2/INPUT'] = np.log2(normed_df.AGO2/(normed_df.INPUT+1e-15))
    normed_df['AGO3/INPUT'] = np.log2(normed_df.AGO3/(normed_df.INPUT+1e-15))
    normed_df = normed_df.sort_values(by='RIP/INPUT', ascending=False)
    
    normed_df.loc[normed_df['FSS/OEWT']==numpy.Infinity,'FSS/OEWT'] = numpy.nan
    normed_df.loc[normed_df['OEMut/OEWT']==numpy.Infinity,'OEMut/OEWT'] = numpy.nan
    normed_df.loc[normed_df['FSS/OEWT']==-numpy.Infinity,'FSS/OEWT'] = numpy.nan
    normed_df.loc[normed_df['OEMut/OEWT']==-numpy.Infinity,'OEMut/OEWT'] = numpy.nan
    normed_df.loc[normed_df['AGO2/INPUT']==-numpy.Infinity,'AGO2/INPUT'] = numpy.nan
    normed_df.loc[normed_df['AGO3/INPUT']==-numpy.Infinity,'AGO3/INPUT'] = numpy.nan
    
    return normed_df

normed_df_detected = get_normed_rnaseq_table(); print(normed_df_detected.shape)
normed_df_detected = normed_df_detected.loc[ set(shape_1000.keys())&set(normed_df_detected.index) ,:].sort_values(by='RIP/INPUT', ascending=False)
print(normed_df_detected.shape)

TransList = PredictStructure(normed_df_detected, seq_1000, shape_1000)
renameTrans(TransList, Gaper)

TransWithStemLoop = stemLoopGetSHAPE(TransList, shape_1000)
stemLoopGetAgo(TransWithStemLoop, AGO2, AGO3)

len(normed_df_detected)
len(TransList)
len(TransWithStemLoop)

trans_dict = { trans.name:trans for trans in TransWithStemLoop }

#########################
## Prepare the data
#########################

def show_rna_ratio(tid_list):
    Len = len(tid_list)
    tRNA  = len([ name for name in tid_list if 'tRNA' in name ]) / Len * 100
    SNORA = len([ name for name in tid_list if 'SNORA' in name ]) / Len * 100
    SNORD = len([ name for name in tid_list if 'SNORD' in name ]) / Len * 100
    miRNA = len([ name for name in tid_list if 'miRNA' in name ]) / Len * 100
    snRNA = len([ name for name in tid_list if 'snRNA' in name ]) / Len * 100
    other = 100-tRNA-SNORA-SNORD-miRNA-snRNA
    print("tRNA: %.2f%%; SNORA: %.2f%%; SNORD: %.2f%%; miRNA: %.2f%%; snRNA: %.2f%%; other: %.2f%%; Number: %d" % (tRNA, SNORA, SNORD, miRNA, snRNA, other, Len))

def dfshape(trans_list, leave=['tRNA','miRNA','SNORA','SNORD']):
    shape_list = []
    tid_list = []
    type_list = []
    enrich_list = []
    cleavage_list = []
    for trans in trans_list:
        mytype = trans.rnatype
        if mytype not in ('tRNA','miRNA','SNORA','SNORD','mRNA-intergenic-lncRNA'):
            mytype = 'other'
        type_list.append(mytype)
        shape_array = trans.stemloop_shape[:]
        tid_list.append(trans.name)
        shape_list.append( [ it if it!='NULL' else np.nan for it in shape_array ] )
        enrich_list.append(trans.enrich_score)
        cleavage_list.append(trans.cleavage_score)
    df = pd.DataFrame(shape_list, index=tid_list)
    
    for i in range(df.shape[1]):
        df.loc[df.iloc[:,i].isna(),i] = df.mean()[i]
    
    Filter = [it in leave for it in type_list]
    df = df.loc[Filter, :]
    type_list = [ it for it in type_list if it in leave ]
    enrich_list = [ enrich_list[i] for i,it in enumerate(type_list) if it in leave ]
    cleavage_list = [ cleavage_list[i] for i,it in enumerate(type_list) if it in leave ]
    
    return df, type_list, enrich_list, cleavage_list

df, type_list, enrich_list, cleavage_list = dfshape(TransWithStemLoop, leave=['tRNA','miRNA','SNORA','SNORD','mRNA-intergenic-lncRNA'])


