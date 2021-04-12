
importCommon()
import DicerMiRNA

#################################
##### Prepare miRNA
#################################

def miRNA_nameConversion():
    nameConversion = {}
    
    miRNA_ref_fa = "/Share/home/zhangqf7/lipan/reference/miRNA_tRNA_snRNA_snoRNA_rRNA_mtRNA_YRNA/miRNA/hsa.fa"
    for line in open(miRNA_ref_fa):
        if line[0] == '>':
            raw_id,new_id = re.findall("Alias=(\\w+);Name=([\\w\\-]+)\\(", line)[0]
            nameConversion[raw_id] = new_id
    
    return nameConversion

def trim_MiRBase_dot(miRNA_shape, miRNA_shape_seq, MiRBase):
    import re
    
    trimmed_shape = {}
    common_id = set(miRNA_shape)&set(MiRBase)
    print("Common id: "+str(len(common_id)))
    for miRID in common_id:
        mir_seq = MiRBase[miRID][0]
        shape_list = miRNA_shape[miRID]
        ## 1. Trim mir_seq
        lindex = re.search("[AUCG]", mir_seq).start()
        rindex = re.search("[AUCG][aucg]*$", mir_seq).start()
        mature_miRNA = mir_seq[lindex:rindex+1].upper()
        ## 2. locate shape seq
        index = miRNA_shape_seq[miRID].find(mature_miRNA)
        if index == -1:
            print("Warning: {} not found".format(miRID))
            continue
        assert len(miRNA_shape[miRID]) == len(miRNA_shape_seq[miRID])
        mature_miRNA_shape = miRNA_shape[miRID][index:index+len(mature_miRNA)]
        
        trimmed_shape[miRID] = ( mature_miRNA_shape, mature_miRNA )
    
    return trimmed_shape

def predict_miRNA_structure(trimmed_shape, OUT=sys.stdout):
    import Colors
    miRStruct = {}
    for miRID in trimmed_shape:
        shape_list, seq = trimmed_shape[miRID][:2]
        if len(seq)<35:
            continue
        curshape = shape_list[:]
        for i in range(len(curshape)):
            if curshape[i]!='NULL':
                curshape[i] = max(min(curshape[i],1),0)
        unconstraint = Structure.predict_structure(seq.upper())
        constraint = Structure.predict_structure(seq.upper(), curshape)
        
        auc = General.calc_AUC_v2(constraint, curshape)
        
        print(">"+miRID, round(auc,3), sep="\t", file=OUT)
        print(seq, file=OUT)
        print(unconstraint,"unconstraint",sep="\t", file=OUT)
        print(constraint,"constraint",sep="\t", file=OUT)
        print(Colors.color_SHAPE(curshape), file=OUT)
        OUT.flush()
        
        miRStruct[miRID] = ( curshape, seq, constraint, unconstraint, auc )
    return miRStruct

def write_miRNA_dot(miRStruct, outFolder):
    outFolder = outFolder.rstrip('/')+'/'
    
    constraint = { k:(miRStruct[k][1],miRStruct[k][2]) for k in miRStruct }
    unconstraint = { k:(miRStruct[k][1],miRStruct[k][3]) for k in miRStruct }
    auc = [ (k,miRStruct[k][4]) for k in miRStruct ]
    open(outFolder+'auc.txt', 'w').writelines( "\n".join([ "{}\t{}".format(it[0],it[1]) for it in auc]) )
    General.write_dot(constraint, outFolder+'constraint.dot')
    General.write_dot(unconstraint, outFolder+'unconstraint.dot')
    ## split version
    current_tid_list = list(constraint.keys())
    
    print ("finish...")

miRNAmapName = miRNA_nameConversion()

inFn = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000"
shape_1000, sequence_1000 = DicerMiRNA.read_shape(inFn, min_valid_ratio=0.2, min_valid_base_num=30)
miRNA_shape = { miRNAmapName[tid.split('_')[1]]:shape_1000[tid] for tid in shape_1000 if tid.startswith('miRNA') }
miRNA_seq = { miRNAmapName[tid.split('_')[1]]:sequence_1000[tid] for tid in sequence_1000 if tid.startswith('miRNA') }

dot_miRNA = General.load_dot("/150T/zhangqf/GenomeAnnotation/miRNA/more/miRNA.dot"); print(len(dot_miRNA))

trimmed_shape =  trim_MiRBase_dot(miRNA_shape, miRNA_seq, dot_miRNA); print(len(trimmed_shape))
out_file = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/structures.txt"
miRStruct = predict_miRNA_structure(trimmed_shape, OUT=open(out_file, 'w')); print(len(miRStruct))

outFolder = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514"
write_miRNA_dot(miRStruct, outFolder)




