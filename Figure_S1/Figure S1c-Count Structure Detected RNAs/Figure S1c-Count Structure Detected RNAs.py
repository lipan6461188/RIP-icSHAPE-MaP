
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
            #if shape.count('NULL')/len(shape)>0.4: continue
            SHAPE[tid] = shape
            Sequence[tid] = seq
    return SHAPE, Sequence

def filter_wellcovtrans(shapeDict, min_cov_num=30, min_cov_ratio=0.5):
    valid_tid_list = []
    for tid in shapeDict:
        shapemap =shapeDict[tid]
        validnum = len(shapemap)-shapemap.count("NULL")
        if validnum >= min_cov_num and 1.0*validnum/len(shapemap) >= min_cov_ratio:
            valid_tid_list.append(tid)
    return valid_tid_list

def count_rnatype(inputTidList):
    typeCount = {}
    for Tid in inputTidList:
        rna_type = Tid.split('_')[0]
        typeCount[rna_type] = typeCount.get(rna_type, 0) + 1
    return typeCount


#########################
## 1. Read shape and sequence
#########################

shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vitro_SMR_SSII_repX/shape_files/"
SHAPE, Sequence = read_shape(shape_folder); print (len(SHAPE))

#########################
## 2. Filtering
#########################

vfl1 = filter_wellcovtrans(SHAPE, min_cov_num=30, min_cov_ratio=0.5); print( len(vfl1) )
typecount1 = count_rnatype(vfl1)
rna_types = [ it[0] for it in sorted(list(typecount1.items()), key=lambda x: x[1], reverse=True) ]
typecount1 = [ (rt, typecount1[rt]) for rt in rna_types ]
df1 = pd.DataFrame(data=typecount1, columns=['type', 'number']); df1.index = df1.type
df1.to_csv("/Share/home/zhangqf7/figs/vitro.csv", sep="\t")


#########################
## 1. Read shape and sequence
#########################

shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files/"
SHAPE, Sequence = read_shape(shape_folder); print (len(SHAPE))

#########################
## 2. Filtering
#########################

vfl1 = filter_wellcovtrans(SHAPE, min_cov_num=30, min_cov_ratio=0.5); print( len(vfl1) )
typecount1 = count_rnatype(vfl1)
rna_types = [ it[0] for it in sorted(list(typecount1.items()), key=lambda x: x[1], reverse=True) ]
typecount1 = [ (rt, typecount1[rt]) for rt in rna_types ]
df1 = pd.DataFrame(data=typecount1, columns=['type', 'number']); df1.index = df1.type
df1.to_csv("/Share/home/zhangqf7/figs/vivo.csv", sep="\t")








