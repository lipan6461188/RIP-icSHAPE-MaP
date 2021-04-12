
importCommon()
import DicerMiRNA

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

def build_tRNA_structurome(inFile):
    tRNA_structurome = {}
    tRNA_dot = General.load_dot(inFile)
    for tid in tRNA_dot:
        seq = tRNA_dot[tid][0].upper()
        dot = (tRNA_dot[tid][1], tid)
        tRNA_structurome[seq] = dot
    return tRNA_structurome

#########################
## 1. Read shape and sequence
#########################

shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000"
SHAPE, Sequence = read_shape(shape_folder); print (len(SHAPE))
tRNA_structurome = build_tRNA_structurome("/150T/zhangqf/GenomeAnnotation/tRNA/GtRNAdb/hg19/human_tRNA.dot")

auc_list = []
for tid in Sequence:
    seq = Sequence[tid].replace('U','T')
    if seq in tRNA_structurome:
        dot = tRNA_structurome[seq][0]
        shape = SHAPE[tid]
        if 1-shape.count('NULL')/len(shape)<0.6:
            continue
        if len(dot)!=len(shape):
            print (dot)
            print (shape)
            break
        auc = General.calc_AUC_v2(dot, shape)
        auc_list.append(auc)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
ax.set_title(str(len(auc_list))+' tRNAs')
ax.set_ylabel('AUC')
data = [ auc_list ]
colors = [Colors.RGB['blue']]
Figures.violinPlot(ax, data, ['tRNA'],colors=colors)
fig.tight_layout()
fig.savefig("figs/tRNA.pdf")
fig.show()


# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_2d.txt'), 'w')

auc_list = []
for tid in Sequence:
    seq = Sequence[tid].replace('U','T')
    if seq in tRNA_structurome:
        dot = tRNA_structurome[seq][0]
        shape = SHAPE[tid]
        if 1-shape.count('NULL')/len(shape)<0.6:
            continue
        if len(dot)!=len(shape):
            print (dot)
            print (shape)
            break
        auc = General.calc_AUC_v2(dot, shape)
        print(f"{tid}\t{auc}", file=OUT)

OUT.close()






