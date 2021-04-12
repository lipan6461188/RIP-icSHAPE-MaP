
importCommon()

high_exp_RNAs = General.load_fasta("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa")

############################
######  读取Coverage
############################

def readCoverage(depthFn, Sequence):
    depth = {}
    for line in open(depthFn):
        data = line.strip().split()
        pos = int(data[1])
        nucdepth = int(data[2])
        if data[0] not in depth:
            depth[data[0]] = [0]*len(Sequence[data[0]])
        depth[data[0]][pos-1] = nucdepth
    return depth

inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/how_many_premiRNA/INPUT_repX_20190514.bedGraph"
DMSO_SMR_SSII_repX = readCoverage(inFn, high_exp_RNAs); print(len(DMSO_SMR_SSII_repX))

def readReadsNum(countFn):
    Count = {}
    for line in open(countFn):
        count,tid = line.strip().split()
        Count[tid] = int(count)
    return Count

inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/how_many_premiRNA/INPUT_repX_20190514.count"
DMSO_SMR_SSII_repX_count = readReadsNum(inFn); print(len(DMSO_SMR_SSII_repX_count))

############################
######  判断miRNA时候被完全覆盖
############################

def get_pre_miRNA(Depth):
    effect_pre_miRNA = []
    for tid in set(Depth):
        depth = Depth[tid]
        detc_start = 0
        detc_end = 0
        for i in range(len(depth)):
            if depth[i]>2:
                detc_start = i+1
                break
        for i in range(len(depth)-1,0,-1):
            if depth[i]>2:
                detc_end = i
                break
        if detc_start<detc_end:
            effective_region = depth[detc_start:detc_end]
            mid = len(effective_region)//2
            mid_region = effective_region[mid-20:mid+20]
            left_region = effective_region[mid-20:mid]
            right_region = effective_region[mid:mid+20]
            if len(effective_region)>40 and 0.5<np.mean(left_region)/np.mean(right_region)<2 and np.mean(mid_region)/np.mean(effective_region)>0.8 and min(mid_region)>0:
                effect_pre_miRNA.append(tid)
    return effect_pre_miRNA

DMSO_SMR_SSII_repX_miRNAs = get_pre_miRNA(DMSO_SMR_SSII_repX); print(len(DMSO_SMR_SSII_repX_miRNAs))

############################
######  Get miRNAs with Structures
############################

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

SHAPE, Sequence = read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000/")
detected_miRNAs = [ it for it in SHAPE if it.startswith('miRNA') ]; print(len(detected_miRNAs))

############################
######  Overlap
############################

len( set(detected_miRNAs)&set(DMSO_SMR_SSII_repX_miRNAs) )

total_cov = 0
detected_cov = 0
for tid in DMSO_SMR_SSII_repX_miRNAs:
    total_cov += DMSO_SMR_SSII_repX_count[tid]
    if tid in detected_miRNAs:
        detected_cov += DMSO_SMR_SSII_repX_count[tid]

print(detected_cov/total_cov)

############################
######  计算相对表达量
############################

relative_exp = []
for tid in DMSO_SMR_SSII_repX_miRNAs:
    exp = DMSO_SMR_SSII_repX_count[tid] / len(high_exp_RNAs[tid])
    relative_exp.append([tid, exp])

relative_exp.sort(key=lambda x: x[1],reverse=True)

colors = []
for tid,exp in relative_exp:
    if tid in detected_miRNAs:
        colors.append('green')
    else:
        print(tid)
        colors.append('blue')

plt.figure(figsize=(15,5))
plt.bar(range(len(relative_exp)), [it[1] for it in relative_exp], color=colors)
plt.ylim(0, 100)
plt.savefig("figs/RIP_DMSO_repX_20190514.pdf")
plt.show()


# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_2c.txt'), 'w')
for name,exp in relative_exp[:150]: 
    if name in detected_miRNAs:
        print(f"{name}\t{exp}\tYes", file=OUT)
    else:
        print(f"{name}\t{exp}\tNo", file=OUT)

OUT.close()




