
import General
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_mutnum(inMutParserFn, sequence, minDepth):
    baseMutCount = {'A':[], 'T':[], 'C':[], 'G':[]}
    sequence = sequence.upper()
    
    lineCount = 0
    for line in open(inMutParserFn):
        lineCount += 1
        if lineCount == 1:
            continue
        data = line.strip().split()
        delnum = sum( [int(it) for it in data[:4]] )
        insnum = sum( [int(it) for it in data[4:9]] )
        mutnum = sum( [int(it) for it in data[9:21]] )
        muldel = int(data[21])
        mulins = int(data[22])
        mulmut = int(data[23])
        comcha = int(data[24]) + int(data[25])
        effdep = int(data[27])
        
        base = sequence[lineCount-2]
        ## Mismatch, Deletion, Insertion, Mismatch(multi), Deletion(multi), Insertion(multi), Complex
        if effdep >= minDepth:
            baseMutCount[base].append( [mutnum,delnum,insnum,mulmut,muldel,mulins,comcha,effdep] )
    
    if len(sequence) != lineCount-1:
        sys.stderr.writelines("Error: Different length")
    
    return baseMutCount

def calc_bulk_mutratio(baseMutCount, max_mutratio=0.2):
    baseMutRatio = {}
    discard = 0
    for base in baseMutCount:
        baseMutRatio[base] = []
        for mutnum,delnum,insnum,mulmut,muldel,mulins,comcha,effdep in baseMutCount[base]:
            mutratio = round(1.0*mutnum/effdep,5)
            delratio = round(1.0*delnum/effdep,5)
            insratio = round(1.0*insnum/effdep,5)
            mulmutratio = round(1.0*mulmut/effdep,5)
            muldelratio = round(1.0*muldel/effdep,5)
            mulinsratio = round(1.0*mulins/effdep,5)
            comcharatio = round(1.0*comcha/effdep,5)
            total_ratio = round(1.0*(mutnum+delnum+insnum+mulmut+muldel+mulins+comcha)/effdep,5)
            if total_ratio >= max_mutratio:
                discard += 1
                continue
            baseMutRatio[base].append( [total_ratio, mutratio, delratio, insratio, mulmutratio, muldelratio, mulinsratio, comcharatio, effdep] )
    
    print ("filter too large mutation ratio -- ", discard)
    return baseMutRatio

def combine_mutnum(baseMutCount_list):
    import copy
    combinedBaseMutCount = copy.deepcopy(baseMutCount_list[0])
    for base in combinedBaseMutCount:
        i = 1
        while i < len(baseMutCount_list):
            combinedBaseMutCount[base] += baseMutCount_list[i][base]
            i += 1
    return combinedBaseMutCount

def isValidMutnum(baseMutCount, minCount):
    isValid = True
    for base in ('A','T','C','G'):
        if len(baseMutCount[base])<minCount:
            isValid = False
    return isValid

ref_seq = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa")
Mut_root = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/5.collect_mutation/"
my_samples = ['DMSO_SMR_SSII_repX', 'NAI_100mm_vivo_SMR_SSII_repX', 'NAI_100mm_vitro_SMR_SSII_repX']

####################################
#### 统计方法: 每一个位点一个值
####################################

only_miRNA = False
mutRatioDict = {}
for libname in my_samples:
    print (libname)
    
    if only_miRNA:
        libFiles = os.popen("ls "+Mut_root+"miRNA*/"+libname).readlines()
        rRNA_5S = os.popen("ls "+Mut_root+"rRNA_human_5S/"+libname).readlines()
        libFiles = [ it.rstrip() for it in libFiles+rRNA_5S ]
    else:
        libFiles = [ it.rstrip() for it in os.popen("ls "+Mut_root+"*/"+libname).readlines() ]
    
    mutCount_list = []
    for parseFile in libFiles:
        tid = parseFile.split('/')[-2]
        sequence = ref_seq[tid].upper()
        if 'N' in sequence: continue 
        baseMutCount = read_mutnum(parseFile, sequence, minDepth=5000)
        if isValidMutnum(baseMutCount, minCount=10):
            mutCount_list.append( baseMutCount )
    
    combinedBaseMutCount = combine_mutnum(mutCount_list)
    if 'DMSO' in libname:
        mutratio = calc_bulk_mutratio(combinedBaseMutCount, max_mutratio=0.05)
    else:
        mutratio = calc_bulk_mutratio(combinedBaseMutCount, max_mutratio=0.20)
    
    mutRatioDict[libname] = mutratio
    datasize = "datasize: %s\t%s\t%s\t%s" % (len(mutratio['A']), len(mutratio['T']), len(mutratio['C']), len(mutratio['G']))
    print(datasize)

Box_list = []
for libname in my_samples:
    Box_list.append( [it[0] for it in mutRatioDict[libname]['A'] ] )
    Box_list.append( [it[0] for it in mutRatioDict[libname]['T'] ] )
    Box_list.append( [it[0] for it in mutRatioDict[libname]['C'] ] )
    Box_list.append( [it[0] for it in mutRatioDict[libname]['G'] ] )

import Figures
fig = plt.figure(1, figsize=(5, 6))
ax = fig.add_subplot(111)
labels = ['A','T','C','G']*3
facecolors = [Colors.RGB['red'],Colors.RGB['purple'],Colors.RGB['blue'],Colors.RGB['green']]*3
Figures.boxPlot(ax, Box_list, labels=labels, facecolors=facecolors, showOutliers=False)
fig.savefig("/Share2/home/zhangqf7/figs/fig2.pdf")
plt.show()


# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_1d.txt'), 'w')

for libname in my_samples:
    print(libname+"\tA\t", file=OUT, end="")
    print( "\t".join([str(it[0]) for it in mutRatioDict[libname]['A'] ]), file=OUT )
    print(libname+"\tT\t", file=OUT, end="")
    print( "\t".join([str(it[0]) for it in mutRatioDict[libname]['T'] ]), file=OUT )
    print(libname+"\tC\t", file=OUT, end="")
    print( "\t".join([str(it[0]) for it in mutRatioDict[libname]['C'] ]), file=OUT )
    print(libname+"\tG\t", file=OUT, end="")
    print( "\t".join([str(it[0]) for it in mutRatioDict[libname]['G'] ]), file=OUT )

OUT.close()









