################################
####  Figure S2E
####  Dicer结合区域片段的长度分布
################################
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import Figures, Colors
import pandas as pd
import numpy as np
import os

inFn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/which_rna_enriched/5.foldChange/foldChange.txt"
df = pd.read_csv(inFn, sep="\t", index_col=0)
enriched_RNAs_2 = df.loc[(df['log2FC']>0),:].index.tolist()

def peak_boundary(a_list, gap=1):
    max_v = sum(a_list[:30])
    last_v = max_v
    start = 0
    i = 1
    while i<len(a_list)-30:
        v = last_v - a_list[i-1] + a_list[i+29]
        last_v = v
        if max_v < v:
            start = i
            max_v = v
        i += 1
    end = start + 31 - gap
    while start > gap-1:
        to_stop = False
        for g in range(1, gap+1):
            #print( start, a_list[start], start-g, a_list[start-g], a_list[start]/a_list[start-g] )
            if a_list[start]/(a_list[start-g]+1) > 2.0:
                to_stop = True
                break
        if to_stop:
            break
        start -= 1
    while end < len(a_list):
        to_stop = False
        for g in range(1, min(gap+1,len(a_list)-end+1) ):
            if a_list[end-1]/(a_list[end-1+g]+1) > 2.0:
                to_stop = True
                break
        if to_stop:
            break
        end += 1
    return start, end

def read_cov():
    Cov = {}
    root='/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/'
    IN = open("%sDMSO_repX.depth" %(root,))
    
    last_tid = ""
    tmp_cov = []
    for line in IN:
        tid,pos,cov = line.strip().split()
        #print(tid,pos,cov)
        if last_tid!=tid:
            if last_tid!="":
                if len(tmp_cov) < 500:
                    Cov[last_tid] = tmp_cov
                #print(len(Cov))
                #if len(Cov)==20: return Cov
            last_tid = tid
            tmp_cov = [ int(cov) ]
        else:
            tmp_cov.append( int(cov) )
    Cov[tid] = tmp_cov
    return Cov

Cov = read_cov()

gap = {'miRNA': [],
        'tRNA': [],
        'snRNA': [],
        'snoRNA': [],
        'miscRNA': [],
        'mRNA':[],
        'others': []}

for tid in enriched_RNAs_2:
    if tid not in Cov:
        print(tid)
        continue
    a,b = peak_boundary(Cov[tid], gap=3)
    # gap.append( b-a+1 )
    if tid.startswith('miRNA'):
        gap['miRNA'].append(b-a+1)
    elif tid.startswith('tRNA'):
        gap['tRNA'].append(b-a+1)
    elif tid.startswith('snRNA'):
        gap['snRNA'].append(b-a+1)
    elif tid.startswith('snoRNA'):
        gap['snoRNA'].append(b-a+1)
    elif tid.startswith('miscRNA'):
        gap['miscRNA'].append(b-a+1)
    elif tid.startswith('mRNA'):
        gap['mRNA'].append(b-a+1)
    else:
        gap['others'].append(b-a+1)

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 6), sharey=True)
axs.set_title('Bind region')
axs.set_ylabel('Length of region')
data = [ gap['miRNA'], gap['tRNA'], gap['snRNA'], gap['snoRNA'], gap['miscRNA'], gap['mRNA'], gap['others']]
colors = [Colors.RGB['green'], Colors.RGB['red'], Colors.RGB['brown'], Colors.RGB['blue'], Colors.RGB['yellow'], Colors.RGB['pink'], Colors.RGB['gray']]
Figures.violinPlot(axs, data, ['miRNA','tRNA', 'snRNA', 'snoRNA', 'miscRNA', 'mRNA', 'others'],colors=colors,rem_ext=0.05)
# axs.legend()
fig.tight_layout()
outputDir = '/Share2/home/zhangqf7/jinsong_zhang/mirna_biogenesis/Results_rebutal'
fig.savefig(os.path.join(outputDir, 'length_dis.pdf'), dpi=120)
plt.show()
counts = list(map(len, data))
medians = list(map(np.median, data))
print(counts)
print(medians)

# print(np.median(np.array([gap['miRNA'], gap['tRNA'], gap['snRNA'], gap['snoRNA'], gap['miscRNA'], gap['others']]), axis=0))