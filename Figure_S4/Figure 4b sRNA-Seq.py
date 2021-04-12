
importCommon()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Figures

class Transcript():
    def __init__(self,cid,name,Type,enrich,cleavage,heatmapseq,heatmapdot,heatmapshape,fullseq,fulldot,fullshape):
        self.cid = int(cid)
        self.name = name
        self.Type = Type
        self.enrich = float(enrich)
        self.cleavage = float(cleavage)
        self.heatmapseq = heatmapseq
        self.heatmapdot = heatmapdot
        self.heatmapshape = [ float(i) if i!='NULL' else i for i in heatmapshape.split(',') ]
        self.fullseq = fullseq
        self.fulldot = fulldot
        self.fullshape = [ float(i) if i!='NULL' else i for i in fullshape.split(',') ]


### Read file
Trans = []
for line in open("Figure_S3B_Data_renamed.txt"):
    data = line.strip().split()
    if len(data)!=12:
        print(line,end="")
        print(len(data))
        break
    del data[3]
    trans = Transcript(*data)
    Trans.append(trans)


table = pd.read_csv("table.csv", sep=",")
table.head()

table_data = {}
for line in table.values:
    tid = line[0]
    i = tid.find('tRNA')
    if i!=-1:
        tid = tid[i:]
    table_data[tid] = (line[1], line[2], line[3])



leave_enrich = True
clusterNames = [[],[],[]]
clusterNames[0] = [ x.name for x in Trans if ((leave_enrich and x.enrich>0) or (not leave_enrich)) and (x.cid==0) ]
clusterNames[1] = [ x.name for x in Trans if ((leave_enrich and x.enrich>0) or (not leave_enrich)) and (x.cid==1) ]
clusterNames[2] = [ x.name for x in Trans if ((leave_enrich and x.enrich>0) or (not leave_enrich)) and (x.cid==2) ]


Foldchange_NoDice425_Vs_293T_total = []
Foldchange_NoDice425RISC_Vs_NoDice425 = []
for cid in (0,1,2):
    for tid in clusterNames[cid]:
        if tid not in table_data:
            print(tid)
            continue
        if table_data[tid][1]!='#VALUE!':
            Foldchange_NoDice425_Vs_293T_total.append( [ np.log2(float(table_data[tid][1])), cid ] )
        if table_data[tid][2]!='#VALUE!':
            Foldchange_NoDice425RISC_Vs_NoDice425.append( [ np.log2(float(table_data[tid][2])), cid ] )


Foldchange_NoDice425_Vs_293T_total = pd.DataFrame(Foldchange_NoDice425_Vs_293T_total, columns=['value','cid']); 
Foldchange_NoDice425RISC_Vs_NoDice425 = pd.DataFrame(Foldchange_NoDice425RISC_Vs_NoDice425, columns=['value','cid']); 
print("Foldchange_NoDice425_Vs_293T_total\n", Foldchange_NoDice425_Vs_293T_total.cid.value_counts())
print("Foldchange_NoDice425RISC_Vs_NoDice425\n", Foldchange_NoDice425RISC_Vs_NoDice425.cid.value_counts())


sns.violinplot(data=Foldchange_NoDice425_Vs_293T_total, x='cid', y='value'); 
plt.ylim(-15, 10)
plt.ylabel("log2foldchange")
plt.title("Foldchange_NoDice425_Vs_293T_total")
plt.show()

sns.violinplot(data=Foldchange_NoDice425RISC_Vs_NoDice425, x='cid', y='value'); 
plt.ylim(-15, 10)
plt.ylabel("log2foldchange")
plt.title("Foldchange_NoDice425RISC_Vs_NoDice425")
plt.show()


df = Foldchange_NoDice425_Vs_293T_total
c1 = df.value[df.cid==0].tolist()
c2 = df.value[df.cid==1].tolist()
c3 = df.value[df.cid==2].tolist()
p12 = scipy.stats.ttest_ind(c1, c2)[1]
p13 = scipy.stats.ttest_ind(c1, c3)[1]
p23 = scipy.stats.ttest_ind(c2, c3)[1]
print("p12=%s;p13=%s;p23=%s" % (p12,p13,p23))

df = Foldchange_NoDice425RISC_Vs_NoDice425
c1 = df.value[df.cid==0].tolist()
c2 = df.value[df.cid==1].tolist()
c3 = df.value[df.cid==2].tolist()
p12 = scipy.stats.ttest_ind(c1, c2)[1]
p13 = scipy.stats.ttest_ind(c1, c3)[1]
p23 = scipy.stats.ttest_ind(c2, c3)[1]
print("p12=%s;p13=%s;p23=%s" % (p12,p13,p23))


#### Source

OUT = open('/private/tmp/source_data/Figure_S4_left.txt', 'w')
cluster_dict = { 0: 'I', 1: 'II', 2: 'III' }
for cid in (0,1,2):
    for tid in clusterNames[cid]:
        if tid not in table_data:
            print(tid)
            continue
        if table_data[tid][1]!='#VALUE!':
            print( tid, cluster_dict[cid], np.log2(float(table_data[tid][1])), sep="\t", file=OUT)

OUT.close()


OUT = open('/private/tmp/source_data/Figure_S4_right.txt', 'w')
cluster_dict = { 0: 'I', 1: 'II', 2: 'III' }
for cid in (0,1,2):
    for tid in clusterNames[cid]:
        if tid not in table_data:
            print(tid)
            continue
        if table_data[tid][2]!='#VALUE!':
            print( tid, cluster_dict[cid], np.log2(float(table_data[tid][2])), sep="\t", file=OUT)

OUT.close()




