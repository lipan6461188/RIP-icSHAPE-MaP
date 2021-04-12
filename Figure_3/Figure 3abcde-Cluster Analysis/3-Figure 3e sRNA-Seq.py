
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


Foldchange_293TRISC_Vs_293T = []
for cid in (0,1,2):
    for tid in clusterNames[cid]:
        if tid not in table_data:
            print(tid)
            continue
        if table_data[tid][0]!='#VALUE!':
            Foldchange_293TRISC_Vs_293T.append( [ np.log2(float(table_data[tid][0])), cid ] )

Foldchange_293TRISC_Vs_293T = pd.DataFrame(Foldchange_293TRISC_Vs_293T, columns=['value','cid']); 
print("Foldchange_293TRISC_Vs_293T\n", Foldchange_293TRISC_Vs_293T.cid.value_counts())


sns.violinplot(data=Foldchange_293TRISC_Vs_293T, x='cid', y='value'); 
plt.ylim(-15, 10)
plt.ylabel("log2foldchange")
plt.title("Foldchange_293TRISC_Vs_293T")
plt.show()


df = Foldchange_293TRISC_Vs_293T
c1 = df.value[df.cid==0].tolist()
c2 = df.value[df.cid==1].tolist()
c3 = df.value[df.cid==2].tolist()
p12 = scipy.stats.ttest_ind(c1, c2)[1]
p13 = scipy.stats.ttest_ind(c1, c3)[1]
p23 = scipy.stats.ttest_ind(c2, c3)[1]
print("p12=%s;p13=%s;p23=%s" % (p12,p13,p23))



#### Source

OUT = open('/private/tmp/source_data/Figure_3e_right.txt', 'w')

cluster_dict = { 0: 'I', 1: 'II', 2: 'III' }
for cid in (0,1,2):
    for tid in clusterNames[cid]:
        if tid not in table_data:
            print(tid)
            continue
        if table_data[tid][0]!='#VALUE!':
            print( tid, cluster_dict[cid], np.log2(float(table_data[tid][0])), sep="\t", file=OUT)
            #Foldchange_293TRISC_Vs_293T.append( [ np.log2(float(table_data[tid][0])), cid ] )

OUT.close()






