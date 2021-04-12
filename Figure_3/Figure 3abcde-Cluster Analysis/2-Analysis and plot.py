
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap


################ Define PCA model

pca = PCA(n_components=2)
pca.fit(df)
d2_df = pca.transform(df)
d2_df = pd.DataFrame(d2_df, index=df.index, columns=['PC1','PC2'])


################################
######    Figure 3b
######    Plot the weight for each components
################################

comp1_variance = round(pca.explained_variance_ratio_[0] * 100, 1)
comp2_variance = round(pca.explained_variance_ratio_[1] * 100, 1)
sns.heatmap(pca.components_, vmax=0.25, vmin=-0.25, center=0, cmap=sns.diverging_palette(220, 20, n=20))
plt.xticks(range(0, 61, 10), range(-30, 31, 10))
plt.yticks([0, 1], [f'Component_1 ({comp1_variance}%)', f'Component_2 ({comp2_variance}%)'])
plt.savefig("figs/pca_kmeans.pdf")
plt.show()

# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_3b.txt'), 'w')
print( "\t".join([str(d) for d in pca.components_[0]]), file=OUT )
print( "\t".join([str(d) for d in pca.components_[1]]), file=OUT )
OUT.close()

################ Define KMeans model

kmeans = KMeans(n_clusters=3, verbose=False).fit(d2_df)

################ Classification

classes_old = kmeans.predict(d2_df)

classes = np.zeros_like(classes_old)
class_nums = [sum(classes_old==0), sum(classes_old==1), sum(classes_old==2)]
class_nums_order = np.argsort(class_nums)
classes[ classes_old == class_nums_order[0] ] = 0
classes[ classes_old == class_nums_order[1] ] = 2
classes[ classes_old == class_nums_order[2] ] = 1

color_map = {'miRNA':Colors.RGB['blue'], 'tRNA':Colors.RGB['red'], 'snoRNA':Colors.RGB['green'], \
            'mRNAExon':Colors.RGB['orange'], 'mRNAIntron':Colors.RGB['orange'], 'intergenic':Colors.RGB['orange'], 'miscRNA':Colors.RGB['orange']}
rnatype_colors = [color_map[it.split("_")[0]] for it in d2_df.index]
linewidths = []
for i,tid in enumerate(d2_df.index):
   if trans_dict[tid].enrich_score<0:
      rnatype_colors[i] += '50'
      linewidths.append(0)
   else:
      linewidths.append(0.5)


########## 原始类中的RNA比例

rna_class = [[],[],[]]
for i, tid in enumerate(d2_df.index):
    #if trans_dict[tid].enrich_score>=0:
    rna_class[classes[i]].append(trans_dict[tid])

show_rna_ratio([trans.rnatype for trans in rna_class[0]]) # Cluster I
show_rna_ratio([trans.rnatype for trans in rna_class[1]]) # Cluster II
show_rna_ratio([trans.rnatype for trans in rna_class[2]]) # Cluster III


################################
######    Figure 3a
######    Plot PCA/KMeans
################################

mesh_step_size = 0.01
x_min = -2.0
x_max = 2.0
y_min = -1.5
y_max = 2.0
plot_symbol_size = 50

xx, yy = np.meshgrid(np.arange(x_min, x_max, mesh_step_size), np.arange(y_min, y_max, mesh_step_size))
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)
plt.figure(figsize=(7,6))
cmap_light = ListedColormap([Colors.RGB['pale_green'], Colors.RGB['pale_blue'], Colors.RGB['pale_red']])
plt.pcolormesh(xx, yy, Z, cmap=cmap_light)
ax = plt.scatter(d2_df['PC1'], d2_df['PC2'], s=plot_symbol_size, c=rnatype_colors, edgecolor='black', linewidths=linewidths)
plt.xlim(xx.min(), xx.max())
plt.ylim(yy.min(), yy.max())
patch0 = mpatches.Patch(color=Colors.RGB['red'], label='tRNA')
patch1 = mpatches.Patch(color=Colors.RGB['blue'], label='pre-miRNA')
patch2 = mpatches.Patch(color=Colors.RGB['green'], label='snoRNA')
patch3 = mpatches.Patch(color=Colors.RGB['orange'], label='mRNA/intergenic/miscRNA')
plt.legend(handles=[patch0, patch1, patch2, patch3])
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title("PCA and K-means")
plt.savefig("figs/pca_kmeans.pdf")
plt.show()

# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_3a.txt'), 'w')

for i,tid in enumerate(d2_df.index):
  if classes[i] == 1: 
    cluster = "I"
  elif classes[i] == 2: 
    cluster = "II"
  else:
    cluster = "III"
  if trans_dict[tid].enrich_score<0:
      print( tid, d2_df.loc[tid, 'PC1'], d2_df.loc[tid, 'PC2'], cluster, "Not Enriched", file=OUT, sep="\t" )
  else:
      print( tid, d2_df.loc[tid, 'PC1'], d2_df.loc[tid, 'PC2'], cluster, "Enriched", file=OUT, sep="\t" )

OUT.close()


################################
######    Figure 3c
######    Plot SHAPE Tracks
################################


color_map = {'miRNA':'blue', 'tRNA':'red', 'snoRNA':'green', 'mRNAExon':'lightyellow', 'mRNAIntron':'lightyellow', 'intergenic':'lightyellow', 'miscRNA':'lightyellow'}

LEFT_ENRICH = True

rna_class = [[], [], []]
if LEFT_ENRICH:
   for classid,tid in zip(classes, df.index):
      if trans_dict[tid].enrich_score>=0:
         rna_class[ classid ].append( trans_dict[tid] )
   
   #rna_class.sort(key=lambda x:len(x), reverse=True)
else:
   for classid,tid in zip(classes, df.index):
      rna_class[ classid ].append( trans_dict[tid] )
   
   #rna_class.sort(key=lambda x: len(x))
   #rna_class[1],rna_class[2] = rna_class[2],rna_class[1]

rna_class[0].sort(key=lambda trans: trans.enrich_score, reverse=True) # Cluster III
rna_class[1].sort(key=lambda trans: trans.enrich_score, reverse=True) # Cluster I
rna_class[2].sort(key=lambda trans: trans.enrich_score, reverse=True) # Cluster II
print( len(rna_class[0]), len(rna_class[1]), len(rna_class[2]) )

show_rna_ratio([trans.rnatype for trans in rna_class[0]])
show_rna_ratio([trans.rnatype for trans in rna_class[1]])
show_rna_ratio([trans.rnatype for trans in rna_class[2]])

OUT = open("/tmp/color_cluster", 'w')
cc = {0:'red', 1:'blue', 2:'green'}
for classid in (0,1,2):
   class_colors = Colors.f("  ", bc=cc[classid])
   for trans in rna_class[classid]:
      colorBlock = Colors.f('  ', bc=color_map[trans.name.split('_')[0]])
      print( Colors.color_SHAPE(trans.stemloop_shape), class_colors, colorBlock, sep=" ", file=OUT)

OUT.close()


## Source data
OUT = open(join(HOME, 'figs/source_data/Figure_3c.txt'), 'w')
cluster_dict = { 0:'I', 1:'II', 2:'III' }

for classid in (0,1,2):
   cluster = cluster_dict[classid]
   for trans in rna_class[classid]:
      colorBlock = Colors.f('  ', bc=color_map[trans.name.split('_')[0]])
      print( trans.name, trans.rnatype, cluster, trans.stemloop_seq, *trans.stemloop_shape, sep="\t", file=OUT)

OUT.close()



##############################
############### Figure_S3B_Data_renamed.txt
############### 首先运行上面的代码，并且：LEFT_ENRICH = False
##############################

OUT = open("/tmp/cluster_data", 'w')
cc = {0:'red', 1:'blue', 2:'green'}
for classid in (0,1,2):
   #class_colors = Colors.f("  ", bc=cc[classid])
   for trans in rna_class[classid]:
        data = []
        sl_str = "{}-{}".format(trans.stem_loop[0], trans.stem_loop[3])
        data.append(str(classid))
        data.append(trans.name)
        data.append(trans.rnatype)
        data.append(sl_str)
        data.append(trans.enrich_score)
        data.append(trans.cleavage_score)
        data.append(round(trans.ago2,3))
        data.append(round(trans.ago3,3))
        data.append(trans.stemloop_seq)
        data.append(trans.stemloop_dot )
        data.append( ",".join( [ str(d) for d in trans.stemloop_shape ] ) )
        data.append(trans.full_seq)
        data.append( trans.full_dot )
        data.append( ",".join( [ str(d) for d in trans.full_shape ] ) )
        
        print( *data, sep="\t", file=OUT)

OUT.close()

nameConversion = DicerMiRNA.miRNA_nameConversion()

IN = open("/tmp/cluster_data")
OUT = open("/tmp/cluster_data_renamed", 'w')
for line in IN:
    data = line.strip().split()
    if data[1].startswith('miRNA'):
        raw_name = data[1].split('_')[1]
        miRNA_id = nameConversion[raw_name]
        data[1] = miRNA_id
    print("\t".join(data), file=OUT)

OUT.close()


################ Plot Enrichment and Cleavage

# Run with LEFT_ENRICH = False

score_df = []
score_df += [(trans.name, trans.enrich_score, trans.cleavage_score, trans.ago2, trans.ago3) for trans in rna_class[0]]
score_df += [(trans.name, trans.enrich_score, trans.cleavage_score, trans.ago2, trans.ago3) for trans in rna_class[1]]
score_df += [(trans.name, trans.enrich_score, trans.cleavage_score, trans.ago2, trans.ago3) for trans in rna_class[2]]

score_df = pd.DataFrame(score_df, columns=['tid', 'enrich', 'cleavage', 'ago2', 'ago3'])
score_df.index = score_df['tid']
score_df = score_df.drop(axis=1, labels=['tid'])

sns.heatmap(data=score_df, vmin=-2.0, vmax=2.0, center=0.0, cmap=sns.color_palette("coolwarm", 50))
plt.tight_layout()
plt.savefig("figs/heatmap.pdf")
plt.show()

################ Plot icSHAPE-Map score profile

def calc_ave_shape_profile(trans_list):
   profile = []
   for i in range(len(trans_list[0].stemloop_shape)):
      profile.append([])
   for trans in trans_list:
      #if trans.enrich_score<1: continue
      for i,shape in enumerate(trans.stemloop_shape):
         if shape!='NULL':
            profile[i].append(shape)
   for i in range(len(profile)):
      profile[i] = round(np.mean(profile[i]), 3)
   return profile

p1 = calc_ave_shape_profile(rna_class[0])
p2 = calc_ave_shape_profile(rna_class[1])
p3 = calc_ave_shape_profile(rna_class[2])

plt.plot(range(len(p1)),p1,color=Colors.RGB['red'],label="Cluster 1")
plt.plot(range(len(p2)),p2,color=Colors.RGB['green'],label="Cluster 2")
plt.plot(range(len(p3)),p3,color=Colors.RGB['blue'],label="Cluster 3")
plt.legend()
plt.savefig("figs/profile.pdf")
plt.show()

## Source data
OUT = open(join(HOME, 'figs/source_data/Figure_3d.txt'), 'w')

print( "I", *p1, sep="\t", file=OUT )
print( "II", *p2, sep="\t", file=OUT )
print( "III", *p3, sep="\t", file=OUT )

OUT.close()

##############################
############### Figure 3e first two figures: Plot Enrichment and Cleavage for 3 clusters
##############################


# 1->Cluster I; 2->Cluster II; 0->Cluster III
rna_class_all = [[], [], []]
for classid,tid in zip(classes, df.index):
   rna_class_all[ classid ].append( trans_dict[tid] )


enrich = [[],[],[]]
cleavage = [[],[],[]]
for clusterid in (0,1,2):
   for trans in rna_class_all[clusterid]:
      if trans.enrich_score>0:
        enrich[clusterid].append(trans.enrich_score)
        cleavage[clusterid].append(trans.cleavage_score)

print([ len(x) for x in enrich ])
print([ len(x) for x in cleavage ])

for i in (0,1,2):
   for j in range(i+1,3):
      p = scipy.stats.ttest_ind(enrich[i], enrich[j])[1]
      print("%d <==> %d: %s" % (i,j,p))

for i in (0,1,2):
   for j in range(i+1,3):
      p = scipy.stats.ttest_ind(cleavage[i], cleavage[j])[1]
      print("%d <==> %d: %s" % (i,j,p))


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
colors = [Colors.RGB['blue'], Colors.RGB['green'], Colors.RGB['red']]
Figures.violinPlot(ax, enrich, ['Cluster1','Cluster2','Cluster3'],colors=colors,rem_ext=0.01)
fig.tight_layout()
fig.savefig("figs/enrich.pdf")
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
colors = [Colors.RGB['blue'], Colors.RGB['green'], Colors.RGB['red']]
Figures.violinPlot(ax, cleavage, ['Cluster1','Cluster2','Cluster3'],colors=colors,rem_ext=0.01)
fig.tight_layout()
fig.savefig("figs/cleavage.pdf")
plt.show()



## Source data
OUT = open(join(HOME, 'figs/source_data/Figure_3e_left.txt'), 'w')
print('RNA\tCluster\tEnrichment score\tCleavage score', file=OUT)

cluster_dict = { 0:'I', 1:'II', 2:'III' }
for clusterid in (0,1,2):
   for trans in rna_class_all[clusterid]:
      if trans.enrich_score>0:
        print( trans.name, cluster_dict[clusterid], trans.enrich_score, trans.cleavage_score, file=OUT, sep='\t')

OUT.close()











