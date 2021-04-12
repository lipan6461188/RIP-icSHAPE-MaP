import Visual

##########################
########   5S rRNA 2019-07-22
##########################

seq,dot = General.load_dot('/150T/zhangqf/GenomeAnnotation/Known_Structures/CRW/human_5S.dot')['human_5S']
human_5S = General.load_SHAPEMap("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/5.collect_mutation/rRNA_human_5S/NAI_100mm_vivo_SMR_SSII_repX.shape")
shape = human_5S['shape_pro_list']

General.calc_AUC_v2(dot, shape)
print(Visual.Plot_RNAStructure_Shape(seq, dot, shape))


# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_1c.txt'), 'w')
for d1,d2,d3 in zip(seq, dot, shape): print(f"{d1}\t{d2}\t{d3}", file=OUT)

OUT.close()

##########################
########   Step plot
##########################

seq,dot = General.load_dot('/150T/zhangqf/GenomeAnnotation/Known_Structures/CRW/human_5S.dot')['human_5S']
human_5S = General.load_SHAPEMap("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/5.collect_mutation/rRNA_human_5S/NAI_100mm_vivo_SMR_SSII_repX.shape", loadAll=True)
shape = human_5S['shape_pro_list']

nai_mut_ratio = [ mut/cov for mut,cov in zip(human_5S['mod_list'],human_5S['mod_cov_list']) ]
dmso_mut_ratio = [ mut/cov for mut,cov in zip(human_5S['dmso_list'],human_5S['dmso_cov_list']) ]

fig,ax = plt.subplots(1, figsize=(12,5))
for i in range(len(dot)):
    if dot[i] == '.':
        rect = plt.Rectangle([i-1,0.0],1,0.1,color=Colors.RGB['yellow'])
        ax.add_patch(rect)

ax.step(range(len(dmso_mut_ratio)), dmso_mut_ratio, c=Colors.RGB['indigo'])
ax.step(range(len(nai_mut_ratio)), nai_mut_ratio, c=Colors.RGB['red'])
fig.savefig(join(HOME, 'figs/step.pdf'))
fig.show()

# Source data
OUT = open(join(HOME, 'figs/source_data/Figure_1b.txt'), 'w')
for d1,d2,d3 in zip(seq, nai_mut_ratio, dmso_mut_ratio): print(f"{d1}\t{d2}\t{d3}", file=OUT)

OUT.close()



