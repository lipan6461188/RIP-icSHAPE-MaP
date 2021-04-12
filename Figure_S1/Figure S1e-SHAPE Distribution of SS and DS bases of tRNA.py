
importCommon()
import DicerMiRNA


#########################
## 1. Read shape and sequence
#########################

shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files_1000"
#shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files"
#shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files"
SHAPE, Sequence = DicerMiRNA.read_shape(shape_folder, postfix='.shape', min_valid_ratio=0.8, min_valid_base_num=0, relocate=False, loadAll=False); print (len(SHAPE))
#tRNA_structurome = build_tRNA_structurome("/150T/zhangqf/GenomeAnnotation/tRNA/GtRNAdb/hg19/human_tRNA.dot")

tRNA_dot, tRNA_annot = General.load_dot("/150T/zhangqf/GenomeAnnotation/tRNA/GtRNAdb/hg19/human_tRNA-2.dot", rem_tVersion=False, load_annotation=True)

score = {}
HMM_score = {}
Structure_score = {}
for tid in tRNA_annot:
    pattern = r'score=(-?[\d\.]+);HMM_score=(-?[\d\.]+);Structure_score=(-?[\d\.]+)'
    s1,s2,s3 = re.findall(pattern, tRNA_annot[tid])[0]
    score[tid],HMM_score[tid],Structure_score[tid] = float(s1),float(s2),float(s3)

highly_confidence_tRNAs = [k for k,v in tRNA_annot.items() if "possible_pseudogene=0;possible_intron=0" in v]

print(f"Number: {len(highly_confidence_tRNAs)}")


########################
## 2. AUC 0.7以上的单双链分布
########################

ds_values = []
ss_values = []
for tid in Sequence:
    if tid.startswith('tRNA_'):
        dot_tid = "tRNA-"+"_".join(tid.split('_')[1:])
        if dot_tid in tRNA_dot and "possible_pseudogene=0;possible_intron=0" in tRNA_annot[dot_tid]:
            tRNA_seq, dot = tRNA_dot[dot_tid]
            assert tRNA_seq == Sequence[tid].replace('U', 'T')
            shape = SHAPE[tid]
            if 1-shape.count('NULL')/len(shape)<0.8:
                continue
            auc = General.calc_AUC_v2(dot, shape)
            if auc>0.70:
                for d,s in zip(dot, shape):
                    if s!='NULL':
                        if d=='.':
                            ss_values.append(np.clip(float(s),0,1))
                        else:
                            ds_values.append(np.clip(float(s),0,1))


print(np.mean(ss_values), np.mean(ds_values))
_, p = scipy.stats.ttest_ind(ss_values, ds_values)


########  Violin plot

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4,5), sharey=True)
axs.set_title(f"tRNA (AUC>0.7)\nn={len(ss_values)}, {len(ds_values)}\np={p:.5f}")
axs.set_ylabel('Reactivity')
data = [ss_values, ds_values]
colors = [ Colors.RGB['red'], Colors.RGB['green'] ]
Figures.violinPlot(axs, data, ['Single-stranded','Double-stranded'],colors=colors)
fig.tight_layout()
plt.savefig(join(HOME, 'figs/tRNA_violinplot.pdf'))
fig.show()


########################
## 3. Double, Single-stranded percentage
########################


stackedBars = [ ]
barLabels = [ ]
for d_min in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
    if d_min == 0.9:
        ss_num = len([ d for d in ss_values if d_min<=d ])
        ds_num = len([ d for d in ds_values if d_min<=d ])
    else:
        ss_num = len([ d for d in ss_values if d_min<=d<d_min+0.1 ])
        ds_num = len([ d for d in ds_values if d_min<=d<d_min+0.1 ])
    total = ss_num + ds_num
    stackedBars.append( [ss_num/total, ds_num/total] )
    barLabels.append(f"{d_min:.1f}-{(d_min+0.1):.1f} \n(n={total})")

stackedLabels = ['Single-stranded', 'Double-stranded']
stackedColors = [Colors.RGB['blue'], Colors.RGB['green']]

Figures.stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors)
plt.tight_layout()
plt.savefig(os.path.join(HOME, 'figs/ratio.pdf'))
plt.show()



















