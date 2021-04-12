
importCommon()
import DicerMiRNA

vivo_dir = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files"
vitro_dir = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vitro_SMR_SSII_repX/shape_files"

vivo_shape, vivo_seq = DicerMiRNA.read_shape(vivo_dir,min_valid_ratio=0.5)
vitro_shape, vitro_seq = DicerMiRNA.read_shape(vitro_dir,min_valid_ratio=0.5)

def clip_shape(shape, minV=0, maxV=1):
    for key in shape:
        shape_list = shape[key]
        for i in range(len(shape_list)):
            if shape_list[i]!='NULL':
                shape_list[i] = max(min(1, float(shape_list[i])), 0)
        shape[key] = shape_list
    return shape

vivo_shape = clip_shape(vivo_shape, minV=0, maxV=1)
vitro_shape = clip_shape(vitro_shape, minV=0, maxV=1)


#############################
#####  window识别的差异
#############################

def windows(diff_arr, winWidth=5):
    i = 3
    window_data = []
    while i<len(diff_arr):
        windata = diff_arr[i:i+winWidth]
        if np.count_nonzero(~np.isnan(windata))>winWidth//2:
            window_data.append(np.nanmean(windata))
        i += 3
    return window_data

def window_mean(diff_arr):
    window_data_list = windows(diff_arr)
    if len(window_data_list)>1:
        return np.mean(window_data_list)
    else:
        return 999

diff_data = []
for key in set(vivo_shape)&set(vitro_shape):
    shape1 = vivo_shape[key]
    shape2 = vitro_shape[key]
    shape1_arr = np.array([ np.nan if d=='NULL' else float(d) for d in shape1 ])
    shape2_arr = np.array([ np.nan if d=='NULL' else float(d) for d in shape2 ])
    diff_arr = np.abs(shape1_arr - shape2_arr)
    trans_diff_mean = np.nanmean(diff_arr)
    trans_diff_median = np.nanmedian(diff_arr)
    win_mean = window_mean(diff_arr)
    diff_data.append([ key,trans_diff_mean,trans_diff_median,np.nan if win_mean==999 else win_mean ])

diff_tabel = pd.DataFrame(diff_data, columns=['tid', 'trans_mean', 'trans_median', 'window_mean'])
diff_tabel = diff_tabel.sort_values(by='trans_mean', ascending=False)

color_map = {
    'lncRNAIntron': Colors.RGB['yellow'],
    'rRNA': Colors.RGB['blue'],
    'snoRNA': Colors.RGB['green'],
    'Y': Colors.RGB['brown'],
    'intergenic': Colors.RGB['yellow'],
    'tRNA': Colors.RGB['red'],
    'miRNA': Colors.RGB['khaki'],
    'snRNA' : Colors.RGB['amber']
}

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

RNA_types = [ d.split('_')[0] for d in diff_tabel.loc[:, 'tid'].tolist() ]
colors = [ color_map[d] for d in RNA_types ]
x = range(diff_tabel.shape[0])
y = [1 for _ in x]
#plt.bar(x, y, color=colors)

plt.figure(figsize=(12, 6))
plt.bar(x, diff_tabel.trans_mean, color=colors)

legend_elements = [Patch(facecolor=color_map[d], label=d) for d in color_map ]
plt.legend(handles=legend_elements, loc='upper right')
plt.xlabel("RNAs sorted by difference")
plt.ylabel("difference[ mean(abs(invivo - invitro)) ]")
plt.savefig(join(HOME, 'figs/VTD-bar.pdf'))
plt.close()


## Source

diff_tabel.loc[:, ['tid','trans_mean']].to_csv(join(HOME, 'figs/source_data/Figure_S2.txt'))


