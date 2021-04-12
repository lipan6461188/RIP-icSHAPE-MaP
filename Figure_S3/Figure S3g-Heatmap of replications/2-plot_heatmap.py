
importCommon()

def list_dir_shape(folder):
    files = os.listdir(folder)
    files = [ it for it in files if it.endswith('.shape') ]
    return files

def calc_mutate_rate(shape):
    nai_rate_list = []
    dmso_rate_list = []
    for nai_num,nai_cov,dmso_num,dmso_cov in zip(shape['mod_list'],shape['mod_cov_list'],shape['dmso_list'],shape['dmso_cov_list']):
        if nai_cov > 1000:
            nai_rate = 1.0*nai_num/nai_cov
            assert nai_rate<=1.0
        else:
            nai_rate = 'NULL'
        if dmso_cov > 1000:
            dmso_rate = 1.0*dmso_num/dmso_cov
            assert dmso_rate<=1.0
        else:
            dmso_rate = 'NULL'
        nai_rate_list.append(nai_rate)
        dmso_rate_list.append(dmso_rate)
    return nai_rate_list, dmso_rate_list

def load_base_mutrate(Rep1, Rep2, Rep3):
    common_trans = set(Rep1) & set(Rep2) & set(Rep3)
    mutrate = []
    TransCount = 0
    
    for tid in common_trans:
        #print ("Loading "+tid+"....")
        shape_1 = General.load_SHAPEMap(ROOT1+"/"+tid,relocate=True,loadAll=True)
        shape_2 = General.load_SHAPEMap(ROOT2+"/"+tid,relocate=True,loadAll=True)
        shape_3 = General.load_SHAPEMap(ROOT3+"/"+tid,relocate=True,loadAll=True)
        nai_rate_1, dmso_rate_1 = calc_mutate_rate(shape_1)
        nai_rate_2, dmso_rate_2 = calc_mutate_rate(shape_2)
        nai_rate_3, dmso_rate_3 = calc_mutate_rate(shape_3)
        for n1,n2,n3,d1,d2 in zip(nai_rate_1,nai_rate_2,nai_rate_3,dmso_rate_1,dmso_rate_2):
            if 'NULL' not in (n1,n2,n3,d1,d2) and max([n1,n2,n3])<0.2 and max([d1,d2])<0.05:
                mutrate.append([n1,n2,n3,d1,d2])
    
    return mutrate

ROOT1="/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Replicates/Out_rep1/shape_files"
ROOT2="/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Replicates/Out_rep2/shape_files"
ROOT3="/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Replicates/Out_rep3/shape_files"

Rep1 = list_dir_shape(ROOT1)
Rep2 = list_dir_shape(ROOT2)
Rep3 = list_dir_shape(ROOT3)

mutrate_list = load_base_mutrate(Rep1, Rep2, Rep3); print(len(mutrate_list))

def calculate_correlation_matrix(score_list, title_list):
    df = General.init_pd_rect(len(score_list), len(score_list), title_list, title_list)
    for i in range(len(score_list)):
        df.iloc[i, i] = 1
        for j in range(i+1,len(score_list)):
            df.iloc[j, i] = df.iloc[i, j] = scipy.stats.pearsonr(score_list[i], score_list[j])[0]
    return df

score_list = [  [d[0] for d in mutrate_list], \
                [d[1] for d in mutrate_list], \
                [d[2] for d in mutrate_list], \
                [d[3] for d in mutrate_list], \
                [d[4] for d in mutrate_list] ]

title_list = ['nai_1','nai_2','nai_3','dmso_1','dmso_2']
df = calculate_correlation_matrix(score_list, title_list)

sns.heatmap(df, annot=True, fmt=".3f", cmap=sns.light_palette(('green'), n_colors=50))
plt.savefig("figs/Figure2_heatmap.pdf")
plt.show()

