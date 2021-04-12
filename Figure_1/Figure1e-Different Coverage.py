
import numpy as np
import General
import scipy
import scipy.stats
import matplotlib.pyplot as plt 
from subprocess import getoutput

assert getoutput("which shapemapper_mutation_counter") != ""

###############################
###  Functions
###############################

def load_parsercounter(inFn):
    mutnum_list = []
    depth_list = []
    
    IN = open(inFn)
    line = IN.readline()
    for line in IN:
        data = line.strip().split()
        mutnum = sum( [ int(it) for it in data[:26] ] )
        effdepth = int(data[27])
        mutnum_list.append(mutnum)
        depth_list.append(effdepth)
    
    return mutnum_list, depth_list

def getFilesList(inFolder, inFileName):
    full_file_path = []
    
    inFolder = inFolder.rstrip('/') + '/'
    middle_folders = os.listdir(inFolder)
    for middle_folder in middle_folders:        
        full_folder = inFolder + middle_folder
        if os.path.isdir(full_folder):
            full_folder += '/'
            filenames = os.listdir(full_folder)
            if inFileName in filenames:
                full_file_path.append( full_folder+inFileName )
    
    return full_file_path

def sample_lines(inFn, outFn, ratio, mode='w'):
    import random
    
    OUT = open(outFn, mode)
    for line in open(inFn):
        if random.random() < ratio:
            if line[-1] != '\n':
                line += '\n'
            OUT.writelines(line)
    
    OUT.close()

def filter_filelist(filelist, mincov=1000, minnum=50, minratio=0.5):
    filtered_filelist = []
    for file in filelist:
        mutnum_list, depth_list = load_parsercounter(file)
        depth = [it for it in depth_list if it>=mincov]
        if len(depth)>=minnum and 1.0*len(depth)/len(depth_list)>=minratio:
            filtered_filelist.append(file)
    return filtered_filelist

def sample_candidates(IN_Folder, OUT_Folder, sequence):
    """
        平衡数据
    """
    OUT_Folder = OUT_Folder.rstrip('/') + '/'
    
    tid = IN_Folder.split('/')[-1]
    IN_PATH = IN_Folder.replace("5.collect_mutation", "4.split_targets")
    dmso_1 = IN_PATH + "/DMSO_SMR_SSII_rep1"
    dmso_2 = IN_PATH + "/DMSO_SMR_SSII_rep2"
    nai_1 = IN_PATH + "/NAI_100mm_vivo_SMR_SSII_rep1"
    nai_2 = IN_PATH + "/NAI_100mm_vivo_SMR_SSII_rep2"
    
    dmso_1_lc = int(getoutput("wc -l "+dmso_1).split()[0])
    dmso_2_lc = int(getoutput("wc -l "+dmso_2).split()[0])
    nai_1_lc = int(getoutput("wc -l "+nai_1).split()[0])
    nai_2_lc = int(getoutput("wc -l "+nai_2).split()[0])
    dmso_lc = dmso_1_lc + dmso_2_lc
    
    min_line = min(dmso_lc, nai_1_lc, nai_2_lc)
    for code in ('nai_1', 'nai_2'):
        lc = eval(code+"_lc")
        in_file = eval(code)
        out_file = OUT_Folder + code
        count_file = out_file + ".count"
        if code == lc:
            os.system("cp %s %s" % (in_file, out_file))
        else:
            ratio = 1.0 * min_line / lc
            sample_lines(in_file, out_file, ratio)
        os.system("shapemapper_mutation_counter -i %s -c %s -s -n %s" % (out_file, count_file, len(sequence[tid])))
    
    if min_line == dmso_lc:
        sample_lines(dmso_1, OUT_Folder+"dmso", 1.0, mode='w')
        sample_lines(dmso_2, OUT_Folder+"dmso", 1.0, mode='a')
    else:
        sample_lines(dmso_1, OUT_Folder+"dmso", 1.0 * min_line / dmso_lc, mode='w')
        sample_lines(dmso_2, OUT_Folder+"dmso", 1.0 * min_line / dmso_lc, mode='a')
    
    os.system("shapemapper_mutation_counter -i %s -c %s -n %s" % (OUT_Folder+"dmso", OUT_Folder+"dmso.count", len(sequence[tid])))

def sample_all(inFolderList, OUT_DIR, sequence):
    for folder in inFolderList:
        tid = folder.split('/')[-1]
        print(tid)
        outfolder = OUT_DIR.rstrip('/')+'/'+ tid
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)
        sample_candidates(folder, outfolder, sequence)

def calc_slinding_correlation(list1,list2,wsize=100,wstep=10):
    """
    list1 and list2:
        [ (coverage, mut_ratio or mut_num),(coverage, mut_ratio or mut_num),(coverage, mut_ratio or mut_num)... ]
    
    Return:
        [ (correlation, ave_cov),(correlation, ave_cov),(correlation, ave_cov)... ]
    """
    import scipy
    import scipy.stats
    
    corr_list = []
    i = 0
    while i+wsize<len(list1):
        w1=list1[i:i+wsize]
        w2=list2[i:i+wsize]
        c1=[it[0] for it in w1]
        c2=[it[0] for it in w2]
        r1=[it[1] for it in w1]
        r2=[it[1] for it in w2]
        
        cor = scipy.stats.pearsonr(r1,r2)[0]
        cov1 = np.mean(c1)
        cov2 = np.mean(c2)
        corr_list.append( [ cor, (cov1+cov2)/2 ] )
        i += wstep
    
    return corr_list

def cdf(data_list, color='red', topdown=False, label=None, plotMedian=True):
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    
    data_list = np.sort(data_list)
    if topdown:
        p = 1 - 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    else:
        p = 1.0 * np.arange(len(data_list))/(len(data_list) - 1)
    plt.plot(data_list, p, color=color, label=label)
    if plotMedian:
        median_x = data_list[ len(data_list)/2 ]
        median_y = p[ len(p)/2 ]
        plt.plot([median_x], [median_y], 'bo')
        plt.axvline(x=median_x, ymin=0, ymax=1, linewidth=2, color='r')


###############################
###  1. 准备数据
###############################

sequence = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa")

inFolder = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/5.collect_mutation"
inFileName = "NAI_100mm_vivo_SMR_SSII_rep1"
nai_1_list = getFilesList(inFolder, inFileName); print(len(nai_1_list))
inFileName = "NAI_100mm_vivo_SMR_SSII_rep2"
nai_2_list = getFilesList(inFolder, inFileName); print(len(nai_2_list))
inFileName = "DMSO_SMR_SSII_rep1"
dmso_1_list = getFilesList(inFolder, inFileName); print(len(dmso_1_list))
inFileName = "DMSO_SMR_SSII_rep2"
dmso_2_list = getFilesList(inFolder, inFileName); print(len(dmso_2_list))

cand_nai_1_list = filter_filelist(nai_1_list, mincov=1000, minnum=30, minratio=0.5); print(len(cand_nai_1_list))
cand_nai_2_list = filter_filelist(nai_2_list, mincov=1000, minnum=30, minratio=0.5); print(len(cand_nai_2_list))
cand_dmso_1_list = filter_filelist(dmso_1_list, mincov=1000, minnum=30, minratio=0.5); print(len(cand_dmso_1_list))
cand_dmso_2_list = filter_filelist(dmso_2_list, mincov=1000, minnum=30, minratio=0.5); print(len(cand_dmso_2_list))

path_in_list = list(set(["/".join(it.split('/')[:-1]) for it in cand_nai_1_list]) & set(["/".join(it.split('/')[:-1]) for it in cand_nai_2_list]) & 
                    set(["/".join(it.split('/')[:-1]) for it in cand_dmso_1_list]) & set(["/".join(it.split('/')[:-1]) for it in cand_dmso_2_list]))

if not os.path.exists('/tmp/Rep_Cov'):
    os.mkdir("/tmp/Rep_Cov")

if not os.path.exists('/tmp/Rep_Cov/NAI'):
    os.mkdir("/tmp/Rep_Cov/NAI")

sample_all(path_in_list, "/tmp/Rep_Cov/NAI/", sequence)

###############################
###  2. 计算相关性
###############################

inFolder = "/tmp/Rep_Cov/NAI/"
folders = os.listdir(inFolder)
random.shuffle(folders)

sequence = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa")

nai_mod_list = []

for folder in folders:
    seq = sequence[folder]
    d_mutnum_list, d_depth_list = load_parsercounter(inFolder+"/"+folder+"/dmso.count")
    n_mutnum_list_1, n_depth_list_1 = load_parsercounter(inFolder+"/"+folder+"/nai_1.count")
    n_mutnum_list_2, n_depth_list_2 = load_parsercounter(inFolder+"/"+folder+"/nai_2.count")
    for base,dn,dd,nn1,nd1,nn2,nd2 in zip(seq,d_mutnum_list,d_depth_list,n_mutnum_list_1,n_depth_list_1,n_mutnum_list_2,n_depth_list_2):
        if dd<100 or nd1<100 or nd2<100 or 1.0*dn/dd>0.05:
            continue
        nai_mod_list.append((nd1,1.0*nn1/nd1,nd2,1.0*nn2/nd2))

print(len(nai_mod_list))
nai_mod_list.sort(key=lambda x: x[0])

nai_mod_list_rep1 = [ it[:2] for it in nai_mod_list ]
nai_mod_list_rep2 = [ it[2:] for it in nai_mod_list ]

corr_list = calc_slinding_correlation(nai_mod_list_rep1,nai_mod_list_rep2,wsize=50,wstep=10)

rl0 = [ it[0] for it in corr_list if 500<it[1]  ] + [0]; print(len(rl0))
rl1 = [ it[0] for it in corr_list if 1000<it[1] ] + [0]; print(len(rl1))
rl2 = [ it[0] for it in corr_list if 2000<it[1] ] + [0]; print(len(rl2))
rl3 = [ it[0] for it in corr_list if 3000<it[1] ] + [0]; print(len(rl3))
rl4 = [ it[0] for it in corr_list if 4000<it[1] ] + [0]; print(len(rl4))
rl5 = [ it[0] for it in corr_list if 5000<it[1] ] + [0]; print(len(rl5))

cdf(rl0, color='#ff5722', topdown=False, label=">500",  plotMedian=False)
cdf(rl1, color='#e91e63', topdown=False, label=">1000", plotMedian=False)
cdf(rl2, color='#4caf50', topdown=False, label=">2000", plotMedian=False)
cdf(rl3, color='#ff9800', topdown=False, label=">3000", plotMedian=False)
cdf(rl4, color='#00bcd4', topdown=False, label=">4000", plotMedian=False)
cdf(rl5, color='#9c27b0', topdown=False, label=">5000", plotMedian=False)

plt.xlim(0.80, 1.0)
plt.legend()
plt.savefig("/Share2/home/zhangqf7/figs/NAI-AC.pdf")
plt.close()

