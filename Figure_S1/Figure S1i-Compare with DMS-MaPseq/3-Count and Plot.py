
import numpy as np
import General
import scipy
import scipy.stats
import matplotlib.pyplot as plt 
from subprocess import getoutput
importCommon()

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
        #print(file)
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
    IN_PATH = IN_Folder.replace("6.shapemapper", "5.split_mutation")
    nai_1 = IN_PATH + "/REP1"
    nai_2 = IN_PATH + "/REP2"
    
    nai_1_lc = int(getoutput("wc -l "+nai_1).split()[0])
    nai_2_lc = int(getoutput("wc -l "+nai_2).split()[0])
    
    min_line = min(nai_1_lc, nai_2_lc)
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
###  1. 准备数据 - Balance data
###############################

sequence = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref/yeast.uniq.fa")

inFolder = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/6.shapemapper"
inFileName = "REP1"
REP1 = getFilesList(inFolder, inFileName); print(len(REP1))
inFileName = "REP2"
REP2 = getFilesList(inFolder, inFileName); print(len(REP2))

REP1_list = filter_filelist(REP1, mincov=1000, minnum=30, minratio=0.5); print(len(REP1_list))
REP2_list = filter_filelist(REP2, mincov=1000, minnum=30, minratio=0.5); print(len(REP2_list))

path_in_list = list(set(["/".join(it.split('/')[:-1]) for it in REP1_list]) & set(["/".join(it.split('/')[:-1]) for it in REP2_list]))

if not os.path.exists('/tmp/Rep_Cov'):
    os.mkdir("/tmp/Rep_Cov")

if not os.path.exists('/tmp/Rep_Cov/DMS'):
    os.mkdir("/tmp/Rep_Cov/DMS")

sample_all(path_in_list, "/tmp/Rep_Cov/DMS/", sequence)

###############################
###  2. 计算相关性
###############################

import numpy as np
import matplotlib.pyplot as plt

sequence = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref/yeast.uniq.fa")

inFolder = "/tmp/Rep_Cov/DMS/"
folders = os.listdir(inFolder)
random.shuffle(folders)

nai_mod_list = []

for folder in folders:
    seq = sequence[folder]
    n_mutnum_list_1, n_depth_list_1 = load_parsercounter(inFolder+"/"+folder+"/nai_1.count")
    n_mutnum_list_2, n_depth_list_2 = load_parsercounter(inFolder+"/"+folder+"/nai_2.count")
    for base,nn1,nd1,nn2,nd2 in zip(seq,n_mutnum_list_1,n_depth_list_1,n_mutnum_list_2,n_depth_list_2):
        
        ######   这里记得去掉TG碱基
        if nd1<100 or nd2<100 or 1.0*nn1/nd1>0.20 or 1.0*nn2/nd2>0.20 or base in ('T', 'G'):
            continue
        nai_mod_list.append((nd1,1.0*nn1/nd1,nd2,1.0*nn2/nd2))

print (len(nai_mod_list))
nai_mod_list.sort(key=lambda x: x[0])

nai_mod_list_rep1 = [ it[:2] for it in nai_mod_list ]
nai_mod_list_rep2 = [ it[2:] for it in nai_mod_list ]

corr_list = calc_slinding_correlation(nai_mod_list_rep1,nai_mod_list_rep2,wsize=50,wstep=10)

rl0 = [ it[0] for it in corr_list if 500<it[1] ] + [0]; print (len(rl0))
rl1 = [ it[0] for it in corr_list if 1000<it[1] ] + [0]; print (len(rl1))
rl2 = [ it[0] for it in corr_list if 2000<it[1] ] + [0]; print (len(rl2))
rl3 = [ it[0] for it in corr_list if 3000<it[1] ] + [0]; print (len(rl3))
rl4 = [ it[0] for it in corr_list if 4000<it[1] ] + [0]; print (len(rl4))
rl5 = [ it[0] for it in corr_list if 5000<it[1] ] + [0]; print (len(rl5))

cdf(rl0, color='#ff5722', topdown=False, label=">500", plotMedian=False)
cdf(rl1, color='#e91e63', topdown=False, label=">1000", plotMedian=False)
cdf(rl2, color='#4caf50', topdown=False, label=">2000", plotMedian=False)
cdf(rl3, color='#ff9800', topdown=False, label=">3000", plotMedian=False)
cdf(rl4, color='#00bcd4', topdown=False, label=">4000", plotMedian=False)
cdf(rl5, color='#9c27b0', topdown=False, label=">5000", plotMedian=False)

plt.xlim(0.80, 1.0)
plt.legend()
plt.savefig("figs/DMS-AC.pdf")
plt.show()

###############################
###  3. 使用Mismatch coverage
###############################

import numpy as np
import matplotlib.pyplot as plt

sequence = General.load_fasta("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref/yeast.uniq.fa")

inFolder = "/tmp/Rep_Cov/DMS/"
folders = os.listdir(inFolder)
random.shuffle(folders)

nai_mod_list = []

for folder in folders:
    seq = sequence[folder]
    n_mutnum_list_1, n_depth_list_1 = load_parsercounter(inFolder+"/"+folder+"/nai_1.count")
    n_mutnum_list_2, n_depth_list_2 = load_parsercounter(inFolder+"/"+folder+"/nai_2.count")
    for base,nn1,nd1,nn2,nd2 in zip(seq,n_mutnum_list_1,n_depth_list_1,n_mutnum_list_2,n_depth_list_2):
        
        ######   这里记得去掉TG碱基
        if nd1<100 or nd2<100 or 1.0*nn1/nd1>0.20 or 1.0*nn2/nd2>0.20 or base in ('T', 'G'):
            continue
        nai_mod_list.append(( nn1, 1.0*nn1/nd1, nn2, 1.0*nn2/nd2))

print (len(nai_mod_list))
nai_mod_list.sort(key=lambda x: x[0])

dl0 = [ it for it in nai_mod_list if 5<(it[0]+it[2])/2<15 ]; print (len(dl0))
dl1 = [ it for it in nai_mod_list if 15<(it[0]+it[2])/2<25 ]; print (len(dl1))
dl2 = [ it for it in nai_mod_list if 25<(it[0]+it[2])/2<35 ]; print (len(dl2))
dl3 = [ it for it in nai_mod_list if 35<(it[0]+it[2])/2<45 ]; print (len(dl3))
dl4 = [ it for it in nai_mod_list if 45<(it[0]+it[2])/2<55 ]; print (len(dl4))

nai_mod_list_rep1_0 = [ it[:2] for it in dl0 ]
nai_mod_list_rep2_0 = [ it[2:] for it in dl0 ]
nai_mod_list_rep1_1 = [ it[:2] for it in dl1 ]
nai_mod_list_rep2_1 = [ it[2:] for it in dl1 ]
nai_mod_list_rep1_2 = [ it[:2] for it in dl2 ]
nai_mod_list_rep2_2 = [ it[2:] for it in dl2 ]
nai_mod_list_rep1_3 = [ it[:2] for it in dl3 ]
nai_mod_list_rep2_3 = [ it[2:] for it in dl3 ]
nai_mod_list_rep1_4 = [ it[:2] for it in dl4 ]
nai_mod_list_rep2_4 = [ it[2:] for it in dl4 ]

corr_list_0 = calc_slinding_correlation(nai_mod_list_rep1_0,nai_mod_list_rep2_0,wsize=50,wstep=10)
corr_list_1 = calc_slinding_correlation(nai_mod_list_rep1_1,nai_mod_list_rep2_1,wsize=50,wstep=10)
corr_list_2 = calc_slinding_correlation(nai_mod_list_rep1_2,nai_mod_list_rep2_2,wsize=50,wstep=10)
corr_list_3 = calc_slinding_correlation(nai_mod_list_rep1_3,nai_mod_list_rep2_3,wsize=50,wstep=10)
corr_list_4 = calc_slinding_correlation(nai_mod_list_rep1_4,nai_mod_list_rep2_4,wsize=50,wstep=10)

rl0 = [ it[0] for it in corr_list_0 ] + [0]; print (len(rl0))
rl1 = [ it[0] for it in corr_list_1 ] + [0]; print (len(rl1))
rl2 = [ it[0] for it in corr_list_2 ] + [0]; print (len(rl2))
rl3 = [ it[0] for it in corr_list_3 ] + [0]; print (len(rl3))
rl4 = [ it[0] for it in corr_list_4 ] + [0]; print (len(rl4))

cdf(rl0, color='#ff5722', topdown=False, label=">10X", plotMedian=False)
cdf(rl1, color='#e91e63', topdown=False, label=">20X", plotMedian=False)
cdf(rl2, color='#4caf50', topdown=False, label=">30X", plotMedian=False)
cdf(rl3, color='#ff9800', topdown=False, label=">40X", plotMedian=False)
cdf(rl4, color='#00bcd4', topdown=False, label=">50X", plotMedian=False)

plt.xlim(0.40, 1.0)
plt.legend()
plt.savefig("figs/DMS-AC-modCov.pdf")
plt.show()





