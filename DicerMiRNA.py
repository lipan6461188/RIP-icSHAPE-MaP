
import General, Colors
import Structure
import seaborn  as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import scipy.stats
import os, sys, re, math

assert sys.version[0] == '3'

################################
#### read SHAPE
################################

def read_shape(inFolder, postfix=".shape", min_valid_ratio=0, min_valid_base_num=0, relocate=False, loadAll=False):
    """
    Return SHAPE, Sequence
    """
    SHAPE = {}; Sequence = {}
    import os
    files = os.listdir(inFolder)
    for file in files:
        if file.endswith(postfix):
            tid = file.rstrip(postfix)
            try:
                data = General.load_SHAPEMap(os.path.join(inFolder,file), relocate=relocate, loadAll=loadAll)
            except ValueError:
                print(f"Warning: {tid} cannot be read: {os.path.join(inFolder,file)}")
                continue
            shape = data['shape_pro_list']
            if len(shape)==0: continue
            seq = data['seq']
            if 1-shape.count('NULL')/len(shape)<min_valid_ratio:
                continue
            if len(shape)-shape.count('NULL')<min_valid_base_num:
                continue
            SHAPE[tid] = shape
            Sequence[tid] = seq
    return SHAPE, Sequence

################################
#### Format
################################

def period_number(Len):
    num_str = ""
    for i in range(1, Len+1):
        num_str += str(i%10)
    return num_str

def period_tick(Len):
    tick_str = ""
    for i in range(1,Len+1):
        if i%10 == 0:
            tick_str += "|"
        elif i%5 == 0:
            tick_str += "*"
        else:
            tick_str += " "
    return tick_str

################################
#### read PDBs
################################

def read_pdb(pdbFn, mcsym=False):
    Seq = ""
    Atom_3d = {}
    Model = {}
    mod_id = 1
    for line in open(pdbFn):
        if line.startswith('MODEL'):
            if Model:
                Atom_3d[mod_id] = Model
                Seq = ""
                Model = {}
                mod_id += 1
                continue
        Record = line[:6].strip()
        Atom = line[12:16].strip()
        Char = line[17:20].strip()
        resSeq = line[22:26].strip()
        X = line[30:38].strip()
        Y = line[38:46].strip()
        Z = line[46:54].strip()
        if Record=='ATOM':# and Atom=='O3\'':
            if mcsym:
                if Atom == r'O3*':
                    Seq += Char
            else:
                if Atom == 'O3\'':
                    Seq += Char
            resSeq = int(resSeq)
            try:
                Model[resSeq][Atom] = (float(X),float(Y),float(Z))
            except KeyError:
                Model[resSeq] = {}
                Model[resSeq][Atom] = (float(X),float(Y),float(Z))
    
    if Model:
        Atom_3d[mod_id] = Model
    
    return Seq, Atom_3d

def distance(Pos1, Pos2):
    return round(math.sqrt((Pos1[0]-Pos2[0])**2+(Pos1[1]-Pos2[1])**2+(Pos1[2]-Pos2[2])**2),2)

def calculate_PDBdistance_matrix(seq, PDBModel, mcsym=False):
    start = 1
    end = len(PDBModel)
    
    df = []
    for i in range(end):
        df.append([0]*end)
    
    for i in range(1,end+1):
        for j in range(i+1,end+1):
            if mcsym:
                df[i-1][j-1] = distance(PDBModel[i]["P"], PDBModel[j]["O3*"])
            else:
                df[i-1][j-1] = distance(PDBModel[i]["C5'"], PDBModel[j]["O3'"])
    
    index_ticksLabels = [ str(index+1)+base+str(len(seq)-index) for index,base in enumerate(seq) ]
    column_ticksLabels = [ str(index+1)+base for index,base in enumerate(seq) ]
    df = pd.DataFrame(df, columns=column_ticksLabels, index=index_ticksLabels)
    return df

def PlotMatrix(df, vmin=50, vmax=70, center=60, figsize=(35,25)):
    """
    seq, pdb = read_pdb("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/Figure4c-example/Rosetta/constraint/miR-16-1-P10DU/final.out.23.pdb")
    df = calculate_PDBdistance_matrix(seq, pdb[1])
    PlotMatrix(df, vmin=50, vmax=70, center=60)
    plt.savefig("figs/haha.pdf")
    plt.close()
    """
    fig,ax = plt.subplots(1,1,figsize=figsize) #plt.figure(figsize=figsize)
    sns.heatmap(data=df, ax=ax, vmin=vmin, vmax=vmax, center=center, annot=True, fmt=".1f", cmap='YlGnBu', cbar=False)
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    return fig, [ax]
    #plt.savefig("figs/haha.pdf")
    #plt.close()

def PlotTwoSmallMatrix(df, vmin=50, vmax=70, center=60, figsize=(8,6)):
    Len = df.shape[0]
    df_left = df.iloc[0:5, 18:26]
    df_right = df.iloc[Len-26:Len-18, Len-5:Len].T
    
    fig, axs = plt.subplots(2,1,figsize=figsize)
    sns.heatmap(data=df_left, vmin=vmin, vmax=vmax, center=center, annot=True, fmt=".1f", cmap=sns.color_palette("RdBu_r", 50), cbar=False, ax=axs[0])
    axs[0].xaxis.tick_top()
    axs[0].yaxis.tick_right()
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=90)
    axs[0].set_yticklabels(axs[0].get_yticklabels(), rotation=0)
    axs[0].xaxis.tick_top()
    axs[0].yaxis.tick_right()
    
    sns.heatmap(data=df_right, vmin=vmin, vmax=vmax, center=center, annot=True, fmt=".1f", cmap=sns.color_palette("RdBu_r", 50), cbar=False, ax=axs[1])
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=90)
    axs[1].set_yticklabels(axs[1].get_yticklabels(), rotation=0)
    axs[1].xaxis.tick_top()
    axs[1].yaxis.tick_right()
    
    return fig, axs

def locate_rosetta_output(inFolder):
    constraint = {}
    
    inFolder = inFolder.rstrip('/')
    folders = os.listdir(inFolder)
    for folder in folders:
        if os.path.isdir(inFolder+'/'+folder):
            files = os.listdir(inFolder+'/'+folder)
            constraint[folder] = []
            for file in files:
                if re.match(r'.*\.out\.\d+\.pdb', file):
                    full_path = inFolder+'/'+folder+'/'+file
                    constraint[folder].append(full_path)
    
    return constraint

def locate_mcsym_output(inFolder):
    constraint = {}
    
    inFolder = inFolder.rstrip('/')
    folders = os.listdir(inFolder)
    for folder in folders:
        if os.path.isdir(inFolder+'/'+folder):
            files = os.listdir(inFolder+'/'+folder)
            constraint[folder] = []
            for file in files:
                if re.match(r'.*\.pdb', file):
                    full_path = inFolder+'/'+folder+'/'+file
                    constraint[folder].append(full_path)
    
    return constraint


def Rosseta_check_str(PDB_file_list, DotDict):
    for tid in PDB_file_list:
        if not isinstance(DotDict[tid], str):
            raise RuntimeError(r"DotDict should be {tid:dot1, tid2:dot2...}")
        
        print("Check",tid,sep=" ",end="\t")
        if tid not in DotDict:
            print("Failed. No 2nd structure in DotDict")
            continue
        if len(PDB_file_list[tid]) == 0:
            print("Failed. No PDB files in PDB_file_list")
            continue
        path = "/".join(PDB_file_list[tid][0].split('/')[:-1])
        dot = open(path+"/dot.str").readline().strip()
        if dot!=DotDict[tid]:
            print("Failed. Different structure")
            return
        print("Success.")

def MCSymLocal_check_str(PDB_file_list, DotDict):
    for tid in PDB_file_list:
        if not isinstance(DotDict[tid], str):
            raise RuntimeError(r"DotDict should be {tid:dot1, tid2:dot2...}")
        
        print("Check",tid,sep=" ",end="\t")
        if tid not in DotDict:
            print("Failed. No 2nd structure in DotDict")
            continue
        if len(PDB_file_list[tid]) == 0:
            print("Failed. No PDB files in PDB_file_list")
            continue
        path = "/".join(PDB_file_list[tid][0].split('/')[:-1])
        contents = open(path+"/"+tid+".mcc").readlines()
        i = 0
        found = False
        while i<len(contents):
            if contents[i].startswith("sequence("):
                if DotDict[tid]!=contents[i+1][2:].strip():
                    print("Failed. Different structure")
                    return
                found = True
            i += 1
        if not found:
            print("Failed. Wrong mcc file")
            return
        print("Success.")

def locate_mature(annot_seq):
    """
    annot_seq               -- Raw annotated seq
    
        such as: GUGGUUAUCCCUGUCCUGUUCGuuuugcucaugucgAAUCGUACAGGGUCAUCCACUU
    """
    lindex = re.search("[aucg]", annot_seq).start()
    rindex = re.search("[aucg][AUCG]*$", annot_seq).start()
    return lindex, rindex

def construct_distance_df(location_table, Mature_Annot_miRNA_dot):
    dist_df = []
    
    for miRID in set(location_table)&set(Mature_Annot_miRNA_dot):
        print("load "+miRID+'...')
        Mature_seq = Mature_Annot_miRNA_dot[miRID][0]
        lindex, rindex = locate_mature(Mature_seq)
        lfrag = Mature_seq[:lindex]
        rfrag = Mature_seq[rindex+1:]
        llen = len(lfrag)  # left length
        rlen = len(rfrag)  # right length
        mlen = rindex-lindex+1  # middle length
        assert llen + rlen + mlen == len(Mature_seq)
        
        rep = 1
        for file in location_table[miRID]:
            Seq, Atom_3ds = read_pdb(file)
            if Seq.upper() != Mature_seq.upper():
                print(Seq.upper()+"\t"+Mature_seq.upper())
                raise RuntimeError("Different sequence")
            Atom_3d = Atom_3ds[1]
            spatial_5p = distance(Atom_3d[1]["C5'"], Atom_3d[lindex]["O3'"])
            spatial_3p = distance(Atom_3d[rindex+2]["C5'"], Atom_3d[len(Mature_seq)]["O3'"])
            dist_df.append( [miRID,rep,llen,rlen,spatial_5p,spatial_3p] )
            rep += 1
    
    dist_df = pd.DataFrame(dist_df, columns=['miRID','rep','llen','rlen','spatial_5p','spatial_3p'])
    return dist_df

def read_auc(inFn):
    data_list = [ (line.strip().split()[0],float(line.strip().split()[1])) for line in open(inFn).readlines() ]
    data_list.sort(key=lambda x: x[1])
    return data_list

def read_pdb_score(pdb_fn):
    for line in open(pdb_fn):
        if line.startswith('silent_score '):
            return float(line.strip().split()[1])
    return None

def batch_readPDB(PDB_file_list, mcsym=False, return_score=False):
    last_seq = ""
    matrix_list = []
    score_list = []
    for pdbFn in PDB_file_list:
        seq, Atom_3d = read_pdb(pdbFn, mcsym=mcsym)
        if last_seq == "":
            last_seq = seq
        else:
            assert last_seq == seq
        matrix = calculate_PDBdistance_matrix(seq, Atom_3d[1], mcsym=mcsym)
        matrix_list.append(matrix)
        score_list.append( read_pdb_score(pdbFn) )
    if return_score:
        return seq, matrix_list, score_list
    else:
        return seq, matrix_list

def matrix_list_statics(matrix_list, func=np.median):
    dim = matrix_list[0].shape[0]
    count_list = General.init_list_rect(dim,dim,0)
    
    for i in range(dim):
        for j in range(dim):
            value_list = []
            for matrix in matrix_list:
                value_list.append( matrix.iloc[i,j] )
            count_list[i][j] = func(value_list)
    
    count_df = pd.DataFrame(count_list, columns=matrix_list[0].columns, index=matrix_list[0].index)
    return count_df

"""
Rosseta_dir = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/Figure4c-example/Rosetta/constraint"
location_table = locate_rosetta_output(Rosseta_dir)
seq, matrix_list = batch_readPDB(PDB_file_list)
count_df = matrix_list_statics(matrix_list, func=np.median)
PlotMatrix(count_df, vmin=50, vmax=70, center=60)
plt.savefig("figs/haha.pdf")
plt.close()
"""

def miRNA_nameConversion():
    nameConversion = {}
    
    miRNA_ref_fa = "/Share/home/zhangqf7/lipan/reference/miRNA_tRNA_snRNA_snoRNA_rRNA_mtRNA_YRNA/miRNA/hsa.fa"
    for line in open(miRNA_ref_fa):
        if line[0] == '>':
            raw_id,new_id = re.findall("Alias=(\\w+);Name=([\\w\\-]+)\\(", line)[0]
            nameConversion[raw_id] = new_id
    
    return nameConversion


################################
###  Load data
################################

def trim_MiRBase_miRNA(MiRBase_miRNA, rePred=False):
    """
    MiRBase_miRNA           -- MiRBase miRNA dot dictionary
    rePred                  -- Re-predict structure
    """
    import re
    
    trimed_miRNA = {}
    for miRNAID in MiRBase_miRNA:
        seq, dot = MiRBase_miRNA[miRNAID]
        seq = seq.replace('T','U').replace('t','u')
        lindex = re.search("[AUCG]", seq).start()
        rindex = re.search("[AUCG][aucg]*$", seq).start()
        mature_miRNA = seq[lindex:rindex+1]
        if re.search("[aucg]", mature_miRNA) == None:
            continue
        
        if rePred:
            mature_dot = Structure.predict_structure(mature_miRNA.upper())
            if mature_dot.find(")") < mature_dot.rfind("("): continue
        else:
            mature_dot = list(dot[lindex:rindex+1])
            
            ## Trim dot
            if mature_dot.count('(') > mature_dot.count(')'):
                diff = mature_dot.count('(') - mature_dot.count(')')
                i = 0
                while diff > 0:
                    if mature_dot[i] == '(':
                        mature_dot[i] = '.'
                        diff -= 1
                    i += 1
            elif mature_dot.count(')') > mature_dot.count('('):
                diff = mature_dot.count(')') - mature_dot.count('(')
                i = len(mature_dot)-1
                while diff > 0:
                    if mature_dot[i] == ')':
                        mature_dot[i] = '.'
                        diff -= 1
                    i -= 1
            else:
                pass
            
            mature_dot = "".join(mature_dot)
            if mature_dot.count('(')==0 or mature_dot.count(')')==0:
                continue
        
        trimed_miRNA[miRNAID] = (mature_miRNA, mature_dot)
    
    return trimed_miRNA

miRBase = None
trimed_miRBase = None
shape_vivo = None
shape_vitro = None
shape_RIP = None
sequence = None

def load_data():
    global miRBase
    global trimed_miRBase
    global shape_vivo
    global shape_vitro
    global shape_RIP
    global sequence
    
    print(Colors.f("Start to load SHAPE-MaP data...", 'yellow'))
    dir_ = '/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files'
    shape_vivo, _ = read_shape(dir_, ".shape", relocate=True)
    print(Colors.f("shape_vivo loaded", 'yellow'))
    dir_ = '/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vitro_SMR_SSII_repX/shape_files'
    shape_vitro, _ = read_shape(dir_, ".shape", relocate=True)
    print(Colors.f("shape_vitro loaded", 'yellow'))
    dir_ = '/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/RIP_NAIN3_repX_20190514/shape_files'
    shape_RIP, _ = read_shape(dir_, ".shape", relocate=True)
    print(Colors.f("shape_RIP loaded", 'yellow'))
    
    sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vitro_SMR_SSII_repX/index/Genome.fa')
    print(Colors.f("sequence loaded", 'yellow'))
    
    miRBase = General.load_dot("/150T/zhangqf/GenomeAnnotation/miRNA/more/miRNA.dot")
    miRBase = {k:miRBase[k] for k in miRBase if k.startswith('hsa')}
    trimed_miRBase = trim_MiRBase_miRNA(miRBase, rePred=False)


