

##### Run lpd lpdir envpy3 first

import subprocess
import threading
import os

def get_pairs():
    pairs = []
    pairs.append( [ 'NAI_100mm_vivo_CIRL_SSII',                     'DMSO_CIRL_SSII',                       'DC_CIRL_SSII' ] )
    pairs.append( [ 'NAI_100mm_vitro_CIRL_SSII',                    'DMSO_CIRL_SSII',                       'DC_CIRL_SSII' ] )
    pairs.append( [ 'NAI_100mm_vivo_SMR_SSII_repX',                 'DMSO_SMR_SSII_repX',                   'DC_SMR_SSII_repX' ] )
    pairs.append( [ 'NAI_100mm_vitro_SMR_SSII_repX',                'DMSO_SMR_SSII_repX',                   'DC_SMR_SSII_repX' ] )
    pairs.append( [ 'NAI_200mm_vitro_SMR_BENR_SSII',                'DMSO_SMR_BENR_SSII' ] )
    pairs.append( [ 'NAI_100mm_vitro_SMR_CENR_SSII',                'DMSO_SMR_BENR_SSII' ] )
    pairs.append( [ 'NAI_100mm_vitro_SMR_BENR_SSII',                'DMSO_SMR_BENR_SSII' ] )
    pairs.append( [ 'NAI_50mm_exvivo_dicer_1_CIRL_CENR_SSII_repX',  'DMSO_dicer_1_CIRL_BENR_SSII_repX',     'DC_dicer_CIRL_CENR_SSII_repX' ] )
    pairs.append( [ 'NAI_50mm_exvivo_dicer_2_CIRL_CENR_SSII_repX',  'DMSO_dicer_2_CIRL_BENR_SSII_repX',     'DC_dicer_CIRL_CENR_SSII_repX' ] )
    pairs.append( [ 'NAI_100mm_vivo_CIRL_CENR_SSII_repX',           'DMSO_CIRL_BENR_SSII',                  'DC_CIRL_CENR_SSII' ] )
    pairs.append( [ 'NAI_100mm_vivo_CIRL_TGIII',                    'DMSO_CIRL_TGIII',                      'DC_CIRL_TGIII' ] )
    pairs.append( [ 'NAI_100mm_vitro_CIRL_TGIII',                   'DMSO_CIRL_TGIII',                      'DC_CIRL_TGIII' ] )
    pairs.append( [ 'NAI_100mm_vitro_CIRL_CENR_SSII_repX',          'DMSO_CIRL_BENR_SSII',                  'DC_CIRL_CENR_SSII' ] )
    pairs.append( [ 'NAI_100mm_vitro_SMR_TGIII',                    'DMSO_SMR_TGIII' ] )
    pairs.append( [ 'NAI_200mm_vitro_SMR_TGIII',                    'DMSO_SMR_TGIII' ] )
    
    opt = {}
    for triplet in pairs:
        opt[ triplet[0] ] = { 'DMSO':triplet[1], 'DC':None }
        if len(triplet) == 3:
            opt[ triplet[0] ][ 'DC' ] = triplet[2]
    
    return opt

def format_pairs(file_list):
    opt = get_pairs()
    
    construct_pairs = []
    
    NAI_samples = [ file for file in file_list if file.startswith("NAI") ]
    for NAI_file in NAI_samples:
        necfiles = opt[NAI_file]
        if necfiles['DMSO'] in file_list:
            if necfiles['DC'] in file_list:
                construct_pairs.append( [necfiles['DMSO'], NAI_file, necfiles['DC']] )
            else:
                construct_pairs.append( [necfiles['DMSO'], NAI_file] )
    
    return construct_pairs

def load_fasta(seqFn, rem_tVersion=False):
    Fasta = {}
    cur_tid = ''
    
    for line in open(seqFn):
        if line[0] == '>':
            cur_tid = line[1:].split()[0]
            if rem_tVersion and '.' in cur_tid:
                cur_tid = ".".join(cur_tid.split(".")[:-1])
            Fasta[ cur_tid ] = ''
        else:
            Fasta[ cur_tid ] += line.strip()
    
    return Fasta

class MakeReacProfClass(threading.Thread):
    def __init__(self, inDir, outDir, seqFn, folderList, file_filter=None, mindepth=5000, maxbg=0.05, steps_to_run=[1,2,3]):
        threading.Thread.__init__(self)
        self.inDir = inDir
        self.outDir = outDir
        self.seqFn = seqFn
        self.folderList = folderList[:]
        self.file_filter = set(file_filter[:])
        self.mindepth = mindepth
        self.maxbg = maxbg
        self.steps_to_run = steps_to_run[:]
    
    def run(self):
        inDir = self.inDir.rstrip('/') + '/'
        outDir = self.outDir.rstrip('/') + '/'
        
        sequence = load_fasta(self.seqFn)
        
        cmd1 = "shapemapper_mutation_counter -i %s -s -c %s -n %s 2>/dev/null 1>/dev/null"
        cmd2 = "make_reactivity_profiles.py --fa %s --rna %s --counts %s --out %s --mindepth %s --maxbg %s 2>/dev/null 1>/dev/null"
        cmd3 = "normalize_profiles.py --disable-end3 --tonorm %s 1>/dev/null"
        #inputFolders = os.listdir(inDir)
        for folder_name in self.folderList:
            print("Run "+folder_name)
            rna_name = folder_name
            tLen = len(sequence[rna_name])
            oFolder = outDir + folder_name
            if not os.path.isdir(inDir+folder_name):
                continue
            if not os.path.exists(oFolder):
                os.mkdir(oFolder)
            inputFiles = os.listdir(inDir+folder_name)
            if self.file_filter:
                inputFiles = set(inputFiles) & self.file_filter
            
            if 1 in self.steps_to_run:
                for iFile in inputFiles:
                    inFile_full = inDir+folder_name+"/"+iFile
                    outFile_full = oFolder + "/" + iFile
                    
                    command = cmd1 % (inFile_full, outFile_full, tLen)
                    output = subprocess.getoutput(command)
            if 2 in self.steps_to_run or 3 in self.steps_to_run:
                con_pairs = format_pairs(inputFiles)
                for pair in con_pairs:
                    dmso_full = oFolder + "/" + pair[0]
                    nai_full = oFolder + "/" + pair[1]
                    #if len(pair) == 3:
                    #    dc_full = oFolder + "/" + pair[2]
                    #else:
                    #    dc_full = ""
                    dc_full = ""
                    outSHAPEfull = oFolder + "/" + pair[1] + ".shape"
                    
                    if 2 in self.steps_to_run:
                        command = cmd2 % (seqFn, rna_name, nai_full+" "+dmso_full+" "+dc_full, outSHAPEfull, self.mindepth, self.maxbg)
                        output = subprocess.getoutput(command)
                    
                    if 3 in self.steps_to_run:
                        command = cmd3 % (outSHAPEfull, )
                        output = subprocess.getoutput(command)
                        if "NormError" in output:
                            print("==>filter "+outSHAPEfull)
                            os.remove(outSHAPEfull)
                        else:
                            print("==>Success "+outSHAPEfull)

def batch_collect_mutations(inDir, outDir, seqFn, nproc=1, file_filter=None, mindepth=5000, maxbg=0.05, steps_to_run=[1,2,3]):
    import subprocess, math
    import _thread
    import time
    import random
    
    inputFolders = os.listdir(inDir)
    random.shuffle(inputFolders)
    if 'rRNA_human_5S' in inputFolders:
        inputFolders.remove('rRNA_human_5S')
        inputFolders.insert(0, 'rRNA_human_5S')
    if 'error' in inputFolders:
        inputFolders.remove('error')
    N_forEachProc = math.ceil(len(inputFolders) / nproc)
    print("Number for each process: "+str(N_forEachProc))
    
    Folders_Lists = []
    i = 0
    while i < len(inputFolders):
        Folders_Lists.append( inputFolders[i:i+N_forEachProc] )
        i += N_forEachProc
    
    thred_list = []
    for folder_list in Folders_Lists:
        thread = MakeReacProfClass(inDir, outDir, seqFn, folder_list, file_filter, mindepth, maxbg, steps_to_run)
        thred_list.append( thread )
    
    print("Number of threads list: "+str(len(thred_list)))
    for thread in thred_list:
        thread.start()
    
    for thread in thred_list:
        thread.join()

inDir = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/5.split_mutation"
outDir = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/6.shapemapper"
seqFn = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/test/dms-mapseq-rep/ref/yeast.uniq.fa"

file_filter = ['REP1', 'REP2']
batch_collect_mutations(inDir, outDir, seqFn, nproc=20, file_filter=file_filter, mindepth=1000, maxbg=0.05, steps_to_run=[1])


