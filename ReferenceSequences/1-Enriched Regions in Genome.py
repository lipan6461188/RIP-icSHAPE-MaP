
importCommon()



################################
#####  1. 从bam文件中取出Peaks
################################


DMSO_BAM = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/DMSO_repX.bam"
NAIN3_BAM = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/NAIN3_repX.bam"


def Bam_cov(bamFn):
    
    temp_pos = "samtools view -bh -F 0x10 "+bamFn
    temp_neg = "samtools view -bh -f 0x10 "+bamFn
    temp_depth = " | samtools depth -m 9999999 -"
    
    print("Start to load positive strand...")
    Peak_regions_pos = GetCov(temp_pos+temp_depth)
    print("Start to load negative strand...")
    Peak_regions_neg = GetCov(temp_neg+temp_depth)
    
    return Peak_regions_pos, Peak_regions_neg

def GetCov(CMD, minCov=1500, minPeakWidth=30):
    #print(CMD)
    #return
    Peak_Regions = []
    
    curChrID = ""
    chrPosStart = 0
    maxDepth = 0
    chrPosLast = 0
    IN = os.popen(CMD)
    for line in IN:
        chrID, chrPos, depth = line.strip().split()
        chrPos = int(chrPos)
        depth = int(depth)
        if depth > minCov:
            if curChrID != chrID:
                if chrPosLast-chrPosStart>minPeakWidth and maxDepth>2200:
                    print("{}:{}-{}".format(curChrID,chrPosStart,chrPosLast))
                    Peak_Regions.append([curChrID,chrPosStart,chrPosLast])
                curChrID = chrID
                chrPosStart = chrPosLast = chrPos
                maxDepth = depth
            else:
                if chrPos-chrPosLast>10:
                    if chrPosLast-chrPosStart>minPeakWidth and maxDepth>2200:
                        print("{}:{}-{}".format(curChrID,chrPosStart,chrPosLast))
                        Peak_Regions.append([curChrID,chrPosStart,chrPosLast])
                    curChrID = chrID
                    chrPosStart = chrPosLast = chrPos
                    maxDepth = depth
                chrPosLast = chrPos
                maxDepth = max(maxDepth, depth)
    if chrPosLast-chrPosStart>minPeakWidth and maxDepth>2200:
        print("{}:{}-{}".format(curChrID,chrPosStart,chrPosLast))
        Peak_Regions.append([curChrID,chrPosStart,chrPosLast])
    return Peak_Regions

DMSO_regions_pos, DMSO_regions_neg = Bam_cov(DMSO_BAM)
NAI_regions_pos, NAI_regions_neg = Bam_cov(NAIN3_BAM)

### Save

nai_pos_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/nai_pos.txt"
nai_neg_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/nai_neg.txt"

open(nai_pos_ofile, "w").writelines( "\n".join([ "\t".join([it[0], str(it[1]), str(it[2])]) for it in  NAI_regions_pos ]) )
open(nai_neg_ofile, "w").writelines( "\n".join([ "\t".join([it[0], str(it[1]), str(it[2])]) for it in  NAI_regions_neg ]) )

dmso_pos_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/dmso_pos.txt"
dmso_neg_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/dmso_neg.txt"

open(dmso_pos_ofile, "w").writelines( "\n".join([ "\t".join([it[0], str(it[1]), str(it[2])]) for it in  DMSO_regions_pos ]) )
open(dmso_neg_ofile, "w").writelines( "\n".join([ "\t".join([it[0], str(it[1]), str(it[2])]) for it in  DMSO_regions_neg ]) )


################################
#####  2. 合并DMSO和NAI Peaks
################################

nai_pos_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/nai_pos.txt"
nai_neg_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/nai_neg.txt"
dmso_pos_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/dmso_pos.txt"
dmso_neg_ofile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/combine_20190514/dmso_neg.txt"

### Read

NAI_regions_pos = [ line.strip().split() for line in open(nai_pos_ofile).readlines()]
NAI_regions_pos = [ [it[0], int(it[1]), int(it[2])] for it in NAI_regions_pos ]
NAI_regions_neg = [ line.strip().split() for line in open(nai_neg_ofile).readlines()]
NAI_regions_neg = [ [it[0], int(it[1]), int(it[2])] for it in NAI_regions_neg ]

DMSO_regions_pos = [ line.strip().split() for line in open(dmso_pos_ofile).readlines()]
DMSO_regions_pos = [ [it[0], int(it[1]), int(it[2])] for it in DMSO_regions_pos ]
DMSO_regions_neg = [ line.strip().split() for line in open(dmso_neg_ofile).readlines()]
DMSO_regions_neg = [ [it[0], int(it[1]), int(it[2])] for it in DMSO_regions_neg ]


def combine_DMSO_NAI_Peaks(DMSO_Peaks, NAI_Peaks):
    DMSO_Peaks.sort(key=lambda x: (x[0], x[1]))
    NAI_Peaks.sort(key=lambda x: (x[0], x[1]))
    combined_Peaks = {}
    for d_peak in DMSO_Peaks:
        for n_peak in NAI_Peaks:
            if d_peak[0] == n_peak[0] and d_peak[1]<=n_peak[2] and n_peak[1]<=d_peak[2]:
                peak_start = min(d_peak[1], n_peak[1])
                peak_end = max(d_peak[2], n_peak[2])
                try:
                    combined_Peaks[ n_peak[0] ].append( (peak_start, peak_end) )
                except KeyError:
                    combined_Peaks[ n_peak[0] ] = [ (peak_start, peak_end) ]
                break
    return combined_Peaks

combined_Peaks_pos = combine_DMSO_NAI_Peaks(DMSO_regions_pos, NAI_regions_pos)
combined_Peaks_neg = combine_DMSO_NAI_Peaks(DMSO_regions_neg, NAI_regions_neg)

################################
#####  3. Read Repeat
################################

def read_chrID_scaffoldID(inFile):
    import re
    
    IDconversion = {}
    for line in open(inFile):
        if line[0] == '#':
            continue
        data = line.strip().split('\t')
        if data[2] == 'region' and data[0].startswith('NC_'):
            #print data[8]
            chrID = re.findall(";chromosome=(\\S+);gbk", data[8])
            if chrID:
                chrID = 'chr' + chrID[0]
            else:
                if "Name=MT" in data[8]:
                    chrID = 'chrM'
                else:
                    continue
            
            IDconversion[ data[0] ] = chrID
    
    return IDconversion

def load_repeat():
    Repeat = {}
    IDconversion = read_chrID_scaffoldID(inFile="/150T/zhangqf/GenomeAnnotation/NCBI/hg38.gff3")
    for line in open("/150T/zhangqf/GenomeAnnotation/RepeatMask/NCBI/hg38_mask"):
        try:
            chrID, start, end, Type = line.strip().split('\t')
        except ValueError:
            print(line.strip())
            continue
        if not re.match(r"\(\S+\)n", Type):
            continue
        #print (Type)
        if chrID not in IDconversion:
            continue
        normal_chrID = IDconversion[chrID]
        try:
            Repeat[normal_chrID].append( (int(start), int(end)) )
        except KeyError:
            Repeat[normal_chrID] = [ (int(start), int(end)) ]
    for chrID in Repeat:
        Repeat[chrID].sort(key=lambda x: x[0])
    return  Repeat

repeat = load_repeat()


################################
#####  4. 把序列取出来，作为参考基因组
################################

Seqer = Seq.seqClass("/150T/zhangqf/GenomeAnnotation/genome/hg38.fa")
gaper_hg38 = GAP.init("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/GTF/hg38.genomeCoor.bed", "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/GTF/hg38_transcriptome.fa")

def WhatIsIt(chrID,chrStart,chrEnd,strand,GAPer):
    transList = GAPer.genomeCoor2transCoor(chrID,chrStart,chrEnd,strand)
    transIDList = [ it[3] for it in transList ]
    ID = ""
    if transList:
        for tid in transIDList:
            ft = GAPer.getTransFeature(tid)
            if ft['gene_type'] in ('guide_RNA','snRNA','snoRNA','tRNA','precursor_RNA'):
                ID = ft['gene_type']+"_"+ft['gene_name']
                break
        if ID == "":
            for tid in transIDList:
                ft = GAPer.getTransFeature(tid)
                if ft['gene_type'] == 'mRNA':
                    ID = ft['gene_type']+"_exon_"+ft['gene_name']
                    break
        if ID == "":
            for tid in transIDList:
                ft = GAPer.getTransFeature(tid)
                if ft['gene_type'] == 'lnc_RNA':
                    ID = ft['gene_type']+"_exon_"+ft['gene_name']
                    break
        if ID == "":
            tid = transIDList[0]
            ft = GAPer.getTransFeature(tid)
            ID = ft['gene_type']+"_exon_"+ft['gene_name']
    else:
        geneList = GAPer.genomeCoor2geneCoor(chrID,chrStart,chrEnd,strand)
        geneIDList = [ it[3] for it in geneList ]
        if geneList:
            for gid in geneIDList:
                ft = GAPer.getGeneParser()[gid]
                if 'mRNA' in ft['gene_type']:
                    ID = "mRNA_intron_"+ft['gene_name']
            if ID == "":
                for gid in geneIDList:
                    ft = GAPer.getGeneParser()[gid]
                    if 'lnc_RNA' in ft['gene_type']:
                        ID = "lnc_RNA_intron_"+ft['gene_name']
            if ID == "":
                gid = geneIDList[0]
                ft = GAPer.getGeneParser()[gid]
                ID = ft['gene_type'][0]+"_intron_"+ft['gene_name']
        else:
            ID = "intergenic"
    
    return ID

def inRepeat(repeats, chrID, Start, End):
    for r_s, r_e in repeats[chrID]:
        if r_s<End and Start<r_e:
            return True
    return False

OUT = open("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/ref_genome.fa", "w")
for chrID in combined_Peaks_pos:
    for start, end in combined_Peaks_pos[chrID]:
        start -= 2
        end += 2
        in_repeat = inRepeat(repeat, chrID, start, end)
        seq = Seqer.fetch(chrID, start, end, '+')
        ID = WhatIsIt(chrID,start,end,'+',gaper_hg38)
        if in_repeat: 
            print ("In repeat: "+seq+"\t"+ID)
            continue
        
        header = ">{}_{}_{}_pos".format(chrID, start, end)+"_"+ID
        
        OUT.writelines(header+"\n")
        OUT.writelines(Seq.flat_seq(seq) +"\n")

for chrID in combined_Peaks_neg:
    for start, end in combined_Peaks_neg[chrID]:
        start -= 2
        end += 2
        in_repeat = inRepeat(repeat, chrID, start, end)
        seq = Seqer.fetch(chrID, start, end, '-')
        ID = WhatIsIt(chrID,start,end,'-',gaper_hg38)
        if in_repeat: 
            print ("In repeat: "+seq+"\t"+ID)
            continue
        
        header = ">{}_{}_{}_neg".format(chrID, start, end)+"_"+ID
        
        OUT.writelines(header+"\n")
        OUT.writelines(Seq.flat_seq(seq) +"\n")

OUT.close()

################################
#####  5. Overlap Map to genome and Older reference
################################

older_ref_fn = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/1.STAR_map/ref/high_exp_RNAs.fa"
older_ref = General.load_fasta(older_ref_fn)

combined_Peaks = {}
for chrID in combined_Peaks_pos:
    combined_Peaks[chrID] = []
    for start, end in combined_Peaks_pos[chrID]:
        combined_Peaks[chrID].append([start, end, '+'])

for chrID in combined_Peaks_neg:
    combined_Peaks[chrID] = combined_Peaks.get(chrID, [])
    for start, end in combined_Peaks_neg[chrID]:
        combined_Peaks[chrID].append([start, end, '-'])

OUT = open("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/losted_ref_seq.fa", 'w')
inTransNum = 0
total = 0
count = 0

for chrID in combined_Peaks:
    for start,end,strand in combined_Peaks[chrID]:
        trans_list = gaper_hg38.genomeCoor2transCoor(chrID,start,end,strand)
        total += 1
        if trans_list:
            inTransNum += 1
            trans_id_list = [ it[3] for it in trans_list ]
            inOneTrans = False; whoInOneTrans = ""
            for tid in trans_id_list:
                ft = gaper_hg38.getTransFeature(tid)
                gene_type = ft['gene_type']
                if gene_type in ('mRNA', 'protein_coding'):
                    continue
                seq = gaper_hg38.getTransSeq(tid)
                mid_seq = seq[5:-5]
                found = False
                for ref_seq in older_ref.values():
                    if mid_seq in ref_seq:
                        count += 1
                        whoInOneTrans = tid
                        found = True
                        break
                if found:
                    inOneTrans = True
                    break
            if not inOneTrans:
                for tid in trans_id_list:
                    ft = gaper_hg38.getTransFeature(tid)
                    if ft['gene_type'] not in ('mRNA', 'protein_coding') and ft['trans_len']<200:
                        OUT.writelines(">{}_{}_{}_ncbi\n".format(ft['gene_type'],tid,ft['gene_name']))
                        OUT.writelines(Seq.flat_seq(gaper_hg38.getTransSeq(tid))+"\n")
                        break

OUT.close()


################################
#####  6. Create IGV batch script
################################

def draw_2nd_struct(inFolder, OUT=sys.stdout, OUTSCRIPT=sys.stdout):
    inFolder = inFolder.rstrip('/')
    file_list = os.listdir(inFolder)
    
    OUTSCRIPT.writelines("snapshotDirectory /tmp/snapshot\n")
    for file in file_list:
        if not file.endswith('.shape'):
            continue
        shape = General.load_SHAPEMap(inFolder+"/"+file)
        seq = shape['seq']
        shape = shape['shape_pro_list']
        dot = Structure.predict_structure(seq, shape)
        OUT.writelines(Visual.Plot_RNAStructure_Shape(seq,dot,shape,title=file.rstrip('.shape'))+"\n")
        #print (file)
        chrID,start,end = file.split('_')[:3]
        start,end = int(start),int(end)
        OUTSCRIPT.writelines("goto {}:{}-{}\n".format(chrID,start-100,end+100))
        OUTSCRIPT.writelines("snapshot "+file.split('.')[0]+".png\n")


inFolder = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/Output2/shape_files"
OUT = open('/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/Output2/structure.cmd', 'w')
OUTSCRIPT = open('/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/3.GenomeSHAPE/Output2/IGV-script.cmd', 'w')
draw_2nd_struct(inFolder, OUT, OUTSCRIPT)
OUT.close()
OUTSCRIPT.close()




