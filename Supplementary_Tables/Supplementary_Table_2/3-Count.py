
importCommon()

GAPer = GAP.init("/150T/zhangqf/GenomeAnnotation/NCBI/hg38.genomeCoor.bed")

def count_Genome_bam(inBamFn, GAPer):
    intergene = 0
    intron = 0
    RNAtype = {}
    
    IN = os.popen("samtools view %s" % (inBamFn, ))
    for line in IN:
        if line[0] == '@': continue
        data = line.split('\t')
        strand = "+" if data[1]=="0" else "-"
        refID = data[2]
        map_start = int(data[3])
        gene_regions = GAPer.genomeCoor2geneCoor(refID, map_start, map_start+25, strand)
        trans_regions = GAPer.genomeCoor2transCoor(refID, map_start, map_start+25, strand)
        
        if len(gene_regions) == 0:
            intergene += 1
        else:
            if len(trans_regions) != 0:
                trans_set = set([ it[3] for it in trans_regions ])
                gene_types = set([ GAPer.getTransFeature(tid)['gene_type'] for tid in trans_set ])
                if 'mRNA' in gene_types or 'protein_coding' in gene_types:
                    RNAtype['mRNA'] = RNAtype.get('mRNA', 0) + 1
                elif 'tRNA' in gene_types:
                    RNAtype['tRNA'] = RNAtype.get('tRNA', 0) + 1
                elif 'precursor_RNA' in gene_types or 'miRNA' in gene_types:
                    RNAtype['miRNA'] = RNAtype.get('miRNA', 0) + 1
                elif 'SRP_RNA' in gene_types:
                    RNAtype['SRP_RNA'] = RNAtype.get('SRP_RNA', 0) + 1
                elif 'snoRNA' in gene_types:
                    RNAtype['snoRNA'] = RNAtype.get('snoRNA', 0) + 1
                elif 'guide_RNA' in gene_types:
                    RNAtype['guide_RNA'] = RNAtype.get('guide_RNA', 0) + 1
                elif 'snRNA' in gene_types:
                    RNAtype['snRNA'] = RNAtype.get('snRNA', 0) + 1
                elif 'vault_RNA' in gene_types:
                    RNAtype['vault_RNA'] = RNAtype.get('vault_RNA', 0) + 1
                elif 'lnc_RNA' in gene_types:
                    RNAtype['lnc_RNA'] = RNAtype.get('lnc_RNA', 0) + 1
                else:
                    if len(gene_types)>=2: print(gene_types)
                    for gene_type in gene_types:
                        RNAtype[gene_type] = RNAtype.get(gene_type, 0) + 1.0/len(gene_types)
            else:
                intron += 1
    
    RNAtype['intron'] = intron
    RNAtype['intergene'] = intergene
    return RNAtype

def read_STAR_Logfinalout(Log_final_out):
    """
    Log_final_out           -- Folder or Log.final.out produced by STAR
    """
    if os.path.isdir(Log_final_out):
        Log_final_out = Log_final_out.rstrip('/')
        dirfiles = os.listdir(Log_final_out)
        if 'Log.final.out' in dirfiles:
            Log_final_out += '/Log.final.out'
        else:
            raise RuntimeError("No Log.final.out file in directory")
    
    Mapping = {}
    for line in open(Log_final_out):
        if '|' in line:
            left, right = line.split('|')
            left = left.strip()
            right = right.strip()
            if 'Number of input reads' == left:
                Mapping['input'] = int(right)
            elif 'Average input read length' == left:
                Mapping['length'] = int(right)
            elif 'Uniquely mapped reads number' == left:
                Mapping['uniq_num'] = int(right)
            elif 'Number of reads mapped to multiple loci' == left:
                Mapping['mult_num'] = int(right)
            elif 'Number of reads mapped to too many loci' == left:
                Mapping['toomany_num'] = int(right)
    
    unmapped = Mapping['input'] - Mapping['uniq_num'] - Mapping['mult_num'] - Mapping['toomany_num']
    Mapping['unmap'] = unmapped
    return Mapping

def count_smallRNA_Bam(inBamFn):
    SmallRNATypes = {}
    IN = os.popen("samtools view %s" % (inBamFn, ))
    for line in IN:
        data = line.split()
        rna_type = data[2].split('_')[0]
        try:
            SmallRNATypes[rna_type] += 1
        except KeyError:
            SmallRNATypes[rna_type] = 1
    return SmallRNATypes

def calculate_reads_distribution(smallRNA_mapping_Logfinalout, genome_mapping_Logfinalout, smallSTAR, genomeSTAR):
    ## Check
    assert smallRNA_mapping_Logfinalout['unmap'] + smallRNA_mapping_Logfinalout['toomany_num'] == genome_mapping_Logfinalout['input']
    assert smallRNA_mapping_Logfinalout['uniq_num'] + smallRNA_mapping_Logfinalout['mult_num'] == sum(smallSTAR.values())
    assert genome_mapping_Logfinalout['uniq_num'] + genome_mapping_Logfinalout['mult_num'] == sum(genomeSTAR.values())
    
    ## smallRNA statistics
    smallRNA = copy.deepcopy(smallSTAR)
    genomeRNA = copy.deepcopy(genomeSTAR)
    for rna_type in smallRNA:
        grna_type = rna_type+"_RNA" if rna_type in ('Y','vault') else rna_type
        if grna_type in genomeRNA:
            smallRNA[rna_type] += genomeRNA[grna_type]
            del genomeRNA[grna_type]
    
    ## genome
    map_to_genome = {}
    if 'lnc_RNA' in genomeRNA:
        genomeRNA['lncRNA'] = genomeRNA.get('lncRNA', 0) + genomeRNA['lnc_RNA']
        del genomeRNA['lnc_RNA']
    for rna_type in genomeRNA:
        if rna_type in ('intergene','lncRNA','mRNA','intron') :
            map_to_genome[rna_type] = map_to_genome.get(rna_type, 0) + genomeRNA[rna_type]
        else:
            map_to_genome['other'] = map_to_genome.get('other', 0) + genomeRNA[rna_type]
    
    smallRNA.update(map_to_genome)
    smallRNA['unmap'] = genome_mapping_Logfinalout['unmap'] + genome_mapping_Logfinalout['toomany_num']
    
    assert sum(smallRNA.values()) == smallRNA_mapping_Logfinalout['input']
    
    return smallRNA

def combine_list(statis, rna_type_list):
    combineList = []
    for rna_type in rna_type_list:
        combineList.append( statis[rna_type] )
    return combineList

def get_statistics(sample_name):
    smallBam = count_smallRNA_Bam(ROOT_1+sample_name+"/Aligned.out.bam")
    smallSTAR = read_STAR_Logfinalout(ROOT_1+sample_name)
    GenomeBam = count_Genome_bam(ROOT_2+sample_name+"/Aligned.out.bam", GAPer)
    GenomeSTAR = read_STAR_Logfinalout(ROOT_2+sample_name)
    statis = calculate_reads_distribution(smallSTAR, GenomeSTAR, smallBam, GenomeBam)
    combList = combine_list(statis, rna_type_list)
    return combList

rna_type_list = ['tRNA', 'chrM', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'other', 'Y', 'intergene', 'lncRNA', 'vault', 'mRNA', 'intron', 'unmap']

ROOT_1 = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/1.STAR_mapping/"
ROOT_2 = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/2.mapGenome/"


###############
#### 20190514 data
###############

INPUT_rep1_20190514 = get_statistics("INPUT_rep1_20190514")
INPUT_rep2_20190514 = get_statistics("INPUT_rep2_20190514")
INPUT_rep3_20190514 = get_statistics("INPUT_rep3_20190514")
RIP_DMSO_rep1_20190514 = get_statistics("RIP_DMSO_rep1_20190514")
RIP_DMSO_rep2_20190514 = get_statistics("RIP_DMSO_rep2_20190514")
RIP_DMSO_rep3_20190514 = get_statistics("RIP_DMSO_rep3_20190514")
RIP_NAIN3_rep1_20190514 = get_statistics("RIP_NAIN3_rep1_20190514")
RIP_NAIN3_rep2_20190514 = get_statistics("RIP_NAIN3_rep2_20190514")
RIP_NAIN3_rep3_20190514 = get_statistics("RIP_NAIN3_rep3_20190514")

samples = ['INPUT_rep1_20190514', 'INPUT_rep2_20190514', 'INPUT_rep3_20190514', 
            'RIP_DMSO_rep1_20190514', 'RIP_DMSO_rep2_20190514', 'RIP_DMSO_rep3_20190514', 
            'RIP_NAIN3_rep1_20190514','RIP_NAIN3_rep2_20190514','RIP_NAIN3_rep3_20190514']

df = pd.DataFrame([eval(file) for file in samples], index=samples, columns=rna_type_list)

df.to_csv("/Share/home/zhangqf7/figs/Pan.csv", sep="\t")

df_ratio = df.divide(df.sum(axis=1), axis='rows') 
df_ratio.to_csv("/Share/home/zhangqf7/figs/Pan_ratio.csv", sep="\t")


































































































###############
#### More More Data
###############


NAI_100mm_vivo_CIRL_CENR_SSII_rep1 = get_statistics("NAI_100mm_vivo_CIRL_CENR_SSII_rep1")
NAI_100mm_vivo_CIRL_CENR_SSII_rep2 = get_statistics("NAI_100mm_vivo_CIRL_CENR_SSII_rep2")
NAI_100mm_vivo_CIRL_SSII = get_statistics("NAI_100mm_vivo_CIRL_SSII")
NAI_100mm_vivo_CIRL_TGIII = get_statistics("NAI_100mm_vivo_CIRL_TGIII")


colnames = ['NAI_100mm_vivo_CIRL_CENR_SSII_rep1', 'NAI_100mm_vivo_CIRL_CENR_SSII_rep2', 'NAI_100mm_vivo_CIRL_SSII', 'NAI_100mm_vivo_CIRL_TGIII']
df = pd.DataFrame([NAI_100mm_vivo_CIRL_CENR_SSII_rep1,NAI_100mm_vivo_CIRL_CENR_SSII_rep2,NAI_100mm_vivo_CIRL_SSII,NAI_100mm_vivo_CIRL_TGIII], index=colnames, columns=rna_type_list)

df.to_csv("/Share/home/zhangqf7/figs/Pan.csv", sep="\t")

df_ratio = df.divide(df.sum(axis=1), axis='rows') 
df_ratio.to_csv("/Share/home/zhangqf7/figs/Pan_ratio.csv", sep="\t")


###############
#### 20190412 data
###############

DMSO_20190412 = get_statistics("DMSO_20190412")
NAI_100mm_10min_rep1_20190412 = get_statistics("NAI_100mm_10min_rep1_20190412")
NAI_100mm_10min_rep2_20190412 = get_statistics("NAI_100mm_10min_rep2_20190412")
NAI_100mm_5min_rep1_20190412 = get_statistics("NAI_100mm_5min_rep1_20190412")
NAI_100mm_5min_rep2_20190412 = get_statistics("NAI_100mm_5min_rep2_20190412")
NAI_50mm_10min_rep1_20190412 = get_statistics("NAI_50mm_10min_rep1_20190412")
NAI_50mm_10min_rep2_20190412 = get_statistics("NAI_50mm_10min_rep2_20190412")
NAI_50mm_5min_rep1_20190412 = get_statistics("NAI_50mm_5min_rep1_20190412")
NAI_50mm_5min_rep2_20190412 = get_statistics("NAI_50mm_5min_rep2_20190412")


colnames = ['DMSO_20190412', 'NAI_100mm_10min_rep1_20190412', 'NAI_100mm_10min_rep2_20190412', 'NAI_100mm_5min_rep1_20190412', 'NAI_100mm_5min_rep2_20190412', 'NAI_50mm_10min_rep1_20190412', 'NAI_50mm_10min_rep2_20190412','NAI_50mm_5min_rep1_20190412','NAI_50mm_5min_rep2_20190412']
df = pd.DataFrame([DMSO_20190412,NAI_100mm_10min_rep1_20190412,NAI_100mm_10min_rep2_20190412,NAI_100mm_5min_rep1_20190412,NAI_100mm_5min_rep2_20190412,NAI_50mm_10min_rep1_20190412,NAI_50mm_10min_rep2_20190412,NAI_50mm_5min_rep1_20190412,NAI_50mm_5min_rep2_20190412], index=colnames, columns=rna_type_list)

df.to_csv("/Share/home/zhangqf7/figs/Pan.csv", sep="\t")

df_ratio = df.divide(df.sum(axis=1), axis='rows') 
df_ratio.to_csv("/Share/home/zhangqf7/figs/Pan_ratio.csv", sep="\t")




###############
#### 20190422 data
###############

DMSO_20190422 = get_statistics("DMSO_20190422")
NAI_100mm_10min_rep1_20190422 = get_statistics("NAI_100mm_10min_rep1_20190422")
NAI_100mm_10min_rep2_20190422 = get_statistics("NAI_100mm_10min_rep2_20190422")
NAI_100mm_5min_rep1_20190422 = get_statistics("NAI_100mm_5min_rep1_20190422")
NAI_100mm_5min_rep2_20190422 = get_statistics("NAI_100mm_5min_rep2_20190422")
NAI_50mm_10min_rep1_20190422 = get_statistics("NAI_50mm_10min_rep1_20190422")
NAI_50mm_10min_rep2_20190422 = get_statistics("NAI_50mm_10min_rep2_20190422")
NAI_50mm_5min_rep1_20190422 = get_statistics("NAI_50mm_5min_rep1_20190422")
NAI_50mm_5min_rep2_20190422 = get_statistics("NAI_50mm_5min_rep2_20190422")

samples = ['DMSO_20190422', 'NAI_100mm_10min_rep1_20190422', 'NAI_100mm_10min_rep2_20190422', 
            'NAI_100mm_5min_rep1_20190422', 'NAI_100mm_5min_rep2_20190422', 'NAI_50mm_10min_rep1_20190422', 
            'NAI_50mm_10min_rep2_20190422','NAI_50mm_5min_rep1_20190422','NAI_50mm_5min_rep2_20190422']

df = pd.DataFrame([eval(file) for file in samples], index=samples, columns=rna_type_list)

df.to_csv("/Share/home/zhangqf7/figs/Pan.csv", sep="\t")

df_ratio = df.divide(df.sum(axis=1), axis='rows') 
df_ratio.to_csv("/Share/home/zhangqf7/figs/Pan_ratio.csv", sep="\t")


###############
#### 20190606 data
###############

INPUT_rep1_20190606 = get_statistics("INPUT_rep1_20190606")
INPUT_rep2_20190606 = get_statistics("INPUT_rep2_20190606")
RIP_rep1_20190606 = get_statistics("RIP_rep1_20190606")
RIP_rep2_20190606 = get_statistics("RIP_rep2_20190606")

samples = ['INPUT_rep1_20190606', 'INPUT_rep2_20190606', 'RIP_rep1_20190606', 'RIP_rep1_20190606'] 
df = pd.DataFrame([eval(file) for file in samples], index=samples, columns=rna_type_list)
df.to_csv("/Share/home/zhangqf7/figs/Pan.csv", sep="\t")
df_ratio = df.divide(df.sum(axis=1), axis='rows') 
df_ratio.to_csv("/Share/home/zhangqf7/figs/Pan_ratio.csv", sep="\t")





