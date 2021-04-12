
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def ref_profile(MD_tag):
    """
    MD_tag          -- 59A11/51^CT18/6G4C20G1A5C5A1^C3A15G1G15
    
    Return the profile corresponding to the raw sequence
    """
    profile = ""
    current = ""
    
    deletion_state = False
    for alpha in list(MD_tag):
        if '0'<=alpha<='9':
            current += alpha
            deletion_state = False
        else:
            if current:
                profile += "0"*int(current)
                current = ""
            if alpha == '^':
                deletion_state = True
            elif alpha in ('A','T','C','G'):
                if deletion_state:
                    profile += alpha
                else:
                    profile += '1'
    if current:
        profile += "0"*int(current)
    
    print(profile)

def mutate_count(Cigar, MD_tag):
    counts = 0
    
    deletion_state = False
    for alpha in list(MD_tag):
        if '0'<=alpha<='9':
            deletion_state = False
        else:
            if alpha == '^':
                deletion_state = True
                counts += 1
            elif alpha in ('A','T','C','G'):
                if not deletion_state:
                    counts += 1
    
    counts += Cigar.count('I')
    
    return counts

##### miRNA and 5SrRNA
def count_mutevent_num(inBam, ratio=0.3):
    mutevent_num = []
    
    if ratio == None:
        commands = "samtools view %s" % (inBam, )
    else:
        commands = "samtools view -s %s %s" % (ratio, inBam)
    
    for line in os.popen(commands):
        data = line.strip().split()
        ref_id, Cigar, MD_tag = data[2], data[5], data[16]
        if ref_id.startswith('miRNA') or ref_id.startswith('rRNA_human_5S'):
            MD_tag = MD_tag.lstrip("MD:Z:")
            mut_count = mutate_count(Cigar, MD_tag)
            mutevent_num.append(mut_count)
    
    return sorted(mutevent_num)

def get_ratio_list(input_list):
    length = len(input_list)
    n0 = input_list.count(0)
    n1 = input_list.count(1)
    n2 = input_list.count(2)
    n3 = input_list.count(3)
    n4 = input_list.count(4)
    n5 = input_list.count(5)
    nl6 = length-n0-n1-n2-n3-n4-n5
    ratio = [ 1.0*n0/length, 1.0*n1/length, 1.0*n2/length, 1.0*n3/length, 1.0*n4/length, 1.0*n5/length, 1.0*nl6/length ]
    numlist = [ n0,            n1,            n2,            n3,            n4,            n5,            nl6 ]
    ratio = [ round(it, 5) for it in ratio ]
    return ratio, numlist



ROOT = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/"
DMSO_SMR_SSII_rep1 = ROOT+"NAI_100mm_vivo_SMR_SSII_repX/mapping/DMSO_1.clean.Aligned.sortedByCoord.out.bam"
DMSO_SMR_SSII_rep2 = ROOT+"NAI_100mm_vivo_SMR_SSII_repX/mapping/DMSO_2.clean.Aligned.sortedByCoord.out.bam"
NAI_100mm_vivo_SMR_SSII_rep1 = ROOT+"NAI_100mm_vivo_SMR_SSII_repX/mapping/INVIVO_1.clean.Aligned.sortedByCoord.out.bam"
NAI_100mm_vivo_SMR_SSII_rep2 = ROOT+"NAI_100mm_vivo_SMR_SSII_repX/mapping/INVIVO_2.clean.Aligned.sortedByCoord.out.bam"
NAI_100mm_vitro_SMR_SSII_rep1 = ROOT+"NAI_100mm_vitro_SMR_SSII_repX/mapping/INVITRO_1.clean.Aligned.sortedByCoord.out.bam"
NAI_100mm_vitro_SMR_SSII_rep2 = ROOT+"NAI_100mm_vitro_SMR_SSII_repX/mapping/INVITRO_2.clean.Aligned.sortedByCoord.out.bam"

sample_names = [ 'DMSO_SMR_SSII_rep1','DMSO_SMR_SSII_rep2',
                'NAI_100mm_vivo_SMR_SSII_rep1','NAI_100mm_vivo_SMR_SSII_rep2',
                'NAI_100mm_vitro_SMR_SSII_rep1','NAI_100mm_vitro_SMR_SSII_rep2' ]

##########################
###   Read data
##########################

MutEventCount = {}
for sample_name in sample_names:
    print(sample_name)
    MutEventCount[sample_name] = count_mutevent_num(eval(sample_name), ratio=None)

##########################
###   Count
##########################

sample_rl_list = []
sample_num_list = []

for sample_name in sample_names:
    rl, num = get_ratio_list(MutEventCount[sample_name])
    sample_rl_list.append(rl)
    sample_num_list.append(num)

rl_df = pd.DataFrame(sample_rl_list, columns=['0','1','2','3','4','5','>=6'])
rl_df.index = sample_names
rl_df.to_csv("/Share/home/zhangqf7/figs/rl_df.csv", sep="\t")

num_df = pd.DataFrame(sample_num_list, columns=['0','1','2','3','4','5','>=6'])
num_df.index = sample_names
num_df.to_csv("/Share/home/zhangqf7/figs/num_df.csv", sep="\t")


##########################
###   Plot 2 -- barplot  : Use all mutations, not 1-2 mutations
##########################

dmso_0_rep1 = rl_df.loc['DMSO_SMR_SSII_rep1','0']
dmso_12_rep1 = 1-dmso_0_rep1

dmso_0_rep2 = rl_df.loc['DMSO_SMR_SSII_rep2','0']
dmso_12_rep2 = 1 - dmso_0_rep2

vivo_0_rep1 = rl_df.loc['NAI_100mm_vivo_SMR_SSII_rep1','0']
vivo_12_rep1 = 1 - vivo_0_rep1

vivo_0_rep2 = rl_df.loc['NAI_100mm_vivo_SMR_SSII_rep2','0']
vivo_12_rep2 = 1 - vivo_0_rep2

vitro_0_rep1 = rl_df.loc['NAI_100mm_vitro_SMR_SSII_rep1','0']
vitro_12_rep1 = 1 - vitro_0_rep1

vitro_0_rep2 = rl_df.loc['NAI_100mm_vitro_SMR_SSII_rep2','0']
vitro_12_rep2 = 1 - vitro_0_rep2


data = []
data.append( (dmso_0_rep1, "dmso", "rep1", "Without mutates") )
data.append( (dmso_12_rep1, "dmso", "rep1", "With mutates") )

data.append( (dmso_0_rep2, "dmso", "rep2", "Without mutates") )
data.append( (dmso_12_rep2, "dmso", "rep2", "With mutates") )

data.append( (vivo_0_rep1, "vivo", "rep1", "Without mutates") )
data.append( (vivo_12_rep1, "vivo", "rep1", "With mutates") )

data.append( (vivo_0_rep2, "vivo", "rep2", "Without mutates") )
data.append( (vivo_12_rep2, "vivo", "rep2", "With mutates") )

data.append( (vitro_0_rep1, "vitro", "rep1", "Without mutates") )
data.append( (vitro_12_rep1, "vitro", "rep1", "With mutates") )

data.append( (vitro_0_rep2, "vitro", "rep2", "Without mutates") )
data.append( (vitro_12_rep2, "vitro", "rep2", "With mutates") )

df_data = pd.DataFrame(data=data, columns=['ratio','type','rep','mutnum'])

sns.barplot(data=df_data, x='type', y='ratio', hue='mutnum')
plt.savefig("figs/xxx.pdf")
plt.show()





