
import DicerMiRNA

root = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/"


#############################
#####  NAI_100mm_vitro_SMR_SSII_repX
#############################

shape,sequence = DicerMiRNA.read_shape(root+"NAI_100mm_vitro_SMR_SSII_repX/shape_files/", 
        min_valid_ratio=0, min_valid_base_num=5, relocate=True, loadAll=False)

OUT = open("/tmp/icSHAPE-MaP/sRNA_icSHAPE-MaP_invitro.txt", 'w')
for tid in sequence.keys():
    i = 1
    for base,score in zip(sequence[tid], shape[tid]):
        print("%s\t%d\t%c\t%s"%(tid, i, base, score), file=OUT)
        i += 1

OUT.close()

#############################
#####  NAI_100mm_vivo_SMR_SSII_repX
#############################

shape,sequence = DicerMiRNA.read_shape(root+"NAI_100mm_vivo_SMR_SSII_repX/shape_files/", 
        min_valid_ratio=0, min_valid_base_num=5, relocate=True, loadAll=False)

OUT = open("/tmp/icSHAPE-MaP/sRNA_icSHAPE-MaP_invivo.txt", 'w')
for tid in sequence.keys():
    i = 1
    for base,score in zip(sequence[tid], shape[tid]):
        print("%s\t%d\t%c\t%s"%(tid, i, base, score), file=OUT)
        i += 1

OUT.close()

#############################
#####  RIP_NAIN3_repX_20190514
#############################

shape,sequence = DicerMiRNA.read_shape(root+"RIP_NAIN3_repX_20190514/shape_files/", 
        min_valid_ratio=0, min_valid_base_num=5, relocate=True, loadAll=False)

OUT = open("/tmp/icSHAPE-MaP/RIP_icSHAPE-MaP.txt", 'w')
for tid in sequence.keys():
    i = 1
    for base,score in zip(sequence[tid], shape[tid]):
        print("%s\t%d\t%c\t%s"%(tid, i, base, score), file=OUT)
        i += 1

OUT.close()




