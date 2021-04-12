import General
import numpy as np
import Structure
import matplotlib.pyplot as plt
import Colors, Visual
import GAP

import DicerMiRNA


RNA_shape,RNA_seq = DicerMiRNA.read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files")


###### RNU7

seq = RNA_seq['snRNA_ENST00000458811.1'][30:]
shape = RNA_shape['snRNA_ENST00000458811.1'][30:]

dot = "...((((((((((.......))))))))))..."
AUC = General.calc_AUC_v2(dot, shape)

print( Visual.Plot_RNAStructure_Shape(seq, dot, shape, mode='label', title='U7 (AUC=%.3f)'%(AUC, )) )



###### SNORA3B

seq = RNA_seq["tRNA_Gln-TTG-2-1"]
shape = RNA_shape["tRNA_Gln-TTG-2-1"]

dot = '((((((...((((........)))).(((((.......)))))....(((((.......))))).)))))).'
AUC = General.calc_AUC_v2(dot, shape)

print( Visual.Plot_RNAStructure_Shape(seq, dot, shape, mode='label', title='Gln-TTG-2-1 (AUC=%.3f)'%(AUC, )) )





#################################
##### 在数据中寻找例子
#################################




import DicerMiRNA
RNA_shape,RNA_seq = DicerMiRNA.read_shape("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/NAI_100mm_vivo_SMR_SSII_repX/shape_files")
human_dot = General.load_dot("/150T/zhangqf/GenomeAnnotation/Rfam/Parsed_Structure/human.dot")
seq_dot = { seq.replace('T','U'):dot for seq,dot in human_dot.values() }

##### 选出有序列的 tid

trans_with_structures = []
count = 0
for tid in RNA_seq:
    for seq in seq_dot:
        if RNA_seq[tid] == seq:
            count += 1
            trans_with_structures.append(tid)
            break

print(count)

##### 比较序列和shape

for tid in trans_with_structures:
    seq = RNA_seq[tid]
    shape = RNA_shape[tid]
    dot = seq_dot[seq]
    if shape.count('NULL')<20:
        auc = General.calc_AUC_v2(dot, shape)
        if auc > 0.8:
            print(tid, auc, dot)

###### snoRNA_ENST00000579879.1

seq = RNA_seq["snoRNA_ENST00000408314.1"]
shape = RNA_shape["snoRNA_ENST00000408314.1"]

dot = seq_dot[seq]; print(dot)
AUC = General.calc_AUC_v2(dot, shape)

print( Visual.Plot_RNAStructure_Shape(seq, dot, shape, mode='label', title='SNORA3B (AUC=%.3f)'%(AUC, )) )






















