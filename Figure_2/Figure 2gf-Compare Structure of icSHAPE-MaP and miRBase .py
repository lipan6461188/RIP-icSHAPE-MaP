importCommon()

miRBase = General.load_dot("/150T/zhangqf/GenomeAnnotation/miRNA/more/miRNA.dot")
miRBase = {k:miRBase[k] for k in miRBase if k.startswith('hsa')}

def read_shape(inFolder):
    SHAPE = {}; Sequence = {}
    inFolder = inFolder.rstrip('/')+'/'
    files = os.listdir(inFolder)
    for file in files:
        if file.endswith(".shape"):
            tid = file.rstrip('.shape')
            data = General.load_SHAPEMap(inFolder+file)
            shape = data['shape_pro_list']
            seq = data['seq']
            #if shape.count('NULL')/len(shape)>0.4: continue
            SHAPE[tid] = shape
            Sequence[tid] = seq
    return SHAPE, Sequence

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

trimed_miRNA = trim_MiRBase_miRNA(miRBase, rePred=False)

seq_based_trimed_miRNA = { trimed_miRNA[k][0].upper():(k, trimed_miRNA[k][1]) for k in trimed_miRNA }

shape_folder = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/Final-Dicer-Run/Out/shape_files"
SHAPE, Sequence = read_shape(shape_folder); print (len(SHAPE))

Data = []
for tid in Sequence:
    if not tid.startswith('miRNA'):
        continue
    seq = Sequence[tid]
    for spliced_seq in seq_based_trimed_miRNA:
        miRBase_name,miRBase_dot = seq_based_trimed_miRNA[spliced_seq]
        if spliced_seq in seq:
            index = seq.find(spliced_seq)
            spliced_shape = SHAPE[tid][index:index+len(spliced_seq)]
            pred_dot = Structure.predict_structure(spliced_seq, spliced_shape, si=-0.4, sm=2.0) ######### 注意这里的参数
            print (">",miRBase_name,tid)
            print (spliced_seq)
            print (miRBase_dot)
            print ( Colors.color_SHAPE(spliced_shape) )
            print (pred_dot)
            Data.append( [tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot] )
            if miRBase_name in ('hsa-mir-24-2', 'hsa-mir-29c'):
                auc = General.calc_AUC_v2(miRBase_dot, spliced_shape)
                title = f"{miRBase_name} AUC={auc:.3f} miRBase"
                cmd_miRBase = Visual.Plot_RNAStructure_Shape(spliced_seq, miRBase_dot, spliced_shape, title=title)
                auc = General.calc_AUC_v2(pred_dot, spliced_shape)
                title = f"{miRBase_name} AUC={auc:.3f} Pred"
                cmd_Pred = Visual.Plot_RNAStructure_Shape(spliced_seq, pred_dot, spliced_shape, title=title)
                print(cmd_miRBase)
                print(cmd_Pred)
                print("")


####################################
# Figure 2e and S3h
####################################

targets = ['miRNA_MI0000440', 'miRNA_MI0000469', 'miRNA_MI0000073', 'miRNA_MI0005559']
for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    if tid in targets:
        pred_auc = General.calc_AUC_v2(pred_dot, spliced_shape)
        miRBase_auc = General.calc_AUC_v2(miRBase_dot, spliced_shape)
        pred_energy = Structure.estimate_energy(spliced_seq, pred_dot, spliced_shape)
        miRBase_energy = Structure.estimate_energy(spliced_seq, miRBase_dot, spliced_shape)
        title = "%s %s pred_auc=%.3f miRBase_auc=%.3f pred_energy=%.1f miRBase_energy=%.1f" % (tid, miRBase_name, pred_auc, miRBase_auc, pred_energy, miRBase_energy)
        print( Visual.Plot_RNAStructure_Shape(spliced_seq, pred_dot, spliced_shape, title=title+" pred", wait=False) )
        print( Visual.Plot_RNAStructure_Shape(spliced_seq, miRBase_dot, spliced_shape, title=title+" miRBase", wait=False) )
        print("")

"""
nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN AGAGCUUAGCUGAUUGGUGAACAGUGAUUGGUUUCCGCUUUGUUCACAGUGGCUAAGUUCUGC -structureDBN "((((((((((((....(((((((((((......).)))..)))))))..)))))))))))).." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "12,14,16,27,29,30,31,32,39,40" -applyBasesStyle2on "13,15,17,25,28" -applyBasesStyle3on "26,33,34,35,37,41,49,50" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,18,19,20,21,22,23,24,36,38,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,59,60,61,62,63"  -spaceBetweenBases "0.8" -title "miRNA_MI0000440 hsa-mir-27b pred_auc=0.859 miRBase_auc=0.806 pred_energy=-59.8 miRBase_energy=-55.9 pred" &
nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN AGAGCUUAGCUGAUUGGUGAACAGUGAUUGGUUUCCGCUUUGUUCACAGUGGCUAAGUUCUGC -structureDBN "((((((((((((....((((((((....(((...)))..))))))))..)))))))))))).." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "12,14,16,27,29,30,31,32,39,40" -applyBasesStyle2on "13,15,17,25,28" -applyBasesStyle3on "26,33,34,35,37,41,49,50" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,18,19,20,21,22,23,24,36,38,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,59,60,61,62,63"  -spaceBetweenBases "0.8" -title "miRNA_MI0000440 hsa-mir-27b pred_auc=0.859 miRBase_auc=0.806 pred_energy=-59.8 miRBase_energy=-55.9 miRBase" &

nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN UCCCUGAGACCCUUUAACCUGUGAGGACAUCCAGGGUCACAGGUGAGGUUCUUGGGAGCC -structureDBN "((((.((((.((((..((((((((............)))))))))))).))))))))..." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "16,25,26,27,29,31,33,34,35" -applyBasesStyle2on "15,17,32" -applyBasesStyle3on "24,28,30,36" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60"  -spaceBetweenBases "0.8" -title "miRNA_MI0000469 hsa-mir-125a pred_auc=0.855 miRBase_auc=0.669 pred_energy=-57.6 miRBase_energy=-44.6 pred" &
nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN UCCCUGAGACCCUUUAACCUGUGAGGACAUCCAGGGUCACAGGUGAGGUUCUUGGGAGCC -structureDBN "..((..(((.((((..((((((((((....))....)))))))))))).)))..))...." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "16,25,26,27,29,31,33,34,35" -applyBasesStyle2on "15,17,32" -applyBasesStyle3on "24,28,30,36" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60"  -spaceBetweenBases "0.8" -title "miRNA_MI0000469 hsa-mir-125a pred_auc=0.855 miRBase_auc=0.669 pred_energy=-57.6 miRBase_energy=-44.6 miRBase" &

nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN AGUUUUGCAUAGUUGCACUACAAGAAGAAUGUAGUUGUGCAAAUCUAUGCAAAACUGA -structureDBN "((((((((((((((((((.((............)).))))))..)))))))))))).." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "24,27,29,30,31,32,33" -applyBasesStyle2on "23,42" -applyBasesStyle3on "25" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,26,28,34,35,36,38,39,40,41,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58" -applyBasesStyle5on "12,37"  -spaceBetweenBases "0.8" -title "miRNA_MI0000073 hsa-mir-19a pred_auc=0.878 miRBase_auc=0.776 pred_energy=-56.7 miRBase_energy=-52.6 pred" &
nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN AGUUUUGCAUAGUUGCACUACAAGAAGAAUGUAGUUGUGCAAAUCUAUGCAAAACUGA -structureDBN "((((((((((((((((((((((.......))))...))))))..)))))))))))).." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "24,27,29,30,31,32,33" -applyBasesStyle2on "23,42" -applyBasesStyle3on "25" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,26,28,34,35,36,38,39,40,41,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58" -applyBasesStyle5on "12,37"  -spaceBetweenBases "0.8" -title "miRNA_MI0000073 hsa-mir-19a pred_auc=0.878 miRBase_auc=0.776 pred_energy=-56.7 miRBase_energy=-52.6 miRBase" &

nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN UGCGGGGCUAGGGCUAACAGCAGUCUUACUGAAGGUUUCCUGGAAACCACGCACAUGCUGUUGCCACUAACCUCAACCU -structureDBN "...((((.((((((.(((((((.((.....)).((((((...)))))).......)))))))))).))).))))....." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "28,31,33,41,42,43,51,53,55,66" -applyBasesStyle2on "15,26,27,32,40,49" -applyBasesStyle3on "29,30,50" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,34,35,36,37,38,39,44,45,46,47,48,52,54,56,57,58,59,60,61,62,63,64,65,67,68,69,70,71,72,73,74,75,76,77,78" -applyBasesStyle5on "79"  -spaceBetweenBases "0.8" -title "miRNA_MI0005559 hsa-mir-744 pred_auc=0.758 miRBase_auc=0.748 pred_energy=-64.4 miRBase_energy=-45.9 pred" &
nohup java -cp /Users/lee/Documents/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN UGCGGGGCUAGGGCUAACAGCAGUCUUACUGAAGGUUUCCUGGAAACCACGCACAUGCUGUUGCCACUAACCUCAACCU -structureDBN "(..((((.((((((.(((((((.................................)))))))))).))).)))).)..." -drawBackbone false -bpStyle simple -basesStyle1 "label=#B61D22" -basesStyle2 "label=#ED9616" -basesStyle3 "label=#194399" -basesStyle4 "label=#040000" -basesStyle5 "label=#828282" -applyBasesStyle1on "28,31,33,41,42,43,51,53,55,66" -applyBasesStyle2on "15,26,27,32,40,49" -applyBasesStyle3on "29,30,50" -applyBasesStyle4on "1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,34,35,36,37,38,39,44,45,46,47,48,52,54,56,57,58,59,60,61,62,63,64,65,67,68,69,70,71,72,73,74,75,76,77,78" -applyBasesStyle5on "79"  -spaceBetweenBases "0.8" -title "miRNA_MI0005559 hsa-mir-744 pred_auc=0.758 miRBase_auc=0.748 pred_energy=-64.4 miRBase_energy=-45.9 miRBase" &
"""

#############################
#####   Plot Heatmap： 每一个miRNA一行，比较预测的结构和miRBase的结构之间的差别
#############################

heatmap_content = []
for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    pred_bpmap = Structure.dot2bpmap(pred_dot)
    miRBase_bpmap = Structure.dot2bpmap(miRBase_dot)
    tmp_list = []
    for i in range(1, len(pred_dot)+1):
        if pred_bpmap.get(i,0)==0 and miRBase_bpmap.get(i,0)!=0:
            # 单链 -> 双链
            tmp_list.append(1)
        elif pred_bpmap.get(i,0)!=0 and miRBase_bpmap.get(i,0)==0:
            # 双链 -> 单链
            tmp_list.append(2)
        elif pred_bpmap.get(i,0)!=miRBase_bpmap.get(i,0):
            # 双链 -> 双链
            tmp_list.append(3)
        else:
            tmp_list.append(0)
    stems = Structure.find_stem(pred_dot, max_stem_gap=3, min_stem_len=2)
    #lm = pred_dot.rfind("(")
    #rm = pred_dot.find(")")
    #center = (lm+rm)//2
    center = (stems[-1][1]+stems[-1][2])//2
    tmp_list = [0]*10+tmp_list+[0]*10
    auc1 = General.calc_AUC_v2(pred_dot, spliced_shape)
    auc2 = General.calc_AUC_v2(miRBase_dot, spliced_shape)
    heatmap_content.append(tmp_list[center+10-30:center+10+30] + [auc1, auc2])

heatmap_df = pd.DataFrame(heatmap_content, index=[d[1] for d in Data])
heatmap_df['sum'] = np.sum(heatmap_df.iloc[:, :60]!=0,axis=1)
heatmap_df = heatmap_df.sort_values(by=['sum'], ascending=False)
diff = heatmap_df['sum'].values
heatmap_df = heatmap_df.drop('sum', axis=1)

plt.figure(figsize=(8,25))
sns.heatmap(heatmap_df.iloc[:,:60])
plt.xticks(range(0,61,10),range(-30,31,10))
plt.savefig(join(HOME, 'figs/miRNA_heatmap.pdf'))
plt.close()

#### 显示有结构差异的碱基数
diff_num = np.array([sum(d!=0) for d in heatmap_df.iloc[:,:60].values])
np.sum(diff_num>0)
np.sum(diff_num>=5)


miR_loopLen = []
Pred_loopLen = []
for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    miR_stems = Structure.find_stem(miRBase_dot, max_stem_gap=3, min_stem_len=2)
    Pred_stems = Structure.find_stem(pred_dot, max_stem_gap=3, min_stem_len=2)
    miR_loopLen.append( miR_stems[-1][2]-miR_stems[-1][1]-1 )
    Pred_loopLen.append( Pred_stems[-1][2]-Pred_stems[-1][1]-1 )


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
data = [ miR_loopLen, Pred_loopLen ]
colors = [Colors.RGB['blue'], Colors.RGB['green']]
Figures.violinPlot(axs, data, ['miRBase', 'Preded'],colors=colors)
_, p = scipy.stats.ttest_rel(miR_loopLen, Pred_loopLen)
axs.set_title(f"paired t-test={p:.4f}")
fig.tight_layout()
plt.savefig(join(HOME, 'figs/miRNA_violin.pdf'))
fig.show()


#### 显示有结构差异的碱基数
diff_num = np.array([sum(d!=0) for d in heatmap_df.iloc[:,:60].values])
np.sum(diff_num>0)
np.sum(diff_num>=5)

### Figure 2g

miR_loopLen = []
Pred_loopLen = []
for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    miR_stems = Structure.find_stem(miRBase_dot, max_stem_gap=3, min_stem_len=2)
    Pred_stems = Structure.find_stem(pred_dot, max_stem_gap=3, min_stem_len=2)
    miR_loopLen.append( miR_stems[-1][2]-miR_stems[-1][1]-1 )
    Pred_loopLen.append( Pred_stems[-1][2]-Pred_stems[-1][1]-1 )


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
data = [ miR_loopLen, Pred_loopLen ]
colors = [Colors.RGB['blue'], Colors.RGB['green']]
Figures.violinPlot(axs, data, ['miRBase', 'Preded'],colors=colors)
_, p = scipy.stats.ttest_rel(miR_loopLen, Pred_loopLen)
axs.set_title(f"paired t-test={p:.4f}")
fig.tight_layout()
plt.savefig(join(HOME, 'figs/miRNA_violin.pdf'))
fig.show()

p((((..((((((.(((((((((.............)))))))))))).)))..))))..

##### Source data Figure 2f

OUT = open(join(HOME, 'figs/source_data/Figure_2f.txt'), 'w')

for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    #pred_bpmap = Structure.dot2bpmap(pred_dot)
    #miRBase_bpmap = Structure.dot2bpmap(miRBase_dot)
    stems = Structure.find_stem(pred_dot, max_stem_gap=3, min_stem_len=2)
    #lm = pred_dot.rfind("(")
    #rm = pred_dot.find(")")
    #center = (lm+rm)//2
    center = (stems[-1][1]+stems[-1][2])//2
    padding_spliced_seq = "p"*10+spliced_seq+"p"*10
    padding_pred_dot = "p"*10+pred_dot+"p"*10
    padding_miRBase_dot = "p"*10+miRBase_dot+"p"*10
    print( miRBase_name, padding_spliced_seq[center+10-30:center+10+30], padding_pred_dot[center+10-30:center+10+30], padding_miRBase_dot[center+10-30:center+10+30], sep="\t", file=OUT )

OUT.close()

##### Source data Figure 2g

OUT = open(join(HOME, 'figs/source_data/Figure_2g.txt'), 'w')

for data in Data:
    tid, miRBase_name, spliced_seq, spliced_shape, pred_dot, miRBase_dot = data
    miR_stems = Structure.find_stem(miRBase_dot, max_stem_gap=3, min_stem_len=2)
    Pred_stems = Structure.find_stem(pred_dot, max_stem_gap=3, min_stem_len=2)
    miR_len = miR_stems[-1][2]-miR_stems[-1][1]-1
    Pred_len = Pred_stems[-1][2]-Pred_stems[-1][1]-1
    print(miRBase_name, Pred_len, miR_len, sep="\t", file=OUT)

OUT.close()



