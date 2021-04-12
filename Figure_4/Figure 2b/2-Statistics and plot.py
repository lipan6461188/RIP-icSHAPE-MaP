
###############################
### 计算三级结构中的距离
###############################

importCommon()
import DicerMiRNA

root = "/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/"

constraint_dot = General.load_dot(root+"constraint.dot")
unconstraint_dot = General.load_dot(root+"unconstraint.dot")

constraint_loc = DicerMiRNA.locate_rosetta_output(root+"Rosetta-constraint")
#unconstraint_loc = DicerMiRNA.locate_rosetta_output(root+"Rosetta-unconstraint")

DicerMiRNA.Rosseta_check_str(constraint_loc, {tid:constraint_dot[tid][1] for tid in constraint_dot})
#DicerMiRNA.Rosseta_check_str(unconstraint_loc, {tid:unconstraint_dot[tid][1] for tid in unconstraint_dot})


def trim_MiRBase_dot(MiRBase):
    import re
    
    miRNA_splite_site = {}
    for miRID in MiRBase:
        if not miRID.startswith('hsa'):
            continue
        mir_seq = MiRBase[miRID][0]
        
        glindex = re.search("[AUCG]", mir_seq).start()+1
        grindex = re.search("[AUCG][aucg]*$", mir_seq).start()+1
        premiRNA = mir_seq[glindex-1:grindex]
        if grindex-glindex<50:
            continue
        
        m5r_index = re.search("[aucg]", premiRNA).start()
        m3l_index = re.search("[aucg][AUCG]*$", premiRNA).start()+2
        mature_5miRNA = premiRNA[:m5r_index]
        mature_3miRNA = premiRNA[m3l_index-1:]
        
        #print(">"+mir_seq)
        #print(premiRNA)
        #print(mature_5miRNA, mature_3miRNA)
        miRNA_splite_site[miRID] = [1, m5r_index, m3l_index, len(premiRNA)]
    
    return miRNA_splite_site


MiRBase = General.load_dot("/150T/zhangqf/GenomeAnnotation/miRNA/more/miRNA.dot"); print(len(MiRBase))
miRNA_splite_site = trim_MiRBase_dot(MiRBase)

distance_df = []
for tid in constraint_loc:
    PDB_file_list = constraint_loc[tid]
    if len(PDB_file_list)==0:
        continue
    print ("Process", tid, "...")
    l5,r5,l3,r3 = miRNA_splite_site[tid]
    seq, matrix_list = DicerMiRNA.batch_readPDB(PDB_file_list)
    Len_5 = r5
    Len_3 = r3-l3+1
    for i,matrix in enumerate(matrix_list):
        d5 = matrix.iloc[0,r5-1]
        d5_minus1 = matrix.iloc[0,r5-2]
        d5_plus1 = matrix.iloc[0,r5]
        d3 = matrix.iloc[l3-1,r3-1]
        d3_minus1 = matrix.iloc[l3,r3-1]
        d3_plus1 = matrix.iloc[l3-2,r3-1]
        distance_df.append( [tid, i+1, Len_5, '5p', 'center', d5 ] )
        distance_df.append( [tid, i+1, Len_5, '5p', 'minus1', d5_minus1 ] )
        distance_df.append( [tid, i+1, Len_5, '5p', 'plus1', d5_plus1 ] )
        distance_df.append( [tid, i+1, Len_3, '3p', 'center', d3 ] )
        distance_df.append( [tid, i+1, Len_3, '3p', 'minus1', d3_minus1 ] )
        distance_df.append( [tid, i+1, Len_3, '3p', 'plus1', d3_plus1 ] )

distance_df = pd.DataFrame(distance_df, columns=['miRNA','rep','CenterLength','5or3','CenterType','distance'])
distance_df.to_csv("/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/miRNA_data.csv")

#############################
##### 读取上面的数据
#############################

distance_df = pd.read_csv('/Share2/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/miRNA_data.csv', index_col=0)

#############################
##### 选出variation不太大的miRNAs
#############################

all_tids = {21:[], 22:[], 23:[]}
approved_tids = {21:[], 22:[], 23:[]}
for nucnum in range(21,24):
    discard = 0
    distance_df_2x = distance_df.loc[(distance_df['CenterLength']==nucnum)&(distance_df['CenterType']=="center")&(distance_df['5or3']=="3p"), ]
    distlist_2x = np.array(distance_df_2x['distance'].tolist())
    lower_qtile = np.quantile(distlist_2x, 0.25)
    upper_qtile = np.quantile(distlist_2x, 0.75)
    M = np.median(distlist_2x)
    lower_bound = M-(M-lower_qtile)*2
    upper_bound = M+(upper_qtile-M)*2
    
    gq = np.quantile(distlist_2x, 0.75) - np.quantile(distlist_2x, 0.25)
    for tid in set(distance_df_2x['miRNA'].tolist()):
        distance = distance_df_2x.loc[distance_df_2x['miRNA']==tid, 'distance'].tolist()
        all_tids[nucnum].append(tid)
        if lower_bound<np.median(distance)<upper_bound:
            sq = np.quantile(distance, 0.75) - np.quantile(distance, 0.25)
            if gq>sq:
                approved_tids[nucnum].append(tid)
            else:
                print('discard ', tid, 'self large variation')
                discard += 1
        else:
            print('discard ', tid, 'global large variation')
            discard += 1
    
    print("nucnum %d, discard %d, left %d" % (nucnum, discard, len(approved_tids[nucnum])))

boxDataCenter_dict = {}
Num = {}
for nucnum in range(21,24):
    sub_df = distance_df.loc[(distance_df.CenterLength==nucnum)&(distance_df['5or3']=="3p"),:]
    Num[nucnum] = len(approved_tids[nucnum])
    boxDataCenter_dict[nucnum] = []
    for miRNA in approved_tids[nucnum]:
        center = sub_df.loc[(sub_df.miRNA==miRNA)&(sub_df.CenterType=="center"),"distance"].tolist()
        boxDataCenter_dict[nucnum].append( [center, miRNA] )

for nucnum in boxDataCenter_dict:
    boxDataCenter_dict[nucnum].sort(key=lambda x: np.median(x[0]))

boxDataCenter =     [it[0] for it in boxDataCenter_dict[21]]+[it[0] for it in boxDataCenter_dict[22]]+[it[0] for it in boxDataCenter_dict[23]]
boxLabelsCenter =   [it[1] for it in boxDataCenter_dict[21]]+[it[1] for it in boxDataCenter_dict[22]]+[it[1] for it in boxDataCenter_dict[23]]


###############################
### 数出那些成熟miRNA比多切一位或者少切一位的miRNA距离更加接近60的miRNA
###############################


miRNA_info_5 = {}
miRNA_info_3 = {}
for row_data in distance_df.values:
    miRID, pos_s, pos_e, p5or3, aligntype, distance = row_data
    if p5or3 == '5p':
        if miRID not in miRNA_info_5:
            miRNA_info_5[miRID] = [0,0,0]
        if aligntype=='center':
            miRNA_info_5[miRID][1] = distance
        elif aligntype=='minus1':
            miRNA_info_5[miRID][0] = distance
        else:
            miRNA_info_5[miRID][2] = distance
    elif p5or3 == '3p':
        if miRID not in miRNA_info_3:
            miRNA_info_3[miRID] = [0,0,0]
        if aligntype=='center':
            miRNA_info_3[miRID][1] = distance
        elif aligntype=='minus1':
            miRNA_info_3[miRID][0] = distance
        else:
            miRNA_info_3[miRID][2] = distance

approved_tids_list = approved_tids[21] + approved_tids[22] + approved_tids[23]
print(len(approved_tids_list))

count_5 = 0
for tid in approved_tids_list:
    m1,c,p1 = miRNA_info_5[tid]
    if abs(c-60)<abs(m1-60) and abs(c-60)<abs(p1-60):
        #print(abs(c-60), abs(m1-60), abs(p1-60))
        count_5 += 1

count_3 = 0
for tid in approved_tids_list:
    m1,c,p1 = miRNA_info_3[tid]
    if abs(c-60)<abs(m1-60) and abs(c-60)<abs(p1-60):
        count_3 += 1

print(count_5, count_3)



############# Center

medians_21 = np.median([ np.median(boxDataCenter[i]) for i in range(Num[21]) ])
medians_22 = np.median([ np.median(boxDataCenter[i+Num[21]]) for i in range(Num[22]) ])
medians_23 = np.median([ np.median(boxDataCenter[i+Num[21]+Num[22]]) for i in range(Num[23]) ])

plt.rcParams['font.family'] = "sans-serif"
fig = plt.figure(1, figsize=(20, 6))
ax = fig.add_subplot(111)
facecolors = [Colors.RGB['pink']]*Num[21]+[Colors.RGB['blue']]*Num[22]+[Colors.RGB['green']]*Num[23]
Figures.boxPlot(ax, boxDataCenter, labels=boxLabelsCenter, facecolors=facecolors)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

r1 = Num[21]/(Num[21]+Num[22]+Num[23])
r2 = (Num[21]+Num[22])/(Num[21]+Num[22]+Num[23])
ax.axhline(y=medians_21, xmin=0, xmax=r1, linewidth=1, color=Colors.RGB['pink'])
ax.axhline(y=medians_22, xmin=r1, xmax=r2, linewidth=1, color=Colors.RGB['blue'])
ax.axhline(y=medians_23, xmin=r2, xmax=1, linewidth=1, color=Colors.RGB['green'])

ax.set_ylim(0, 75)
ax.set_title("Center")
fig.tight_layout()
fig.savefig("figs/distance_dir/distance_center_2.pdf")
plt.close()

## 显著性

medians_21s = [ np.median(boxDataCenter[i]) for i in range(Num[21]) ]
medians_22s = [ np.median(boxDataCenter[i+Num[21]]) for i in range(Num[22]) ]
medians_23s = [ np.median(boxDataCenter[i+Num[21]+Num[22]]) for i in range(Num[23]) ]

p12 = scipy.stats.ttest_ind(medians_21s, medians_22s)[1]
p13 = scipy.stats.ttest_ind(medians_21s, medians_23s)[1]
p23 = scipy.stats.ttest_ind(medians_22s, medians_23s)[1]

rejected, corected_pvalue = multicomp.fdrcorrection0([p12, p23, p13])


##### Source 
OUT = open(join(HOME, 'figs/source_data/Figure_4b.txt'), 'w')

length_list = [21]*Num[21] + [22]*Num[22] + [23]*Num[23]
for tid,dist_array, length in zip(boxLabelsCenter, boxDataCenter, length_list):
    for dist in dist_array:
        print( tid, length, dist, file=OUT, sep='\t' )

OUT.close()



