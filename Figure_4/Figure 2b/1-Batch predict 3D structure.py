
from D3 import Rosetta
import General, os

miRNA_dot = General.load_dot("/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/constraint.dot")
miRNAs = list(miRNA_dot.keys())

job_list = []
Folder = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/6.miRNA_Structure/20190514/Rosetta-constraint/"
for miRID in miRNAs:
    print(miRID)
    seq, dot = miRNA_dot[miRID]
    outFolder = Folder+miRID+'/'
    os.system("mkdir -p "+outFolder)
    helix_setup_job,rna_denovo_job,extract_job = Rosetta.pred_3D_rosetta(seq, dot, miRID, outFolder, gen_modellimit=500, topnum=50, verbose=False, queue="Z-ZQF")
    job_list.append( (miRID, helix_setup_job,rna_denovo_job,extract_job) )

for miRID,j1,j2,j3 in job_list:
    j3.wait()
    print(miRID+" has finished...")

## nohup python run.py &


