
import General

def format_sequence_name(inFile):
    namealias = {}
    readSeq = {}
    cur_tid = ""
    for line in open(inFile):
        if line.startswith(">"):
            if line.startswith(">ID="):
                pattern = "ID=(\\w+);Alias=(\\w+);Name=([\\w\\-\\d]+)"
                results = re.findall(pattern, line)
                if len(results) != 1:
                    print "Error: "+line.strip()
                    return
                id1,id2,id3 = results[0]
                cur_tid = "miRNA_%s" % (id1,)
            
            elif line.count('|') == 3:
                tid,gid,gname,gtype = line[1:-1].split('|')
                cur_tid = "%s_%s" % (gtype, tid)
            
            elif line.startswith(">Homo_"):
                tid = line[1:-1].split()[0].replace("Homo_sapiens_tRNA-", "")
                cur_tid = "tRNA_%s" % (tid, )
            
            else:
                tid = line[1:-1].split()[0]
                if tid == "chrM":
                    cur_tid = tid
                else:
                    cur_tid = "rRNA_" + line[1:-1]
            
            namealias[cur_tid] = line[1:-1]
            if cur_tid in readSeq:
                    print "Error: %s duplicate" % (cur_tid, )
                    return
            readSeq[cur_tid] = ""
        else:
            readSeq[cur_tid] += line.strip().upper()
    
    return readSeq, namealias


short_rna, namealias = format_sequence_name("/Share/home/zhangqf7/lipan/reference/miRNA_tRNA_snRNA_snoRNA_rRNA_mtRNA_YRNA/combineRNA/miRNA_tRNA_snRNA_snoRNA_rRNA_mtRNA.uniq.fa")

outFile = "/Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/ref/smallRNA.fa"
General.write_fasta(short_rna, outFile)

"""
cd /Share/home/zhangqf7/lipan/precursor_SHAPEMAP/statisticalReadsDist/ref
icSHAPE-pipe starbuild -i smallRNA.fa -o ./ 
"""

