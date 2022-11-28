from codecs import EncodedFile
from lib2to3.pgen2.pgen import generate_grammar
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from matplotlib.pyplot import plot
import pandas as pd
import numpy as np
import os
import collections
from Bio.Phylo.PAML import codeml

inf0 =  "cytb_Leucocytozoon.fst"
ali0 =  "aligned.fasta"
ali1 = "aligned_reverse.fasta"
dist0 = 'dist.txt'
dist1 = 'upper_dist.txt'
 
def fasDict(file):
    seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(file, "fasta")}
    return seq_dict
 
def translate(seq): 
    #Standard Translation table
    dna_table = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L' # correct
                ,'TCT':'S','TCC':'S','TCA':'S','TCG':'S' # correct
                ,'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*' # correct
                ,'TGT':'C','TGC':'C','TGA':'*','TGG':'W' # correct
                ,'CTT':'L','CTC':'L','CTA':'L','CTG':'L' # correct
                ,'CCT':'P','CCC':'P','CCA':'P','CCG':'P' # correct
                ,'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q' # correct
                ,'CGT':'R','CGC':'R','CGA':'R','CGG':'R' # correct
                ,'ATT':'I','ATC':'I','ATA':'I','ATG':'M' # correct
                ,'ACT':'T','ACC':'T','ACA':'T','ACG':'T' # correct
                ,'AAT':'N','AAC':'N','AAA':'K','AAG':'K' # correct
                ,'AGT':'S','AGC':'S','AGA':'R','AGG':'R' # correct
                ,'GTT':'V','GTC':'V','GTA':'V','GTG':'V' # correct
                ,'GCT':'A','GCC':'A','GCA':'A','GCG':'A' # correct
                ,'GAT':'D','GAC':'D','GAA':'E','GAG':'E' # correct
                ,'GGT':'G','GGC':'G','GGA':'G','GGG':'G'} # correct
 
    amino_acid = ''
    first = 0
 
    while first < len(seq):
        third = first + 3
        codon = seq[first:third]
        if 'N' in codon or 'X' in codon:
            amino_acid += '?'
            first += 3
        elif codon in dna_table:
            amino_acid += dna_table[codon]
            first += 3
        else:
            amino_acid += '?'
            first += 3
 
    return amino_acid
 
def align(file):  
    in_file = file
    mafft_cline = MafftCommandline('mafft', input=in_file)
    stdout, stderr = mafft_cline()
    with open("aligned.fasta", "w") as handle:
        handle.write(stdout)
    align = AlignIO.read("aligned.fasta", "fasta")

def outputTranslate(inf0):  
    inf = open(inf0, "r")
    seqs = fasDict(inf)
    faa_file0 = "cytb_Leucocytozoon.faa"
    faa_file = open(faa_file0, "w")
    for header, seq in seqs.items():
        aa = translate(seq)
        #print (">" + header)
        print (">" + header, file = faa_file)
        count = 0
        while count <= len(aa):
            #print (aa[count:count + 50])
            print (aa[count:count + 50], file= faa_file)
            count += 50
    faa_file.close()
    inf.close()
    align(faa_file0)

def outputRevFile(ali1):
    outputTranslate(inf0)
    nucl_seqs = fasDict(inf0)
    aa_seqs = fasDict(ali0)
    oP = open(ali1, "w")
    for key, values in nucl_seqs.items():
        n_seq = nucl_seqs[key]
        if key in aa_seqs.keys():
            aligned = str(aa_seqs[key])
            raw_aa = aligned.replace('-', '')
            count = 0
            nucl_aligned = ''
        for letter in aligned:
            if letter == '-':
                nucl_aligned += '---'
            elif letter == '?':
                nucl_aligned += '---'
                count += 1
            elif letter == raw_aa[count]:
                nucl = count * 3
                nucl_aligned += n_seq[nucl:nucl + 3]
                count += 1
        print('>' + key, file = oP)
        print(nucl_aligned, file = oP)
        
def tree(ali1):
    mafftAlignment = AlignIO.read(ali1, "fasta")
    AlignIO.write(mafftAlignment,ali1+".phy","phylip-relaxed")
    phymlcmd = PhymlCommandline("phyml",input=ali1+".phy", model = "GTR")
    phymltextout,phymltexterr = phymlcmd()   
    tree_filename = ali1 + ".phy_phyml_tree.txt"
    tree = Phylo.read(tree_filename, "newick")
    print("Ascii tree:\n")
    Phylo.draw_ascii(tree)

def dist(ali1, dist0):
    aln = AlignIO.read(open(ali1+".phy"), 'phylip')
    calculator = DistanceCalculator('blastn')
    dm = calculator.get_distance(aln)
    f = open(dist0, "w")
    print(dm, file = f)
    #print (np.max(dm))
    f.close()
    
def readDist(dist0, dist1):
    cmd = "awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }' " + dist0 + ' > ' + dist1
    os.system(cmd)
    myfile = pd.read_table(dist1, index_col=0)
    myrevfile = myfile.iloc[::-1]
    myrevfile.replace(0, np.NaN, inplace=True)
    a, b = myrevfile.stack().idxmin()
    y, z = myrevfile.stack().idxmax()
    matrix = myrevfile.to_numpy()
    max = np.nanmax(matrix)
    min = np.nanmin(matrix)
    mean = np.nanmean(matrix, dtype=np.float64)
    #print (np.nanmax(matrix))
    #print(myrevfile.loc[[a], [b]])
    #print(myrevfile.loc[[y], [z]])
    return a,b,y,z, max, min, mean

def gettingFeatures():
    minFile = open("min_dist.fa", "w")
    maxFile = open("max_dist.fa", "w")
    a,b,y,z, max, min , mean = readDist(dist0, dist1)
    seqs = fasDict(ali1)
    print("the average distance among the given sequences is", mean)
    print ("the minimun distance is", min,  "between sequences" , a, "and", b) 
    print ("the maximum distance is", max,  "between sequences" , y, "and", z) 
    for key, value in  seqs.items():
        if a in key:
            print(">" + str(key), file = minFile)
            print(value, file = minFile)
        elif b in key:
            print(">" + str(key), file = minFile)
            print(value, file = minFile)
        elif y in key:
            print(">" + str(key), file = maxFile)
            print(value, file = maxFile)
        elif z in key:
            print(">" + str(key), file = maxFile)
            print(value, file = maxFile)
    minFile.close()
    maxFile.close()   

def ENc(ali1):
    print("gene,Nc,number_of_codons")
    with open(ali1, 'r') as f:
        head = ''
        mean = 0
        count = 0
        for line in f:
            if line[0] == ">":
                if len(head) != 0 and len(seq) != 0:
                    k,n,fcf,tc = prep_counts_for_Nc(seq)
                    Nc = calc_Nc(k, n, fcf)
                    mean += Nc
                    count += 1
                    print(head[1:],Nc,tc,sep=',',flush=True)
                head = line.rstrip()
                seq = ''
            else:
                seq += line.rstrip().upper()
        if len(head) != 0 and len(seq) != 0:
            k,n,fcf,tc = prep_counts_for_Nc(seq)
            Nc = calc_Nc(k, n, fcf)
            mean += Nc
            count += 1
            print(head[1:],Nc,tc,sep=',',flush=True)
        print(mean/count)

def prep_counts_for_Nc(codseq):
    codon2CF = {
    # 1 fold
    'ATG' : 'M10', # Met
    'TGG' : 'W10', # Trp

    # 2 fold
    'TTT' : 'F21', 'TTC' : 'F22', # Phe
    'TAT' : 'Y21', 'TAC' : 'Y22', # Tyr
    'CAT' : 'H21', 'CAC' : 'H22', # His
    'CAA' : 'Q21', 'CAG' : 'Q22', # Gln
    'AAT' : 'N21', 'AAC' : 'N22', # Asn
    'AAA' : 'K21', 'AAG' : 'K22', # Lys
    'GAT' : 'D21', 'GAC' : 'D22', # Asp
    'GAA' : 'E21', 'GAG' : 'E22', # Glu
    'TGT' : 'C21', 'TGC' : 'C22', # Cys
    'TTA' : 'L21', 'TTG' : 'L22', # Leu
    'AGT' : 'S21', 'AGC' : 'S22', # Ser
    'AGA' : 'R21', 'AGG' : 'R22', # Arg

    #3 fold
    'ATT' : 'I31', 'ATC' : 'I32', 'ATA' : 'I33', # Ile

    #4 fold
    'GTT' : 'V41', 'GTC' : 'V42', 'GTA' : 'V43', 'GTG' : 'V44', # Val
    'CCT' : 'P41', 'CCC' : 'P42', 'CCA' : 'P43', 'CCG' : 'P44', # Pro
    'ACT' : 'T41', 'ACC' : 'T42', 'ACA' : 'T43', 'ACG' : 'T44', # Thr
    'GCT' : 'A41', 'GCC' : 'A42', 'GCA' : 'A43', 'GCG' : 'A44', # Ala
    'GGT' : 'G41', 'GGC' : 'G42', 'GGA' : 'G43', 'GGG' : 'G44', # Gly
    'CTT' : 'L41', 'CTC' : 'L42', 'CTA' : 'L43', 'CTG' : 'L44', # Leu
    'TCT' : 'S41', 'TCC' : 'S42', 'TCA' : 'S43', 'TCG' : 'S44', # Ser
    'CGT' : 'R41', 'CGC' : 'R42', 'CGA' : 'R43', 'CGG' : 'R44', # Arg
    }
    count_codons = collections.Counter([ codon2CF[codseq[i:(i+3)]] for i in range(0,len(codseq),3) if codon2CF.get(codseq[i:(i+3)]) ])
    m = collections.Counter([ cf[0:-1] for cf in count_codons.keys()])
    n = { cf : 0 for cf in m.keys() }
    for cf in count_codons.keys():
        n[cf[0:-1]] += count_codons[cf]
    total_codons = sum([ n[cf] for cf in n.keys() ])
    fcf = { cf : 0 for cf in n.keys() }
    for cf in count_codons.keys():
        fcf[cf[0:-1]] += ((count_codons[cf] + 1)/(n[cf[0:-1]] + m[cf[0:-1]])) ** 2

    k = collections.Counter([ cf[1] for cf in m.keys()])

    k_fcf_count = { fold : 0 for fold in k.keys() }
    k_n_count = { fold : 0 for fold in k.keys() }

    for cf in m.keys():
        k_n_count[cf[1]] += n[cf]
        k_fcf_count[cf[1]] += (n[cf]*fcf[cf])

    return k, k_n_count, k_fcf_count, total_codons

def calc_Nc(k, k_n_count, k_fcf_count):
    Nc = k['1']
    for fold in k.keys():
        if fold != '1':
            Nc += (k[fold] * k_n_count[fold])/k_fcf_count[fold]
    return Nc

def gaps(ali1):
    with open("nogaps.txt", "w") as o:
        for record in AlignIO.read(ali1, "fasta"):
            record.seq = record.seq.ungap("-")
            SeqIO.write(record, o, "fasta")

def cut_gap_in_blocks(ali1):
    align = AlignIO.read(ali1, "fasta")
    edited = align[:, :11] #+ align[:, 12:]
    print(edited)


# outputRevFile(ali1)
# tree(ali1)
# dist(ali1, dist0)
# gettingFeatures()
# ENc(ali1)
cut_gap_in_blocks(ali1)
# 1. dn/ds ratio per each aligment
# 2. hyphy?
# 3. Codon bias plot
# 4. ENc
# 5. distribution codon vs dn/ds ratio (genes)
# 6. Gene Class (function, aligment, revTrans, tree build, dn/ds, )
# 7. Tree bulding
# 8. alimengt
