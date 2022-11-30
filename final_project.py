from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.SeqUtils import GC, GC_skew
from itertools import combinations
from dnds import dnds, pnps, substitutions, dnds_codon, dnds_codon_pair, syn_sum, translate
import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np
import os
import itertools
import pandas as pd
import remove_gaps
import subprocess
import collections


EXTENSIONS_FASTA = ['.fna']
EXTENSIONS_FAA_aligned   = ['_aligned.fasta']
EXTENSIONS_FST_reversed = ['.fst']
EXTENSIONS_FST_trimmed = ['_reverse_trimmed.fst']


class aligment():
    def get_file_list_fasta(self, root_dir, E):
        file_list = []
        for root, directories, filenames in os.walk(root_dir):
            for filename in filenames:
                if any (ext in filename for ext in E):
                    file_list.append(os.path.join(root, filename))
        file_list.sort()
        return file_list
       
    def fasDict(self, file):
        if file.endswith(".fna"):
            fasta_file = self.get_file_list_fasta("test", EXTENSIONS_FASTA)
            list_dict = []
            for i in fasta_file:
                seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(i, "fasta")}
                list_dict.append(seq_dict)
            return list_dict
            
        elif file.endswith('_aligned.fasta'):
            fasta_file = self.get_file_list_fasta("test", EXTENSIONS_FAA_aligned)
            list_dict = []
            for i in fasta_file:
                seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(i, "fasta")}
                list_dict.append(seq_dict)
            return list_dict
                
        elif file.endswith('_reverse_trimmed.fst'):
            fasta_file = self.get_file_list_fasta("test", EXTENSIONS_FST_trimmed)
            list_dict = []
            for i in fasta_file:
                seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(i, "fasta")}
                list_dict.append(seq_dict)
            return list_dict   
     
    def translation(self,seq): 
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
                    ,'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
                    ,'---' : "-"} # correct
     
        amino_acid = ''
        first = 0
     
        while first < len(seq):
            third = first + 3
            codon0 = seq[first:third]
            if 'N' in codon0 or 'X' in codon0:
                amino_acid += '?'
                first += 3
            elif codon0 in dna_table:
                amino_acid += dna_table[codon0]
                first += 3
            else:
                amino_acid += '?'
                first += 3
     
        return amino_acid
     
    def align(self, file):  
        in_file = file
        mafft_cline = MafftCommandline('mafft', input=in_file)
        stdout, stderr = mafft_cline()
        name = in_file.split(".")[0]
        with open(name + "_aligned.fasta", "w") as handle:
            handle.write(stdout)
        align = AlignIO.read("aligned.fasta", "fasta")

    def outputTranslate(self,file):
        fasta_name = self.get_file_list_fasta("test", EXTENSIONS_FASTA)
        for file in fasta_name:
            data = self.fasDict(file)
        list_name = []
        for i in fasta_name:
            name = i.split("/")[1]
            faa_file_name = name.split(".")[0] + ".faa"
            list_name.append(faa_file_name)

        for seqs, name in zip(data, list_name):
            faa_file = open(r'test/'+ name, "w")
            for header, seq in seqs.items():
                aa = self.translation(seq)
                print (">" + header, file = faa_file)
                count = 0
                while count <= len(aa):
                    print (aa[count:count + 50], file= faa_file)
                    count += 50   
        faa_file.close()
        for i in list_name:
            self.align(r'test/' + i)

    def outputRevFile(self):
        
        fasta_name = self.get_file_list_fasta("test", EXTENSIONS_FASTA)
        
        for i in fasta_name:
            self.outputTranslate(i)
            
        faa_name = self.get_file_list_fasta("test", EXTENSIONS_FAA_aligned)  

        for a, b in zip (fasta_name, faa_name):
            self.outputTranslate(a)
            nucl_seqs = self.fasDict(a)
            aa_seqs = self.fasDict(b)

        for i, name, j in zip(nucl_seqs, fasta_name, aa_seqs):
            file_name = name.split("/")[1]
            file_name1 = file_name.split(".")[0]
            oP = open(r'test/'+ file_name1 + "_reverse.fst", "w")
            for key, values in i.items():
                n_seq = i[key]
                if key in i.keys():
                    aligned = str(j[key])
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
            
class genomicFeatures():

    def gettingFeatures(self):
        x = aligment()
        fasta_name = x.get_file_list_fasta("test", EXTENSIONS_FASTA)
        f = open("test/genomicFeatures.txt", "w")
        #print("ID_seq", "\tGC_content", "\tGC_skew")
        print("ID_seq", "GC_content", "GC_skew", file = f)
        for i in fasta_name:
            alignment = x.fasDict(i)
        
        name_seq = []
        for j in alignment:
            for key, values in sorted(j.items()):
                #print(key, "\t" + str(GC(values)),"\t" + str(GC_skew(values, window=100000)[0]))
                print(key, "\t" + str(GC(values)),"\t" + str(GC_skew(values, window=100000)[0]), file = f)
                name_seq.append(key)
        f.close()
        
        df = pd.read_csv("test/genomicFeatures.txt", sep=" ")
        ids = [*set(name_seq)]
        count = 0
        colors = ["aliceblue", "antiquewhite",  "aqua", "aquamarine",
                    "blueviolet" , "brown", "burlywood", "cadetblue",
                   "chartreuse",  "chocolate", "cornflowerblue", "crimson"]
        
        lists = []
        for i in range(len(ids)):
            count += 1
            name = "df" + str(count)
            lists.append(name)
        # first = lists[0]
        first = df[df["ID_seq"] == ids[0]] 
        # fig = plt.figure(figsize=(10,5))
        ax = first.plot(x='GC_skew', y='GC_content', kind='scatter', c=colors[0], label=ids[0])
        
        del lists[0]
        del ids[0]
        del colors[0]
        
        for i, j, k in zip(lists, ids, colors):
            i = df[df['ID_seq'] == j]
            i.plot(x='GC_skew', y='GC_content', kind='scatter', ax=ax, c=k, label=j)
        
        plot_name="test/GC_content_skew_scatterplot.jpg"
        plt.legend(bbox_to_anchor=(1,1), loc='upper left', prop={'size': 4})
        # plt.legend(bbox_to_anchor=(1,1), loc='upper left')

        
        # l0 = df.ID_seq
        # x0 = df.GC_skew
        # y0 = df.GC_content
        
        # fig = plt.figure(figsize=(10,5))
        # for x, y, z in zip(x0, y0, l0):
            # plt.scatter(x, y, label = z) 
        # plt.xlabel('GC Skew')
        # plt.ylabel('GC Content') 
        # plt.title('hola mundo')
        # plt.legend(bbox_to_anchor=(1,1), loc='upper left')
        # fig.savefig('line plot.jpg', bbox_inches='tight', dpi=150)
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)
        
class codon_usage():
    def ENc(self):
        x = aligment()
        fasta_name = x.get_file_list_fasta("test", EXTENSIONS_FASTA)
        g = open("test/ENc.txt", "w")
        # print("ID_seq", "\tENc", "\tcodon_number")
        print("ID_seq", "ENc", "codon_number", file = g)
        for i in fasta_name:
            with open(i, 'r') as f:
                head = ''
                mean = 0
                count = 0
                for line in f:
                    if line[0] == ">":
                        if len(head) != 0 and len(seq) != 0:
                            k,n,fcf,tc = self.prep_counts_for_Nc(seq)
                            Nc = self.calc_Nc(k, n, fcf)
                            mean += Nc
                            count += 1
                            # print(head[1:],Nc,tc,sep='\t',flush=True)
                            print(head[1:],Nc,tc,sep=' ',flush=True, file = g)
                        head = line.rstrip()
                        seq = ''
                    else:
                        seq += line.rstrip().upper()
                if len(head) != 0 and len(seq) != 0:
                    k,n,fcf,tc = self.prep_counts_for_Nc(seq)
                    Nc = self.calc_Nc(k, n, fcf)
                    mean += Nc
                    count += 1
                    # print(head[1:],Nc,tc,sep='\t',flush=True)
                    print(head[1:],Nc,tc,sep=' ',flush=True, file = g)
                # print(mean/count)
        g.close()

    def prep_counts_for_Nc(self, codseq):
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

    def calc_Nc(self, k, k_n_count, k_fcf_count):
        Nc = k['1']
        for fold in k.keys():
            if fold != '1':
                Nc += (k[fold] * k_n_count[fold])/k_fcf_count[fold]
        return Nc

    def plottingENc_CN(self):
        self.ENc()
        df = pd.read_csv("test/ENc.txt", sep=" ")
        name_seq = df['ID_seq'].tolist()
        ids = [*set(name_seq)]
        
        count = 0
        colors = ["aliceblue", "antiquewhite",  "aqua", "aquamarine",
                    "blueviolet" , "brown", "burlywood", "cadetblue",
                   "chartreuse",  "chocolate", "cornflowerblue", "crimson"] 
        lists = []
        for i in range(len(ids)):
            count += 1
            name = "df" + str(count)
            lists.append(name)
        first = df[df["ID_seq"] == ids[0]] 
        ax = first.plot(y='codon_number', x='ENc', kind='scatter', c=colors[0], label=ids[0])
        del lists[0]
        del ids[0]
        del colors[0]
        
        for i, j, k in zip(lists, ids, colors):
            i = df[df['ID_seq'] == j]
            i.plot(y='codon_number', x='ENc', kind='scatter', ax=ax, c=k, label=j)
        
        plot_name="ENc_codon_number_scatterplot.jpg"
        plt.legend(bbox_to_anchor=(1,1), loc='upper left', prop={'size': 4})
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)

    def plotENc_box(self):
        self.ENc()
        df = pd.read_csv("test/ENc.txt", sep=" ")
        del df['codon_number']
        df1 = df.pivot(columns = "ID_seq", values = "ENc")
        ax = df1.plot.box(rot=90, fontsize = 5)
        plot_name="test/ENc_boxplot.jpg"
        plt.tight_layout()
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)

    def distributionENc(self):
        self.ENc()
        df = pd.read_csv("test/ENc.txt", sep=" ")
        del df['codon_number']
        ax = df.plot.hist(figsize=[10, 8], bins=40, alpha=0.5)
        plot_name="test/ENc_histogram.jpg"
        plt.tight_layout()
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)

class dn_ds_pn_ps():
    def cut_gap_in_blocks(self):
        x = aligment()
        fasta_name = x.get_file_list_fasta("test", EXTENSIONS_FST_reversed)
        for i in fasta_name:
            cmd = f'python remove_gaps.py -c "{i}" --trim_gappy 0'
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            out, err = p.communicate() 

    def dn_ds(self):
        x = aligment()
        fasta_name = x.get_file_list_fasta("test", EXTENSIONS_FST_trimmed)

        h = open("test/dnds.txt", "w")
        print("pairwaise_comparison", "dn/ds", "pn/ps", file = h)
        for a in fasta_name:
            alignment0 = x.fasDict(a)
        
        for i in alignment0:
            keys, values = zip(*i.items())
            xList0 = list(values)
            xList1 = itertools.combinations(xList0, 2)

            yList0 = list(keys)
            yList1 = itertools.combinations(yList0, 2)
            
            seq = []
            name_seq = []
            for j, k in zip(list(xList1),list(yList1)):
                name_seq.append(k)
                seq.append(j)
                
            for l, m in zip(seq, name_seq):
                n = '-'.join(m)
                try:
                    xxx1 = float(pnps(str(l[0]), str(l[1])))
                    xxx2 = dnds(str(l[0]), str(l[1]))
                except ValueError:
                    xxx1 = 0
                    xxx2 = 0
                except ZeroDivisionError:
                    xxx1 = 0
                    xxx2 = 0
                if xxx2 == -0.0:
                    xxx3 = 0
                    print(n, xxx3, int(xxx1), file = h)
                else:
                    print(n, xxx2, xxx1, file = h)
        h.close()
            
    def plotting_dnds(self):
        df = pd.read_csv("test/dnds.txt", sep=" ")
        del df['pn/ps']
        df1 = df.pivot(columns = "pairwaise_comparison", values = "dn/ds")
        ax = df1.plot.box(rot=90, fontsize = 5)
        plot_name="test/dn_ds_boxplot.jpg"
        plt.tight_layout()
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)

    def plotting_pnps(self):
        df = pd.read_csv("test/dnds.txt", sep=" ")
        del df['dn/ds']
        df1 = df.pivot(columns = "pairwaise_comparison", values = "pn/ps")
        ax = df1.plot.box(rot=90, fontsize = 5)
        plot_name="test/pn_ps_boxplot.jpg"
        plt.tight_layout()
        ax.figure.savefig(plot_name,format='jpeg',dpi=150)

def main():
    class_ali = aligment()
    fasta_file = class_ali.get_file_list_fasta("test", EXTENSIONS_FASTA)

    class_ali.outputRevFile()

    class_genFea = genomicFeatures()
    class_genFea.gettingFeatures()
    
    class_ENc = codon_usage()
    class_ENc.distributionENc()
    class_ENc.plottingENc_CN()
    class_ENc.plotENc_box()
    
    class_dnds = dn_ds_pn_ps()
    class_dnds.cut_gap_in_blocks()
    class_dnds.dn_ds()
    class_dnds.plotting_dnds()
    class_dnds.plotting_pnps()
    

if __name__ == "__main__":
    main()
