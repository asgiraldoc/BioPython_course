from typing import Counter
import re
import argparse
import os
from dnds import dnds, pnps, substitutions, dnds_codon, dnds_codon_pair, syn_sum, translate
from Bio import SeqIO
import itertools


parser = argparse.ArgumentParser()
parser.add_argument("-f", dest="fa", required=True, type=str, help="fasta folder")
args = parser.parse_args()

EXTENSIONS_FASTA = ['.fasta']

def get_file_list_fasta(root_dir, E):
    file_list = []
    for root, directories, filenames in os.walk(root_dir):
        for filename in filenames:
            if any (ext in filename for ext in E):
                file_list.append(os.path.join(root, filename))
    file_list.sort()
    return file_list

fasta_file = get_file_list_fasta(args.fa, EXTENSIONS_FASTA)
list_files = []
for i in fasta_file:
    list_files.append(i)

sequence_1 = ""
sequence_2 = ""

for i in list_files:
    with open(i, "r") as f:
        seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(i, "fasta")}
        res = list(itertools.combinations(seq_dict, 2))
        for i, j in res:
            if i in seq_dict and j in seq_dict:
                sequence_1 += str((seq_dict[i]))
                sequence_2 += str((seq_dict[j]))
                try:
                    xxx1 = round(float(pnps(sequence_1, sequence_2)), 3)
                    xxx2 = round(dnds(sequence_1, sequence_2), 3)
                except ZeroDivisionError:
                    xxx1 = 0
                    xxx2 = 0
                print(i, j,'\t', xxx1,'\t',xxx2, '\t')
