import collections
import sys

def ENc():
    input_file = sys.argv[1]
    print("gene,Nc,number_of_codons")
    with open(input_file, 'r') as f:
        head = ''
        for line in f:
            if line[0] == ">":
                if len(head) != 0 and len(seq) != 0:
                    k,n,fcf,tc = prep_counts_for_Nc(seq)
                    Nc = calc_Nc(k, n, fcf)
                    print(head[1:],Nc,tc,sep=',',flush=True)
                head = line.rstrip()
                seq = ''
            else:
                seq += line.rstrip().upper()
        if len(head) != 0 and len(seq) != 0:
            k,n,fcf,tc = prep_counts_for_Nc(seq)
            Nc = calc_Nc(k, n, fcf)
            print(head[1:],Nc,tc,sep=',',flush=True)

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

if __name__ == "__main__":
    ENc()