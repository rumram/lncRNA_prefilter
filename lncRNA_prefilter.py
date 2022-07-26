import sys
from Bio import SeqIO
import re


# The script extracts FASTA sequences from multifasta file created with StringTie (hence 'MSTRG' IDs).
# Standard codon table is used.
# Sequences with ORF longer than 100 aa are discarded as well as sequences with name starting with 'ENS' (Ensembl IDs).
# The script searches ORF in 3 possible translation frames.
# Reverse frames not required since StringTie provides only forward direction.


# Provide FASTA file input as first argument.
# Provide FASTA file name as second argument.

def translation(nucleotide_seq):
    # Translate nucleotide codon to amino acid
    # TTG, CTG always treated as start codons!
    threes = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'M', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'M',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

    protein = ""
    if len(nucleotide_seq) % 3 != 0:
        z = len(nucleotide_seq) % 3
        rem = len(nucleotide_seq[:-z])
    else:
        rem = len(nucleotide_seq)
    for i in range(0, rem, 3):
        codon = nucleotide_seq[i:i + 3]
        protein += threes[codon]
    return protein


names = []

with open(sys.argv[1]) as hand:
    for seqid in SeqIO.parse(hand, "fasta"):
        if seqid.name.startswith("MSTRG"):
            for i in range(3):
                if re.findall(r'M[^_]{98,}_', translation(seqid.seq[i:])):
                    names.append(seqid.name)
    dset = set(names)
    dlist = list(dset)

filtered_seqs = []
with open(sys.argv[1]) as hand:
    for seqid in SeqIO.parse(hand, "fasta"):
        if seqid.name not in dlist and not seqid.name.startswith("ENS"):
            filtered_seqs.append(seqid)

SeqIO.write(filtered_seqs, sys.argv[2], "fasta")
