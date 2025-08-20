# DNA Toolkit file
import collections
from DNA_Structures import *


# Check the sequence to make sure it is a valid DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotide:
            return False
    return tmpseq

def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def transcription(seq):
    """DNA to RNA Transcription, Replacing Replacing Thymine with Uracil"""
    return seq.replace("T", "U")

def reverse_complement(seq):
    """Swapping Adenine with Thymine and  Guanine with Cytosine, Reversing newly generated strings"""
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])

def gc_content(seq):
    """GC Content present in a DNA/RNA Sequence"""
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def gc_content_subsec(seq, k=20):
    """GC Content for DNA/RNA of sub-length K, which is 20 by default, edit to your choice"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res