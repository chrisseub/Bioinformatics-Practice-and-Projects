# DNA Toolset/Code testing file
from DNAToolkit import *
from utilities import colored
import random

# Create a random DNA seq for testing
rndDNAStr = "".join(random.choice(Nucleotide) for _ in range(50))

DNAStr = (validateSeq(rndDNAStr))


#Firtst Tutorial Execution, for NA Count
print(f'\nSequence: {colored(DNAStr)}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')


#Second Tutorial Execution, for transcription
print(f'[3] + DNA to RNA Transcription: {colored(transcription(DNAStr))}\n')
print(f"[4] + DNA String + Reverse Complement:\n5' {colored(DNAStr)} 3'")

print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {colored(reverse_complement(DNAStr))}  5'\n")


#Third Tutorial Execution, for GC content determination
print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n')