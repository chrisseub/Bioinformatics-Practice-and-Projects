# DNA Toolset/Code testing file
from DNAToolkit import *
import random

# Create a random DNA seq for testing
rndDNAStr = "".join(random.choice(Nucleotide) for _ in range(50))

DNAStr = (validateSeq(rndDNAStr))

#Firtst Tutorial Execution
print(f'\nSequence: {DNAStr}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')

#Second Tutorial Execution
print(f'[3] + DNA to RNA Transcription: {transcription(DNAStr)}\n')