# DNA Toolset/Code testing file
from DNAToolkit import *
import random

# Create a random DNA seq for testing
rndDNAStr = "".join(random.choice(Nucleotide) for _ in range(50))

DNAStr = (validateSeq(rndDNAStr))
print(countNucFrequency(DNAStr))