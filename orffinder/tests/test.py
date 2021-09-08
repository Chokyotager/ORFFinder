from Bio import SeqIO
from Bio.Seq import Seq

import orffinder

sequence = list(SeqIO.parse("gene.fasta", "fasta"))[0]

out = orffinder.getORFProteins(sequence, minimum_length=75, remove_nested=True)
print(out)
