# ORFFinder
ORFFinder in Python. Inspired by NCBI's version: https://www.ncbi.nlm.nih.gov/orffinder/

### Installation:
`pip3 install orffinder`

### Usage
Import the package

```py
from Bio import SeqIO
import orffinder

sequence = SeqIO.read("gene.fasta", "fasta")
orffinder.getORFs(sequence, minimum_length=75, remove_nested=True)
```

### Documentation
**getORFs()**
```
Returns the loci of discovered ORFs in a dictionary format.

sequence: sequence in Biopython Seq or String format.
minimum_length: minimum size of ORF in nucleotides. Default: 75
start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
remove_nested: remove all ORFs completely encased in another. Default: False
trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
```

**getORFNucleotides()**
```
Returns the loci of discovered ORFs in a dictionary format.

sequence: sequence in Biopython Seq or String format.
return_loci: return the loci together with the nucleotide sequences. Default: False
minimum_length: minimum size of ORF in nucleotides. Default: 75
start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
remove_nested: remove all ORFs completely encased in another. Default: False
trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
```

**getORFProteins()**
```
Returns the loci of discovered ORFs in a dictionary format.

sequence: sequence in Biopython Seq or String format.
translation_table: translation table as per BioPython. Default: 1
return_loci: return the loci together with the protein sequences. Default: False
minimum_length: minimum size of ORF in nucleotides. Default: 75
start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
remove_nested: remove all ORFs completely encased in another. Default: False
trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
```
