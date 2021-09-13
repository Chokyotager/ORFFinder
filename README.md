<div align="center">
  <br />
  <p>
    <a href="https://github.com/Chokyotager/ORFFinder"><img src="https://github.com/Chokyotager/ORFFinder/blob/main/images/ORFFinder.png?raw=true" alt="banner" /></a>
  </p>
  <br />
  <p>
    <a href="https://pepy.tech/project/orffinder"><img src="https://pepy.tech/badge/orffinder" alt="Downloads" /></a>
  </p>
</div>

# ORFFinder
ORFFinder in Python. Inspired by NCBI's version: https://www.ncbi.nlm.nih.gov/orffinder/

Finds the open reading frame (6-frame scan) on a given 5' to 3' nucleotide.

### Installation:
`pip3 install orffinder`

### Terminal Usage
Two command-line executable commands are available: `orffinder-to-gtf` `orffinder-to-sequence`.

Documentation for these commands can be retrieved by specifying `<command> -h`.

### API Usage
Import the package

**IMPORTANT: Your DNA/RNA strand should always be from the 5' to 3' direction when input!**
![Transcription direction](https://cdn.kastatic.org/ka-perseus-images/1da89713b9aa8067742244d916749e72561bb3cc.png)
(Image credit: Khan Academy)

```py
from Bio import SeqIO
from orffinder import orffinder

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
Returns a list of Biopython Seq objects or loci of discovered ORFs with Biopython Seq objects in a dictionary format.

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
Returns a list of Biopython Seq objects or loci of discovered ORFs with Biopython Seq objects in a dictionary format.

sequence: sequence in Biopython Seq or String format.
translation_table: translation table as per BioPython. Default: 1
return_loci: return the loci together with the protein sequences. Default: False
minimum_length: minimum size of ORF in nucleotides. Default: 75
start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
remove_nested: remove all ORFs completely encased in another. Default: False
trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
```

### Dependencies
Biopython (https://biopython.org/)
