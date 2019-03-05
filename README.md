[![LICENSE](https://img.shields.io/badge/license-LGPL-blue.svg)](https://github.com/akshayparopkari/BiasAway/blob/master/LICENSE)

BiasAway
---                            


We provide here the README of the BiasAway software developed in Python. The
software provides a user with four approaches for generating a background
useful to enrichment analyses:

1) mononucleotide shuffled target sequence to preserve the mononucleotide composition of the target sequences,
2) dinucleotide shuffled target sequence to preserve the dinucleotide composition of the target sequences,
3) genomic sequences matched to the mononucleotide composition of each target sequence to preserve the non-random association of nucleotides,
4) sliding window of mononucleotide shuffled target sequence,
5) sliding window of dinucleotide shuffled target sequence, 
6) genomic sequences matched in windows of internal mononucleotide composition for each target sequence.

The latter three backgrounds (BiasAway 4-6) are variants of the former three
backgrounds (BiasAway 1-3), in which we utilized a sliding window over the
ChIP-Seq sequences to determine a distribution for local regions of
composition. The background sequence set is then generated (mono- or
di-nucleotide shuffle) or selected from a pool of genomic sequences (genomic
composition match) to match the distribution of window compositions for each
target sequence. These latter backgrounds were considered because due to
evolutionary changes such as insertion of repetitive sequences, local
rearrangements, or biochemical missteps, the target sequences may have
sub-regions of distinct nucleotide composition.

You can get more information by reading [Improving detection and analysis of transcription
factor binding sites within ChIP-Seq data](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-472).

---

System requirements:

* Python should be installed on your machine (only version 2.7 has been
        tested).
* Biopython should be installed. See
    http://biopython.org for instructions on how to install it.


---

Download:
This is a fork of the original [BiasAway](https://github.com/wassermanlab/BiasAway) tool.

---

Usage:

To get the help on how to use the software, you can run:
- `$ python BiasAway -h`

Help for each specific way of generating background can be obtain by typing:
- `$ python BiasAway m -h` or
- `$ python BiasAway f -h` or
- `$ python BiasAway d -h` or
- `$ python BiasAway w -h` or
- `$ python BiasAway g -h` or
- `$ python BiasAway c -h`

Note that `BiasAway g` and `BiasAway c` ask for a background repository, a
background file can also be giving optionally. If the background repository
contains files corresponding to %GC bins, the software will use these files and
ignore any other background file. If the background repository is empty,
sequences from the background files will be used and %GC bins files will be
constructed and stored in the background directory for further usage.
Using pre-computed %GC bins files help in speeding up the background generation
process.
Preprocessed background repositories of _human_ and _mouse_ can be found at:
http://cisreg.cmmt.ubc.ca/BiasAway_background/

---

Reference:

The BiasAways software have been described in [Improving detection and analysis of transcription
factor binding sites within ChIP-Seq data](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-472).
Please cite this publication when using BiasAway.
