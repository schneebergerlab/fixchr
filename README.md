# fixchr
## Introduction
This package selects homologous chromosomes between two genomes by comparing whole-genome alignments between them. Additionally, it generates dotplots for quick checking of the output.


## Installation:
The easiest method to install fixchr is using anaconda:
```
conda install -c bioconda fixchr    # TO BE DONE
```
For manual installation the pre-requisites are:
1. Python >= 3.8
2. Python libraries. These can be installed in a conda environment using:
```
conda install numpy=1.21.2 pandas=1.2.4 matplotlib=3.3.4 setuptools pysam=0.19.0
```
Then download fixchr and install:
```
git clone https://github.com/schneebergerlab/fixchr.git
cd fixchr
python setup.py install
```

After this fixchr should be installed and in your environment. Test it by printing the help message:
```
fixchr -h
```

The package provides two binaries:
1. `fixchr` for filtering chromosomes
2. `dotplot` for visualising alignments using SAM/BAM/PAF/coord(from mummer) files

## Inputs requirements
1. Chromosome-level assemblies for the genomes to be compared
2. Whole-genome alignment between them

## Example command
```
fixchr -c alignment.paf -r ref.fa.gz -q qry.fa.gz -F P
dotplot -c alignment.paf -r ref.fa.gz -q qry.fa.gz -F P
```

