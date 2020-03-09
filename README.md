# motif-mark

Motif Mark is a simple application which allows us to visualize 
certain motifs of interest in sequences from multi-line fasta files.

## Install

`git clone https://github.com/jgallowa07/motif-mark.git`

and be sure to install all listed in `requirements.txt`

## Usage

```bash
usage: motif_mark.py [-h] [-fasta FASTA] [-motifs MOTIFS] [-cmap CMAP]
                     [-out OUT]

This program will take as input a fasta file, and a file with motif sequences
of interest found in transcripts. and produce a plot with the exons and motifs
labled. A subplot containing a linear cartoon will be produced in one pdf
where each sequence in the fasta gets a single line.

optional arguments:
  -h, --help      show this help message and exit
  -fasta FASTA    fasta file containing transcripts where we are looking for
                  motifs This file can be
  -motifs MOTIFS  properly formatted list of valid motifs to be labled on
                  cartoon. This should be a text file where each motif is
                  seperated by a new line
  -cmap CMAP      A color map of your choosing. The default is plasma, but
                  viridis, magma, inferno, and cividis are all acceptable
                  inputs
  -out OUT        The name of the output pdf cartoon, default is
                  motif_mark.pdf
```

## Example

running

`cd scripts && python3 motif_mark.py -fasta ../test_files/Figure_1.fasta -motifs ../test_files/Fig_1_motifs.txt -out example.pdf -cmap plasma`

will produce 

![alt text](./test_files/example.pdf)














