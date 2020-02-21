#!/usr/bin/env python3
import argparse
from collections import defaultdict
import sys
#import gzip
import numpy as np
from helpers import *
from classes import *

# TODO This is a rough draft.
parser = argparse.ArgumentParser(description='This \
    program will take as input a fasta file, \
    and a file with motif sequences of interest found in transcripts. \
    and produce a plot with the exons and motifs labled. \
    A subplot containing a linear cartoon will be produced in one pdf \
    where each sequence in the fasta gets a single line.')

parser.add_argument('-fasta', type=str,
    help='fasta file containing transcripts where we are looking for motifs \
    This file can be ')

parser.add_argument('-motifs', type=str, 
    help=' properly formatted \
    list of valid motifs to be labled on cartoon. \
    This should be a text file where each motif is seperated \
    by a new line')

parser.add_argument('-out', type=str, 
    help='The name of the output cartoon')

# TODO
# could add color palette
# could add output pref, svg of pdf.


if __name__ == "__main__":
    
    args = parser.parse_args()
    
    motif_controller = MotifController(args.motifs)
    print('\n'.join(motif_controller.patterns))
    #for color in motif_controller.colors:
    #    print(color)
    sequence_iterator =  fasta_record_iter(open(args.fasta))
    
    first = next(sequence_iterator)
    print('\n'.join(first))

    seq_info = SequenceInformation(first[0],first[1],motif_controller)
    #print(f"motif: {seq_info.motifs}")
    for motif in seq_info.motifs:
        print(f"start: {motif.start_position}, end: {motif.end_position}, m_idx: \
{motif.motif_index}, sequence: {first[1][motif.start_position:motif.end_position]}")
    print(f"exon start: {seq_info.exon_start}, exon_end: {seq_info.exon_end}")

    


    

    
        
