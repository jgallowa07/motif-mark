#!/usr/bin/env python3
import argparse
from collections import defaultdict
import sys
from helpers import *
from classes import *


if __name__ == "__main__":

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

    parser.add_argument('-cmap', type=str, default='plasma', \
        help=' A color map of your choosing. The default is \
        plasma, but viridis, magma, inferno, and cividis are all \
        acceptable inputs')

    parser.add_argument('-out', type=str, default='motif_mark.pdf', \
        help='The name of the output pdf cartoon, default is \
        motif_mark.pdf')
    
    args = parser.parse_args()
    
    motif_controller = MotifController(args.motifs, args.cmap)
    sequence_iterator = fasta_record_iter(open(args.fasta))
    sequence_info_list = [SequenceInformation(seq[0],seq[1],motif_controller) for seq in sequence_iterator]
    SequenceDrawer(motif_controller, sequence_info_list, args.out)
