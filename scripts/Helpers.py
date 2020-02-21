import os
import sys


def fasta_record_iter(fasta_fp):
    """
    This iterator takes a fasta file pointer,
    and yeild tuples where the first element is the string header,
    and the second element is the string containing the full sequence.
    iterates until EOF is reached
    """
    while True:
        header = fasta_fp.readline().strip()
        if header == '':
            break
        sequence = ''
        while True:
            cur_pos = fasta_fp.tell()
            sequence_line = fasta_fp.readline().strip()
            if sequence_line == '':
                yield (header, sequence)
                break
            if sequence_line[0] != '>':
                sequence += sequence_line.strip()
            else:
                fasta_fp.seek(cur_pos)
                yield (header, sequence)
                break
    return None


def reverse_compliment(sequence):
    '''
    This function takes in a string representing a 
    sequence of DNA, and returns it's reverse compliment.
    '''
    map_nucleotide_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
    compliment = ""    
    for nucleotide in sequence:
        compliment += map_nucleotide_dict[nucleotide]
    return ''.join(list(reversed(compliment)))


