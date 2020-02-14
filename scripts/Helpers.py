import os
import sys


def convert_fasta_to_phylip(fasta_file):
    """
    This function will take in a fasta file and
    produce another file with the same base name and the extention
    .phylip which 
    """
    fafp = open(fasta_file,"r")
    fyfp = open(fasta_file+".phylip","w")
    fa_iter = fasta_record_iter(fafp)
    for header, sequence in fa_iter:
        fyfp.write(f"{header}\t{sequence}\n")
    return None

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

def fastq_records_iterator(fast1_fp):
    '''
    This function will take in a file pointer 
    and returns an iterator which yeilds four lines of a file
    (stripped of \n) in a list of strings. 
    '''
    while True:
        record = []
        for i in range(4):
            record.append(file_pointer.readline().strip())
        if record[0] == '':
            break
        else:
            yield tuple(record)
    return None

def average_quality_score(phred_string):
    '''
    Take in a string of quality scores and compute the average.  
    '''
    average_score = 0
    for score in phred_string:
        average_score += convert_phred(score)
    return average_score/len(phred_string)

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33

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


