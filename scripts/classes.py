"""
Classes.py

Author: Jared Galloway
Date: 2/20/2020

This file will contain all class definitions
to be used when drawing motifs.
"""

#IMPORTS
import io
import sys
import re
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

class MotifController(object):
    """
    This class will take the motif file path and return 
    return an onject which will describe the colors and regex patterns 
    for each motif. 

    Along with storing patterns and colors
    """


    def __init__(self, motif_file_path):

        self.motifs = []
        self.colors = []
        self.patterns = []
        self.num_motifs = 0

        fp = open(motif_file_path)
        for i,line in enumerate(fp):
            motif = line.strip()
            # TODO Error checking here.
            self.motifs.append(motif)
            self.patterns.append(self._create_regex_pattern_from_motif(motif))

        self.num_motifs = i+1
        color_map = cm.get_cmap('viridis', self.num_motifs)
        self.colors = color_map.colors
        

    def _create_regex_pattern_from_motif(self, motif):
        """
        string -> string

        This private function will take in a motif and return the regex
        pattern to account for sequence ambiguity.
        """

        pattern_dict = {
            "U":"T", "u":"t",
            "Y":"[CT]", "R":"[AG]", "W":"[AT]", "S":"[GC]", "K":"[TG]", "M":"[CA]",
            "D":"[AGT]", "V":"[ACG]", "H":"[ACT]", "B":"[CGT]", "N":"[ATCG]", "X":"[ATCG]",
            "y":"[ct]", "r":"[ag]", "w":"[at]", "s":"[gc]", "k":"[tg]", "m":"[ca]",
            "d":"[agt]", "v":"[acg]", "h":"[act]", "b":"[cgt]", "n":"[atcg]", "x":"[atcg]"
            }
        
        pattern = ""
        for base in motif:
            if base in 'atcgATCG':
                pattern += base
            else:
                pattern += pattern_dict[base]

        return pattern


class MotifPattern(object):
    """
    this class is essentially just a dictionary to hold 
    information about motifs as they're found.
    """

    def __init__(self,
        start_position,
        end_position,
        motif_index
        ):
        
        self.start_position = start_position
        self.end_position = end_position
        self.motif_index = motif_index


class SequenceInformation(object):
    """
    """
    

    def __init__(self,
        header,
        sequence,
        motif_controller
        ):

        self.length = len(sequence)
        self.start_position = None
        self.end_position = None
        self.exon_start = None
        self.exon_end = None
        self.motifs = []

        self._extract_motif_information(sequence, motif_controller)
        # TODO start/end position from header
        # TODO exons 

    def _extract_motif_information(self, sequence, mot_con):
        """
        string, string -> dict.

        This function will take in a single entry from a fasta file,
        presumably the transcript containing motifs and EXONS. 
        Exons in the string should be UPPER-CASE. The function will 
        then populate and return a dictionary which contains all relevent 
        information necessary to draw a picture of that transcript.

        As far as I can tell, the necessary information is the 
            
            - start/end position of the transcript.
            - [(motif position / motif type_index)]
            
        """
        
        for motif_idx, (pattern, motif) in enumerate(zip(mot_con.patterns,mot_con.motifs)):
            window_size = len(motif)
            for seq_idx in range(self.length - window_size + 1):
                # TODO check for exon start/end
                
                match = re.match(pattern, sequence[seq_idx:seq_idx + window_size])
                if match == None:
                    continue
                else:
                    self.motifs.append(MotifPattern(
                        start_position = seq_idx,
                        end_position = seq_idx + window_size,
                        motif_index = motif_idx
                    ))

        return None
            
            
            
class SequenceDrawer(object):
    """
    This class will handle drawing sequences
    given 
    """

    pass









