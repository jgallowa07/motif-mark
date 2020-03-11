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
import cairo
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

    def __init__(self, motif_file_path, cmap):

        self.motifs = []
        self.colors = []
        self.patterns = []
        self.num_motifs = 0

        fp = open(motif_file_path)
        for i,line in enumerate(fp):
            motif = line.strip()
            self.motifs.append(motif)
            self.patterns.append(self._create_regex_pattern_from_motif(motif))

        self.num_motifs = i+1
        color_map = cm.get_cmap(cmap, self.num_motifs)
        self.colors = color_map.colors
        
    def _create_regex_pattern_from_motif(self, motif):
        """
        string -> string

        This private function will take in a motif and return the regex
        pattern to account for sequence ambiguity.
        """

        pattern_dict = {
            "u":"[ut]","t":"[ut]",
            "y":"[ct]", "r":"[ag]", "w":"[at]", "s":"[gc]", "k":"[tg]", "m":"[ca]",
            "d":"[agt]", "v":"[acg]", "h":"[act]", "b":"[cgt]", "n":"[atcg]", "x":"[atcg]"
            }
        
        pattern = ""
        for base in motif.lower():
            if base in 'atcg':
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

        self.header = header 
        self.num_nuc = len(sequence)
        self.exon_start = None
        self.exon_end = None
        self.motifs = []

        self._extract_motif_information(sequence, motif_controller)

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
        exon_flag = False
        for motif_idx, (pattern, motif) in enumerate(zip(mot_con.patterns,mot_con.motifs)):
            window_size = len(motif)
            for seq_idx in range(self.num_nuc - window_size + 1):

                # TODO question : is there only one exon per seq?
                if str.isupper(sequence[seq_idx]) and not exon_flag:
                    self.exon_start = seq_idx
                    exon_flag = True
                if str.islower(sequence[seq_idx]) and exon_flag:
                    self.exon_end = seq_idx
                    exon_flag = False
                
                match = re.match(pattern, sequence[seq_idx:seq_idx + window_size].lower())
                if match != None:
                    self.motifs.append(MotifPattern(
                        start_position = seq_idx,
                        end_position = seq_idx + window_size,
                        motif_index = motif_idx
                    ))

        return None
            
            
# TODO            
class SequenceDrawer(object):
    """
    This class will handle drawing sequences
    given all sequence data objects in a list,
    a respective motif controller,
    and an output pdf filepath
    """


    def __init__(self, motif_controller, sequence_info_list, output):
        # feed a list of sequence info object, and the respective 
        # motif controller to draw sequences.
        self.spread = 128
        self.top_space = 128

        self.motif_controller = motif_controller
        self.sequence_info_list = sequence_info_list
        self.output = output
        self.width = 1024
        self.height = self.spread * len(sequence_info_list) + self.top_space
        self.surface = cairo.PDFSurface(self.output, self.width, self.height)        
        self.ctx = cairo.Context(self.surface)
        self.ctx.scale(self.width, self.height)
        self.ctx.rectangle(0, 0, 1, 1)  # Rectangle(x0, y0, x1, y1)
        self.ctx.set_source_rgb(1, 1, 1)
        self.ctx.fill()

        self.draw_all_sequences()

    def draw_all_sequences(self):
        
        
        # set some constants
        right_margin = 0.2
        left_margin = 0.1
        motif_height = 0.01
        exon_height = 0.025
        header_margins = 40
        legend_margins = 20
        box_text_margin = 10

        longest = max([seq_data.num_nuc for seq_data in self.sequence_info_list])

        # make legend first
        top_left_x = 1 - right_margin + (legend_margins / self.height)
        top_left_y = (legend_margins / self.height)
        
        legend_width = (1 - legend_margins / self.width) - top_left_x
        legend_height = len(self.sequence_info_list) * (legend_margins / self.height)
        
        self.ctx.rectangle(top_left_x, top_left_y, legend_width, legend_height)
        self.ctx.set_source_rgba(0, 0, 0, 0.9)
        self.ctx.set_line_width(0.001)
        self.ctx.stroke()

        for i in range(self.motif_controller.num_motifs):

            y = top_left_y + (i * (legend_margins / self.height)) + ((legend_margins / self.height) / 2)
            width = box_text_margin / self.width

            self.ctx.rectangle(top_left_x + (legend_margins / self.width), y - width / 2, width, width)
            r = self.motif_controller.colors[i][0]
            g = self.motif_controller.colors[i][1]
            b = self.motif_controller.colors[i][2]
            self.ctx.set_source_rgba(r, g, b, 0.7)
            self.ctx.fill()
            
            self.ctx.set_source_rgb(0, 0, 0)            
            self.ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.ctx.set_font_size(0.01)
            self.ctx.move_to((2*width) + top_left_x + (legend_margins / self.width), y)
            self.ctx.show_text(self.motif_controller.motifs[i])

        # draw each sequence one at a time
        for i,seq_info in enumerate(self.sequence_info_list):

            # where is the "center" height for each sequence
            y = ((i * self.spread) + (self.spread // 2) + self.top_space) / self.height
            
            # draw header
            self.ctx.set_source_rgb(0, 0, 0)            
            self.ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.ctx.set_font_size(0.015)
            self.ctx.move_to(left_margin, y - (header_margins/self.height))
            self.ctx.show_text(seq_info.header)
            
            # draw sequence line
            right_x = seq_info.num_nuc / longest - right_margin
            total_length = right_x - left_margin

            self.ctx.move_to(left_margin, y)
            self.ctx.line_to(right_x, y)
            self.ctx.set_source_rgb(0, 0, 0)
            self.ctx.set_line_width(0.005)
            self.ctx.stroke()

            # draw exon
            start_prop = seq_info.exon_start / seq_info.num_nuc
            end_prop = seq_info.exon_end / seq_info.num_nuc
            
            top_left_x = (start_prop * total_length) + left_margin
            top_left_y = y - exon_height
            width = ((end_prop * total_length) + left_margin) - top_left_x
            height = exon_height * 2
            
            self.ctx.rectangle(top_left_x, top_left_y, width, height)
            self.ctx.set_source_rgba(0, 0, 0, 0.9)
            self.ctx.set_line_width(0.001)
            self.ctx.stroke()
                    
            # draw each motif
            for motif in seq_info.motifs:
                
                start_prop = motif.start_position / seq_info.num_nuc
                end_prop = motif.end_position / seq_info.num_nuc
                
                top_left_x = (start_prop * total_length) + left_margin
                top_left_y = y - motif_height
                width = ((end_prop * total_length) + left_margin) - top_left_x
                height = motif_height * 2
                
                self.ctx.rectangle(top_left_x, top_left_y, width, height)
                r = self.motif_controller.colors[motif.motif_index][0]
                g = self.motif_controller.colors[motif.motif_index][1]
                b = self.motif_controller.colors[motif.motif_index][2]
                self.ctx.set_source_rgba(r, g, b, 0.7)
                self.ctx.fill()

        return None


