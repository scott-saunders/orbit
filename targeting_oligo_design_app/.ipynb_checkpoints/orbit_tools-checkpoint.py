import pandas as pd
import numpy as np
import holoviews as hv
import panel as pn

#from functions import *
from Bio.Seq import Seq
from Bio.SeqIO import parse
from bokeh.models import BasicTickFormatter
from holoviews import opts

def get_replichore(pos, ori = 3923882.5, ter = 1590250.5 ):
    
    """
    Determine the replichore of a bacterial chromosome for a certain position. Requires origin and terminus positions. Assumes E. coli like organization.
    
    pos : int
        Genomic coordinate of interest.
    ori : float
        Genomic coordinate of the origin of replication. 
    ter : float
        Genomic coordinate of the replication terminus.
    """
    
    pos = int(pos)
    
    if((pos<0)| (pos>4641652)):
        raise TypeError("position must be within genome.")
    
    if((pos > ori) | (pos<ter)):
       rep = 1
    elif((pos<ori) & (pos>ter)):
       rep = 2
    
    return rep

def get_target_oligo(left_pos, right_pos, genome, homology = 90, attB_dir = '+', attB_fwd_seq = 'ggcttgtcgacgacggcggtctccgtcgtcaggatcat',  verbose = False):
    """
    Given a set of parameters, get an ORBIT oligo that targets the lagging strand. 
    Left and right positions are absolute genomic coordinates that specify the final nucleotides to keep unmodified in the genome, 
    everything in between will be replaced by attB. In other words the left position nucleotide is the final nt before attB in the oligo.
    The right position nt is the first nt after attB in the oligo.
    
    This function determines the lagging strand by calling `get_replichore()` on the left_pos.
    Typically attB_dir should be set to the same direction as the gene of interest, such that the integrating plasmid will insert with payload facing downstream.
    attB_fwd_seq can be modified, and the total homology can be modified, but should be an even number since homology arms are symmetric. 
    
    Verbose prints helpful statements for testing functionality.
    
    Parameters
    -----------------
    left_pos : int
        Left genomic coordinate of desired attB insertion. attB is added immediately after this nt.
    right_pos : int
        Right genomic coordinate of desired attB insertion. attB is added immediately before this nt.
    genome : str
        Genome as a string.
    homology : int (even)
        Total homology length desired for oligo. Arm length = homology / 2.
    attB_dir : chr ('+' or '-')
        Desired direction of attB  based on genomic strand. Typically same direction as gene.
    attB_fwd_seq : str
        Sequence of attB to insert between homology arms.
    verbose : bool
        If true, prints details about genomic positions and replichore.
    Returns
    ---------------
    oligo : str
        Targeting oligo against lagging strand, including the attB sequence in the correct orientation.
    """
    
    left_pos = int(left_pos)
    
    right_pos = int(right_pos)
    
    # Arm length is 1/2 total homology. Arms are symmetric
    arm_len = int(homology / 2)
    
    # Arms from genome string. Note 0 indexing of string vs. 1 indexing of genomic coordinates.
    # As written, should be inclusive.
    left_arm = genome[(left_pos - arm_len):left_pos]
    
    right_arm = genome[(right_pos - 1):(right_pos - 1 + arm_len)]

    # Generate attB reverse sequence
    seq_attB = Seq(attB_fwd_seq)
    attB_rev_seq = str(seq_attB.reverse_complement())
    
    # Replichore 1
    if get_replichore(left_pos) == 1:
        
        rep = 1
        
        # Reverse complement replichore 1 sequences.
        left_arm_seq = Seq(left_arm)
        left_arm_prime = str(left_arm_seq.reverse_complement())
        
        right_arm_seq = Seq(right_arm)
        right_arm_prime = str(right_arm_seq.reverse_complement())
        
        # Determine attB direction and paste fwd/rev seq accordingly
        if attB_dir == '+':
            
            oligo = right_arm_prime + attB_rev_seq + left_arm_prime
            
        elif attB_dir == '-':
            
            oligo = right_arm_prime + attB_fwd_seq + left_arm_prime
    
    # Replichore 2
    elif get_replichore(left_pos) == 2:
        
        rep = 2
        
        # '+' arm sequence used. Determine attB direction and paste accordingly.
        if attB_dir == '+':
            
            oligo = left_arm + attB_fwd_seq + right_arm
        
        elif attB_dir == '-':
            
            oligo = left_arm + attB_rev_seq + right_arm    
            
    # Verbose print statements
    if verbose:
        
        print('left_arm_coord = ', left_pos - arm_len,' : ', left_pos)
        print('right_arm_coord = ', right_pos - 1, ' : ', right_pos -1 + arm_len)
        print('Replichore = ', rep)
    
    return oligo


def get_pos_details(left_pos, right_pos, homology, direction):
    
    left_pos = int(left_pos)
    right_pos = int(right_pos)
    
    replichore = get_replichore(left_pos)
    
    arm_len = int(homology) / 2
    
    rep_dir = str(replichore) + direction
    
    rep_dir_dict = {
        "1+": "`5' |-- Right_arm (Downstream) --|-- attB_rev --|-- Left_arm (Upstream) --| 3'`",
        "1-": "`5' |-- Right_arm (Upstream) --|-- attB_fwd --|-- Left_arm (Downstream) --| 3'`",
        "2+": "`5' |-- Left_arm (Upstream) --|-- attB_fwd --|-- Right_arm (Downstream) --| 3'`",
        "2-": "`5' |-- Left_arm (Downstream) --|-- attB_rev --|-- Right_arm (Upstream) --| 3'`"
    }
    
    #rep_dir_dict.get(rep_dir, "No info available")
    
    left_arm_str = '\n\n**Left arm:** `(' + str(int(left_pos - arm_len)) + ' - ' + str(left_pos) + ') nt`'
    right_arm_str = '\n\n**Right arm:** `(' + str(right_pos) + ' - ' + str(int(right_pos + arm_len)) + ') nt`'
    rep_str = '\n\n**Replichore:** `' + str(replichore) + '`'
    dir_str = '\n\n**attB direction:** `' + direction + '`'
    oligo_len = '\n\n**Oligo length:** `' + str(int(homology + 38)) + ' nt`'
    
    md =   left_arm_str + right_arm_str + rep_str + dir_str + oligo_len + '\n\n**Oligo structure:** ' + rep_dir_dict.get(rep_dir, "No info available")

    
    return pn.pane.Markdown(md, width = 1000)

def plot_nearby(left_pos, right_pos, homology, df_genes):
    
    left_pos = int(left_pos)
    right_pos = int(right_pos)
    
    arm_len = homology / 2
    
    arms = {'start': [left_pos-arm_len, right_pos], 'stop': [left_pos, right_pos+arm_len],'arm':['left','right'],'target_oligo':['target_oligo','target_oligo']}
    
    left_line = hv.VLine(left_pos).opts(color = 'black', line_width = 1)
    right_line = hv.VLine(right_pos).opts(color = 'black', line_width = 1)
    
    arms = hv.Segments(arms, kdims = ['start', 'target_oligo','stop','target_oligo'])
    
    genome_segments = hv.Segments(df_genes, kdims = ['left_pos','Direction', 'right_pos','Direction']).opts(tools = ['hover'])
    
    genome_points = hv.Scatter(df_genes, 'left_pos','Direction') * hv.Scatter(df_genes, 'right_pos','Direction')
    
    genome_labels = hv.Labels(df_genes,kdims = ['center_pos','Direction'], vdims = 'gene_label' ).opts(text_font_size='8pt', text_color='gray', xoffset = 0)
    
    genome_plot = genome_segments * genome_points* genome_labels *left_line * right_line * arms

    return genome_plot.opts(xlim = (left_pos - 1000, right_pos + 1000), width = 1000,xformatter = BasicTickFormatter(use_scientific = False))
    