#!/usr/bin/env python
# coding: utf-8

#------------------------------------------------------------------------------------#
# Function seq2logo adapted from Logomaker python package
#------------------------------------------------------------------------------------#

import pandas as pd
import logomaker as lm

def SeqPlot(weight_matrix):
    """
    Creates a simple sequence logo from a weight matrix
    Weight matrix should be a dict in a list, [{}]
    """
    SeqPlot = lm.Logo(df = pd.DataFrame(weight_matrix),
                      fade_below=0.5,
                      shade_below=0.5,
                      figsize=(15,7))
    
    # set axis labels
    SeqPlot.ax.set_xlabel('Amino acid position', fontsize=14)
    SeqPlot.ax.set_ylabel('bits', fontsize=14)
    
    # Highlighting specific or range of amino acids
    #SeqPlot.highlight_position(2,color='lightgray',alpha=0.5)
    #SeqPlot.highlight_position_range(2,6,alpha=0.5,color='lightgray')
    return SeqPlot

#------------------------------------------------------------------------------------#
# Citation in BibTex format
#------------------------------------------------------------------------------------#
  
"""
@article {Tareen635029,
	author = {Tareen, Ammar and Kinney, Justin B.},
	title = {Logomaker: Beautiful sequence logos in python},
	elocation-id = {635029},
	year = {2019},
	doi = {10.1101/635029},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Sequence logos are visually compelling ways of illustrating the biological properties of DNA, RNA, and protein sequences, yet it is currently difficult to generate such logos within the Python programming environment. Here we introduce Logomaker, a Python API for creating publication-quality sequence logos. Logomaker can produce both standard and highly customized logos from any matrix-like array of numbers. Logos are rendered as vector graphics that are easy to stylize using standard matplotlib functions. Methods for creating logos from multiple-sequence alignments are also included.Availability and Implementation Logomaker can be installed using the pip package manager and is compatible with both Python 2.7 and Python 3.6. Source code is available at http://github.com/jbkinney/logomaker.Supplemental Information Documentation is provided at http://logomaker.readthedocs.io.Contact jkinney{at}cshl.edu.},
	URL = {https://www.biorxiv.org/content/early/2019/05/13/635029},
	eprint = {https://www.biorxiv.org/content/early/2019/05/13/635029.full.pdf},
	journal = {bioRxiv}
}
"""