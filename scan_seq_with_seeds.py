from Bio import SeqIO

from Bio.Seq import Seq

from Bio.Alphabet import generic_dna

from Bio.Alphabet import generic_rna

import subprocess

import os

import re

from Bio.SeqUtils import GC

import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import interactive

import random

from matplotlib.backends.backend_pdf import PdfPages

from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

##### parameters

utr_Database="circHIPK3.fa" 

mirna="mir_of_int_ctrl.fa" 

seed_size=6

nb_max_mismatch=0

output_file='miRofINT_scan_circHIPK3.bed'

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def window(fseq, window_size=5):
    for i in range(len(fseq) - window_size + 1):
        yield fseq[i:i+window_size]
		
		
f1 = open(output_file, 'w')

slen=seed_size+1

features=[]

for miR in SeqIO.parse(mirna, "fasta"):
	mir_seed=Seq(str(miR.seq)[1:slen],generic_rna)  #"AAGGCAC"
	print(miR.name)
	print(str(mir_seed))
	for utr in SeqIO.parse(utr_Database, "fasta"):
		pos = 0
		for seq in window(str(utr.seq), len(str(mir_seed))):
			if(hamming2(str(seq.upper()), str(mir_seed.back_transcribe().reverse_complement()))<=nb_max_mismatch):
				f1.write(utr.id+"\t"+str(pos)+"\t"+str(pos+len(str(mir_seed)))+"\t"+miR.id+"\t"+str(hamming2(str(seq.upper()), str(mir_seed.back_transcribe().reverse_complement())))+"\t"+ "+"+"\t"+str(seq.upper())+"\n")
				features.append(GraphicFeature(start=pos, end=pos+len(str(mir_seed)), strand=+1, color="#ccccff",label=re.sub(r'mmu-', '', miR.id)))
			pos=pos+1
			#print(pos)
			
f1.close()

record = CircularGraphicRecord(sequence_length=1100, features=features) #
record.plot(figure_width=2)
plt.show()
