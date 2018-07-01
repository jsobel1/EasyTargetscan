#!/usr/bin/python

# Easy TargetScan.py

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
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import random

from matplotlib.backends.backend_pdf import PdfPages

from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

import argparse
parser = argparse.ArgumentParser(description='Allows to use TargetScan from FASTA files')
parser.add_argument("-utr", type=str, help="a FASTA file conataining sequences to scan",required=True)
parser.add_argument("-mir", type=str, help="a FASTA file conataining mature miRNA sequences",required=True)
parser.add_argument("-output", type=str, help="prefix of the output file",default="test_EasyTargetScan")
#parser.add_argument("-seed", type=int, help="seed lenght: by default 8", default=8)
parser.add_argument("-Taxon", type=int, help="mouse by default:10090", default=10090)
args = parser.parse_args()

utr_Database=args.utr 

mirna=args.mir 

seed_size=7

taxon=args.Taxon 

nb_max_mismatch=0

output_file=args.output

path_to_targetscan='targetscan_70.pl'

def get_color(seed_match_type):
	return {
        '6mer': "#00ff99",
        '7mer-1a': "#9999ff",
		'7mer-m8':"#ff66cc",
		'8mer-1a': "#ff0000"
    }.get(seed_match_type, "#ccccff") 

f1 = open("miR_seeds_temp.txt", 'w')
slen=seed_size+1

for miR in SeqIO.parse(mirna, "fasta"):
	mir_seed=Seq(str(miR.seq)[1:slen],generic_rna) 
	f1.write(miR.id+"\t"+str(mir_seed)+"\t"+str(taxon)+"\n")	
f1.close()

for utr in SeqIO.parse(utr_Database, "fasta"):
	f2 = open("UTRs_temp.txt", 'w')
	f2.write(utr.id+"\t"+str(taxon)+"\t"+str(utr.seq)+"\n")
	f2.close()
	output_file_utr='_'.join([utr.id,output_file])
	p = subprocess.Popen(' '.join(["perl", path_to_targetscan, "miR_seeds_temp.txt", "UTRs_temp.txt", output_file_utr]), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	p.wait()
	features=[]

	with open(output_file_utr) as fp:
		next(fp)
		for line in fp:
			print(line)
			content=line.split("\t")
			features.append(GraphicFeature(start=int(content[3]), end=int(content[4]), strand=+1, color=get_color(content[8]),label=re.sub(r'mmu-', '', content[1])))

	record = GraphicRecord(sequence_length=len(str(utr.seq)), features=features) # Circular
	record.plot(figure_width=12)
	plt.title(' '.join([utr.id,'sequence']))
	patch1 = mpatches.Patch(color="#00ff99", label='6mer')
	patch2 = mpatches.Patch(color="#9999ff", label='7mer-1a')
	patch3 = mpatches.Patch(color="#ff66cc", label='7mer-m8')
	patch4 = mpatches.Patch(color= "#ff0000", label='8mer-1a')
	plt.legend(handles=[patch1,patch2,patch3,patch4])
	plt.show()




