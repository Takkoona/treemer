# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 08:51:59 2017

@author: Chengyang
"""

from Bio import Phylo, SeqIO, AlignIO
from Trinity import Trinity
import argparse
#%%

parser = argparse.ArgumentParser(description="Trim seq by tree.")
parser.add_argument('seqFile', type=str,
                    help='Sequence file in fasta format')
parser.add_argument('alignFile', type=str,
                    help='Alignment file from the sequence file')
parser.add_argument('treeFile', type=str,
                    help='Tree file from the alignment file')

args = parser.parse_args()

seq_file = args.seqFile
align_file = args.alignFile
align_fmt = "clustal"
tree_file = args.treeFile

seqs = SeqIO.index(seq_file, 'fasta')
aligns = AlignIO.read(align_file, align_fmt)
tree = Phylo.read(tree_file, 'newick')

x = Trinity(seqs, aligns, tree)

x.set_similarity(0.9)
#x.set_level(3)
x.set_blast()

clstr = x.trim_by_tree()

for cluster in clstr:
    for trichord in cluster:
        if trichord.prsrv:
            print str(trichord) + '*'
        else:
            print trichord
    print ''