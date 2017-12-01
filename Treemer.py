# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 08:51:59 2017

@author: Chengyang
"""

from Bio import Phylo, SeqIO, AlignIO
from Trinity import Trinity
import argparse

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
tree_file = args.treeFile

try:
    seqs = SeqIO.index(seq_file, 'fasta')
except IOError:
    print "Make sure the sequence in fasta format"

align_fmts = ["clustal", "phylip", "fasta",
              "emboss", "stockholm"]
for fmt in align_fmts:
    try:
        aligns = AlignIO.read(align_file, fmt)
    except IOError:
        print "The alignment format is not {}".format(fmt)
    except ValueError:
        print "Read {} format failed".format(fmt)
try:
    aligns
except NameError:
    print "{} is not any supported format".format(align_file)

tree_fmts = ["newick", "nexus", "phyloxml", "nexml"]
for fmt in tree_fmts:
    try:
        tree = Phylo.read(tree_file, 'newick')
    except IOError:
        print "The tree format is not {}".format(fmt)

try:
    tree
except NameError:
    print "{} is not any supported format".format(tree_file)

x = Trinity(seqs, aligns, tree)

x.set_similarity(0.97)
#x.set_level(1)
#print x.check_num()

clstr = x.trim_by_tree()
#
#y = x.aligned
#
for cluster in clstr:
    for trichord in cluster:
        if trichord.prsrv:
            print '{}*'.format(trichord)
        else:
            print trichord
    print ''
    
#for i in y:
#    if 'NL/118/2001' in i or 'BI/16190/1968' in i:
#        print i, y[i]