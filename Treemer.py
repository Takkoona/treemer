# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 08:51:59 2017

@author: Chengyang
"""

from Bio import Phylo, SeqIO, AlignIO
from Trinity import Trinity
import argparse, os

parser = argparse.ArgumentParser(description="Trim seq by tree.")
parser.add_argument('seqFile', type=str,
                    help="Sequence file in fasta format")
parser.add_argument('alignFile', type=str,
                    help="Alignment file from the sequence file")
parser.add_argument('treeFile', type=str,
                    help="Tree file from the alignment file")
parser.add_argument('-t', '--threshold', type=float,
                    help="Set the similarity threshold")
parser.add_argument('-s', '--sites', type=int, nargs='+',
                    help="The input sites are consistent within the cluster")
parser.add_argument('-c', '--clustering', action='store_true',
                    help="Output the clustering in a file")
parser.add_argument('-l', '--level', type=int,
                    help="Set the reduction level")

args = parser.parse_args()

seq_file = args.seqFile
align_file = args.alignFile
tree_file = args.treeFile
threshold = args.threshold
sites = args.sites
clstr_file = args.clustering
level = args.level

try:
    seqs = SeqIO.index(seq_file, 'fasta')
except IOError:
    print "Make sure the sequence in fasta format"

def attempt_read(read_fun, file_path, fmts):
    file_name = os.path.basename(file_path)
    for fmt in fmts:
        try:
            data = read_fun(file_path, fmt)
            break
        except Exception:
            print "{} is not in {} format".format(file_name, fmt)
    try:
        data
    except NameError:
        print "{} is not in any supported format".format(file_name)
    else:
        return data

align_fmts = ["clustal", "phylip", "fasta", "emboss", "stockholm"]  
tree_fmts = ["newick", "nexus", "phyloxml", "nexml"]

aligns = attempt_read(AlignIO.read, align_file, align_fmts)
tree = attempt_read(Phylo.read, tree_file, tree_fmts)

x = Trinity(seqs, aligns, tree)
if threshold is not None:
    x.set_similarity(threshold)
if sites is not None:
    x.set_sites(*sites)
x.set_level(level)
clstr = x.trim_by_tree()

def cluster_out():
    out = ''
    n_clstr = 1
    for cluster in clstr:
        out += 'cluster {}\n'.format(n_clstr)
        n_clstr += 1
        for trichord in cluster:
            out += '{}\n'.format(trichord)
        out += '\n'
    return out

if clstr_file is False:
    print '\n', cluster_out()
else:
    out_file = seq_file + '.clstr'
    with open(out_file, 'w') as f:
        f.write(cluster_out())

#y = x.aligned
#    
#for i in y:
#    if 'NL/118/2001' in i or 'BI/16190/1968' in i:
#        print i, y[i]