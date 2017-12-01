# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 08:41:27 2017

@author: Chengyang
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
from StringIO import StringIO
from collections import defaultdict
import numpy as np
import os

class Trichord(object):
    
    def __init__(self, s_record, a_record, t_path):
        self.s_record = s_record
        self.a_record = a_record
        self.t_path = t_path
        self._path_pos = -1
        self.prsrv = False
        
    def __str__(self):
        s_id = self.s_record.id
        if self.prsrv is False:
            return  s_id
        else:
            return '{}*'.format(s_id)
    
    def proceed(self):
        pro_pos = self._path_pos - 1
        if len(self.t_path) >= abs(pro_pos):
            self._path_pos = pro_pos
    
    def set_prsrv(self):
        self.prsrv = True
    
    @property
    def next_clade(self):
        pro_pos = self._path_pos - 1
        if len(self.t_path) >= abs(pro_pos):
            return self.t_path[pro_pos]
        else:
            return self.clade
        
    @property
    def clade(self):
        return self.t_path[self._path_pos]
    
    @classmethod
    def from_list(cls, assembly):
        s_record, a_record, t_path = assembly
        return cls(s_record, a_record, t_path)


class Trinity(object):
    
    def __init__(self, seqs, aligns, tree):
        self.seqs = seqs
        self.aligns = aligns
        self.tree = tree
        self.aligned = {}
        self.similarity = 0.95
        self.sites = []
        self.level = None
        self.calc_fun = self.__similarity_by_msa
        
    def set_similarity(self, threshold):
        assert 0 <= threshold <= 1
        self.similarity = threshold
        
    def set_sites(self, *args):
        if args:
            for site in args:
                assert isinstance(args, int)
                self.sites.append(site)
        else:
            self.sites = []
    
    def set_level(self, level):
        assert level is None or isinstance(level, int) and level > 0
        self.level = level
        
    def use_msa(self):
        self.calc_fun = self.__similarity_by_msa
    
    def use_blastp(self):
        self.blast = NcbiblastpCommandline
        self.calc_fun = self.__similarity_by_blast
        
    def use_blastn(self):
        self.blast = NcbiblastnCommandline
        self.calc_fun = self.__similarity_by_blast
        
    def check_num(self):
        n_seq = self.seqs.__len__()
        n_align = self.aligns.__len__()
        n_tips = self.tree.count_terminals()
        assert n_seq == n_align, \
        "Different number of records in sequence and alignment"
        assert n_align == n_tips, \
        "Different number of records in alignment and tree"
        return n_seq
        
    def trim_by_tree(self):
        n_seq = self.check_num()
        tips = self.tree.get_terminals()
        trichords = []
        for a_index in range(n_seq):
            a_record = self.aligns[a_index]
            a_id = a_record.id
            assembly = [None, a_record, None]
            for s_id in self.seqs:
                if s_id.startswith(a_id):
                    s_record = self.seqs[s_id]
                    assert s_record is not None, \
                    "Alignment {} not found in sequence".format(a_id)
                    assembly[0] = s_record
            for tip in tips:
                t_id = tip.name
                if a_id.startswith(t_id):
                    t_path = self.tree.get_path(tip)
                    assert t_path is not None, \
                    "Tip {} not found in tree".format(a_id)
                    if isinstance(t_path, list):
                        t_path = tuple([self.tree.root] + t_path)
                    else:
                        t_path = tuple([self.tree.root, t_path])
                    assembly[2] = t_path
            trichord = Trichord.from_list(assembly)
            trichords.append(trichord)
        old_clstr = defaultdict(set)
        for trichord in trichords:
            clade = trichord.clade
            old_clstr[clade].add(trichord)
        level = self.level
        while level > 0 or self.level is None:
            if level is not None:
                level -= 1
            clstr_prvs = old_clstr
            trichords, old_clstr = self.__clustering(trichords, clstr_prvs)
            if old_clstr.values() == clstr_prvs.values():
                break
        final_clstr = defaultdict(set)
        for trichord in trichords:
            clade = trichord.clade
            final_clstr[clade].add(trichord)
        for cluster in final_clstr.values():
            shortest = None
            for trichord in cluster:
                tip = trichord.t_path[-1]
                dist = self.tree.distance(tip)
                if shortest is None:
                    shortest = dist
                    prsrv_record = trichord
                elif dist < shortest:
                    shortest = dist
                    prsrv_record = trichord
            prsrv_record.set_prsrv()
            yield cluster
        
    def __clustering(self, trichords, clstr_prvs):
        clstr_tmp = defaultdict(set)
        for trichord in trichords:
            clade = trichord.next_clade
            clstr_tmp[clade].add(trichord)
        prvs_clustering = clstr_prvs.values()
        for cluster in clstr_tmp.values():
            if cluster not in prvs_clustering:
                similar = self.calc_fun(cluster)
                site_consrv = self.__site_conservation(cluster)
                if all(similar + site_consrv):
                    for trichord in cluster:
                        trichord.proceed()
        return trichords, clstr_tmp
        
    def __site_conservation(self, cluster):
        consrv = []
        for site in self.sites:
            a_sites = [trichord.a_record[site - 1] for trichord in cluster]
            concord = all(x == a_sites[0] for x in a_sites)
            consrv.append(concord)
        return consrv
        
    def __similarity_by_msa(self, cluster):
        similar = []
        for trichord in cluster:
            qry_record = trichord.a_record
            qry_id = qry_record.id
            qry_seq = qry_record.seq
            qry_len = float(len(trichord.s_record))
            qry = np.array(list(qry_seq))
            for trichord in cluster:
                sbjct_record = trichord.a_record
                sbjct_id = sbjct_record.id
                pairing = tuple([qry_id, sbjct_id])
                if pairing in self.aligned:
                    p_match = self.aligned[pairing]
                else:
                    sbjct_seq = sbjct_record.seq
                    sbjct = np.array(list(str(sbjct_seq).replace('-', '_')))
                    p_match = np.sum(qry == sbjct)/qry_len
                    self.aligned[pairing] = p_match
                similar.append(p_match > self.similarity)
        return similar
        
    def __similarity_by_blast(self, cluster):
        directory = "./.Treemer_tmp/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        similar = []
        for trichord in cluster:
            qry_record = trichord.s_record
            qry_id = qry_record.id
            qry_seq = SeqRecord(qry_record.seq, id=qry_id)
            qry_file = directory + "query"
            SeqIO.write(qry_seq, qry_file, "fasta")
            for trichord in cluster:
                sbjct_record = trichord.s_record
                sbjct_id = sbjct_record.id
                pairing = tuple(set([qry_id, sbjct_id]))
                sbjct_seqs = []
                if pairing in self.aligned:
                    p_match = self.aligned[pairing]
                else:
                    sbjct_seq = SeqRecord(sbjct_record.seq,id=sbjct_id)
                    sbjct_seqs.append(sbjct_seq)
                    sbjct_file = directory + "subject"
                    SeqIO.write(sbjct_seq, sbjct_file, "fasta")
                    blast_cline = self.blast(query=qry_file,
                                             subject=sbjct_file,
                                             outfmt=5)
                    stdout, stderr = blast_cline()
                    blast_records = NCBIXML.parse(StringIO(stdout))
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            length = float(alignment.length)
                            for hsp in alignment.hsps:
                                p_match = hsp.identities/length
                                self.aligned[pairing] = p_match
                    similar.append(p_match > self.similarity)
            return similar

if __name__ == '__main__':
    from Bio import Phylo, AlignIO
    seqs = SeqIO.index("./dummy/test.fasta", 'fasta')
    aligns = AlignIO.read("./dummy/test.aln", 'clustal')
    tree = Phylo.read("./dummy/test.dnd", 'newick')
    
    trinity = Trinity(seqs, aligns, tree)
    clstr = trinity.trim_by_tree()