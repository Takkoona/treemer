# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 17:23:24 2017

@author: Chengyang
"""
from Bio.Alphabet import NucleotideAlphabet, ProteinAlphabet
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from collections import defaultdict
import numpy as np

class Dichord(object):
    
    def __init__(self, a_record, t_path):
        self.a_record = a_record
        self.t_path = t_path
        self._path_pos = -1
        self.prsrv = False
        
    def __str__(self):
        a_id = self.a_record.id
        if self.prsrv is False:
            return  a_id
        else:
            return '{}*'.format(a_id)
    
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
        
    @property
    def tip(self):
        return self.t_path[-1]
    
    @classmethod
    def from_list(cls, assembly):
        return cls(*assembly)

     
class Binity(object):
    
    def __init__(self, aligns, tree):
        self.aligns = aligns
        self.tree = tree
        assert isinstance(aligns, MultipleSeqAlignment), \
        "The first input is not an instance of MultipleSeqAlignment"
        assert isinstance(tree, Tree), \
        "The second input is not an instance of Tree"
        assert len(self) is not None
        
    def __len__(self):
        n_align = self.aligns.__len__()
        n_tips = self.tree.count_terminals()
        try:
            assert n_align == n_tips
            return n_align
        except:
            "Different record/tip number: {} in alignment and {} tree".format(n_align, n_tips)
        
    @property
    def seq_type(self):
        alphabet = self.aligns._alphabet
        if alphabet is None:
            return None
        elif isinstance(alphabet, NucleotideAlphabet):
            return 'nucleotide'
        elif isinstance(alphabet, ProteinAlphabet):
            return 'protein'
        
    def get_dichords(self):
        n_seq = len(self)
        for a_index in range(n_seq):
            a_record = self.aligns[a_index]
            a_id = a_record.id
            t_path = self.tree.get_path(a_id)
            assert t_path is not None, \
            "Alignment record {} not found in the tree".format(a_id)
            if isinstance(t_path, list):
                t_path = tuple([self.tree.root] + t_path)
            else:
                t_path = tuple([self.tree.root, t_path])
            assembly = [a_record, t_path]
            dichord = Dichord.from_list(assembly)
            yield dichord


class TraversePaths(object):
        
    def __init__(self, binity):
        self.binity = binity
        self.aligned = {}
        self.similarity = 0.95
        self.sites = []
        self.level = None
        self.select_prsrv = self.__keep_the_similar
        
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
            
    def set_selection(self, method = 'similar'):
        if method is 'nearest':
            self.select_prsrv = self.__keep_the_nearest
        elif method is 'similar':
            self.select_prsrv = self.__keep_the_similar
        else:
            raise Exception("Method only accepts nearest or similar")
    
    def set_level(self, level = None):
        assert level is None or isinstance(level, int) and level > 0
        self.level = level

    def trim_by_tree(self):
        dichords = []
        old_clstr = defaultdict(set)
        for dichord in self.binity.get_dichords():
            clade = dichord.clade
            old_clstr[clade].add(dichord)
            dichords.append(dichord)
        level = self.level
        stage = 0
        while level > 0 or self.level is None:
            if level is not None:
                level -= 1
            stage += 1
            print "Doing level {} reduction".format(stage)
            clstr_prvs = old_clstr
            dichords, old_clstr = self.__clustering(dichords, clstr_prvs)
            if old_clstr.values() == clstr_prvs.values():
                print "\nLevel {} reduction is the same as level {}\n".format(stage, stage - 1)
                break
        else:
            print "\nReduction by level done\n"
        final_clstr = defaultdict(set)
        for dichord in dichords:
            clade = dichord.clade
            final_clstr[clade].add(dichord)
        del dichords
        print "Pruning clusters...\n"
        for cluster in final_clstr.values():
            cluster = self.select_prsrv(cluster)
            for dichord in cluster:
                if not dichord.prsrv:
                    self.binity.tree.prune(dichord.tip)
            yield cluster
        
    def __clustering(self, dichords, clstr_prvs):
        clstr_tmp = defaultdict(set)
        for dichord in dichords:
            clade = dichord.next_clade
            clstr_tmp[clade].add(dichord)
        prvs_clustering = clstr_prvs.values()
        for cluster in clstr_tmp.values():
            if cluster not in prvs_clustering:
                similar = self.__similarity_by_msa(cluster)
                site_consrv = self.__site_conservation(cluster)
                if all(similar + site_consrv):
                    for dichord in cluster:
                        dichord.proceed()
        return dichords, clstr_tmp
        
    def __site_conservation(self, cluster):
        consrv = []
        for site in self.sites:
            a_sites = [dichord.a_record[site - 1] for dichord in cluster]
            concord = all(x == a_sites[0] for x in a_sites)
            consrv.append(concord)
        return consrv
        
    def __similarity_by_msa(self, cluster):
        similar = []
        for dichord in cluster:
            qry_record = dichord.a_record
            qry_id = qry_record.id
            qry_seq = qry_record.seq
            qry_len = float(len(qry_seq.ungap('-')))
            qry = np.array(list(qry_seq))
            for dichord in cluster:
                sbjct_record = dichord.a_record
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
        
    def __keep_the_nearest(self, cluster):
        shortest = None
        for dichord in cluster:
            tip = dichord.tip
            dist = self.binity.tree.distance(tip)
            if shortest is None:
                shortest = dist
                prsrv_record = dichord
            elif dist < shortest:
                shortest = dist
                prsrv_record = dichord
        prsrv_record.set_prsrv()
        return cluster
        
    def __keep_the_similar(self, cluster):
        nearest = None
        for dichord in cluster:
            a_id = dichord.a_record.id
            min_s = None
            for pairing in self.aligned:
                if pairing[0] is a_id:
                    similarity = self.aligned[pairing]
                    if min_s is None:
                        min_s = similarity
                    elif similarity < min_s:
                        min_s = similarity
            if nearest is None:
                nearest = min_s
                prsrv_record = dichord
            elif min_s > nearest:
                nearest = min_s
                prsrv_record = dichord
        prsrv_record.set_prsrv()
        return cluster