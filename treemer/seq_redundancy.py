from collections import defaultdict

from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Phylo.BaseTree import Clade
from Bio.SeqRecord import SeqRecord

from treemer.util import TreeAlignmentMatch, TipSeqLinker
from treemer.errors import MethodNotFoundError, ThresholdError, SiteNotValidError


class Pruner(TreeAlignmentMatch):
    to_keep_methods = ["similar", "nearest", "consensus"]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.similarity = 0.95
        self.sites = []
        self.compared = {}

    def __str__(self) -> str:
        return 'similarity: {}\nsite: {}'.format(
            self.similarity,
            ' '.join(self.sites) if self.sites else None
        )

    def set_similarity(self, threshold: float = 1.0) -> None:
        threshold = float(threshold)
        if threshold <= 0:
            raise ThresholdError(threshold, "equal or less than", 0)
        elif threshold > 1:
            raise ThresholdError(threshold, "higher than", 1)
        self.similarity = threshold

    def set_site(self, sites: iter = None) -> None:
        self.sites = []
        if sites is not None:
            for site in sites:
                if site - 1 not in self.valid_sites:
                    raise SiteNotValidError(site)
                self.sites.append(site)

    def __select_and_prune(self, select: 'function') -> None:
        for clade, cluster in self.clstr.items():
            setattr(select(clade, cluster), 'to_keep', True)
            for dichord in cluster:
                if not dichord.to_keep:
                    self.tree.prune(dichord.tip)

    def __keep_the_similar(self, clade: Clade, cluster: list) -> TipSeqLinker:
        if len(cluster) == 1:
            (dichord,) = cluster
            return dichord
        elif len(cluster) == 2:
            dichords = {
                dichord: self.tree.distance(dichord.tip, clade)
                for dichord in cluster
            }
            return min(dichords, key=dichords.get)
        else:
            dichords = {
                dichord: min(
                    self.compared[pairing] for pairing in self.compared
                    if dichord.a_record.id in pairing
                ) for dichord in cluster
            }
            return max(dichords, key=dichords.get)

    def __keep_the_nearest(self, clade: Clade, cluster: list) -> TipSeqLinker:
        dichords = {
            dichord: self.tree.distance(dichord.tip, clade)
            for dichord in cluster
        }
        return min(dichords, key=dichords.get)

    def __get_the_consensus(self, clade: Clade, cluster: list) -> TipSeqLinker:
        rep_dichord = self.__keep_the_nearest(clade, cluster)
        description = ', '.join(dichord.a_record.id for dichord in cluster)
        self.aligns_as_seqs.update({
            rep_dichord.a_record.id: SeqRecord(
                seq=AlignInfo.SummaryInfo(
                    MultipleSeqAlignment(
                        dichord.a_record for dichord in cluster)
                ), id=clade.name, description=description
            )
        })
        rep_dichord.tip.name = description
        return rep_dichord

    def __qualified(self, cluster: list) -> bool:
        length = len(cluster)
        for i in range(length):
            query = cluster[i]
            query_content = query.site_content(self.sites)
            for j in range(i + 1, length):
                subject = cluster[j]
                subject_content = subject.site_content(self.sites)
                pairing = (query.a_record.id, subject.a_record.id)
                if pairing in self.compared:
                    similarity = self.compared[pairing]
                else:
                    similarity = query.compare(subject)
                    self.compared[pairing] = similarity
                if similarity < self.similarity or query_content != subject_content:
                    return False
        return True

    def trim_tree(self, to_keep: str = "similar") -> None:
        if to_keep is "similar":
            select = self.__keep_the_similar
        elif to_keep is "nearest":
            select = self.__keep_the_nearest
        elif to_keep is "consensus":
            select = self.__get_the_consensus
        else:
            raise MethodNotFoundError(to_keep, self.to_keep_methods)
        while True:
            old_clstr = self.clstr
            self.clstr = defaultdict(list)
            for dichord in self.dichord_list:
                self.clstr[dichord.next_clade].append(dichord)
            if len(self.clstr) == len(old_clstr):
                del old_clstr
                self.clstr = defaultdict(list)
                for dichord in self.dichord_list:
                    self.clstr[dichord.clade].append(dichord)
                break
            previous_clustering = list(old_clstr.values())
            for cluster in self.clstr.values():
                if cluster not in previous_clustering and self.__qualified(cluster):
                    for dichord in cluster:
                        dichord.proceed()
                elif cluster in previous_clustering:
                    previous_clustering.remove(cluster)
        self.__select_and_prune(select)
