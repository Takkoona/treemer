from collections import defaultdict

from Bio import AlignIO, Phylo, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree, Clade

from treemer.errors import TipNotMatchedError


class TipSeqLinker(object):
    def __init__(self, a_record: SeqRecord, t_path: tuple) -> None:
        self.a_record = a_record
        self.t_path = t_path
        self._path_pos = -1
        self.to_keep = False

    def __str__(self) -> str:
        return '{}*'.format(self.a_record.id) if self.to_keep else self.a_record.id

    def __repr__(self) -> str:
        return '{} -> {}'.format(self.a_record.id, abs(self._path_pos))

    def __iter__(self) -> iter:
        return iter(self.a_record)

    def get_pos(self):
        return self._path_pos

    def proceed(self) -> None:
        pro_pos = self._path_pos - 1
        if len(self.t_path) >= abs(pro_pos):
            self._path_pos = pro_pos

    def compare(self, other: 'TipSeqLinker') -> float:
        match, length = 0.0, 0.0
        for i, j in zip(self.a_record, other.a_record):
            if i == j != '-':
                match += 1.0
            elif i == j == '-':
                continue
            length += 1.0
        return match / length

    def site_content(self, sites: iter) -> tuple:
        return tuple(self.a_record[site - 1] for site in sites)

    @property
    def next_clade(self) -> Clade:
        pro_pos = self._path_pos - 1
        return self.t_path[pro_pos] if len(self.t_path) >= abs(pro_pos) else self.clade

    @property
    def clade(self) -> Clade:
        return self.t_path[self._path_pos]

    @property
    def tip(self) -> Clade:
        return self.t_path[-1]


class TreeAlignmentMatch(object):
    def __init__(self, tree: Tree, aligns_as_seqs: dict) -> None:
        self.tree = tree
        self.aligns_as_seqs = aligns_as_seqs
        (
            self.clstr,
            self.dichord_list,
            self.converted_names,
            self.valid_sites
        ) = self.check_matching()

    def __repr__(self) -> str:
        return '\n'.join(tip.name for tip in self.tree.get_terminals())

    def __len__(self) -> int:
        return self.tree.count_terminals()

    @classmethod
    def read(
            cls, aligns_handle: str, aligns_format: str,
            tree_handle: str, tree_format: str
    ) -> 'TreeAlignmentMatch':
        aligns = AlignIO.read(aligns_handle, aligns_format)
        tree = Phylo.read(tree_handle, tree_format)
        if aligns_format is 'fasta':
            aligns_as_seqs = SeqIO.to_dict(
                SeqIO.parse(aligns_handle, 'fasta')
            )
        else:
            from io import StringIO
            aligns_as_seqs = SeqIO.to_dict(
                SeqIO.parse(StringIO(aligns.format('fasta')), 'fasta')
            )
        return cls(tree, aligns_as_seqs)

    def check_matching(self) -> tuple:
        init_clstr = defaultdict(list)
        dichord_list = []
        converted_names = {}
        related_aligns = None
        for n, tip in enumerate(self.tree.get_terminals()):
            tip_name = tip.name
            try:
                seq_record = self.aligns_as_seqs[tip_name]
                dichord = TipSeqLinker(
                    seq_record,
                    (self.tree.root, *self.tree.get_path(tip))
                )
            except KeyError:
                raise TipNotMatchedError(tip)
            init_clstr[tip].append(dichord)
            dichord_list.append(dichord)
            new_seq_id = 'seq{}'.format(n)
            converted_names[tip_name] = new_seq_id
            converted_names[new_seq_id] = tip_name
            if related_aligns is None:
                related_aligns = MultipleSeqAlignment([seq_record])
            else:
                related_aligns.extend([seq_record])
        return (
            init_clstr, dichord_list, converted_names,
            tuple(range(related_aligns.get_alignment_length()))
        )
