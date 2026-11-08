from enum import StrEnum

import numpy as np
from Bio.Align import substitution_matrices
from itertools import combinations
import numpy.typing as npt

from enum import StrEnum

from Bio.Align import PairwiseAligner


class Mode(StrEnum):
    LOCAL = 'local'
    GLOBAL = 'global'
    FOGSAA = 'fogsaa'


def get_aligner(mode: Mode):
    return PairwiseAligner(mode=mode.value, substitution_matrix=substitution_matrices.load('BLOSUM62'))


def self_score(seq: str, identity_scores: substitution_matrices.Array) -> float:
    """Return the similarity score that would be obtained if seq was aligned to itself"""
    return sum(identity_scores[c] for c in seq)


def get_adj_mat(seqs: list[str], aligner: PairwiseAligner, similarity_threshold: float=0.8) -> npt.NDArray[np.bool_]:
    """
    Create boolean adjacency matrix by binarizing the normalized similarity matrix at `similarity_threshold`.
    """
    n = len(seqs)
    mat = np.zeros((n, n))
    identity_scores = np.diag(aligner.substitution_matrix)
    for i in range(n): mat[i, i] = self_score(seqs[i], identity_scores )
    d = np.diag(mat)
    norm = np.sqrt(np.outer(d,d))
    for i, j in combinations(range(n), r=2):
        a = seqs[i]
        b = seqs[j]
        mat[i, j] = mat[j, i] = aligner.score(a, b)
    return (mat / norm) > similarity_threshold
