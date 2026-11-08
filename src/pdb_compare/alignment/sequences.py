from enum import StrEnum
import os
from concurrent.futures import ProcessPoolExecutor

import argparse
from operator import itemgetter
from typing import LiteralString
import string
import numpy as np
from pathlib import Path
from pprint import pprint
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Structure import Structure
from Bio.Align import substitution_matrices
from tqdm import tqdm
from itertools import combinations, islice
from collections import defaultdict, Counter
import numpy.typing as npt

from scipy.sparse import csr_array, csgraph
from enum import StrEnum
from typing import TypedDict
from contextlib import contextmanager
from Bio.PDB.Structure import Structure
from collections.abc import Iterator, Iterable
from pathlib import Path
import argparse
import tempfile
from shutil import copy2
from Bio.PDB.cealign import CEAligner

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
