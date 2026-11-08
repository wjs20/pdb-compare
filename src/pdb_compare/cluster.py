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
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from tqdm import tqdm
from itertools import combinations, islice
from collections import defaultdict, Counter

from scipy.sparse import csr_array, csgraph

parser = PDBParser()
ppbuilder = PPBuilder()

np.set_printoptions(linewidth=180, suppress=True)

AA_CODES_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "UNK": "X",
    "GAP": "-",
}

def extract_sequence(chain):
    return ''.join(AA_CODES_3TO1.get(res.resname, '') for res in chain)


def self_score(seq: str, identity_scores: substitution_matrices.Array) -> float:
    """Return the similarity score that would be obtained if seq was aligned to itself"""
    return sum(identity_scores[c] for c in seq)


def get_adj_mat(seqs: list[str], aligner: PairwiseAligner, similarity_threshold: float=0.8):
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


def seq2id(strucs: dict[str, Structure]) -> dict[str, list[tuple[str, str]]]:
    """
    strucs: dict of Structures
        e.g. {'8ja3.pdb': <Structure id=8ja3.pdb>, '8wrz.pdb': <Structure id=8wrz.pdb>}

    adj:
        [[ True False False False  True False False False]
        [False  True False False False False False False]
        ...
        [False False  True False False False  True False]
        [False False False  True False False False  True]]


    returns:
        {'DIQ...': [('8wrz.pdb', 'L')],
         'DIQ...': [('8ja3.pdb', 'M')],
         'EELRHC': [('8ja3.pdb', 'U')],
         'RVF...': [('8ja3.pdb', 'B')],
         'SCALAK': [('8wrz.pdb', 'V')],
         'TRV...': [('8wrz.pdb', 'A')],
         'VQL...': [('8wrz.pdb', 'H')],
         'VQL...': [('8ja3.pdb', 'I')]}
    """
    struc_seqs = {}
    for astruc in strucs.values():
        chains = list(astruc.get_chains())
        struc_seqs[astruc.id] = {chain.id: extract_sequence(chain) for chain in chains}

    records = defaultdict(list)
    for struc, seqs in struc_seqs.items():
        for chain, seq in seqs.items(): records[seq].append((struc, chain))
    return dict(records)


def get_clusters(adj: np.ndarray) -> list[int]:
    """
    returns: list of length == len(adj) containing numeric chain cluster ids.
    """
    csr = csr_array(adj)
    return csgraph.connected_components(csr)[1].tolist()


def get_canonical_chains(clusters: list[int]) -> list[LiteralString]:
    """
    returns: list of length == len(adj) containing new canonical chain ids.
    ['A', 'B', 'C', 'D', 'A', 'E', 'C', 'D']
    """
    assert (n := len(set(clusters))) < len(string.ascii_uppercase), f'{n} distinct entities found, ran out of chain ids'
    return [string.ascii_uppercase[i] for i in clusters]



def map_clusters2seqs(clusters: list[int], seqs: list[str]) -> dict[int, list[str]]:
    """
    takes list of cluster ids and list of sequences and creates mapping from cluster id to the sequences associated with that cluster id.

    Args:
        clusters e.g. [0, 1, 2, 3, ... 3, 0, 39]
        sequences: e.g. ['TRV...', 'SCA...', 'DIQ...', 'TRV...', ...  'DIQ...', 'TRV...', 'EEYAL']

    Returns:
        map from cluster to list of sequences {0: ['TRV...', 'RVF...', 'TRV...', ... ], 1: ['SCALAK'], 2: ['VQL...', 'VQL...', 'VQL...',... ], 3: ['DIQ...', 'DIQ...', 'DIQ...', ..., ], 4: ['EELRHC'], 5: ['EVQ...', 'VQL...', 'VQL...',... ], ..., 38: ['TRV...', 'VFK...'], 39: ['EEYAL']})
    """
    d = defaultdict(list)
    for clust, seq in zip(clusters, seqs, strict=True): d[clust].append(seq)
    return dict(d)
