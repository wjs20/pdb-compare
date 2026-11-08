#!/usr/bin/env python3
from collections.abc import Iterator, Iterable
from pathlib import Path
import argparse
import tempfile
from shutil import copy2
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.cealign import CEAligner

class ChainSelect(Select):

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

def load_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(Path(pdb_file).stem, str(pdb_file))

def extract_chain_structure(structure, chain_id):
    """
    Create a temporary structure containing only a single chain.
    Used only for determining the alignment.
    """
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp.name, select=ChainSelect(chain_id))
        return load_structure(tmp.name)


def apply_transform(structure, rotation, translation):
    """
    Apply a rigid-body transformation to every atom
    in a structure.
    """
    for atom in structure.get_atoms():
        atom.transform(rotation, translation)


def read_mapping(lines: Iterable[str]) -> list[dict[str, str | Path]]:
    """
    Read chain mapping file in the format <filepath>|<chain id> e.g.
    The chain id should be the file specific reference to the chain that will be used for registration.

        example/input/8wrz.pdb|A
        example/input/8ja3.pdb|B
        ...
        example/input/8j8v.pdb|F
        example/input/7dfa.pdb|A
        example/input/8j97.pdb|A
    """

    records = []
    for line_number, line in enumerate(lines, start=1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            pdb_path, chain_id = line.split('|')
        except ValueError:
            raise ValueError(f'Invalid mapping on line {line_number}: {line}')
        records.append({'pdb_path': Path(pdb_path), 'chain_id': chain_id})
    return records


def align_to_reference(ref_pdb, target_pdb, output_pdb, ref_chain=None, target_chain=None):
    """
    Use CEAligner on selected chains but apply
    the resulting transform to the entire target
    structure.
    """
    reference_full = load_structure(ref_pdb)
    target_full = load_structure(target_pdb)
    reference_alignment = reference_full
    target_alignment = target_full
    if ref_chain:
        reference_alignment = extract_chain_structure(reference_full, ref_chain)
    if target_chain:
        target_alignment = extract_chain_structure(target_full, target_chain)
    aligner = CEAligner()
    aligner.set_reference(reference_alignment)
    aligner.align(target_alignment, transform=False)
    rotation, translation = aligner._rigid_motion
    apply_transform(target_full, rotation, translation)
    io = PDBIO()
    io.set_structure(target_full)
    io.save(str(output_pdb))
    return float(aligner.rms)
