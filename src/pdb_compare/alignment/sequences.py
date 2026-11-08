from typing import TypedDict
from contextlib import contextmanager
from Bio.PDB.Structure import Structure
from collections.abc import Iterator, Iterable
from pathlib import Path
import argparse
import tempfile
from shutil import copy2
from Bio.PDB.cealign import CEAligner


def apply_transform(structure, rotation, translation):
    """
    Apply a rigid-body transformation to every atom
    in a structure.
    """
    for atom in structure.get_atoms():
        atom.transform(rotation, translation)


class ChainRecord(TypedDict):
    pdb_path: Path
    chain_id: str


def parse_mapping(line):
    try:
        pdb_path, chain_id = line.split('|')
    except ValueError:
        raise ValueError(f'Invalid mapping for: {line}')
    return dict(pdb_path=Path(pdb_path), chain_id=chain_id)


def parse_mappings(lines: Iterable[str]) -> list[ChainRecord]:
    """
    Read chain mapping file in the format <filepath>|<chain id> e.g.
    The chain id should be the file specific reference to the chain that will be used for registration.

        example/input/8wrz.pdb|A
        example/input/8ja3.pdb|B
        ...
        example/input/7dfa.pdb|A
        example/input/8j97.pdb|A
    """

    lines = map(str.strip, lines)
    lines = filter(lambda line: line and not line.startswith('#'), lines)
    return list(map(parse_mapping, lines))


@contextmanager
def partial_structure(structure, keep):
    """
    Context manager to return a partial chain view on an existing structure.
    Used in alignment algorithm to align based on subset of chains without having to write out and read a new structure.
    """
    if not keep:
        yield structure
        return

    model = structure[0]
    original = list(model.child_list)

    if keep not in model:
        raise KeyError(f"Chain {keep!r} not found")

    detached = [c for c in original if c.id != keep]

    try:
        for chain in detached:
            model.detach_child(chain.id)
        yield structure
    finally:
        for chain in detached:
            if chain.id not in model:
                model.add(chain)
        # included to retain input chain ordering
        model.child_list[:] = original


def align_to_reference(aligner, target_pdb, target_chain):
    with partial_structure(target_pdb, target_chain):
        aligner.align(target_pdb, transform=False)
        rotation, translation = aligner._rigid_motion
    apply_transform(target_pdb, rotation, translation)
    return float(aligner.rms)


def get_structure_aligner(reference_pdb: Structure, reference_chain: str) -> CEAligner:
    aligner = CEAligner()
    with partial_structure(reference_pdb, reference_chain):
        aligner.set_reference(reference_pdb)
    return aligner


def align_structures(reference: Structure, reference_chain: str, targets: list[tuple[Structure, str]]) -> None:
    aligner = get_structure_aligner(reference, reference_chain)
    print(f"{'Structure':<50}{'Chain':<8}{'RMSD (Å)':>12}")
    print('-' * 72)
    for t, c in targets:
        name = Path(t.id).name
        try:
            rmsd = align_to_reference(aligner, t, c)
            print(f'{name:<50}{str(c):<8}{rmsd:12.3f}')
        except Exception as exc:
            print(f'{name:<50}{str(c):<8}ERROR: {exc}')
