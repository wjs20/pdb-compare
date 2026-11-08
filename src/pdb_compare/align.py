#!/usr/bin/env python3
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

def read_mapping_file(mapping_file):
    """
    Mapping file format:

        path/to/structure1.pdb|A
        path/to/structure2.pdb|B
        path/to/structure3.pdb|A
    """
    records = []
    with open(mapping_file) as fh:
        for line_number, line in enumerate(fh, start=1):
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

def parse_args():
    parser = argparse.ArgumentParser(description='Align PDB structures using Biopython CEAligner.')
    parser.add_argument('-m', '--mapping-file', help='File containing rows:\n<pdb_path>|<chain_id>')
    parser.add_argument('-o', '--outdir', default='aligned', help='Output directory.')
    return parser.parse_args()

def main():
    args = parse_args()
    records = read_mapping_file(args.mapping_file)
    reference_pdb = records[0]['pdb_path']
    reference_chain = records[0]['chain_id']
    targets = records[1:]
    print(f"Using reference: {reference_pdb}{(f' (chain {reference_chain})' if reference_chain else '')}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    copy2(reference_pdb, outdir / reference_pdb.name)
    print(f"{'Structure':<50}{'Chain':<8}{'RMSD (Å)':>12}")
    print('-' * 72)
    for record in targets:
        pdb_path = record['pdb_path']
        chain_id = record['chain_id']
        output_pdb = outdir / pdb_path.name
        try:
            rmsd = align_to_reference(ref_pdb=reference_pdb, target_pdb=pdb_path, output_pdb=output_pdb, ref_chain=reference_chain, target_chain=chain_id)
            print(f'{pdb_path.name:<50}{str(chain_id):<8}{rmsd:12.3f}')
        except Exception as exc:
            print(f'{pdb_path.name:<50}{str(chain_id):<8}ERROR: {exc}')
if __name__ == '__main__':
    main()
