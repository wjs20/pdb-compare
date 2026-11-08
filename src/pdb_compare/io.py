from collections.abc import Iterable
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Structure import Structure
from tqdm import tqdm

from concurrent.futures import ProcessPoolExecutor


def parse_structure(path: Path):
    parser = PDBParser(QUIET=True)
    return str(path), parser.get_structure(str(path), path)


def get_strucs(paths: Iterable[Path], workers=1) -> dict[str, Structure]:
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = pool.map(parse_structure, paths)
        return dict(tqdm(results, total=len(list(paths)), desc="Loading structures"))


def write_structure(structure: Structure, path: Path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path))
