import shutil
import argparse
from pathlib import Path
import sys

from . import clustering
from . import alignment
from . import io



def handle_align(args) -> int:
    if sys.stdin.isatty():
        lines = open(args.mapping_file).read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()
    records = alignment.parse_mappings(lines)
    referece_path = records[0]['pdb_path']
    reference_chain = records[0]['chain_id']
    targets = records[1:]
    print(f"Using reference: {referece_path}{(f' (chain {reference_chain})' if reference_chain else '')}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(referece_path, outdir / referece_path.name)

    _, reference = io.parse_structure(referece_path)

    path2chain = {t['pdb_path']: t['chain_id'] for t in targets}
    strucs = io.get_strucs(path2chain.keys(), args.workers)
    targets = [(struc, path2chain[Path(path)]) for path, struc in strucs.items()]

    alignment.align_structures(reference, reference_chain, targets)

    for t, _ in targets:
        name = Path(t.id).name
        io.write_structure(t, outdir / name)
    return 0


def handle_cluster(args) -> int:
    if sys.stdin.isatty():
        if args.input is None:
            raise argparse.ArgumentError(args.input, message='No input from stdin. Pass input via --input argument')
        lines = args.input.read_text().splitlines()
    else:
        lines = sys.stdin.read().splitlines()
    paths = [Path(line) for line in lines]
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"{path} not found")

    strucs = io.get_strucs(paths, args.workers)
    records = clustering.seq2id(strucs)
    seqs = list(records)
    clusters = clustering.cluster_records(seqs, args.mode, args.threshold)

    if args.display == 'clusters':
        canonicalized = clustering.fmt_canonical(clusters, records)
        for c in canonicalized: print(c)
        return 0
    elif args.display == 'sequences':
        for c, seqs in clusters.items():
            print('cluster: ', c)
            print('-'*30)
            for seq in seqs:
                print(seq)
            print('-'*30)
    return 0
