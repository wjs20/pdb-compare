import shutil
import argparse
from pathlib import Path
import sys

from . import cluster
from . import align
from . import io


def handle_cluster(args) -> int:
    similarity_threshold = args.threshold
    mode = args.mode

    aligner = cluster.get_aligner(cluster.Mode(mode))

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

    records = cluster.seq2id(strucs)
    seqs = list(records)
    adj = cluster.get_adj_mat(seqs, aligner, similarity_threshold)
    clusters = cluster.get_clusters(adj)
    d = cluster.map_clusters2seqs(clusters, seqs)

    if args.display == 'clusters':
        canonicalized = cluster.fmt_canonical(d, records)
        for c in canonicalized: print(c)
        return 0
    elif args.display == 'sequences':
        for c, seqs in d.items():
            print('cluster: ', c)
            print('-'*30)
            for seq in seqs:
                print(seq)
            print('-'*30)
    return 0


def handle_align(args) -> int:
    if sys.stdin.isatty():
        lines = open(args.mapping_file).read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()
    records = align.parse_mappings(lines)
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

    align.align_structures(reference, reference_chain, targets)

    for t, _ in targets:
        name = Path(t.id).name
        io.write_structure(t, outdir / name)
    return 0


def get_user_args():
    parser = argparse.ArgumentParser(description='Align PDB structures using Biopython CEAligner.')

    parser.add_argument('-w', '--workers', type=int, default=1, help='Number of parallel workers to use', required=False)
    subparsers = parser.add_subparsers(required=True)

    cluster_parser = subparsers.add_parser('cluster', help='tools for clustering pdb structure chains')
    cluster_parser.add_argument('-d', '--display', choices=['sequences', 'clusters'])
    cluster_parser.add_argument('-t', '--threshold', type=float, default=0.8, help='Normalized similarity score threshold. Should be less than 1.', required=False)
    cluster_parser.add_argument('-m', '--mode', default='local', choices=['local', 'global'], help='Choose local or global sequence alignment', required=False)
    cluster_parser.add_argument('-i', '--input', type=Path, help='file containing pdb files in the format <pdb file>:<chain id>', required=False)
    cluster_parser.set_defaults(func=handle_cluster)

    align_parser = subparsers.add_parser('align', help='tools for aligning pdb structure chains')
    align_parser.add_argument('-f', '--mapping-file', help='File containing rows:<pdb_path>|<chain_id>')
    align_parser.add_argument('-o', '--outdir', default='aligned', help='Output directory.')
    align_parser.set_defaults(func=handle_align)

    return parser.parse_args()


def main():
    args = get_user_args()
    try:
        status = args.func(args)
    except Exception as e:
        raise SystemExit(str(e))
    else:
        raise SystemExit(status)


if __name__ == '__main__':
    main()
