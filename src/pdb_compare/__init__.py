import argparse
from pathlib import Path

from . import cli


def get_user_args():
    parser = argparse.ArgumentParser(description='Align PDB structures using Biopython CEAligner.')

    parser.add_argument('-w', '--workers', type=int, default=1, help='Number of parallel workers to use', required=False)
    subparsers = parser.add_subparsers(required=True)

    cluster_parser = subparsers.add_parser('cluster', help='tools for clustering pdb structure chains')
    cluster_parser.add_argument('-d', '--display', choices=['sequences', 'clusters'])
    cluster_parser.add_argument('-t', '--threshold', type=float, default=0.8, help='Normalized similarity score threshold. Should be less than 1.', required=False)
    cluster_parser.add_argument('-m', '--mode', default='local', choices=['local', 'global'], help='Choose local or global sequence alignment', required=False)
    cluster_parser.add_argument('-i', '--input', type=Path, help='file containing pdb files in the format <pdb file>:<chain id>', required=False)
    cluster_parser.set_defaults(func=cli.handle_cluster)

    align_parser = subparsers.add_parser('align', help='tools for aligning pdb structure chains')
    align_parser.add_argument('-f', '--mapping-file', help='File containing rows:<pdb_path>|<chain_id>')
    align_parser.add_argument('-o', '--outdir', default='aligned', help='Output directory.')
    align_parser.set_defaults(func=cli.handle_align)

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
