from Bio import AlignIO
from Bio import SeqIO
import argparse
import sys


# This script converts a stockhom formatted file into fasta, which
# is necessary for mfdca retrieval.
# arg1 -> stockhom file path
# arg2 -> path of file output
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--file_path',
        type=str,
        required=True,
        help="stockhom file path to convert to fasta",
    )

    parser.add_argument(
        '--output_path',
        type=str,
        required=True,
        help="path of fasta file output",
    )

    args = parser.parse_args()
    
    file_path = args.file_path
    output_name = args.output_path
    try:
        alignment = AlignIO.read(file_path, 'stockholm')
    except:
        print('could not read file -> not in stockholm format')
    
    # convert to fasta file
    with open(output_name, "w") as handle:
        count = SeqIO.write(alignment, handle, "fasta")
    handle.close()

if __name__ == "__main__":
    main()