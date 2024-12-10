from Bio import AlignIO
from Bio import SeqIO
import random
import sys

# This script takes a random subset of a fasta file and
# creates a new fasta file containing the subset.
# arg1 -> the path of the fasta file to get the subset from
# arg2 -> the number of sequences for the new fasta file
# arg3 -> the path of the new fasta file
def main():
    file_path = sys.argv[1]
    N = int(sys.argv[2])
    output_name = sys.argv[3]
    try:
        sequences = list(SeqIO.parse(file_path, "fasta"))
    except:
        print('could not parse file -> not in fasta format')

    random.shuffle(sequences)
    print('Original number of sequences: ', len(sequences))
    print('New number of sequences: ', N)
    SeqIO.write(sequences[0:N], output_name, "fasta")

if __name__ == "__main__":
    main()