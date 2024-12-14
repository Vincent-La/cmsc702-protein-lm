from Bio import AlignIO
from Bio import SeqIO
import random
import sys
import argparse
import numpy as np
from difflib import SequenceMatcher

# This script takes a random subset of a fasta file and
# creates a new fasta file containing the subset.
# arg1 -> the path of the fasta file to get the subset from
# arg2 -> the number of sequences for the new fasta file
# arg3 -> the left index for bounding the range of sites
# arg4 -> the right index for bounding the range of sites
# arg5 -> the path of the new fasta file
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--file_path',
        type=str,
        required=True,
        help="the path of the fasta file to get the subset from",
    )

    parser.add_argument(
        '--N',
        type=int,
        required=True,
        help="the number of sequences for the new fasta file",
    )

    parser.add_argument(
        '--left',
        type=int,
        required=False,
        help="the number of sequences for the new fasta file",
    )

    parser.add_argument(
        '--right',
        type=int,
        required=False,
        help="the number of sequences for the new fasta file",
    )

    parser.add_argument(
        '--output_path',
        type=str,
        required=True,
        help="the path of the new fasta file",
    )
    
    args = parser.parse_args()
    
    file_path = args.file_path
    N = args.N
    L = args.left
    R = args.right
    output_name = args.output_path
    try:
        seen = set()
        records = []

        flag = True
        if L==None or R==None:
            flag = False
        for record in SeqIO.parse(file_path, 'fasta'): 
            if flag and record[L:R].seq not in seen:
                seen.add(record[L:R].seq)
                records.append(record[L:R])
            elif record.seq not in seen:
                seen.add(record.seq)
                records.append(record)       
                
    except:
        print('could not parse file -> not in fasta format')

    random.shuffle(records)
    SeqIO.write(records[0:N], output_name, "fasta")

if __name__ == "__main__":
    main()