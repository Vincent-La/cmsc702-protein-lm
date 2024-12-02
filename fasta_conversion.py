from dca.dca_class import dca
from Bio import AlignIO
from Bio import SeqIO
import os
import sys

def main():
    file_path = sys.argv[1]
    output_name = sys.argv[2]
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