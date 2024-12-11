from dca.dca_class import dca
import pandas as pd
import os
import sys
import pickle

# This script takes in a fasta file as input (arg 1) and the output file name
# (arg 2) and returns a csv file containing the mfdcas.
#
# NOTE: the input fasta file should be small ~ 2000 sequences.
# You should first take a subset of the entire MSA file (500000 sequences total)
# and create the corresponding fasta file. Then use the new fasta file containing
# 2000 sequences as input for this scrpt to generate the .csv file.
def main():
    # DCA script
    sequence_file = sys.argv[1]
    output_filename = sys.argv[2]
    protein_family = dca(sequence_file)

    df = pd.DataFrame(protein_family.DI)
    df.rename(columns={0: 'first_site', 1: 'second_site', 2: 'mf_dca'}, inplace=True)
    df['first_site'] = df['first_site'].astype(int)
    df['second_site'] = df['second_site'].astype(int)
    df.to_csv(output_filename+'.csv', index=False)

if __name__ == "__main__":
    main()