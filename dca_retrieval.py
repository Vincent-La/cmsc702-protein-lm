from dca.dca_class import dca
import pandas as pd
import os
import sys
import pickle

def main():
    # DCA script
    sequence_file = sys.argv[1]
    output_filename = sys.argv[2]
    protein_family = dca(sequence_file)

    df = pd.DataFrame(protein_family.DI)
    df.rename(columns={0: 'first_site', 1: 'second_site', 2: 'mf_dca'}, inplace=True)
    df.to_csv(output_filename+'.csv')

if __name__ == "__main__":
    main()