from dca.dca_class import dca
import os
import sys
import pickle

def main():
    # DCA script
    sequence_file = sys.argv[1]
    output_filename = sys.argv[2]
    protein_family = dca(sequence_file)
    
    with open(output_filename, 'wb') as f:
    	pickle.dump(protein_family, f)


if __name__ == "__main__":
    main()