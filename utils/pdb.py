''''
Utility functions for working with Biopython structure objects
'''
import numpy as np
from Bio.SeqUtils import seq1
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from typing import List, Dict


'''
    Extract raw sequence string from structure file

    structure: Biopython structure
'''
def seq_from_structure(structure:Structure):

    aa_seq = ''

    for residue in next(structure.get_chains()):
        # remove heteroatoms
        het = residue.get_full_id()[3][0]
        if het == ' ':
            aa_seq += residue.resname

    return seq1(aa_seq)


'''
    Extract list of Biopython Residue objects from pdb structure
'''
def residues_from_structure(structure:Structure):
    
    residues = []
    for residue in next(structure.get_chains()):
        # remove heteroatoms
        het = residue.get_full_id()[3][0]
        if het == ' ':
            residues.append(residue)

    return residues


'''
    return min Euclidean distance b/w residues res1 and res2
'''
def min_residue_dist(res1:Residue, res2:Residue):
    res1_atoms = [a for a in res1.get_atoms()]
    res2_atoms = [a for a in res2.get_atoms()]

    min_dist = float('inf')
    for a_i in res1_atoms:
        for a_j in res2_atoms:
            min_dist = min(min_dist, np.linalg.norm(a_i - a_j))
    
    return min_dist

'''
    get all pairwise distance for every residues in residues
'''
def get_pairwise_dists(residues: List[Residue]):
    pairwise_dists = {}

    for i in range(len(residues)):
        for j in range(i):
            #pairwise_dists[(i,j)] = np.linalg.norm(atom_coords[i] - atom_coords[j])
            pairwise_dists[(i,j)] = min_residue_dist(residues[i], residues[j])

    return pairwise_dists

'''
    Extract all residues in contact from pairwise_dists dict

    structure: pdb structure object
    contact_dist_cutoff: euclidean distance cutoff to consider residues to be "in contact" in angstroms
    seq_separation_cutoff: minimum separation between residue indices to consider
'''
def get_contacts_from_structure(structure:Structure, contact_dist_cutoff=8.0, seq_separation_cutoff=4):
    
    residues = residues_from_structure(structure)
    pairwise_dists = get_pairwise_dists(residues)

    contacts = []

    # NOTE: strictly i > j pairs
    for i in range(len(residues)):
        for j in range(i):
            if abs(i - j) > seq_separation_cutoff and pairwise_dists[(i,j)] <= contact_dist_cutoff:
                contacts.append((i, j))

    return contacts
