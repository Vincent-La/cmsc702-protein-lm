from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.SeqUtils import seq1
import numpy as np
import os
import argparse
from transformers import AutoTokenizer, EsmForProteinFolding, EsmModel
from utils.model_names import ESM_FOLD, ESM_CONTACT_HEAD
from scripts.esm_inference import parse_output
import torch

'''
    Extract raw sequence string from structure file

    pdb_file: path to pdb structure file
'''
def seq_from_pdb(pdb_path):

    parser = PDBParser()
    # idt this structure name key matters
    protein_structure = parser.get_structure('p', pdb_path)

    aa_seq = ''

    for residue in next(protein_structure.get_chains()):
        # remove heteroatoms
        het = residue.get_full_id()[3][0]
        if het == ' ':
            aa_seq += residue.resname

    return seq1(aa_seq)


'''
    ESM inference on a pdb structure

    model_name: huggingface model string
    pdb_path: path to pdb structure file
    device: pytorch device
'''
def inference_from_pdb(model_name, pdb_path, device):
    
    # raw sequence from pdb file
    seq = seq_from_pdb(pdb_path)
    # Init tokenizer
    tokenizer =  AutoTokenizer.from_pretrained(model_name)
    batch = [seq]
    
    print('STARTING INFERENCE...')
    if model_name == ESM_FOLD:
        model = EsmForProteinFolding.from_pretrained(model_name)
        model.to(device)

        with torch.no_grad():
            inputs = tokenizer(
                batch,
                return_tensors='pt',
                padding=True,
                add_special_tokens=False
            )
            inputs.to(device)
            out = model(**inputs)
            contact_matrix = parse_output(out,batch)
            contact_matrix = contact_matrix[0]

            
    elif model_name == ESM_CONTACT_HEAD:
        model = EsmModel.from_pretrained(model_name)
        model.to(device)

        with torch.no_grad():
            inputs = tokenizer(batch, 
                                return_tensors = 'pt', 
                                padding = True,
                                add_special_tokens=True)
            
            inputs.to(device)
            contact_matrix = model.predict_contacts(inputs['input_ids'], inputs['attention_mask']).cpu().numpy()
            contact_matrix = contact_matrix[-1, :len(seq), :len(seq)]

    print('FINISHED')
    return contact_matrix
