from transformers import AutoTokenizer, EsmForProteinFolding, EsmModel
from utils.model_names import ESM_FOLD, ESM_CONTACT_HEAD
from scripts.esm_inference import parse_output
import torch
from utils.pdb import seq_from_structure
from Bio.PDB import PDBParser

'''
    ESM inference on a pdb structure

    model_name: huggingface model string
    pdb_path: path to pdb structure file
    device: pytorch device
'''
def inference_from_pdb(model_name, pdb_path, device):
    
    parser = PDBParser()
    # idt this structure name key matters
    structure = parser.get_structure('p', pdb_path)
    # raw sequence from pdb structure
    seq = seq_from_structure(structure)

    # Init tokenizer
    tokenizer =  AutoTokenizer.from_pretrained(model_name)
    batch = [seq]
    
    print('STARTING INFERENCE...')
    if model_name == ESM_FOLD:
        model = EsmForProteinFolding.from_pretrained(model_name)
        model.to(device)

        with torch.no_grad():
            inputs = tokenizer(batch,
                               return_tensors='pt',
                               padding=True,
                               add_special_tokens=False)
            
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


'''
    Calculate the precision@k given predicted contact matrix

    contact_matrix: square matrix of predicted probabilities 
    true_contacts: list of tuples containing true contacts 
    k: number of top-k scores to consider
    seq_separation_cutoff: minimum distance b/w residues to consider contacts

'''
def precision_at_k(contact_matrix, true_contacts, k, seq_separation_cutoff= 4):

    assert contact_matrix.shape[0] == contact_matrix.shape[1], 'contact matrix should be square!'
    L = contact_matrix.shape[0]

    paired_probs = {}
    for i in range(L):
        for j in range(i):
            if abs(i - j) > seq_separation_cutoff:
                paired_probs[(i,j)] = contact_matrix[i,j]

    # sorted key,val pairs (decending by prob)
    top_k_pairs = [k for k,v in sorted(paired_probs.items(), key = lambda item: -item[1])][:k]

    TP = 0
    for pair in top_k_pairs:
        if pair in true_contacts:
            TP += 1

    precision = TP / k
    return precision

'''
    Calculate the recall@k given predicted contact matrix

    contact_matrix: square matrix of predicted probabilities 
    true_contacts: list of tuples containing true contacts 
    k: number of top-k scores to consider
    seq_separation_cutoff: minimum distance b/w residues to consider contacts

'''
def recall_at_k(contact_matrix, true_contacts, k, seq_separation_cutoff= 4):

    assert contact_matrix.shape[0] == contact_matrix.shape[1], 'contact matrix should be square!'
    L = contact_matrix.shape[0]

    paired_probs = {}
    for i in range(L):
        for j in range(i):
            if abs(i - j) > seq_separation_cutoff:
                paired_probs[(i,j)] = contact_matrix[i,j]

    # sorted key,val pairs (decending by prob)
    pairs = [k for k,v in sorted(paired_probs.items(), key = lambda item: -item[1])]
    

    TP = 0
    for pair in true_contacts[:k]:
        if pair in pairs:
            TP += 1

    recall = TP / k
    return recall
