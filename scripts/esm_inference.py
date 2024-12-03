import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, EsmForProteinFolding
from scipy.special import softmax
import numpy as np

from Bio import AlignIO
import os
import argparse
from tqdm import tqdm


# Custom pytorch dataset for holding raw sequence data
class SequenceDataset(Dataset):
    def __init__(self, msa_file):
        super().__init__()
        self.msa_file = msa_file
        # self.tokenizer = tokenizer
        self.seqs = self._get_sequences()
    
    def __len__(self):
        return len(self.seqs)
    
    def __getitem__(self, idx):
        item = self.seqs[idx]
        return item

    def _get_sequences(self):
        alignment = AlignIO.read(self.msa_file, 'stockholm')
        seqs = []
        for record in alignment:
            seqs.append(str(record.seq).replace('-', ''))
        
        return seqs

# adapted from ColabFold notebook: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb#scrollTo=CcyNpAvhTX6q
def parse_output(output, batch):
    # pae = (output["aligned_confidence_probs"][0].detach().cpu() * np.arange(64)).mean(-1) * 31
    # plddt = output["plddt"][0,:,1].detach().cpu()

    # bins = np.append(0,np.linspace(2.3125,21.6875,63))
    # sm_contacts = softmax(output["distogram_logits"].detach().cpu(),-1)[0]

    # # NOTE: pretty sure this just looks for all distances <8 angstroms?  
    # sm_contacts = sm_contacts[...,bins<8].sum(-1)
    # xyz = output["positions"][-1,0,:,1].detach().cpu()
    # mask = (output["atom37_atom_exists"][0,:,1] == 1).detach().cpu()
    # o = {"pae":pae[mask,:][:,mask],
    #     "plddt":plddt[mask],
    #     "sm_contacts":sm_contacts[mask,:][:,mask],
    #     "xyz":xyz[mask]}

    bins = np.append(0,np.linspace(2.3125,21.6875,63))
    sm_contacts = softmax(output["distogram_logits"].detach().cpu(),-1)
    sm_contacts = sm_contacts[...,bins<8].sum(-1)
    mask = (output["atom37_atom_exists"][:,:,1] == 1).detach().cpu()

    contacts = []

    # over each protein sequence in batch
    for batch_idx in range(sm_contacts.shape[0]):
        contacts_matrix = sm_contacts[batch_idx][mask[batch_idx],:][:,mask[batch_idx]]
        contacts.append(contacts_matrix)

        # assert that if protein sequence is size L, then contacts_matrix is LxL
        assert len(batch[batch_idx]) == contacts_matrix.shape[0]

    return contacts

# Compute contact matrix for each sequence in MSA using PL model 
def plm_inference_per_sequence(model, 
                               tokenizer, 
                               msa_file, 
                               batch_size,
                               device):

    # setup dataset and dataloader
    dataset = SequenceDataset(msa_file)
    eval_dataloader = DataLoader(dataset, batch_size = batch_size)
    
    # compute contact matrix for each sequence and append to a python list
    # b/c numpy does not support jagged matrices
    contact_matrices = []
    with torch.no_grad():
        for step, batch in enumerate(tqdm(eval_dataloader, f"PLM Inference on {os.path.basename(msa_file)}")):
            inputs = tokenizer(batch, 
                               return_tensors = 'pt', 
                               padding = True,
                            #    truncation=True,
                               add_special_tokens=False)

            inputs.to(device)

            out = model(**inputs)
            contacts = parse_output(out, batch)
            contact_matrices.extend(contacts)
    
    return contact_matrices

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--msa_file',
        type=str,
        required=True,
        help="Path of MSA (Stockholm format) file for inference",
    )

    parser.add_argument(
        '--plm_model',
        type = str,
        required=False,
        default = 'facebook/esmfold_v1',
        help = 'Huggingface model string for Protein Langauge Model'
    )

    parser.add_argument(
        '--batch_size',
        type = int,
        required=False,
        default=16,
        help = 'Batch size for PLM inference'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help = 'Output folder to store contact matrices'
    )

    parser.add_argument(
        '--job_name',
        type=str,
        required=True,
        help = 'Name of inference job (to name output files)'
    )

    args = parser.parse_args()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # Init model
    model = EsmForProteinFolding.from_pretrained(args.plm_model)
    model.to(device)

    # Init tokenizer
    tokenizer =  AutoTokenizer.from_pretrained(args.plm_model)

    # Inference
    contact_matrices = plm_inference_per_sequence(model, 
                                                  tokenizer,
                                                  args.msa_file,
                                                  args.batch_size,
                                                  device)
    
    print(f'Writing output to {args.output_dir}')
    # save out contact matrices
    for i,matrix in enumerate(contact_matrices):
        np.save(os.path.join(args.output_dir, f'{i}_{args.job_name}.npy'), matrix)
