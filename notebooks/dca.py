from Bio import AlignIO
import os

# get the cadherin data
cadherin_msa_path = os.path.join('..', 'data', 'cadherin', 'PF00028.alignment.seed')
alignment = AlignIO.read(cadherin_msa_path, 'stockholm')

# obtain the MSA alignment
seqs = []
for record in alignment:
    seqs.append(str(record.seq).replace('-', ''))
    # seqs.append(record.seq)

# DCA script below