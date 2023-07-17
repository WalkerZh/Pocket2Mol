import os

import torch
import torch.nn as nn


class Insert_esm_feature(nn.Module):
    def __init__(self, esm_feature_dir, esm_input_dim=2560):
        super(Insert_esm_feature, self).__init__()
        self.esm_feature_dir = esm_feature_dir
        self.esm_input_dim = esm_input_dim
        self.esm_mlp = Simple_MLP(input_dim=esm_input_dim)

    def forward(self, h_compose, protein_pdb_id, protein_atom_to_res_id, protein_atom_to_chain,
                        protein_length, idx_protein, device):
        # insert the 256 dim feature into `h_compose`
        # processed_esm_feature = torch.zeros(h_compose[0].size(0), 256) # 256 must be passed by config

        h_compose_0 = h_compose[0]
        feature_box = self.load_esm_feature(protein_pdb_id, protein_atom_to_res_id, protein_atom_to_chain,
                                            protein_length, device, ).to(device)
        esm_feature = torch.zeros(h_compose[0].size(0), self.esm_input_dim).to(device)
        # esm_feature.size() = [compose_len, 2560]
        # assert esm_feature_.size() == h_compose_0.size()
        protein_or_ligand = [True if idx in idx_protein else False for idx in range(0, h_compose_0.size(0))]
        
        idx_feature_box = 0
        for idx, p_or_l in enumerate(protein_or_ligand):
            if p_or_l: # idx in idx_protein
                esm_feature[idx] = feature_box[idx_feature_box]
                idx_feature_box += 1

        h_compose_0 += self.esm_mlp(esm_feature)
        h_compose[0] = h_compose_0
        return h_compose
    
    def load_esm_feature(self, protein_pdb_id, protein_atom_to_res_id, protein_atom_to_chain,
                         protein_length, device, ):
        esm_feature_dir = self.esm_feature_dir
        feature_box = torch.zeros(protein_atom_to_res_id.size(0), self.esm_input_dim).to(device)
        next_protein_start = 0
        for protein_len, pdb_id, protein_chain in zip(protein_length, protein_pdb_id, protein_atom_to_chain):
            esm_feature_filename = os.path.join(esm_feature_dir, f"{pdb_id}.pt")
            esm_temp = torch.load(esm_feature_filename)                 
            protein_res_id = protein_atom_to_res_id[next_protein_start: next_protein_start+protein_len]

            '''
            Modified the .pt file format (maybe faster?)

            new format: | concat tensor(esm_feature) | {"chain": tensor_idx}(dict_idx) | {"chain": start}(dict_start) |
            all_start = [int(dict_start[chain]) for chain in protein_chain]
            all_idx = [int(dict_idx[chain]) for chain in protein_chain]
            protein_res_id = protein_res_id - all_start + all_idx
            feature_box[next_protein_start: next_protein_start+protein_len] = torch.index_select(esm_feature, 0, protein_res_id.clone().detach())
            '''

            for idx, res_id_and_chain in enumerate(zip(protein_res_id, protein_chain)):
                res_id, chain = res_id_and_chain
                esm_feature, start = esm_temp[chain][0].to(device), int(esm_temp[chain][1])
                res_id -= start
                feature_box[next_protein_start+idx] = esm_feature[res_id]

            next_protein_start += protein_len
        
        return feature_box

class Simple_MLP(nn.Module):
    def __init__(self, input_dim=2560, hidden_dim=512, output_dim=256): # 15B:5210 / 3B2560
        super(Simple_MLP, self).__init__()
        self.linear1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.ReLU()
        self.linear2 = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        x = self.linear1(x)
        x = self.relu(x)
        x = self.linear2(x)
        return x
