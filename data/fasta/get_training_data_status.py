# new data format for esm feature
# check the status of all crossdocked2020 train and test data

import csv
import json
import os

import torch
from tqdm import tqdm

AA_NAME_SYM = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

crossdocked_data_path = "../crossdocked_pocket10"
pdb_files_path = "./pdb_files"
esm2_36_3B_feature_path_old = "/home/yinxia/blob1/v-weizhan/fasta_feature/esm2_t36_3B_UR50D_v0"
esm2_36_3B_feature_path_new = "/home/yinxia/blob1/v-weizhan/fasta_feature/esm2_t36_3B_UR50D"

# crossdocked(20288) -> pdb file(20287) -> fasta: 4abz / no_fasta / no_start / other extract problem
# pdb file -> esm2 feature: length(fasta) > 800 / no_fasta / no_start

crossdocked_fasta_file = "./crossdocked2fasta_pocket10_v1.csv"


def esm_feature_to_new_format(old_path, new_path):
    '''
    new format:
        | pdb_id | dict{chain: fasta} | start | ... | 
    '''
    dict_old = {}
    all_files_old = os.listdir(old_path)
    for files_old in tqdm(all_files_old):
        split_files_old = files_old.split(".")[0].split("_")
        pdb_id, chain = split_files_old[0], split_files_old[1]

        if pdb_id not in dict_old.keys():
            dict_old[pdb_id] = [chain]
        else:
            dict_old[pdb_id].append(chain)

    for pdb_id, chains in tqdm(dict_old.items()):
        dict_new = {}
        for chain in chains:
            pt_file_name = f"{pdb_id}_{chain}.pt"
            temp_load = torch.load(os.path.join(old_path, pt_file_name))
            dict_new[chain] = (temp_load[0], temp_load[1])
        torch.save(dict_new, os.path.join(new_path, f"{pdb_id}.pt"))
    return

def csv_to_dict(file_path):
    '''
    {pdb_id :{'chain: (fasta, start)}}
    with special requirements (length/start/fasta)
    '''
    ret_dict = {}
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)
        for row in reader:
            pdb_id, chain, fasta, start = row[0], row[1], row[2], row[3]
            # if len(fasta)>800 or fasta=="no_fasta" or start == "no_start": continue # ***
            if pdb_id in ret_dict.keys():
                ret_dict[pdb_id][chain] = (fasta, start)
            else:
                ret_dict[pdb_id] = {}
                ret_dict[pdb_id][chain] = (fasta, start)
    # print(len(ret_dict))
    return ret_dict 

def extract_keys(pdb_file_path): # pocket_file_path, 
    '''
    Get list of aa in pocket 
    Format: (atom_name, res_name, chain, res_id) e.g. ('N', 'TYR', 'A', '188')
    '''
    with open(pdb_file_path, "r", encoding="utf8") as fp:
        pdb_lines = [line.strip() for line in fp]
    pdb_list, chain_list = [], []
    for line in pdb_lines:
        if line[0:6].strip() != 'ATOM': continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21:22].strip()
        res_id = int(line[22:26])
        pdb_list.append((atom_name, res_name, chain, res_id))
        if chain not in chain_list:
            chain_list.append(chain)
    return pdb_list, chain_list

def check_aa(pdb_list, fasta_dict):
    # print(fasta)
    for pdb_ in pdb_list:
        _, res_name, res_chain, res_id = pdb_
        res_name_1 = AA_NAME_SYM[res_name]
        res_id = int(res_id)
        
        try:
            idx = res_id - int(fasta_dict[res_chain][1])
            if res_name_1 != fasta_dict[res_chain][0][idx]:
                return 0 # res not match
        except:
            return -1 # length not match
    return 1

def return_check_status(pocket_file_path, pdb_chain_dict):
    analysis_subdir = pocket_file_path.split("/")[-2].split("_")
    subdir_start, subdir_end = int(analysis_subdir[2]), int(analysis_subdir[3])
    analysis_file_name = pocket_file_path.split("/")[-1].split("_")
    pdb_id, main_chain = analysis_file_name[0].strip(), analysis_file_name[1].strip()
    
    # checking 1
    if pdb_id not in pdb_chain_dict.keys():
        return ("Error11", f"PDB id not exist in pdb dict!")
    pocket_list, chain_list = extract_keys(pocket_file_path)
    for chain in chain_list:
        if chain not in pdb_chain_dict[pdb_id].keys():
            return ("Error12", f"Chain {chain} not exist in pdb dict!")
        pdb_chain_info = pdb_chain_dict[pdb_id][chain]
        if pdb_chain_info[0] == "no_fasta":
            return ("Error13", f"Chain {chain} have no extracted fasta!")
        if pdb_chain_info[1] == "no_start":
            return ("Error14", f"Chain {chain} have no extracted start number!")
        if len(pdb_chain_info[0]) > 800:
            return ("Error15", f"Chain {chain} too long! {len(pdb_chain_info[0])}")
    
    # checking 2
    # pdb_chain_info = pdb_chain_dict[pdb_id][main_chain]
    # if int(pdb_chain_info[1]) != subdir_start:
    #     return ("Error21", f"Start num not match! {subdir_start}(in subdir) | {pdb_chain_info[1]}")
    # if len(pdb_chain_info[0]) != subdir_end-subdir_start+1:
    #     return ("Error22", f"Length of main fasta not match! {subdir_end-subdir_start+1} | {len(pdb_chain_info[0])}")

    # checking 3
    fasta_dict = pdb_chain_dict[pdb_id]
    for pdb_ in pocket_list:
        _, res_name_3, res_chain, res_id = pdb_
        res_name_1 = AA_NAME_SYM[res_name_3]
        res_id = int(res_id)
        chain_fasta = fasta_dict[res_chain][0]
        chain_start = int(fasta_dict[res_chain][1])
        idx = res_id - chain_start
        if idx < 0:
            return ("Error31", f"Chain {res_chain} res_id under the start! {res_id}-{chain_start}<0")
        if idx >= len(chain_fasta):
            return ("Error32", f"Chain {res_chain} res_id exceed the limit! {res_id}-{chain_start}>={len(chain_fasta)}")
        if res_name_1 != chain_fasta[idx]:
            return ("Error33", f"Chain {res_chain} res_name not match! {res_name_1} | {res_id}")

    return ("Pass", "")

def checking_data_status(crossdocked_data_path):
    '''
    For Crossdocked2020 dataset
    | protein file name | ligand file name | subdir | #path | pdb id | main chain | start | end | status |
    '''
    total_pass, total = 0, 0
    pdb_chain_dict = csv_to_dict(crossdocked_fasta_file)
    out_file_name = "crossdocked_data_status.csv"
    with open(out_file_name, 'w', newline='') as csvfile:
        fieldnames = ['protein_file_name', 'ligand_file_name', 'subdir', 'path', 'pdb_id', 'main_chain', 'start', 'end', 'status', 'reason']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        subdir_list = os.listdir(crossdocked_data_path)
        for subdir in tqdm(subdir_list):
            temp_path = os.path.join(crossdocked_data_path, subdir)
            if os.path.isfile(temp_path):
                continue

            file_list = os.listdir(temp_path)
            for file_name in file_list:
                if file_name.endswith(".sdf"):
                    continue
                pocket_file_path = os.path.join(temp_path, file_name)
                status, reason = return_check_status(pocket_file_path, pdb_chain_dict)

                total += 1
                if status == "Pass": total_pass += 1

                analysis_subdir = pocket_file_path.split("/")[-2].split("_")
                subdir_start, subdir_end = int(analysis_subdir[2]), int(analysis_subdir[3])
                analysis_file_name = pocket_file_path.split("/")[-1].split("_")
                pdb_id, main_chain = analysis_file_name[0].strip(), analysis_file_name[1].strip()
                ligand_file_name = "_".join(file_name.split(".")[0].split("_")[:-1]) + ".sdf"

                writer.writerow(
                    {'protein_file_name': file_name, 
                    'ligand_file_name': ligand_file_name, 
                    'subdir': subdir,
                    'path': pocket_file_path,
                    'pdb_id': pdb_id,
                    'main_chain': main_chain,
                    'start': subdir_start,
                    'end': subdir_end,
                    'status': status,
                    'reason': reason,
                    }
                )

    return total_pass, total

if __name__ == "__main__":
    esm_feature_to_new_format(esm2_36_3B_feature_path_old, esm2_36_3B_feature_path_new)
    # total_pass, total = checking_data_status(crossdocked_data_path)
    # print(f"| Total: {total} | Total pass: {total_pass} |")