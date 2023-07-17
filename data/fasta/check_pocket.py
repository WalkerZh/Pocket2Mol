import csv
import os
import sys

from tqdm import tqdm

AA_NAME_SYM = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

def csv_to_dict(file_path):
    '''
    {pdb_id :{'chain: (fasta, start)}}
    '''
    ret_dict = {}
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)
        for row in reader:
            pdb_id, chain, fasta, start = row[0], row[1], row[2], row[3]
            if pdb_id in ret_dict.keys():
                ret_dict[pdb_id][chain] = (fasta, start)
            else:
                ret_dict[pdb_id] = {}
                ret_dict[pdb_id][chain] = (fasta, start)
    # print(len(ret_dict))
    return ret_dict 

def extract_(pdb_file_path): # pocket_file_path, 
    '''
    get list of aa in pocket
    or get list of aa in complete .ent file(pdb)
    ATOM ()
    '''
    with open(pdb_file_path, "r", encoding="utf8") as fp:
        pdb_lines = [line.strip() for line in fp]
    pdb_list = []
    for line in pdb_lines:
        if line[0:6].strip() != 'ATOM': continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21:22].strip()
        res_id = int(line[22:26])
        pdb_list.append((atom_name, res_name, chain, res_id))

    return pdb_list # pocket_list, 

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

def check_pocket(pocket_list, pdb_list):
    '''
    given pocket and complete list
    check whether the pocket is from corresponding pdb
    '''
    print(f"Len(pocket file) = {len(pocket_list)}, Len()=.")
    for pp in pocket_list:
        # print(pp)
        if pp not in pdb_list:
            return False
    return True

if __name__ == "__main__":
    # pocket_file_path = "./1m4n_A_rec_1m7y_ppg_lig_tt_min_0_pocket10.pdb"
    # pdb_file_path = "./pdb_files/pdb1m4n.ent"
    
    # pocket_list, pdb_list = extract_(pocket_file_path), extract_(pdb_file_path)
    # pocket_set, pdb_set = set(pocket_list), set(pdb_list)

    # if pocket_set.issubset(pdb_set):
    #     print("Pass")
    # else: print("Failed")

    file_path = "./crossdocked2fasta_pocket10_v1.csv"
    pdb_dict = csv_to_dict(file_path)
    # print(f"Total data: {len(pdb_dict)}")

    crossdocked_dir = "../crossdocked_pocket10/"
    subdir_list = os.listdir(crossdocked_dir)

    no_pdb_file_list = []
    
    error_times = 0
    error_list = []
    
    no_pass = 0
    no_pass_list = []

    for subdir in tqdm(subdir_list):
        temp_path = os.path.join(crossdocked_dir, subdir)
        if os.path.isfile(temp_path): continue

        file_list = os.listdir(temp_path)
        for file_name in file_list:
            if file_name.endswith(".sdf"): continue

            pocket_file_path = os.path.join(temp_path, file_name)
            # print(pocket_file_path)
            pocket_list = extract_(pocket_file_path)

            pdb_id = pocket_file_path.split("/")[-1].split("_")[0].strip()
            chain = pocket_file_path.split("/")[-1].split("_")[1].strip()

            # if (pdb_id, chain) not in pdb_dict.keys():
            #     no_pdb_file_list.append((pdb_id, chain))
            #     continue
            if pdb_id not in pdb_dict.keys() or chain not in pdb_dict[pdb_id].keys():
                no_pdb_file_list.append((pdb_id, chain))
                continue  

            fasta_dict = pdb_dict[pdb_id]
            check_ = check_aa(pocket_list, fasta_dict)

            if check_ == 0:
                no_pass_list.append(pocket_file_path)
                no_pass += 1
            elif check_ == -1:
                error_list.append(pocket_file_path)
                error_times += 1

    print(f"No pass pockets: {no_pass}")
    print(f"Error times: {error_times}")
    
    no_pass_list = list(set(no_pass_list))
    error_list = list(set(error_list))
    no_pdb_file_list = list(set(no_pdb_file_list))

    with open("no_pass_pocket.txt", "w", encoding="utf8") as fw:
        print("No pass files:", file=fw)
        for no_p in no_pass_list:
            print(no_p, file=fw)

        print("-"*40, file=fw)
        print("Error files:", file=fw)
        for e in error_list:
            print(e, file=fw)
        
        print("-"*40, file=fw)
        print("No pdb files:", file=fw)
        for no_p in no_pdb_file_list:
            print(no_p, file=fw)

