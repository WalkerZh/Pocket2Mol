import csv
import logging
import os
from io import StringIO

import requests
from Bio import SeqIO
from Bio.PDB import PDBList
from tqdm import tqdm

logging.basicConfig(  
    filename='pdb_id_error.log',
    level=logging.ERROR,  
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S', 
)


def get_all_fasta_from_pdb_id(pdb_id):
    pdb_file = f"./pdb_files/pdb{pdb_id}.ent"
    with open(pdb_file, "r") as file:
        pdb_content = file.read()

    pdb_handle = StringIO(pdb_content)
    fasta_handle = StringIO()
    SeqIO.convert(pdb_handle, "pdb-atom", fasta_handle, "fasta")

    fasta_sequence = fasta_handle.getvalue()

    if fasta_sequence:
        chain_fasta_dict = get_all_fasta(fasta_sequence)

        ret_list = []
        for k, v in chain_fasta_dict.items():
            chain, fasta = k, v
            start = get_start_idx(pdb_file, chain)
            if start == 12345678:
                if fasta: ret_list.append((pdb_id, chain, fasta, "no_start"))
                else: ret_list.append((pdb_id, chain, "no_fasta", "no_start"))

            else:
                if fasta: ret_list.append((pdb_id, chain, fasta, str(start)))
                else: ret_list.append((pdb_id, chain, "no_fasta", str(start)))    

        return ret_list

    else:
        print(pdb_id)
        return None



def get_all_fasta(fasta_sequence):
    fasta_list = fasta_sequence.split(">")[1:]
    temp_dict = {}
    for ff in fasta_list:
        ff = ff.strip()
        f_chain = str(ff.split("\n")[0].split(":")[1].strip())
        f_seq = "".join(ff.split("\n")[1:])
        temp_dict[f_chain] = f_seq
    return temp_dict


def get_start_idx(pdb_file, chain):
    with open(pdb_file, "r", encoding="utf8") as f:
        for line in f:
            line = line.strip()
            if line[0:6].strip() == 'ATOM' and line[21:22].strip() == chain:
                return int(line[22:26])
    return 12345678

if __name__ == "__main__":
    pdb_id_file = "./pdb_id.txt"
    with open(pdb_id_file, "r", encoding="utf8") as fr:
        pdb_id_list = [line.strip() for line in fr]
    
    csv_name = "crossdocked2fasta_pocket10.csv"

    error_pdb_list = []
    error_times = 0
    with open(csv_name, 'w', newline='') as csvfile:
        fieldnames = ['pdb_id', 'chain', 'fasta_sequence', 'start'] # 'protein_filename', 'ligand_filename' ,
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for pdb_id in tqdm(pdb_id_list):
            if pdb_id == "4abz": continue
            ret_list = get_all_fasta_from_pdb_id(pdb_id)
            if ret_list == None:
                error_times += 1
                error_pdb_list.append(pdb_id)
            else:
                for ret in ret_list:
                    pdb_id, chain, fasta, start = ret
                    writer.writerow({'pdb_id': pdb_id, 'chain': chain, 'fasta_sequence': str(fasta), 'start': start})
        
    
    print(f"Error times: {error_times}.")

    with open("error_pdb_id_from_pdbid2fasta.txt", "w", encoding="utf8") as fw:
        for e_p in error_pdb_list:
            print(e_p, file=fw)
    
