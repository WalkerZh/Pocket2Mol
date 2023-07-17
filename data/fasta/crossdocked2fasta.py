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


def get_fasta_from_pdb_biopython_no_file(pdb_id, chain):
    pdb_file = f"./pdb_files/pdb{pdb_id}.ent"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        with open(pdb_file, "wb") as f:
            f.write(response.content)
    else:
        pdb_list = PDBList()
        pdb_file = pdb_list.retrieve_pdb_file(pdb_id, file_format="pdb", pdir="pdb_files")
    
    if not os.path.exists(pdb_file):
        return "", -100

    with open(pdb_file, "r") as file:
        pdb_content = file.read()

    pdb_handle = StringIO(pdb_content)
    fasta_handle = StringIO()
    SeqIO.convert(pdb_handle, "pdb-atom", fasta_handle, "fasta")

    fasta_sequence = fasta_handle.getvalue()

    if fasta_sequence:
        start = get_start_idx(pdb_file, chain)
        fasta = get_chain_fasta(fasta_sequence, chain)
        if start:
            return fasta, start
        else:
            return "", -100
    else:
        return "", -100

def get_chain_fasta(fasta_sequence, chain):
    fasta_list = fasta_sequence.split(">")[1:]
    temp_dict = {}
    for ff in fasta_list:
        ff = ff.strip()
        f_chain = str(ff.split("\n")[0].split(":")[1].strip())
        f_seq = "".join(ff.split("\n")[1:])
        temp_dict[f_chain] = f_seq
    return temp_dict[str(chain)]

def get_start_idx(pdb_file, chain):
    with open(pdb_file, "r", encoding="utf8") as f:
        for line in f:
            line = line.strip()
            if line[0:6].strip() == 'ATOM' and line[21:22].strip() == chain:
                return int(line[22:26])
    return None
  
if __name__ == "__main__":
    crossdocked_dir = "../crossdocked_pocket10/"
    subdir_list = os.listdir(crossdocked_dir)

    csv_name = "crossdocked2fasta_pocket10.csv"
    exist_pdb = []

    with open(csv_name, 'w', newline='') as csvfile:
        fieldnames = ['pdb_id', 'chain', 'fasta_sequence', 'start'] # 'protein_filename', 'ligand_filename' ,
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for subdir in tqdm(subdir_list):
            temp_path = os.path.join(crossdocked_dir, subdir)
            if os.path.isfile(temp_path): continue

            file_list = os.listdir(temp_path)
            for file_name in file_list:
                if file_name.endswith(".sdf"): continue
                # ligand_filename = "_".join(file_name.split("_")[:-1]) + ".sdf"
                # if ligand_filename not in file_list:
                #     print(ligand_filename)
                #     continue

                protein_filename = file_name
                pdb_id = str(protein_filename.split("_")[0].strip())
                chain = str(protein_filename.split("_")[1].strip())
                if (pdb_id, chain) in exist_pdb:
                    continue
                else:
                    exist_pdb.append((pdb_id, chain))
                    fasta, start = get_fasta_from_pdb_biopython_no_file(pdb_id, chain)
                    if fasta == "" and start == -100:
                        logging.error(f"PDB id {pdb_id} in {file_name} from {temp_path} is not exist!")
                        continue
                    writer.writerow({'pdb_id': pdb_id, 'chain': chain, 'fasta_sequence': str(fasta), 'start': start})  
