import csv
import os
import pickle
import lmdb
import torch
from torch.utils.data import Dataset
from tqdm.auto import tqdm

from ..protein_ligand import PDBProtein, parse_sdf_file
from ..data import ProteinLigandData, torchify_dict


class PocketLigandPairDataset(Dataset):

    def __init__(self, raw_path, transform=None):
        super().__init__()
        self.raw_path = raw_path.rstrip('/')
        self.index_path = os.path.join(self.raw_path, 'index.pkl')
        self.processed_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_processed.lmdb')
        self.name2id_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_name2id.pt')
        self.transform = transform
        self.db = None

        self.keys = None

        if not (os.path.exists(self.processed_path) and os.path.exists(self.name2id_path)):
            self._process()
            self._precompute_name2id()

        self.name2id = torch.load(self.name2id_path)

    def _connect_db(self):
        """
            Establish read-only database connection
        """
        assert self.db is None, 'A connection has already been opened.'
        self.db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=False,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        with self.db.begin() as txn:
            self.keys = list(txn.cursor().iternext(values=False))

    def _close_db(self):
        self.db.close()
        self.db = None
        self.keys = None
        
    def _process(self):
        db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=True,
            subdir=False,
            readonly=False, # Writable
        )
        with open(self.index_path, 'rb') as f:
            index = pickle.load(f)

        num_skipped = 0

        # open csv file to check the status
        num_status_error, num_key_error = 0, 0
        data_status = {}
        status_path = "./data/fasta/crossdocked_data_status.csv"
        with open(status_path, newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            for row in reader:
                protein_file, ligand_file, subdir, status = row[0], row[1], row[2], row[8]
                pocket_filename = f"{subdir}/{protein_file}"
                ligand_filename = f"{subdir}/{ligand_file}"
                data_status[(pocket_filename, ligand_filename)] = status
        
        with db.begin(write=True, buffers=True) as txn:
            for i, (pocket_fn, ligand_fn, _, rmsd_str) in enumerate(tqdm(index)):
                if pocket_fn is None: continue
                if (pocket_fn, ligand_fn) not in data_status.keys():
                    num_key_error += 1
                    print('Key error, skipping (%d) %s' % (num_key_error, ligand_fn, ))
                    continue
                if data_status[(pocket_fn, ligand_fn)] != "Pass":
                    num_status_error += 1
                    print('Status error, skipping (%d) %s' % (num_status_error, ligand_fn, ))
                    continue
                try:
                    pocket_dict = PDBProtein(os.path.join(self.raw_path, pocket_fn)).to_dict_atom()
                    ligand_dict = parse_sdf_file(os.path.join(self.raw_path, ligand_fn))
                    data = ProteinLigandData.from_protein_ligand_dicts(
                        protein_dict=torchify_dict(pocket_dict),
                        ligand_dict=torchify_dict(ligand_dict),
                    )
                    data.protein_filename = pocket_fn
                    data.protein_pdb_id = pocket_fn.split("/")[1].split("_")[0]
                    data.ligand_filename = ligand_fn
                    # with open("./input_data.txt", "w") as fw:
                    #     for k, v in data:
                    #         if isinstance(v, torch.Tensor):
                    #             print(f"{k}({v.size()}): {v}", file=fw)
                    #         elif isinstance(v, list):
                    #             print(f"{k}({len(v)}): {v}", file=fw)
                    #         else:
                    #             print(f"{k}: {v}", file=fw)
                    txn.put(
                        key = str(i).encode(),
                        value = pickle.dumps(data)
                    )
                except:
                    num_skipped += 1
                    print('Skipping (%d) %s' % (num_skipped, ligand_fn, ))
                    continue
        print(f"| Num_skipped {num_skipped} | Num_key_error {num_key_error} | Num_status_error {num_status_error} |")
        db.close()

    def _precompute_name2id(self):
        name2id = {}
        for i in tqdm(range(self.__len__()), 'Indexing'):
            try:
                data = self.__getitem__(i)
            except AssertionError as e:
                print(i, e)
                continue
            name = (data.protein_filename, data.ligand_filename)
            name2id[name] = i
        torch.save(name2id, self.name2id_path)
    
    def __len__(self):
        if self.db is None:
            self._connect_db()
        return len(self.keys)

    def __getitem__(self, idx):
        if self.db is None:
            self._connect_db()
        key = self.keys[idx]
        data = pickle.loads(self.db.begin().get(key))
        data.id = idx
        assert data.protein_pos.size(0) > 0
        if self.transform is not None:
            data = self.transform(data)
        # fw1 = open("./transformed_data.txt", "w")
        # fw2 = open("./transformed_data_key.txt", "w")
        # for k, v in data:
        #     if isinstance(v, torch.Tensor):
        #         print(f"{k}({v.size()}): {v}", file=fw1)
        #     elif isinstance(v, list):
        #         print(f"{k}({len(v)}): {v}", file=fw1)
        #     else:
        #         print(f"{k}: {v}", file=fw1)
        #     print(f"{k}", file=fw2)
        # fw1.close()
        # fw2.close()
        return data # type(data) = utils.data.ProteinLigandData
        

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    PocketLigandPairDataset(args.path)
