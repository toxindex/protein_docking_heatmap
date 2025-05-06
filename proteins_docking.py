import os
import argparse
import httpx
# import itertools
import pandas as pd  # Required for handling DataFrame objects like proteins and ligands
import protein_ligand_data

import time
    
def submit_jobs(callback_url, proteins, ligands):
    OUTPUT_DIR = "./docking_results"
    for _, protein_row in proteins.iterrows():
        uniprot_id = protein_row["uniprot_id"]
        for _, ligand_row in ligands.iterrows():
            ligand = ligand_row["SMILES"]
            payload = {
                "uniprot_id": uniprot_id,
                "ligand": ligand,
                "callback_url": callback_url
            }
            
            # skip if data already collected for this pair
            fname = protein_ligand_data.make_valid_fname(
                payload["uniprot_id"], payload["ligand"]    
            )
            fpath = os.path.join(OUTPUT_DIR, fname)

            if os.path.exists(fpath):
                print(f"Skipping {uniprot_id} + {ligand[:10]}... (already exists)")
                continue            
            
            try:
                response = httpx.post(
                    "https://diffdock.toxindex.com/start_docking_uniprot",
                    json = payload
                )
                response.raise_for_status()
                print(f"Submitted {uniprot_id} + {ligand}...")
                
                time.sleep(0.5)  # Rate limiting
            except Exception as e:
                print(f"Failed to submit {uniprot_id} + {ligand}...: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Submit batch docking requests with webhook support."
    )
    parser.add_argument(
        "--webhook",
        required = True,
        help = "The callback URL where docking results will be sent."
    )
    parser.add_argument(
        "--protein_set",
        choices = ["thyroid", "autism"],
        required = True,
        help = "Choose which set of proteins to use: 'thyroid' or 'autism'."
    )

    args = parser.parse_args()

    callback_url = args.webhook
    protein_set = args.protein_set
    
    if protein_set == "thyroid":
        proteins = protein_ligand_data.get_thyroid_proteins()
        ligands  = protein_ligand_data.get_thyroid_ligands()
    elif protein_set == "autism":
        proteins = protein_ligand_data.get_autism_proteins()
        ligands  = protein_ligand_data.get_autism_ligands()
    else:
        raise ValueError("Invalid protein set selection.")
        
    submit_jobs(callback_url, proteins, ligands)
