import os
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import protein_ligand_data

def main(proteins, ligands, OUTPUT_DIR, protein_set, use_exp = False):
    # collect data
    docking_score = protein_ligand_data.get_docking_score(proteins, ligands, OUTPUT_DIR)

    '''
    # remove rows and columns with more than threshold NaN values
    nan_threshold = 0.5
    print(f"Initial shape: {docking_score.shape}")
    # print mean number of NaN values in each column
    print(f"Mean number of NaN values in each column: {np.isnan(docking_score).mean(axis = 0)}")
    # docking_score = np.nan_to_num(docking_score, nan = -3.0)
    col_mask = np.isnan(docking_score).mean(axis = 0) > nan_threshold  # get NaN columns as boolean mask
    docking_score = docking_score[:, ~col_mask]
    row_mask = np.isnan(docking_score).mean(axis = 1) > nan_threshold  # get NaN rows as boolean mask
    docking_score = docking_score[~row_mask, :]
    print(f"Shape after removing NaN rows and columns: {docking_score.shape}")
    '''
    mask = np.isnan(docking_score)
    
    if use_exp:
        # convert to e^x scale
        docking_score = np.exp(docking_score)
        # replace any nan values with 0
        docking_score = np.nan_to_num(docking_score, nan = 0.0)
    else:
        # # replace any NaN values with -3 (TODO: more robust way to handle this)
        # docking_score = np.nan_to_num(docking_score, nan = -3.0)

        # docking_score = np.where(mask, np.nanmean(docking_score, axis = 1, keepdims = True), docking_score)
        # docking_score = np.nan_to_num(docking_score, nan = 0.0)
        col_mean = np.nanmean(docking_score, axis = 0)
        inds = np.where(mask)
        docking_score[inds] = np.take(col_mean, inds[1])
        # any remaining NaN values should be replaced with 0
        docking_score = np.nan_to_num(docking_score, nan = 0.0)

    # # Remove rows or columns with any non-finite values
    # finite_rows = np.all(mask, axis=1)
    # finite_cols = np.all(mask, axis=0)
    # docking_score = docking_score[finite_rows][:, finite_cols]
    print(f"Any nan values: {np.isnan(docking_score).any()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Generate heatmap from docking results."
    )
    parser.add_argument(
        "--protein_set",
        choices = ["thyroid", "autism"],
        required = True,
        help = "Choose which set of proteins to use: 'thyroid' or 'autism'."
    )
    parser.add_argument(
        "--output_dir",
        default = "./docking_results",
        help = "Directory containing docking results."
    )

    args = parser.parse_args()

    protein_set = args.protein_set
    proteins, ligands = protein_ligand_data.get_proteins_ligands(protein_set)
    
    main(proteins, ligands, args.output_dir, protein_set, args.use_exp)