import os
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import textwrap
import protein_ligand_data

# def display_heatmap(docking_score, proteins, ligands, protein_set, mask):
#     n_ligand, n_protein = docking_score.shape
#     # generate heatmap
#     fig, ax = plt.subplots(figsize = (10, 8))
#     cax = ax.matshow(docking_score, cmap = 'viridis')
#     # add a label to the colorbar
#     cbar = fig.colorbar(cax, ax = ax, fraction = 0.02, pad = 0.04)
#     cbar.set_label('Docking Score', rotation = 270, labelpad = 15)

#     # Set axis labels
#     ax.set_title(f"Docking Scores for {protein_set.capitalize()} Proteins and Ligands")
#     ax.set_xticks(range(n_protein))
#     ax.set_yticks(range(n_ligand))
#     ax.set_xticklabels(
#         proteins["name"][~mask].values,
#         # rotation = 60,
#         rotation = 90,
#         ha = 'right'
#     )
#     ax.set_yticklabels(ligands["name"])

#     fig.tight_layout()

#     plt.show()  # Display the heatmap

# def display_heatmap(docking_score, proteins, ligands, protein_set, mask):
#     # Filter names if a mask is provided
#     protein_names = proteins["name"][~mask].values
#     ligand_names = ligands["name"].values

#     # Break long protein names across lines
#     def wrap_name(name, width = 20):
#         return '\n'.join([name[i:i+width] for i in range(0, len(name), width)])

#     wrapped_protein_names = [wrap_name(name) for name in protein_names]

#     # Create a DataFrame for seaborn
#     df = pd.DataFrame(docking_score, index = ligand_names, columns = wrapped_protein_names)

#     # Generate clustermap
#     g = sns.clustermap(
#         df,
#         cmap = "viridis",
#         figsize = (14, 6),
#         linewidths = 0.5,
#         xticklabels = True,
#         yticklabels = True,
#         cbar_kws = {'label': 'Docking Score'}
#     )

#     # Adjust font sizes
#     g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize = 10)
#     g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize = 12)

#     # Title
#     plt.title(f"Docking Scores for {protein_set.capitalize()} Proteins and Ligands", fontsize = 14, pad = 20)
#     plt.show()

# def display_heatmap(docking_score, proteins, ligands, protein_set, row_mask, col_mask):
#     # Filter names if a mask is provided
#     protein_names = proteins["name"][~col_mask].values
#     ligand_names = ligands["name"][~row_mask].values

#     # Word-wrap protein names at whitespace
#     def wrap_name(name, width = 20):
#         return textwrap.fill(name, width = width)

#     # wrapped_ligand_names = [wrap_name(name, 15) for name in ligand_names]
#     wrapped_ligand_names = ligand_names
#     wrapped_protein_names = [wrap_name(name) for name in protein_names]

#     # Create DataFrame
#     df = pd.DataFrame(
#         docking_score,
#         index = wrapped_ligand_names,
#         columns = wrapped_protein_names
#     )

#     # Create clustermap with adjusted colorbar
#     g = sns.clustermap(
#         df,
#         cmap = "viridis",
#         # figsize = (16, 6),
#         # linewidths = 0.5,
#         xticklabels = True,
#         yticklabels = False,
#         cbar_kws = {"label": "Docking Score"},
#         # cbar_pos = (0.9, 0.3, 0.03, 0.4),  # (x, y, width, height) in figure coords
#     )

#     # Adjust font sizes
#     g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize = 14)
#     g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize = 10)

#     # Add title
#     g.figure.suptitle(
#         f"Docking Scores for {protein_set.capitalize()} Proteins and Ligands",
#         fontsize = 18,
#         y = 1.05
#     )

#     plt.show()

def display_heatmap(docking_score, proteins, ligands, protein_set, mask):
    # Filter names if a mask is provided
    protein_names = proteins["name"].values
    ligand_names = ligands["name"].values

    # Word-wrap protein names at whitespace
    def wrap_name(name, width = 20):
        return textwrap.fill(name, width = width)

    # wrapped_ligand_names = [wrap_name(name, 15) for name in ligand_names]
    wrapped_ligand_names = ligand_names
    wrapped_protein_names = [wrap_name(name) for name in protein_names]

    # Create DataFrame
    df = pd.DataFrame(
        docking_score,
        index = wrapped_ligand_names,
        columns = wrapped_protein_names
    )

    # Create clustermap with adjusted colorbar
    g = sns.clustermap(
        df,
        cmap = "viridis",
        # figsize = (16, 6),
        # linewidths = 0.5,
        xticklabels = True,
        yticklabels = False,
        cbar_kws = {"label": "Docking Score"},
        mask = mask,
        # cbar_pos = (0.9, 0.3, 0.03, 0.4),  # (x, y, width, height) in figure coords
    )

    # Adjust font sizes
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize = 14)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize = 10)

    # Add title
    g.figure.suptitle(
        f"Docking Scores for {protein_set.capitalize()} Proteins and Ligands",
        fontsize = 18,
        y = 1.05
    )

    # g.ax_heatmap.imshow(mask, cmap='gray_r', alpha=0.5, aspect='auto', zorder=2)

    plt.show()

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
    
    # display data
    # display_heatmap(docking_score, proteins, ligands, protein_set, row_mask, col_mask)
    display_heatmap(docking_score, proteins, ligands, protein_set, mask)

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
    parser.add_argument(
        "--use_exp",
        default = False,
        help = "Use e^x scale for docking scores."
    )

    args = parser.parse_args()

    protein_set = args.protein_set
    
    proteins, ligands = protein_ligand_data.get_proteins_ligands(protein_set)
    
    main(proteins, ligands, args.output_dir, protein_set, args.use_exp)