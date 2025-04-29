import os
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import textwrap
import protein_ligand_data

def get_docking_score(proteins, ligands, OUTPUT_DIR):
    n_protein = len(proteins)
    n_ligand  = len(ligands)
    score_size = (n_ligand, n_protein)

    docking_score = np.full(score_size, np.nan)
    for i, ligand_row in ligands.iterrows():
        ligand = ligand_row["SMILES"]
        for j, protein_row in proteins.iterrows():
            uniprot_id = protein_row["uniprot_id"]
            fname = protein_ligand_data.make_valid_fname(uniprot_id, ligand)
            fpath = os.path.join(OUTPUT_DIR, fname)

            if os.path.exists(fpath):
                with open(fpath) as f:
                    data = json.load(f)
                    docking_score[i, j] = data["result"]["docking_score"]

    return docking_score

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

def display_heatmap(docking_score, proteins, ligands, protein_set, mask):
    # Filter names if a mask is provided
    protein_names = proteins["name"][~mask].values
    ligand_names = ligands["name"].values

    # Word-wrap protein names at whitespace
    def wrap_name(name, width = 20):
        return textwrap.fill(name, width = width)

    wrapped_ligand_names = [wrap_name(name, 15) for name in ligand_names]
    wrapped_protein_names = [wrap_name(name) for name in protein_names]

    # Create DataFrame
    df = pd.DataFrame(docking_score, index = wrapped_ligand_names, columns = wrapped_protein_names)

    # Create clustermap with adjusted colorbar
    g = sns.clustermap(
        df,
        cmap = "viridis",
        figsize = (16, 6),
        linewidths = 0.5,
        xticklabels = True,
        yticklabels = True,
        cbar_kws = {"label": "Docking Score"},
        # cbar_pos = (0.9, 0.3, 0.03, 0.4),  # (x, y, width, height) in figure coords
    )

    # Adjust font sizes
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize = 14)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize = 16)

    # Add title
    g.figure.suptitle(
        f"Docking Scores for {protein_set.capitalize()} Proteins and Ligands",
        fontsize = 18,
        y = 1.05
    )

    plt.show()

def main(proteins, ligands, OUTPUT_DIR, protein_set):
    # collect data
    docking_score = get_docking_score(proteins, ligands, OUTPUT_DIR)
    # remove NaN columns
    mask = np.isnan(docking_score).any(axis = 0)  # get NaN columns as boolean mask
    docking_score = docking_score[:, ~mask]
    # display data
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

    args = parser.parse_args()

    protein_set = args.protein_set
    
    if protein_set == "thyroid":
        proteins = protein_ligand_data.get_thyroid_proteins()
        ligands  = protein_ligand_data.get_thyroid_ligands()
    elif protein_set == "autism":
        proteins = protein_ligand_data.get_autism_proteins()
        ligands  = protein_ligand_data.get_autism_ligands()
    else:
        raise ValueError("Invalid protein set selection.")
    
    main(proteins, ligands, args.output_dir, protein_set)