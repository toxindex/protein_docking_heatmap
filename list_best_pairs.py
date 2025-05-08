import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import pandas as pd
import protein_ligand_data
import textwrap

def plot_ligand_table(best_ligands, docking_scores, protein_names):
    n_rows, n_cols = best_ligands.shape
    fig_width = 1.5 * n_cols
    fig_height = 2 * n_rows
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set color limits: floor at 0 if all scores are positive
    vmin = 0 if np.all(docking_scores >= 0) else None
    im = ax.imshow(docking_scores, cmap="viridis", aspect="auto", vmin=vmin)

    # Estimate character width for wrapping: adjust divisor to tune wrapping tightness
    est_chars_per_cell = int((fig_width * 10) / n_cols)

    wrapper = textwrap.TextWrapper(
        width=est_chars_per_cell,
        break_long_words=True,
        break_on_hyphens=True
    )

    norm = Normalize(vmin=np.nanmin(docking_scores), vmax=np.nanmax(docking_scores))
    cmap = cm.get_cmap("viridis")

    # Add wrapped ligand names with adaptive text color
    for i in range(n_rows):
        for j in range(n_cols):
            name = best_ligands[i, j]
            wrapped_name = wrapper.fill(name)
            score = docking_scores[i, j]

            rgba = cmap(norm(score))
            r, g, b = rgba[:3]
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            text_color = "black" if luminance > 0.5 else "white"

            ax.text(j, i, wrapped_name, ha="center", va="center", color=text_color, fontsize=8)
    
    # Set axis ticks and labels
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(protein_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([f"{i+1}" for i in range(n_rows)], fontsize=9)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Docking Score")

    # Formatting
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)
    ax.set_xlabel("Protein", fontsize=10)
    ax.set_ylabel("Ranked Ligands", fontsize=10)

    plt.tight_layout()
    plt.show()



# def main(proteins, ligands, OUTPUT_DIR, top_n_pairs):
#     # collect data
#     docking_score = protein_ligand_data.get_docking_score(proteins, ligands, OUTPUT_DIR)
#     # replace NaN with -infinity
#     docking_score = np.nan_to_num(docking_score, nan = -np.inf)
#     # get the n greatest scores per protein 
#     best_score_inds = np.flipud(np.argsort(docking_score, axis = 0))[:top_n_pairs]
#     cols = np.arange(docking_score.shape[1])
#     print(docking_score[best_score_inds, cols])

#     # give the names of the top scoring ligands per protein
#     best_ligands = [
#         list(ligands.iloc[best_score_inds[:, i]]["name"].values)
#         for i in range(docking_score.shape[1])
#     ]
#     best_ligands = np.array(best_ligands).T

#     print(best_ligands)

# def main(proteins, ligands, OUTPUT_DIR, top_n_pairs):
#     docking_score = protein_ligand_data.get_docking_score(proteins, ligands, OUTPUT_DIR)
#     docking_score = np.nan_to_num(docking_score, nan = -np.inf)

#     best_score_inds = np.flipud(np.argsort(docking_score, axis = 0))[:top_n_pairs]

#     best_scores = np.array([
#         docking_score[best_score_inds[:, i], i]
#         for i in range(docking_score.shape[1])
#     ]).T

#     best_ligands = np.array([
#         ligands.iloc[best_score_inds[:, i]]["name"].values
#         for i in range(docking_score.shape[1])
#     ]).T

#     protein_names = proteins["name"].values  # Assuming proteins is a DataFrame

#     # if all of the scores are -inf for a column, remove that column

#     plot_ligand_table(best_ligands, best_scores, protein_names)

def main(proteins, ligands, OUTPUT_DIR, top_n_pairs):
    docking_score = protein_ligand_data.get_docking_score(proteins, ligands, OUTPUT_DIR)
    docking_score = np.nan_to_num(docking_score, nan=-np.inf)

    # Get the n greatest scores per protein
    best_score_inds = np.flipud(np.argsort(docking_score, axis=0))[:top_n_pairs]

    best_scores = np.array([
        docking_score[best_score_inds[:, i], i]
        for i in range(docking_score.shape[1])
    ]).T

    best_ligands = np.array([
        ligands.iloc[best_score_inds[:, i]]["name"].values
        for i in range(docking_score.shape[1])
    ]).T

    protein_names = proteins["name"].values  # Assuming proteins is a DataFrame

    # Identify columns where all scores are -inf
    valid_columns = ~np.all(best_scores == -np.inf, axis=0)

    # Filter out invalid columns
    best_score_inds = best_score_inds[:, valid_columns]
    best_scores = best_scores[:, valid_columns]
    best_ligands = best_ligands[:, valid_columns]
    protein_names = protein_names[valid_columns]

    # Plot the ligand table
    plot_ligand_table(best_ligands, best_scores, protein_names)

    

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
        "--top_n_pairs",
        type = int,
        default = 4,
        help = "Display the top n ligands for each protein."
    )
    parser.add_argument(
        "--output_dir",
        default = "./docking_results",
        help = "Directory containing docking results."
    )

    args = parser.parse_args()

    protein_set = args.protein_set
    proteins, ligands = protein_ligand_data.get_proteins_ligands(protein_set)
    
    main(proteins, ligands, args.output_dir, args.top_n_pairs)