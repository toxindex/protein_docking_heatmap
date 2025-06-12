import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import textwrap
import protein_ligand_data


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

def main(proteins, ligands, OUTPUT_DIR, protein_set, use_exp = False, nan_value = None, clip_values = [None, None]):
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
    elif nan_value is not None:
        docking_score = np.nan_to_num(docking_score, nan = nan_value)
        mask = np.zeros(docking_score.shape)
    else:
        col_mean = np.nanmean(docking_score, axis = 0)
        inds = np.where(mask)
        docking_score[inds] = np.take(col_mean, inds[1])
        # any remaining NaN values should be replaced with 0
        docking_score = np.nan_to_num(docking_score, nan = 0.0)

    docking_score = np.clip(docking_score, clip_values[0], clip_values[1])

    print(f"Any nan values: {np.isnan(docking_score).any()}")
    print(f"Min docking score = {np.min(docking_score)}")
    print(f"Max docking score = {np.max(docking_score)}")
    
    # display data
    display_heatmap(docking_score, proteins, ligands, protein_set, mask)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Generate heatmap from docking results."
    )
    parser.add_argument(
        "--protein_set",
        choices = ["thyroid", "autism", "cancer"],
        required = True,
        help = "Choose which set of proteins to use: 'thyroid', 'autism', or 'cancer'."
    )
    parser.add_argument(
        "--output_dir",
        default = "./docking_results",
        help = "Directory containing docking results"
    )
    parser.add_argument(
        "--use_exp",
        default = False,
        help = "Use e^x scale for docking scores."
    )
    parser.add_argument(
        "--nan_value",
        default = None,
        type = float,
        help = "Numerical value to which nan values should be masked"
    )
    parser.add_argument(
        "--clip_values",
        metavar = ('MIN', 'MAX'),
        nargs = 2,
        type = float,
        default = [None, None],
        help = "Restrict docking scores to be between two values, min and max. 'None' if not restricting that value."
    )

    args = parser.parse_args()

    protein_set = args.protein_set
    
    proteins, ligands = protein_ligand_data.get_proteins_ligands(protein_set)
    
    main(proteins, ligands, args.output_dir, protein_set, args.use_exp, args.nan_value, args.clip_values)