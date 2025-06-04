import pandas as pd
from pubchem_query import get_pubchem_data


def main(input_txt, output_pkl, output_csv = None):
    with open(input_txt, 'r') as f:
        substance_names = [line.strip() for line in f if line.strip()]

    # Remove duplicates
    substance_names = list(set(substance_names))
    
    df_output = get_pubchem_data(substance_names)
    
    # Save as Pickle
    df_output.to_pickle(output_pkl)
    print(f"Saved DataFrame to {output_pkl}")
    
    # Optionally also save as CSV
    if output_csv:
        df_output.to_csv(output_csv, index = False)
        print(f"Saved DataFrame to {output_csv}")
    
    return df_output

if __name__ == "__main__":
    input_txt = "./hydrocarbon_chemicals.txt"
    output_pkl = "./hydrocarbons_with_smiles.pkl"
    output_csv = "./hydrocarbons_with_smiles.csv"
    
    df = main(input_txt, output_pkl, output_csv)
    print(df)
