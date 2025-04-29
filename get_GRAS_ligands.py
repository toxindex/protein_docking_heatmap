import pandas as pd
import requests
import time

def get_cid(substance_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{substance_name}/cids/JSON"
    try:
        response = requests.get(url, timeout = 5)
        response.raise_for_status()
        cids = response.json().get('IdentifierList', {}).get('CID', [])
        return cids[0] if cids else None
    except Exception:
        return None

def get_smiles(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url, timeout = 5)
        response.raise_for_status()
        properties = response.json().get('PropertyTable', {}).get('Properties', [])
        if properties and 'CanonicalSMILES' in properties[0]:
            return properties[0]['CanonicalSMILES']
    except Exception:
        return None
    return None

def main(input_csv, output_pkl, output_csv = None):
    # Read the CSV
    df_input = pd.read_csv(
        input_csv, 
        skiprows = 4,          # Skip the garbage
        quotechar = '"',       # Handle quoted fields
        skipinitialspace = True  # Remove extra leading spaces after commas
    )
    
    # Assume the substance names are in the first column
    # substance_names = df_input.iloc[:, 0].dropna().tolist()
    substance_names = df_input['GRAS Substance'].dropna().tolist()
    
    data = []
    for name in substance_names:
        cid = get_cid(name)
        if cid is None:
            continue
        smiles = get_smiles(cid)
        if smiles is None:
            continue
        link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        data.append({'name': name, 'SMILES': smiles, 'PubChem': link})
        time.sleep(0.2)  # polite delay
    
    df_output = pd.DataFrame(data)
    
    # Save as Pickle
    df_output.to_pickle(output_pkl)
    print(f"Saved DataFrame to {output_pkl}")
    
    # Optionally also save as CSV
    if output_csv:
        df_output.to_csv(output_csv, index = False)
        print(f"Saved DataFrame to {output_csv}")
    
    return df_output

if __name__ == "__main__":
    input_csv = "./SCOGS.csv"
    output_pkl = "./scogs_with_smiles.pkl"
    output_csv = "./scogs_with_smiles.csv"
    
    df = main(input_csv, output_pkl, output_csv)
    print(df)
