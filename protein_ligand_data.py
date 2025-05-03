import pandas as pd
import hashlib

def make_valid_fname(uniprot_id, ligand):
    valid_ligand = ligand.replace('/', '_')
    fname = f"{uniprot_id}_{valid_ligand}.json"
    if len(fname) > 255:
        # hash the ligand to make it shorter
        hash_object = hashlib.sha256(ligand.encode())
        hash_ligand = hash_object.hexdigest()
        # create a new filename with the hash
        fname = f"{uniprot_id}_{hash_ligand}.json"
    
    return fname

def get_thyroid_proteins():
    data = {
            "name": [
                "Thyroid hormone receptor alpha",
                "Thyroid hormone receptor beta",
                "Retinoid X receptor α",
                "Peroxisome proliferator‑activated receptor α",
                "Thyroperoxidase",
                "Sodium‑iodide symporter",
                "Thyroglobulin",
                "Thyroid stimulating hormone receptor",
                "Type I deiodinase",
                "Type II deiodinase",
                "Type III deiodinase",
                "UDP‑glucuronosyltransferase 1A1",
                "Sulfotransferase 1A1",
                "Thyroxine‑binding globulin",
                "Transthyretin",
                "Albumin",
                "Monocarboxylate transporter 8",
                "Monocarboxylate transporter 10",
                "Nuclear receptor corepressor 1",
                "Steroid receptor coactivator 1"
            ],
            "gene_symbol": [
                "THRA",
                "THRB",
                "RXRA",
                "PPARA",
                "TPO",
                "SLC5A5",
                "TG",
                "TSHR",
                "DIO1",
                "DIO2",
                "DIO3",
                "UGT1A1",
                "SULT1A1",
                "SERPINA7",
                "TTR",
                "ALB",
                "SLC16A2",
                "SLC16A10",
                "NCOR1",
                "NCOA1"
            ],
            "uniprot_id": [
                "P10827",
                "P10828",
                "P19793",
                "Q07869",
                "P07202",
                "Q92911",
                "P01266",
                "P16473",
                "P49895",
                "Q92813",
                "P55073",
                "P22309",
                "P50225",
                "P05543",
                "P02766",
                "P02768",
                "P36021",
                "Q8TF71",
                "O75376",
                "Q15788"
            ]
    }

    proteins = pd.DataFrame(data)

    # remove entries not in database
    uniprot_to_remove = ["P07202", "Q92911", "P01266", "P49895", "Q92813", "P55073", "P22309", "P36021", "Q8TF71"]
    proteins = proteins[~proteins["uniprot_id"].isin(uniprot_to_remove)]
    proteins.reset_index(drop=True, inplace=True)  # Reset indices

    return proteins

    
def get_thyroid_ligands():
    ligands = pd.read_pickle("./scogs_with_smiles.pkl")
    return ligands
    
'''
    Below is a working short-list of 20 proteins that toxicologists most often
    dock against when they want to know whether a chemical could plausibly push
    a susceptible brain toward an autism-like phenotype.
'''
def get_autism_proteins():
    data = {
        "name": [
            "Aryl-hydrocarbon receptor",
            "Paraoxonase-1",
            "Glutathione S-transferase Mu 1",
            "Cytochrome P450 1A2",
            "SH3- and ankyrin-repeat protein 3",
            "Contactin-associated protein-like 2",
            "Chromodomain-helicase-DNA-binding protein 8",
            "Voltage-gated Na⁺ channel α-subunit 2",
            "Neurexin-1",
            "Serotonin transporter",
            "Oxytocin receptor",
            "Reelin",
            "GABA-A receptor β3 subunit",
            "Dopamine D2 receptor",
            "Phosphatase and tensin homolog",
            "Tuberin",
            "mTOR serine/threonine kinase",
            "Methyl-CpG-binding protein 2",
            "Fragile-X messenger ribonucleoprotein",
            "Brain-derived neurotrophic factor",
        ],
        "gene_symbol": [
            "AHR",
            "PON1",
            "GSTM1",
            "CYP1A2",
            "SHANK3",
            "CNTNAP2",
            "CHD8",
            "SCN2A",
            "NRXN1",
            "SLC6A4",
            "OXTR",
            "RELN",
            "GABRB3",
            "DRD2",
            "PTEN",
            "TSC2",
            "MTOR",
            "MECP2",
            "FMR1",
            "BDNF",
        ],
        "uniprot_id": [
            "P35869",
            "Q15165",
            "P09488",
            "P05177",
            "Q9BYB0",
            "Q9UHC6",
            "Q9HCK8",
            "Q99250",
            "Q9ULB1",
            "P31645",
            "P30559",
            "P78509",
            "P28472",
            "P14416",
            "P60484",
            "P49815",
            "P42345",
            "P51608",
            "Q06787",
            "P23560",
        ],
    }
    proteins = pd.DataFrame(data)

    # remove entries not in database
    uniprot_to_remove = ["Q15165", "Q9BYB0", "Q9UHC6", "Q9ULB1", "P30559", "P49815"]
    proteins = proteins[~proteins["uniprot_id"].isin(uniprot_to_remove)]
    proteins.reset_index(drop=True, inplace=True)  # Reset indices

    return proteins

'''
    GRAS substances
    see:
        https://www.hfpappexternal.fda.gov/scripts/fdcc/index.cfm?set=SCOGS
''' 
def get_autism_ligands():
    # data = {
    #     'name': [
    #         "acetic acid",
    #         # "acetyl distarch adipate",
    #         "aconitic acid",
    #         "adipic acid",
    #         # "acetylated distarch phosphate",
    #         "allyl isothiocyanate",
    #         "aluminum ammonium sulfate",
    #     ],
    #     'SMILES': [
    #         "CC(=O)O",
    #         # "CC1C(C(C(C(O1)CO)OC2C(C(C(C(O2)COC3C(C(C(C(O3)CO)OC)O)O)OC4C(C(C(C(O4)CO)OC)O)OC(=O)CCCCC(=O)OC5C(OC(C(C5O)O)OC)CO)O)O)O)OC(=O)C",
    #         "C(/C(=C/C(=O)O)/C(=O)O)C(=O)O",
    #         "C(CCC(=O)O)CC(=O)O",
    #         # "CC1C(C(C(C(O1)CO)OC2C(C(C(C(O2)COC3C(C(C(C(O3)CO)OC)O)O)OC4C(C(C(C(O4)CO)OC)O)OC(=O)C)O)O)O)OP(=O)(O)OC5C(OC(C(C5O)O)OC)CO",
    #         "C=CCN=C=S",
    #         "[NH4+].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[Al+3]",
    #     ],
    #     'PubChem': [
    #         "https://pubchem.ncbi.nlm.nih.gov/compound/Acetic-Acid",
    #         # "https://pubchem.ncbi.nlm.nih.gov/compound/24847850",
    #         "https://pubchem.ncbi.nlm.nih.gov/compound/Aconitic-Acid",
    #         "https://pubchem.ncbi.nlm.nih.gov/compound/Adipic-Acid",
    #         # "https://pubchem.ncbi.nlm.nih.gov/compound/24832109",
    #         "https://pubchem.ncbi.nlm.nih.gov/compound/Allyl-Isothiocyanate",
    #         "https://pubchem.ncbi.nlm.nih.gov/compound/Aluminum-ammonium-sulfate",
            
    #     ]
    # }
    # ligands = pd.DataFrame(data)
    ligands = pd.read_pickle("./scogs_with_smiles.pkl")
    return ligands