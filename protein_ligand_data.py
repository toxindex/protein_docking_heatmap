import pandas as pd

def make_valid_fname(uniprot_id, ligand):
    valid_ligand = ligand.replace('/', '_')
    fname = f"{uniprot_id}_{valid_ligand}.json"
    
    return fname

def get_thyroid_proteins():
    return [""]
    
def get_thyroid_ligands():
    return [""]
    
'''
    Below is a working short-list of 20 proteins that toxicologists most often
    dock against when they want to know whether a chemical could plausibly push
    a susceptible brain toward an autism-like phenotype.
'''
def get_autism_proteins():
    data = {
        "name": [
            "aryl-hydrocarbon receptor",
            "paraoxonase-1",
            "glutathione s-transferase mu 1",
            "cytochrome p450 1a2",
            "sh3- and ankyrin-repeat protein 3",
            "contactin-associated protein-like 2",
            "chromodomain-helicase-dna-binding protein 8",
            "voltage-gated na⁺ channel α-subunit 2",
            "neurexin-1",
            "serotonin transporter",
            "oxytocin receptor",
            "reelin",
            "gaba-a receptor β3 subunit",
            "dopamine d2 receptor",
            "phosphatase and tensin homolog",
            "tuberin",
            "mtor serine/threonine kinase",
            "methyl-cpg-binding protein 2",
            "fragile-x messenger ribonucleoprotein",
            "brain-derived neurotrophic factor",
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
        ]
    }
    proteins = pd.DataFrame(data)
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