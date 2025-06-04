import os
import json
import numpy as np
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
                "Peroxisome proliferator-activated receptor α",
                "Thyroperoxidase",
                "Sodium-iodide symporter",
                "Thyroglobulin",
                "Thyroid stimulating hormone receptor",
                "Type I deiodinase",
                "Type II deiodinase",
                "Type III deiodinase",
                "UDP-glucuronosyltransferase 1A1",
                "Sulfotransferase 1A1",
                "Thyroxine-binding globulin",
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
    ligands = pd.read_pickle("./scogs_with_smiles2.pkl")
    return ligands

def get_cancer_proteins():
    data = {
        "name": [
            "ABL proto-oncogene 1",
            "PKB alpha",
            "Anaplastic lymphoma kinase",
            "B-cell lymphoma 2",
            "B-Raf kinase",
            "Cyclin D1",
            "Cyclin-dependent kinase 4",
            "Epidermal growth factor receptor",
            "HER2/Neu",
            "Fibroblast growth factor receptor 3",
            "H-Ras GTPase",
            "K-Ras GTPase",
            "N-Ras GTPase",
            "Isocitrate dehydrogenase 1",
            "Janus kinase 2",
            "c-Kit receptor",
            "HGF receptor",
            "c-Myc transcription factor",
            "N-Myc",
            "L-Myc",
            "NOTCH1",
            "NRF2",
            "PDGFRA",
            "PI3K catalytic subunit α",
            "RET receptor tyrosine kinase",
            "Src tyrosine kinase",
            "Telomerase reverse transcriptase",
            "p53 E3 ubiquitin ligase",
            "β-Catenin",
            "MYD88",
            "Calreticulin",
            "TRKA kinase",
            "Enhancer of zeste homolog 2",
            "Homeobox A9",
            "Metastasis-associated lung adenocarcinoma transcript 1",
            "p53",
            "Retinoblastoma protein",
            "Adenomatous polyposis coli",
            "BRCA1",
            "BRCA2",
            "p16^INK4A and p14^ARF",
            "PTEN",
            "Von Hippel-Lindau",
            "LKB1",
            "DPC4",
            "Neurofibromin 1",
            "Merlin",
            "Wilms tumor 1",
            "Hamartin",
            "Tuberin",
            "Menin",
            "Patched 1",
            "E-cadherin",
            "ARID1A",
            "BAF180",
            "BRCA1-associated protein 1",
            "FBW7",
            "KEAP1",
            "ATM",
            "ATR",
            "MLH1",
            "MSH2",
            "MLL2",
            "MLL3",
            "UTX",
            "CHEK2",
            "DICER1",
            "SPOP",
            "ASXL1",
            "CBP",
            "p300",
            "CYLD",
        ],
        "gene_symbol": [
            "ABL1",
            "AKT1",
            "ALK",
            "BCL2",
            "BRAF",
            "CCND1",
            "CDK4",
            "EGFR",
            "ERBB2",
            "FGFR3",
            "HRAS",
            "KRAS",
            "NRAS",
            "IDH1",
            "JAK2",
            "KIT",
            "MET",
            "MYC",
            "MYCN",
            "MYCL",
            "NOTCH1",
            "NFE2L2",
            "PDGFRA",
            "PIK3CA",
            "RET",
            "SRC",
            "TERT",
            "MDM2",
            "CTNNB1",
            "MYD88",
            "CALR",
            "NTRK1",
            "EZH2",
            "HOXA9",
            "MALAT1",
            "TP53",
            "RB1",
            "APC",
            "BRCA1",
            "BRCA2",
            "CDKN2A",
            "PTEN",
            "VHL",
            "STK11",
            "SMAD4",
            "NF1",
            "NF2",
            "WT1",
            "TSC1",
            "TSC2",
            "MEN1",
            "PTCH1",
            "CDH1",
            "ARID1A",
            "PBRM1",
            "BAP1",
            "FBXW7",
            "KEAP1",
            "ATM",
            "ATR",
            "MLH1",
            "MSH2",
            "KMT2D",
            "KMT2C",
            "KDM6A",
            "CHEK2",
            "DICER1",
            "SPOP",
            "ASXL1",
            "CREBBP",
            "EP300",
            "CYLD",
        ],
        "uniprot_id": [
            "P00519",
            "P31749",
            "Q9UM73",
            "P10415",
            "P15056",
            "P24385",
            "P11802",
            "P00533",
            "P04626",
            "P22607",
            "P01112",
            "P01116",
            "P01111",
            "O75874",
            "O60674",
            "P10721",
            "P08581",
            "P01106",
            "P04198",
            "P12524",
            "P46531",
            "Q16236",
            "P16234",
            "P42336",
            "P07949",
            "P12931",
            "O14746",
            "Q00987",
            "P35222",
            "Q99836",
            "P27797",
            "P04629",
            "Q15910",
            "P31269",
            "N/A",
            "P04637",
            "P06400",
            "P25054",
            "P38398",
            "P51587",
            "P42771",
            "P60484",
            "P40337",
            "Q15831",
            "Q13485",
            "P21359",
            "P35240",
            "P19544",
            "Q92574",
            "P49815",
            "O00255",
            "Q13635",
            "P12830",
            "O14497",
            "Q86U86",
            "Q92560",
            "Q969H0",
            "Q14145",
            "Q13315",
            "Q13535",
            "P40692",
            "P43246",
            "O14686",
            "Q8NEZ4",
            "O15550",
            "O96017",
            "Q9UPY3",
            "O43791",
            "Q8IXJ9",
            "Q92793",
            "Q09472",
            "Q9NQC7",
        ],
    }
    proteins = pd.DataFrame(data)

    # remove entries not in database
    uniprot_to_remove = ["N/A", "P04198", "P12524", "P27797", "P31269", "Q15831", "Q92574", "P49815", "O14497", "Q92560", "Q969H0", "Q13315", "Q13535", "Q8IXJ9"]
    proteins = proteins[~proteins["uniprot_id"].isin(uniprot_to_remove)]
    proteins.reset_index(drop=True, inplace=True)  # Reset indices

    return proteins

def get_cancer_ligands():
    ligands = pd.read_pickle("./hydrocarbons_with_smiles.pkl")
    return ligands

def get_proteins_ligands(protein_set):
    if protein_set == "thyroid":
        proteins = get_thyroid_proteins()
        ligands = get_thyroid_ligands()
    elif protein_set == "autism":
        proteins = get_autism_proteins()
        ligands = get_autism_ligands()
    elif protein_set == "cancer":
        proteins = get_cancer_proteins()
        ligands = get_cancer_ligands()
    else:
        raise ValueError("Invalid protein set. Choose 'thyroid', 'autism', or 'cancer'.")

    return proteins, ligands

def get_docking_score(proteins, ligands, OUTPUT_DIR):
    n_protein = len(proteins)
    n_ligand  = len(ligands)
    score_size = (n_ligand, n_protein)

    docking_score = np.full(score_size, np.nan)
    for i, ligand_row in ligands.iterrows():
        ligand = ligand_row["SMILES"]
        for j, protein_row in proteins.iterrows():
            uniprot_id = protein_row["uniprot_id"]
            fname = make_valid_fname(uniprot_id, ligand)
            fpath = os.path.join(OUTPUT_DIR, fname)

            if os.path.exists(fpath):
                with open(fpath) as f:
                    data = json.load(f)
                    docking_score[i, j] = data["result"]["docking_score"]

    return docking_score