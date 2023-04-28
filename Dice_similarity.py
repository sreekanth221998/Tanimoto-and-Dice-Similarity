# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:01:47 2023

@author: Sreekanth
"""
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem

# Read the dataset
df = pd.read_csv("oligomers_50.csv")
df = df.copy()

for index, row in df.iterrows():
    fp_previous = None
    for i in range(1, 51):
        if i == 1:
            mol1 = Chem.MolFromSmiles(row["Monomer"])
            fp1 = AllChem.GetMorganFingerprint(mol1, 4)
            mol2 = Chem.MolFromSmiles(row["RU_W"])
            fp2 = AllChem.GetMorganFingerprint(mol2, 4)
            coefficient = DataStructs.DiceSimilarity(fp1, fp2)
            df.at[index, "Dice_coefficient_Monomer_RU_W"] = coefficient
            fp_previous = fp2
        else:
            mol_current = Chem.MolFromSmiles(row[f"Oligomer_dp{i}"])
            fp_current = AllChem.GetMorganFingerprint(mol_current, 4)
            coefficient = DataStructs.DiceSimilarity(fp_previous, fp_current)
            df.at[index, f"Dice_coefficient_Oligomer_dp{i-1}_Oligomer_dp{i}"] = coefficient
            fp_previous = fp_current

# Save the modified dataframe to a new CSV file
df.to_csv("Dice_similarity.csv", index=False)