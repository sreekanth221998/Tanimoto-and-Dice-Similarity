#@author: Sreekanth

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
            coefficient = DataStructs.TanimotoSimilarity(fp1, fp2)
            df.at[index, "Tanimoto_coefficient_Monomer_RU_W"] = coefficient
            fp_previous = fp2
        else:
            mol_current = Chem.MolFromSmiles(row[f"Oligomer_dp{i}"])
            fp_current = AllChem.GetMorganFingerprint(mol_current, 4)
            coefficient = DataStructs.TanimotoSimilarity(fp_previous, fp_current)
            df.at[index, f"Tanimoto_coefficient_Oligomer_dp{i-1}_Oligomer_dp{i}"] = coefficient
            fp_previous = fp_current

# Save the modified dataframe to a new CSV file
df.to_csv("tanimoto_similarity.csv", index=False)







