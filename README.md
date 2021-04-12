# WWBC DB exposome
Objective: prepare an chemical database for MS non targeted analysis. 


## Dependancies
- python 3.8 (development on Conda - django-env)
- RDKit (> 3.1): http://rdkit.org/docs/index.html
- molVS (> 1): https://molvs.readthedocs.io/en/latest/index.html
- CompDesc ($pip install -i https://test.pypi.org/simple/ CompDesc)
- biotransformer (http://biotransformer.ca/)
- request database.init for DTXSID search

## 0 Dataset
Extract different list of chemicals....

## 1 Step to clean (./../results/prepDB/')

- Classic protocol for cleaning from best practice in the field using molvS (https://molvs.readthedocs.io/en/latest/)
- H removal
- SMILES standardization
- sanitize (pH7) 
- ion stripping

Molecular descriptors are computed using CompDesc including MolWeight


## 2 Filter list of chemical (./../results/prepDB/')
- filter chemical based on MW - [100-1000] - check MolWeight from cleaned SMILES and from MASS_ISOTROPIQUE if not cleaned SMILES

## 3 Compute metabolites
- Run biotransformer from the SMILES_cleaned
- limit metabolite for chemicals list  ["Drug_UCSF_PXYS", "Drug_most comon and haz", "Disinfectant", "FRs", "PFAS", "MC", "MGDev", "ERactive_bin", "E2Up_bin", "P4Up_bin", "pesticidemammarytumors_bin", "nitroPAH_bin", "pellizzari_bin", "phthalate"]
- Clean SMILES for each metabolite
- Define name using DTXSID-precursor--reaction
- cross with the DSSTOX DB local version (lastest 2020)
- filter MW [100-1000] from values of MW from the biotransformer sofware and lipinski violations (< 3)

## 4 Combine chemicals in DB list and metabolites ("./../results/...final.csv")
- Add ID (int) column as unique ID 

# TODO list
- ~~4-07-2021: reload CASRN from DTXSID and add ' in the CASRN col~~
- ~~4-07-2021: check col Phthalate duplicate (P-AM)~~
