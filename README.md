# NTA for woman nurse and firefighters
The goal of this project is to realize Non targeted analysis
## 1. Women Worker Biomonitoring Chemical Database Development (v.2)
### 1.1. Datasets - input files
> - UCSF haz and drugs inventory (<i>UCSF_Haz_Drugs Inventory_FINAL 2_2018_wDTSXIDs.csv</i>)
> - List of drugs from survey on nurses (<i>nurses_otherdrugs_rar_10272021.csv</i>)
> - First version of the WWBC database (<i>WWBC_MS_database_10.28.21.csvv</i>)

### 1.2. Step to prepare chemicals
- Prepare the chemicals structures using the best practices and including:
  - H removal
  - SMILES standardization
  - sanitize (pH7) 
  - ion stripping
- Compute isotopic MW
- Map each chemicals on the DSSTOX database (https://comptox.epa.gov/dashboard/)
- filter chemicals based on MW, excluded chemicals with low MW (< 100Da) and high MW (> 1000Da) and chemicals that do not match any DSSTOX id

### 1.3. Compute metabolites
- Run biotransformer software on the standardize SMILES (http://biotransformer.ca/)
- prepare metabolites using the same protocol that before
- Define name of metabolite: <b>MTX-id</b>
- Filter metabolites based on MW, excluded metabolites with low MW (< 100Da) and high MW (> 1000Da) MW and metabolite with poor absorption Lipinski violations (> 3)

## Step 2: Non-targeted LC-QTOF/MSfull scan analysis procedure (and developing the features list)

<br>
<br>


## Step 3: Review of features and criteria to select features for fragmentation

<br>
<br>

## Step 4: Fragmentation and annotation

<br>
<br>

## Step 5: Chemical selection criteria for targeted confirmation and quantification 

<br>
<br>
<br>

---

# Installation
Develop on Ubuntu-wsl on Windows 10 (wsl2 - ubuntu 20.4)
## Dependencies
- python (3.9)
- biotransformer (http://biotransformer.ca/) => request java install
- request database.init for DTXSID search
- R (4.1.0)

### Python libraries
- RDKit (> 3.1): http://rdkit.org/docs/index.html
- molVS (> 1): https://molvs.readthedocs.io/en/latest/index.html
- CompDesc ($pip install -i https://test.pypi.org/simple/ CompDesc)

### R libraries + dependencies
- ncdf4 ($sudo apt-get install r-cran-ncdf4)
- libmagick ($sudo apt install libmagick++-dev)
- libxml2-dev ($sudo apt-get install libxml2-dev)
- libcurl4-openssl-dev ($sudo apt -y install libcurl4-openssl-dev)
- gsl ($sudo apt install r-cran-gsl)
- all of the libraries can be installed using the script: <b>./R/requirements.R</b>



