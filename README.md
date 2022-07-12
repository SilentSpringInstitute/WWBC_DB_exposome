# NTA for woman nurse and firefighters
The goal of this project is to realize Non targeted analysis for nurse<br>
The pipeline is divided in 5 steps:
- Step 1: Women Worker Biomonitoring Chemical Database Development (v.2)
- Step 2: Non-targeted LC-QTOF/MSfull scan analysis procedure (and developing the features list)
- Step 3: Review of features and criteria to select features for fragmentation
- Step 4: Fragmentation and annotation
- Step 5: Chemical selection criteria for targeted confirmation and quantification 

<br>
Please see the full pipeline in the py note book: xxxx

<br><br>

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



