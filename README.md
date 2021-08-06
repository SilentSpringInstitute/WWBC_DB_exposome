# WWBC DB exposome
The goal of this project is to realize Non targeted analysis
## 1. Step 1: prepare chemical databse knowledge based for annotation
### 1.1. Dataset
Extract different list of chemicals....

### 1.2. Step to clean chemicals (./../results/prepDB/')

- Classic protocol for cleaning from best practice in the field using molvS (https://molvs.readthedocs.io/en/latest/)
- H removal
- SMILES standardization
- sanitize (pH7) 
- ion stripping

Molecular descriptors are computed using CompDesc including MolWeight


### 1.3. Filter list of chemical (./../results/prepDB/')
- filter chemical based on MW - [100-1000] - check MolWeight from cleaned SMILES and from MASS_ISOTROPIQUE if not cleaned SMILES

### 1.4. Compute metabolites
- Run biotransformer from the SMILES_cleaned
- limit metabolite for chemicals list  ["Drug_UCSF_PXYS", "Drug_most comon and haz", "Disinfectant", "FRs", "PFAS", "MC", "MGDev", "ERactive_bin", "E2Up_bin", "P4Up_bin", "pesticidemammarytumors_bin", "nitroPAH_bin", "pellizzari_bin", "phthalate"]
- Clean SMILES for each metabolite
- Define name using DTXSID-precursor--reaction
- cross with the DSSTOX DB local version (lastest 2020)
- filter MW [100-1000] from values of MW from the biotransformer sofware and lipinski violations (< 3)

### 1.5.Combine chemicals in DB list and metabolites ("./../results/...final.csv")
- Add ID (int) column as unique ID 


## Step 2: From NTA make the annotation

## Step 3: Filter results for fragmentation



# Installation
Develop on Ubuntu on Windows 10 (wsl2 - ubuntu 20.4)
## Dependancies
- python (3.9)
- biotransformer (http://biotransformer.ca/) => request java install
- request database.init for DTXSID search
- R (4.1.0)

### Python libraries
- RDKit (> 3.1): http://rdkit.org/docs/index.html
- molVS (> 1): https://molvs.readthedocs.io/en/latest/index.html
- CompDesc ($pip install -i https://test.pypi.org/simple/ CompDesc)

### R libraries + dependancies
- ncdf4 ($sudo apt-get install r-cran-ncdf4)
- libmagick ($sudo apt install libmagick++-dev)
- libxml2-dev ($sudo apt-get install libxml2-dev)
- libcurl4-openssl-dev ($sudo apt -y install libcurl4-openssl-dev)
- gsl ($sudo apt install r-cran-gsl)
- all of the libraries can be installed using the script: ./R/requirements.R



# TODO list
- ~~4-07-2021: reload CASRN from DTXSID and add ' in the CASRN col~~
- ~~4-07-2021: check col Phthalate duplicate (P-AM)~~
- 21-07-2021: reorganize source to integrate all pieces in a pipeline



# diff between database 4-27 and 6-30 for record
  [1] "DTXSID0020103"          "DTXSID00231231"         "DTXSID0024715"         
  [4] "DTXSID0026701"          "DTXSID0040258"          "DTXSID0044525"         
  [7] "DTXSID0045597"          "DTXSID0046270"          "DTXSID0046274"         
 [10] "DTXSID0047117"          "DTXSID0047955"          "DTXSID0048767"         
 [13] "DTXSID0049232"          "DTXSID0049234"          "DTXSID0049286"         
 [16] "DTXSID0049363"          "DTXSID0049367"          "DTXSID0049369"         
 [19] "DTXSID0051863"          "DTXSID1024334"          "DTXSID1040326"         
 [22] "DTXSID1045584"          "DTXSID1045586"          "DTXSID1045798"         
 [25] "DTXSID1047940"          "DTXSID1047942"          "DTXSID1047946"         
 [28] "DTXSID1049275"          "DTXSID1049304"          "DTXSID1049356"         
 [31] "DTXSID1049805"          "DTXSID1051850"          "DTXSID1062253"         
 [34] "DTXSID2020131"          "DTXSID2027836"          "DTXSID2033281"         
 [37] "DTXSID2035013"          "DTXSID2043909"          "DTXSID2047939"         
 [40] "DTXSID2049216"          "DTXSID2049218"          "DTXSID2049391"         
 [43] "DTXSID2063420"          "DTXSID3020338"          "DTXSID3021518"         
 [46] "DTXSID3023005"          "DTXSID3024603"          "DTXSID3029283"         
 [49] "DTXSID3031521"          "DTXSID3033983"          "DTXSID3039240"         
 [52] "DTXSID3045227"          "DTXSID3045562"          "DTXSID3047920"         
 [55] "DTXSID3047922"          "DTXSID3047926"          "DTXSID3047970"         
 [58] "DTXSID3049251"          "DTXSID3049253"          "DTXSID3049388"         
 [61] "DTXSID3051707"          "DTXSID3068209"          "DTXSID3068413"         
 [64] "DTXSID4020115"          "DTXSID4020327"          "DTXSID4024983"         
 [67] "DTXSID4027076"          "DTXSID4029692"          "DTXSID4042123"         
 [70] "DTXSID4044379"          "DTXSID4045555"          "DTXSID4046492"         
 [73] "DTXSID4046604"          "DTXSID4046608"          "DTXSID4047335"         
 [76] "DTXSID4047387"          "DTXSID4047888"          "DTXSID4047913"         
 [79] "DTXSID4047915"          "DTXSID4047917"          "DTXSID4047965"         
 [82] "DTXSID4047967"          "DTXSID4047969"          "DTXSID4049244"         
 [85] "DTXSID4049290"          "DTXSID4049292"          "DTXSID4049298"         
 [88] "DTXSID4049325"          "DTXSID4049402"          "DTXSID4071429"         
 [91] "DTXSID5020235"          "DTXSID5020811"          "DTXSID5021461"         
 [94] "DTXSID5023796"          "DTXSID5023825"          "DTXSID5025948"         
 [97] "DTXSID5027279"          "DTXSID5034981"          "DTXSID5041803"         
[100] "DTXSID5044316"          "DTXSID5045596"          "DTXSID5046013"         
[103] "DTXSID5046308"          "DTXSID5047790"          "DTXSID5047871"         
[106] "DTXSID5047904"          "DTXSID5047954"          "DTXSID5048924"         
[109] "DTXSID5048928"          "DTXSID5049231"          "DTXSID5049233"         
[112] "DTXSID5049285"          "DTXSID5049368"          "DTXSID5052385"         
[115] "DTXSID5065562"          "DTXSID5068603"          "DTXSID6020226"         
[118] "DTXSID6020513"          "DTXSID6024048"          "DTXSID60241694"        
[121] "DTXSID6024331"          "DTXSID6024961"          "DTXSID6025018"         
[124] "DTXSID6026373"          "DTXSID6027214"          "DTXSID6032061"         
[127] "DTXSID6034265"          "DTXSID6034762"          "DTXSID6042313"         
[130] "DTXSID6044723"          "DTXSID6045583"          "DTXSID6047943"         
[133] "DTXSID6047947"          "DTXSID6049274"          "DTXSID6049351"         
[136] "DTXSID6060478"          "DTXSID7020340"          "DTXSID7020506"         
[139] "DTXSID7020899"          "DTXSID7020928"          "DTXSID70227388"        
[142] "DTXSID7026314"          "DTXSID7027097"          "DTXSID7034672"         
[145] "DTXSID7034836"          "DTXSID7040316"          "DTXSID7040449"         
[148] "DTXSID7041126"          "DTXSID7044504"          "DTXSID7044552"         
[151] "DTXSID7044558"          "DTXSID7044928"          "DTXSID7046207"         
[154] "DTXSID7047930"          "DTXSID7067459"          "DTXSID7069283"         
[157] "DTXSID8020337"          "DTXSID8020622"          "DTXSID8020703"         
[160] "DTXSID8021301"          "DTXSID8021351"          "DTXSID8031861"         
[163] "DTXSID80322696"         "DTXSID8034376"          "DTXSID8038770"         
[166] "DTXSID8041909"          "DTXSID8044622"          "DTXSID8045989"         
[169] "DTXSID8046400"          "DTXSID8047183"          "DTXSID8049206"         
[172] "DTXSID8049252"          "DTXSID8049254"          "DTXSID8051704"         
[175] "DTXSID8066842"          "DTXSID9020875"          "DTXSID90233502"        
[178] "DTXSID9026219"          "DTXSID9029695"          "DTXSID9034361"         
[181] "DTXSID9035204"          "DTXSID9041578"          "DTXSID9042124"         
[184] "DTXSID9042415"          "DTXSID9044796"          "DTXSID9047889"         
[187] "DTXSID9047914"          "DTXSID9047916"          "DTXSID9048988"         
[190] "DTXSID9049241"          "DTXSID9049291"          "DTXSID9049297"         
[193] "DTXSID9049299"          "DTXSID9058319"          "DTXSID9058361"         
[196] "DTXSID9062194"          "DTXSID9070531"          "DTXSID5034519"         
[199] "DTXSID0042167-BTMR0168" "DTXSID0042167-BTMR0206" "DTXSID2021151-BTMR0166"
[202] "DTXSID5041726-BTMR0166" "DTXSID7020184-BTMR0166" "DTXSID8021224-BTMR0166"
[205] "DTXSID8031865-BTMR0167" "DTXSID9027073-BTMR0167" "DTXSID9027073-BTMR0244"