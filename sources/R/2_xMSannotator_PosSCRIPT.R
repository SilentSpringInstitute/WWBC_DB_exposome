####R script for annotation of m/z features processed via xcms
## INSTALL xMSannotator
# install.packages(c('BiocManager','data.table','digest'))
# remotes::install_github('omegahat/XMLSchema')
# remotes::install_github('cran/SSOAP')
# BiocManager::install(c("SSOAP","KEGGREST","pcaMethods","Rdisop","GO.db","matrixStats","WGCNA"))
# devtools::install_github("yufree/xMSannotator")

library(xMSannotator)
# library(dplyr)

getwd()

setwd("/Users/jesstro/")
###Upload the xcms-preprocessed files 



dataC <- read.csv("Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processingPOS/20210423POS_xset.csv")
df<-dataC[,c(2,5,10:141)] ###Keep only columns: mz, rt and each samples

colnames(df)[2]<-"time"
df<-unique(df)

max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.7 #correlation threshold between adducts/isotopes
max_isp=5
mass_defect_window=0.01

num_nodes<-2   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption

db_name="Custom" #other options: KEGG, LipidMaps, T3DB, HMDB
status=NA 
num_sets<-3000 #number of sets into which the total number of database entries should be split into;

mode<-"pos" #ionization mode
queryadductlist=c("M+H","M+NH4","M+Na","M+H-H2O","M+H-2H2O") #for positive ionization 
# queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H") #for negative ionization
adduct_weights<-cbind.data.frame(Adduct=c("M+H","M-H"),Weight=c(5,5))

outloc<-"C:/Users/jesstro/Desktop/Pos_annotator_20210708/" ##Folder where to save annotation results

# ###Upload MS library
# dsstox <- read.csv("c:/Users/jesstro/Box/SHE lab/WWBC/Raw data [de-ID]/WWBC NTA Results & Database/CompTox_17March2019_SelectMetaData.csv", header = T)
# dsstox<-dsstox[,c(3,1,6,7)] ###Keep only columns: DTXSID, chemical name, formula and monoisotopic mass
# colnames(dsstox)[1]<-"DTXSID"
# colnames(dsstox)[2]<-"Name"
# colnames(dsstox)[3]<-"Formula"
# colnames(dsstox)[4]<-"MonoisotopicMass"
# dsstox$MonoisotopicMass<-as.numeric(as.character(dsstox$MonoisotopicMass)) ##Should be numeric
# dsstox<-filter(dsstox, MonoisotopicMass>=80 & MonoisotopicMass<=1200) ##Search chemical with mass ranging from 80 to 1,200 DA

wwbcDB <- read.csv("Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/Code/data_files_csv/WWBC_MS_database_6.30.21_prepForAnotation - WWBC_MS_database_6.30.21_prepForAnotation.csv", header = TRUE, stringsAsFactors = FALSE)
wwbcDB <- wwbcDB[, c(3,5,10,11)] ###Keep only columns: DTXSID, chemical name, formula and molecular weight cleaned
colnames(wwbcDB)
colnames(wwbcDB)[3] <- "Formula"
colnames(wwbcDB)[2] <- "Name"
colnames(wwbcDB)[4] <- "MonoisotopicMass"

#Get NA's after converting MonoisotopicMass to numeric. they appear to be blank cells. remove before continuing. 
## confirm with vincent and alyssia if they are indeed blank

wwbcDB$Formula<-as.character(wwbcDB$Formula)

#remove duplicates
 wwbcDB<-wwbcDB[!duplicated(wwbcDB$DTXSID),]

#convert MonoisotopicMass to numeric and keep those >=80 and <=1200
wwbcDB$MonoisotopicMass<-as.numeric(as.character(wwbcDB$MonoisotopicMass)) ##Should be numeric

# sum(wwbcDB$MonoisotopicMass>1200)
# wwbcDB<-filter(wwbcDB, MonoisotopicMass>=80 & MonoisotopicMass<=1200)

# wwbcDB <- filter(wwbcDB, MonoisotopicMass!="")
sum(is.na(wwbcDB$MonoisotopicMass))

#name database customDB for use in the next line of code
customDB<-wwbcDB


##Annotation - Make sure you are setting the correct ionization mode (pos or neg)
annotres<-multilevelannotation(
  df, 
  max.mz.diff = max.mz.diff, 
  max.rt.diff = max.rt.diff, 
                               cormethod = "pearson", 
                               num_nodes = num_nodes, queryadductlist = queryadductlist, mode = mode,
                               outloc=outloc, db_name = db_name, adduct_weights = NA, num_sets = num_sets,
                               allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, customIDs = NA,
                               missing.value = NA, deepsplit = 2, 
                               networktype = "unsigned",
                               minclustsize = 10, module.merge.dissimilarity = 0.2, filter.by = c("M+H"),
                               redundancy_check = TRUE,
                               min_ions_perchem = 1, biofluid.location = NA, origin = NA,
                               status = status, boostIDs = NA, max_isp = 5,
                               MplusH.abundance.ratio.check = FALSE, 
                               customDB = customDB, HMDBselect = "union", mass_defect_window = 0.01,
                               mass_defect_mode = "pos",  pathwaycheckmode = "pm")
 