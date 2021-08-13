#!/usr/bin/env Rscript
library(xMSannotator)





prepXSET = function(p_xset){
  
  dataC <- read.csv(p_xset)
  dataC<-dataC[,c(2,5,10:141)]

  # rename rt -> time
  colnames(dataC)[which(colnames(dataC) == "rt")] = "time"
  
  #df<-dataC[,c(2,5,10:141)] ###Keep only columns: mz, rt and each samples including blank values DIWB
  #colnames(df)[2]<-"time"
  dataC<-unique(dataC)
  
  return(dataC)
}


prepCustomDB = function(p_DB_to_match){
  
  wwbcDB <- read.csv(p_DB_to_match, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  # only save
  l_cols = c("ID",  "name_original", "formula_cleaned", "Molweight_cleaned")
  #"DB_name", "DTXSID",
  wwbcDB <- wwbcDB[, l_cols] ###Keep only columns: DTXSID, chemical name, formula and molecular weight cleaned
  
  # need to rename for annotation
  # rename formula_cleaned -> formula
  colnames(wwbcDB)[which(colnames(wwbcDB) == "formula_cleaned")] = "Formula"
  
  # rename Molweight_cleaned -> MonoisotopicMass
  colnames(wwbcDB)[which(colnames(wwbcDB) == "Molweight_cleaned")] = "MonoisotopicMass"
  
  # rename name_original -> Name
  colnames(wwbcDB)[which(colnames(wwbcDB) == "name_original")] = "Name"
  
  
  # format col
  wwbcDB$Formula<-as.character(wwbcDB$Formula)
  wwbcDB$Name<-as.character(wwbcDB$Name)
  wwbcDB$ID<-as.character(wwbcDB$ID)
  wwbcDB$MonoisotopicMass<-as.numeric(as.character(wwbcDB$MonoisotopicMass)) 
  
  # remove when MonoisotopicMass => NA
  wwbcDB = wwbcDB[!is.na(wwbcDB$MonoisotopicMass),]
  
  return(wwbcDB) 
  
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pr_row_data = args[1]
pr_out = args[2]
mode = args[3]

## linux - wsl2
p_xset = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg/xset.csv"
mode = "neg"
pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg"
p_DB_to_match = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/WWBC_MS_database_4.7.21_prepForAnotation.csv"

##windows native
#p_xset = "c://Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg/xset.csv"
#mode = "neg"
#pr_out = "c://Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg"
#p_DB_to_match = "c://Users/AlexandreBorrel/research/SSI/NTA/results/WWBC_MS_database_4.7.21_prepForAnotation.csv"



# XSET
d_Xset = prepXSET(p_xset)

# DB custom
d_DB = prepCustomDB(p_DB_to_match)



### Parameters ###
##################
max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.7 #correlation threshold between adducts/isotopes
max_isp=5
mass_defect_window=0.01
num_nodes = 2   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption
db_name="Custom" #other options: KEGG, LipidMaps, T3DB, HMDB
status=NA 
num_sets<-3000 #number of sets into which the total number of database entries should be split into;

# ionization criteria depending of mode
if (mode == "neg"){
  queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H")
}else if(mode == "pos"){
  queryadductlist=c("M+H","M+NH4","M+Na","M+H-H2O","M+H-2H2O")
}
#adduct_weights<-cbind.data.frame(Adduct=c("M+H","M-H"),Weight=c(5,5))

##Annotation - Make sure you are setting the correct ionization mode (pos or neg) ##
####################################################################################
annotres<-multilevelannotation(d_Xset,
                               max.mz.diff = max.mz.diff,
                               max.rt.diff = max.rt.diff, 
                               cormethod = "pearson",
                               num_nodes = num_nodes,
                               queryadductlist = queryadductlist,
                               mode = mode,
                               outloc=pr_out,
                               db_name = db_name,
                               adduct_weights = NA,
                               num_sets = num_sets,
                               allsteps = TRUE,
                               corthresh = 0.7,
                               NOPS_check = TRUE,
                               customIDs = NA,
                               missing.value = NA,
                               deepsplit = 2, 
                               networktype = "unsigned",
                               minclustsize = 10,
                               module.merge.dissimilarity = 0.2,
                               filter.by = c("M-H"),
                               redundancy_check = TRUE,
                               min_ions_perchem = 1,
                               biofluid.location = NA,
                               origin = NA,
                               status = status,
                               boostIDs = NA,
                               max_isp = 5,
                               MplusH.abundance.ratio.check = FALSE, 
                               customDB = d_DB, 
                               HMDBselect = "union", 
                               mass_defect_window = 0.01,
                               mass_defect_mode = "pos",  
                               pathwaycheckmode = "pm")


