#!/usr/bin/env Rscript

# From: 
#title: "Negative FF_fragmentation selecton"
#author: "Jessica Trowbridge"
#date: "6/4/2020"
  
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(broom)
library(perm)
library(tidyr)
library(cvcqv)
library(EnvStats)
library(stringr)



# Selecting chemical features for fragmentation in Negative Mode

#This file cleans and selects features for fragmentation. The data imported here has been peak corrected with xcms and then pre-processed urins XMSannotator. The resulting .csv file is the relative peak intensities and features identified through matching to the WWBC library. In this document I will ~~apply a method blank correction criteria~~ and identify chemicals of interest for fragmentation. 


#Step 1—Firefighter related?

#A.	Flag if FB is <50% of Visit 1 average
#B.	Compound detected in at least 3 participants in visit 1, AND
#C.	Different peak area (higher or lower) in visit 1 versus visit 2 or 3. (boxplot of feature by visit and MB), OR
#D.	Different detection frequency (higher or lower) in visit 1 versus visit 2 or 3 (but detected in at least 3 participants)

#Step 2—Other novel chemical of interest?

#A.	Feature detected in at least 3 participants; flag
#B.	Flag feature if FB is <50% of average visit 1 
#C.	Matches to BC-relevant chemicals (eg. give score = 1 if increases E2 and score = 1 if increases P4up, MC, ER active, add up for a BC score) OR
#D.	Matches to PFAS, & FR listed in WWBD OR
#E.	UV protectors (benzophenone-3, oxybenzone, etc.  look at BCPP website)
#F.	Exposure predicted as confidence check
#G.  Not a UCSF drug list chemical


#Step 3-
#Meets FB filter and ubiquitous breast cancer relevant

#Additional details of chemicals included in the fragmentation list: 

#Criteria for selecting chemicals of interest for confirmation:
#1.	Structural criteria
#a.	Acrylates
#b.	PAH, combustion products
#c.	Halogenated compounds (FR, PFAS)
#OR
#2.	Mammary carcinogen (MC) related chemicals
#OR
#3.	Is a metabolite

#later: 

#a)	not been measured in large biomonitoring studies. 
#b)	Exposure predicted


prepXSETForPostAnalysis = function(p_xset){
  
  #read in data. 
  d_xset <- read.csv(p_xset) 
  #d_xset <- d_xset[ , c(2,5,10:137)]
  # colnames(d_xset)
  
  #make feature_id variable and reorganize columns
  d_xset$mz<-round(d_xset$mz,5)
  d_xset$rt <- round(d_xset$rt, 1)
  d_xset$featureid<-as.character(paste(d_xset$mz,d_xset$rt,sep = "-"))
  
  #dataB_neg_xset<-d_xset[,-131]
  #row.names(dataB_neg_xset)<-neg_xset[,131]
  
  return(d_xset)
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pr_row_data = args[1]
pr_out = args[2]

p_xset = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg/xset.csv"
mode = "neg"
p_id = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/FF_raw_data/FF_id_visit.csv"
p_annotation = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg/Stage5.csv"
p_DB = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/WWBC_MS_database_4.7.21_prepForAnotation.csv"

# read xset data
d_xset = prepXSETForPostAnalysis(p_xset)

#read id by visit
d_id = read.csv(p_id, header = T)


###Extract FF id and visit
d_id$dtsc_id<-substr(d_id$id,1,7)
d_id$FF_id<-substr(d_id$id,9,12)
d_id$visit<-substr(d_id$id,14,20)
d_id$dtsc_run <- substr(d_id$id, 30,33)

#make a sample ID to match features data.
d_id$dtsc_fid <- as.character(paste(d_id$dtsc_id,d_id$dtsc_run,sep = "."))  


#Participants 216 and 227 have more than 3 sample collections. They had repeat fires and were labeles as fires 4,5,6. Renumber with 501 and 501 and relabel sample collections. 
#participant 216 did 2 complete sets. 
#participant 227 had a second fire before sample 3. So renumber sample 3 to 1 and so forth. 

#Renumber the 2 participants who gave samples for 2 fires and assign the visits 1,2,3 so they are included in the analysis. 
d_id <- d_id%>% mutate(new_visit= ifelse(FF_id == "F227" & visit == "Visit 3", "Visit 1", 
                                         ifelse(FF_id =="F227" & visit == "Visit 4", "Visit 2", 
                                                ifelse(FF_id =="F227" &visit == "Visit 5", "Visit 3", 
                                                       ifelse( FF_id =="F216" & visit =="Visit 4", "Visit 1", 
                                                               ifelse(FF_id =="F216" & visit == "Visit 5", "Visit 2", 
                                                                      ifelse(FF_id =="F216" & visit == "Visit 6", "Visit 3", visit)))))), 
                       newFF_id = ifelse(FF_id == "F227" & visit == "Visit 3", "F227.2", 
                                         ifelse(FF_id == "F227" & visit == "Visit 4", "F227.2", 
                                                ifelse(FF_id == "F227" & visit == "Visit 5", "F227.2",
                                                       ifelse(FF_id == "F216" & visit == "Visit 4", "F216.2",
                                                              ifelse(FF_id == "F216" & visit == "Visit 5", "F216.2",
                                                                     ifelse(FF_id == "F216" & visit == "Visit 6", "F216.2", FF_id)))))))


# ## Before Matching: 
# * exclude if missing is >50% for visit 1
# * calculate DF of each visit
# * PCA to check for batch effects. 
# 
# 
# ## After matching with xmsannotator output: 
# * calculate a DF for each group/visit
# ** make long-form data with intensity and filenames (sample id?)
# #first
# 
# calculate DF for each feature
# mean intensity by group and SD and CV
# filter out endogenous molecules: CV<30 between the 3 time points becuase not much variation. 
# When comparing between timepoints, log transform. 
# 
# 
# 
# The value under each ID is the relative intensity of each feature. 
# - calculate the CV for the two replicate samples. If large CV (>50) then the signal is not stable and can make NA or 0. This is the technical variation (versus the biological variation which would be the CV across the 3 samples)
# - combine the two runs with average
# - if the relative intensity is below 1000 then is noise (can convert to 0? or NA? or leave as is, but limit the calculations to values above 1000)
# - then can calculate the detection frequency with the cut off of 1000.
# 


# colnames(neg_xset)
#concert to long-format to use tidy to summarize statististics
df_long <- gather(d_xset, "sample_id", "intensity", BD01346.r001:DIWB2.B7 )

#keep relevant vars and rename the last one
d_id <- d_id[, c("dtsc_id","dtsc_fid", "visit", "FF_id", "new_visit", "newFF_id")]


#merge with dtsc data so we have visit # and ID
long_df <- merge(df_long, d_id, by.x = "sample_id", by.y = "dtsc_fid", all.x = TRUE)  #do not drop the Blank and QAQC variables, they do not have a study ID will need to make a new variable with group: visit 1,2,3, MB
#remove the blank samples that are only present in visit 1.

# head(dataA$featureid)

long_df$type<-substr(long_df$sample_id,1,2)

#label field blanks and Water blanks
long_df$newFF_id[is.na(long_df$newFF_id)] <- "FBlank"

long_df2 <- long_df%>%
  mutate(new_visit_blk = ifelse(newFF_id == "F197"|newFF_id== "F223"|newFF_id== "F231"|newFF_id== "F232", "FB", new_visit),
         new_visit_blk = ifelse(is.na(new_visit) & type == "DI", "DIWB", new_visit_blk),
         new_visit_blk = ifelse(is.na(new_visit) & type =="MB", "mb", new_visit_blk),
         group= factor(new_visit_blk))

# table(long_df2$new_visit_blk)

#DO NOT DO THIS: remove the field blank samples.
# long_df2 <- long_df2 %>%
#   filter(newFF_id != "F197"& newFF_id != "F223" & newFF_id != "F231"& newFF_id != "F232")  


## log tansform the intensity to calculate the mean intensity of visits and the Blanks?

## calculate the CV, mean, sd for each feature over for each visit and mb
##not log_transformed. Vincent suggests use log2
long_df2$intensity[long_df2$intensity==0]<-NA

print("Summary for intensity:")
(summary(long_df2$intensity))
# keep_mz <- cv_visit%>%
#  filter(coef_var>30) #206232 samples have cv>30; if CV<30 not much variation 
print("Number of NA for intensity:")
(sum(is.na(long_df2$intensity)))

print("Count visi:")
(table(d_id$new_visit))

#regarding the use of 0 versus NA. WE decided today 8/28/2020 that 
#we would substitute 0 for the non-detects in calculating the summary statistics
#that are used to select the non-matched features for fragmentation. We will 
#also use the median to compare difference between visit 1, 2 and 3 
#(looking at the absolute difference of intensity and detection frequency of the
#visits) and select the 10 features with the largest absolute difference 
#for each mode and for visit 1 vs 2 and 1 vs 3 N= 80?


#calculate summary statistics by for each feature. 
v1vs3_df <- long_df2 %>%
  group_by(mz, rt)%>%
  mutate(intensity1 = ifelse(is.na(intensity),0, intensity),
         log2_intensity = log2(intensity))%>%
  summarize( det_count_v1 = sum(!is.na(intensity[new_visit_blk == "Visit 1"])),
             det_count_v2 = sum(!is.na(intensity[new_visit_blk== "Visit 2"])),
             det_count_v3 = sum(!is.na(intensity[new_visit_blk == "Visit 3"])),
             # det_count_mb = sum(!is.na(intensity[new_visit_blk == "mb"])),
             det_count_FB = sum(!is.na(intensity[new_visit_blk == "FB"])),
             df_v1 = (det_count_v1/40)*100,
             df_v2 = (det_count_v2/38)*100, 
             df_v3 = (det_count_v3/38)*100,
             # df_mb = (det_count_mb/4)*100,
 
             #VISIT1
             mean_visit1 = mean(intensity1[new_visit_blk =="Visit 1"], na.rm = TRUE),
             log_mean_v1= mean(log2(intensity[new_visit_blk=="Visit 1"]), na.rm = TRUE),
             sd_visit1 = sd(intensity1[new_visit_blk =="Visit 1"], na.rm = TRUE),
             median_v1 = median(intensity1[new_visit_blk =="Visit 1"], na.rm = TRUE),
             gm_visit1 = geoMean(intensity[new_visit_blk =="Visit 1"], na.rm = TRUE),
             gsd_visit1 = geoSD(intensity[new_visit_blk =="Visit 1"], na.rm = TRUE),
             cv_visit1 = (sd(intensity1[new_visit_blk == "Visit 1"],na.rm = TRUE)
                          /mean(intensity1[new_visit_blk =="Visit 1"],na.rm = TRUE))*100,
             
             #VISIT2
             mean_visit2 = mean(intensity1[new_visit_blk =="Visit 2"], na.rm = TRUE),
             median_v2 = median(intensity1[new_visit_blk =="Visit 2"], na.rm = TRUE),
             sd_visit2 = sd(intensity1[new_visit_blk =="Visit 2"], na.rm = TRUE),
             gm_visit2 = geoMean(intensity[new_visit_blk =="Visit 2"], na.rm = TRUE),
             gsd_visit2 = geoSD(intensity[new_visit_blk =="Visit 2"], na.rm = TRUE),
             cv_visit2 =  (sd(intensity1[new_visit_blk == "Visit 2"],na.rm = TRUE)
                           /mean(intensity1[new_visit_blk =="Visit 2"],na.rm = TRUE))*100,
             
             #VISIT 3
             mean_visit3 = mean(intensity1[new_visit_blk =="Visit 3"], na.rm = TRUE),
             sd_visit3 = sd(intensity1[new_visit_blk =="Visit 3"], na.rm = TRUE),
             median_v3 = median(intensity1[new_visit_blk =="Visit 3"], na.rm = TRUE),
             gm_visit3 = geoMean(intensity[new_visit_blk =="Visit 3"], na.rm = TRUE),
             gsd_visit3 = geoSD(intensity[new_visit_blk =="Visit 3"], na.rm = TRUE),
             cv_visit3 =  (sd(intensity1[new_visit_blk == "Visit 3"],na.rm = TRUE)
                           /mean(intensity1[new_visit_blk =="Visit 3"],na.rm = TRUE))*100,
             
             #ALLVISITS
             mean_allvisit = mean(intensity1, na.rm = TRUE), 
             gm_allvisit = geoMean(intensity, na.rm = TRUE),
             gsd_allvisit = geoSD(intensity, na.rm = TRUE),
             cv_allvisit = (sd(intensity1, na.rm = TRUE)/mean(intensity1, na.rm = TRUE))*100, 
             # mean_MB = mean(intensity1[new_visit_blk == "mb"], na.rm = TRUE),
             # sd_MB = sd(intensity1[new_visit_blk == "mb"], na.rm = TRUE),
             # cv_MB = sd(intensity1[new_visit_blk == "mb"],na.rm = TRUE)/mean(intensity1[new_visit_blk =="mb"],na.rm = TRUE),
             #FIELDBLANK
             mean_FB = mean(intensity1[new_visit_blk=="FB"], na.rm = TRUE),
             log_mean_FB= mean(log2(intensity[new_visit_blk=="FB"])),
             
             mean_DIWB = mean(intensity1[new_visit_blk=="DIWB"], na.rm = TRUE),
             PvalV1_v3 = wilcox.test(intensity1[group=="Visit 1"], intensity1[group=="Visit 3"],exact = F)$p.value, 
             PvalV1_V2 = wilcox.test(intensity1[group=="Visit 1"], intensity1[group=="Visit 2"],exact = F)$p.value,
             # max_mb = max(intensity1[group == "mb"])
  )


# sigp_score = ifelse(conc_p < 0.1, 1, 0)
#if TRUE then mean at time 1 is greater than mean at time 2. 

v1vs3_df <- v1vs3_df%>%
  mutate(DFv1_min = ifelse(det_count_v1>=6 , 1, 0),
         mean_FB = ifelse(mean_FB == "NaN",NA, mean_FB),
         # log_mean_v1 = ifelse(log_mean_v1=="NaN", NA, mean_FB),
         percent_FBvV1= (mean_FB/mean_visit1)*100,
         logFBvV1 = (log_mean_FB/log_mean_v1)*100, 
         meanv1_gt_v3 = ifelse(mean_visit1 > mean_visit3, 1, 0),
         meanv1_gt_v2 = ifelse(mean_visit1 > mean_visit2, 1, 0),
         geom_diff1v3 = ifelse(gm_visit1 > gm_visit3, 1, 0),
         geom_diff1v2 = ifelse(gm_visit1 > gm_visit2, 1, 0),
         # sigp_score = ifelse(conc_p < 0.1, 1, 0), 
         # sig_diff_mean_v1_v3 = ifelse(sigp_score == 1 & v1_gt_v3 == 1, 1, 0 ),
         #detected in at least 3 in V1
         df1v3_difference = abs(df_v1 - df_v3),
         df1v2_difference = abs(df_v1 - df_v2)
         # dif10per2 = det_count_v2/ det_count_v1, 
         # dif10per3 = det_count_v3/ det_count_v1,
         # ave_int3 = mean_visit1/mean_visit3, 
         # ave_int2 = mean_visit1/mean_visit2
  )

#make a feature Id with a combo of mz and rt. We will merge on this later
v1vs3_df$mz<-round(v1vs3_df$mz,5)
v1vs3_df$rt <- round(v1vs3_df$rt, 1)
v1vs3_df$featureid<-as.character(paste(v1vs3_df$mz,v1vs3_df$rt,sep = "-"))



#05/21 we realized that we should use the field blank for filtering our samples rather than the MB. 
###Check for percent intensity in MB compared to Visit 1
#exclude blank intensity >10% compared to visit 1 intensity: excludemb = 1 
#exclude if DF is <20%
# table(v1vs3_df$cv_allvisit >50)
# 
# table(v1vs3_df$percent_MBallvisit <=10)
# table(v1vs3_df$percent_MBv1 <= 10)
#exclude features if the average intensity of the blank is >10% of the average intensity of all visits. 
# df10per<-filter(v1vs3_df,percent_MBallvisit<=10) # keep 1218 features after correcting for MB


#8/13/2020 we decided to use a filter 5x greater for features versus the blank. <20%
# df2<-filter(v1vs3_df,percent_MBallvisit<=20)
# df2 <- v1vs3_df


## merge with annotated chemicals and SSI database ##
#####################################################
df2 <- v1vs3_df



#import gss data to merge with 
d_annotation <- read.csv(file = p_annotation)

#make the feature id variable on which to merge. 
colnames(d_annotation)[which(colnames(d_annotation) == "time")] = "rt"
d_annotation$mz<-round(d_annotation$mz,5)
d_annotation$rt <- round(d_annotation$rt, 1)
d_annotation$featureid<-as.character(paste(d_annotation$mz,d_annotation$rt,sep = "-"))

# rename chemical_id to ID for the merge with the DB
names(d_annotation)[names(d_annotation)=="chemical_ID"] <- "ID"

#### ======= check
#keep relevant variables
#d_annotation = d_annotation[,c("featureid","chemical_ID","Confidence","score","mz","rt","MatchCategory","theoretical.mz", "delta_ppm","Name","Formula","MonoisotopicMass","Adduct")]

#rename to match previous datasets
#names(d_annotation)[names(d_annotation)=="chemical_ID"] <- "DTXSID"
#names(d_annotation)[names(d_annotation)=="Name"]<-"Compound"


#FILTER/Remove confidence <1. 

d_annotation <- d_annotation%>% 
  filter(Confidence >0)


###Merge with source and tox data
#import SSI features database with tox data: 

d_DB =  read.csv(p_DB, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# rename ID to 


#merge with SSI tox database by the DB_name and DTXSID, use the DTXSID and DB_name with labeling for the parent compound starting with MTX

df_annotation<-merge(d_annotation, d_DB, by= "ID", all= TRUE)
df_annotation$mz<-round(df_annotation$mz,5)
df_annotation$rt <- round(df_annotation$rt, 1)
df_annotation$featureid<-as.character(paste(df_annotation$mz,df_annotation$rt,sep = "-"))

#merge with annotation data
df4_neg <- merge(df2, df_annotation, by= "featureid", all.x = TRUE)

#identify features that matched and did not match the database
df4_neg$Adduct<-as.character(df4_neg$Adduct) ##We can do it at the adduct type level  

df4_neg[["Adduct"]][is.na(df4_neg[["Adduct"]])]<-"No match"

#count the number of features with a match
#unique features
df4_neg %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ###942features match and 11114 no match

# unique DTXSIDs
df4_neg %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(DTXSID)) ###1425 unique matches

# 
# df4_neg%>%
#   group_by(Adduct)%>%
#   summarize(count= length(DTXSID)) #298 total matches.

#keep those that match
df5_neg <- df4_neg%>%
  filter( Adduct == "M-H")

# df5_neg%>%
#   group_by(Adduct)%>%
#   summarise(count=n_distinct(featureid)) #178 unique features, #298 matches. 
# 
# df5_neg %>%
#   group_by(Adduct)%>%
#   summarise(count=n_distinct(featureid)) ###179 features matched
# 
# df5_neg %>%
#   group_by(Adduct)%>%
#   summarise(count=n_distinct(DTXSID)) ###245 unique matches




#filter by 
# Confidence >= 1
# diff_v1.3 ==TRUE
# test <- full_df%>%
#   filter(diff_v1.3 =="TRUE", 
#          Confidence>0
#          )

# write_csv(df5_neg, file = "data_files_csv/20210513_df5_Neg_fulllist.csv")

# make a variable with sums of MC chemicals, If fragmented in 2020 and if is ubuiquitous


#if need to convert variables to numeric rather than characters. 
# df5_pos <- df5_pos%>%
#   mutate(PFAS = ifelse(is.na(PFAS), 0, ifelse(PFAS==1, 1, 0)),
#          FRs = ifelse(is.na(FRs), 0, ifelse(FRs==1, 1, 0)),
#          # nitroPAH_bin = ifelse(nitroPAH_bin=="1", 1, 0),
#          exposurepred_bin = ifelse(is.na(exposurepred_bin), 0,ifelse(exposurepred_bin==1, 1, 0)),
#          ERactive_bin = ifelse(is.na(ERactive_bin), 0, ifelse(ERactive_bin==1, 1, 0)),
#          E2Up_bin = ifelse(is.na(E2Up_bin), 0, ifelse(E2Up_bin==1, 1, 0)),
#          P4Up_bin = ifelse(is.na(P4Up_bin), 0, ifelse(P4Up_bin==1, 1, 0)),
#          MC = ifelse(is.na(PFAS), 0, ifelse(MC==1, 1, 0)),
#          pesticidemammarytumors_bin = ifelse(is.na(PFAS), 0,ifelse(pesticidemammarytumors_bin==1, 1, 0)))


#create a sum of mammary carcinogens
df5_neg$MC_sum <- rowSums(df5_neg[,c("ERactive_bin", "E2Up_bin","P4Up_bin", "MC", "pesticidemammarytumors_bin")], na.rm = TRUE)



#find the number of lim_df5_neg that were sent for fragmentation. 
neg_frag <- read_csv("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2020Fragmentation list/list_DTXSID_fragmentation_nurses_FF_neg.csv")
# small_merge <- merge(lim_df5_neg, neg_frag, by = "DTXSID") #n= 70 match were sent for Fragmentation. 

#merge the matches between fragmentation data and the limited suspect list
confirmed_neg <- read.csv(file = "/Users/trowbridgej/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_neg/20210513_confirmedfeaturesfrag_neg.csv")

confirmed_dtxsid <- confirmed_neg$DTXSID[confirmed_neg$confirmed_match>=1] 


table(confirmed_neg$confirmed_match>=1) #13 features from the present list were confirmed in the prior fragmentation

#read the xls file that we annotated.
priority_df <- read_csv("data_files_csv/20210511Neg_features_FF_prioritized.csv")

dtxsidof_interest <- priority_df$DTXSID[priority_df$of_interest>=1]


#create new variables of interest. Filter on these to select chemicals for fragmentation
#make a new variable for if they were sent for fragmentation
#create a new variable for those that are acrylates

#frag_yn = was fragmented in the first round of analyses. 
#confirmed_frag = was confirmed in the prior round of fragmentation
#of_interest= selected as a chemical of interest from visual inspection of chemical names
#acrylate= is an acrylate from Alex's formula calculations
#is_ucsfdrug = was on the UCSF list of toxic drugs, we are excluding these chemicals for firefighters
#sig_dif1vs2or3= annotates 1 if the pvalue was <0.1 between visits. 
#ubiqui_MC = denotes ubiquitous mammary carcinogens. Mammary carcinogens that were found in most/all participants. MC_sum>1 and the significnt differnce was not significant. 

df5_neg <- df5_neg%>%
  mutate(new_id= ifelse(str_detect(DTXSID, "MTX"), parent_id, DTXSID),
         frag_yn = ifelse(DTXSID %in% neg_frag$DTXSID, "yes", "no"), 
         confirmed_frag= ifelse(DTXSID %in% confirmed_dtxsid, "confirmed", "notconf"),
         of_interest_prev = ifelse(DTXSID %in% dtxsidof_interest, "yes", NA),
         acrylate = ifelse(grepl("acrylate", name_original),1, 0), 
         is_ucsfdrug = ifelse(grepl('Drug', source), 1, 0),
         sig_dif1v2or3 = ifelse(PvalV1_v3 <= 0.1 | PvalV1_V2<= 0.1, 1, 0), 
         ubiqui_MC = ifelse(MC_sum >=1 & sig_dif1v2or3 ==0, 1, 0)
  )



#Add in chemical data based on structure from Alex Borrel

#read in negative data: 
nitro_df <- read_csv("data_files_csv/NEG_match_PAH.csv")
colnames(nitro_df)
#keep dtxsid and new columns
nitro_df <- nitro_df[ , c(1, 4:8)]


nitro_df$nitro_sum<- rowSums(nitro_df[,c(2:6)], na.rm = TRUE)

table(nitro_df$nitro_sum)

#merge with the nitro data 
test <- merge(df5_neg, nitro_df, by.x = "new_id", by.y = "DTXSID", all.x = TRUE) #new_id brings in the parent ID of the metabolite)

#import PAH and Drug data from Alex
drug_PAH_df <- read_csv("data_files_csv/NEG_drug_PAH.csv") 
colnames(drug_PAH_df)
drug_PAH_df <- drug_PAH_df[ , c(3, 6:13)]

table(drug_PAH_df$PAHLIST)
drug_PAH_df <- drug_PAH_df%>%
  mutate(drug_any = ifelse(SWGDRUG =="Y", 1, ifelse(STATINS =="Y", 1, ifelse(UOATARGPHARMA =="Y", 1, ifelse(VETDRUGS=="Y", 1, ifelse(DRUGBANK =="Y", 1, ifelse(ANTIBIOTICS =="Y", 1, ifelse(ZINC15PHARMA =="Y", 1, 0))))))))

drug_PAH_df <- drug_PAH_df%>%
  filter(DTXSID != "-")

drug_PAH_df <- distinct(drug_PAH_df)

table(drug_PAH_df$PAHLIST) #no PAH in the negative mode. 

# drug_PAH_df %>%
#   summarise(count=n_distinct(DTXSID)) ###334 unique matches. 
# 
df5_neg <- merge(test, drug_PAH_df,by.x = "new_id", by.y = "DTXSID", all.x = TRUE)




#remove those with FB filter <=50 and those detected in fewer than 3 participants of visit 1 
#filter the data by FB filter, minimum detection frequency, and significant difference between visit 1 and 2 or 3
# filter out those that are UCSF drugs in the variable source
lim_df5_neg <- df5_neg%>% 
  filter(percent_FBvV1<=50,  DFv1_min ==1)


#filter based on the following criteria: 

# Matches BC relevant chemicals, PFAS, FR, UV protectors, Nitro compound, Acrylate, PAH, Not a drug from the UCSF drug list. Unbiquitous mammary carcinogen, is a metabolite (based on )

#for sunscreens how to select in the database? what list to use? 

#is chemical of interest

lim_df5_test <- lim_df5_neg%>%
  rowwise(featureid)%>%
  mutate(of_interest_all = sum(c(PFAS, FRs, nitro_sum, acrylate, MC_sum, ubiqui_MC ), na.rm =TRUE), 
         of_interest_cat = ifelse(Drug_UCSF_PXYS!=1 & of_interest_all>0, 1, 0))


table(lim_df5_test$of_interest_cat)

lim_df5_test %>% 
  filter(ubiqui_MC==1)%>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ###19

selected_neg <- distinct(lim_df5_neg, featureid, .keep_all = TRUE)

#save the limited dataframe for analysis and feature selection
#use lim_df5_neg for making the boxplots
write.csv(lim_df5_test, "20210716Neg_features_FF_limited.csv")

# write.csv(selected_neg, "20210521_neg_Features_ff_selected.csv")
# write.csv(df5_neg, "20210501list_matched_NEGfeaturesFF.csv")


## ##Annotation ends here##

# futher data exploration below. 





#compare the number of features that are newly matched in this second annotation. 
#read in first annotation
old_annotation_sansMB <- read.csv("2020_NegAnnotation_noMBfilter.csv", stringsAsFactors = FALSE)

# merge with df5_neg by DTXSID
test_annotation <- merge(df5_neg, old_annotation_sansMB, by = "DTXSID", all.x = TRUE)

test_annotation <- test_annotation %>%
  mutate(new_feature = as.factor(ifelse(DTXSID == "-", 1, 0)))

# table(is.na(test_annotation$DTXSID))

# table(test_annotation$new_feature) #132 new features

test_annotation$DB_name[test_annotation$new_feature==1]

#first make original data usable by averaging r001 and r002
#  take average of .r001 and .r002

new_df <- neg_xset[, 1:2]


#make for loop
# dim(pos_xset)

positions = seq(3,130, by = 2)
sample_id = substr(colnames(neg_xset)[positions], 1, 7)

for (aux in seq_along(positions)) {
  pos1 = positions[aux]
  pos2 = pos1 + 1
  a = rowMeans(neg_xset[, pos1:pos2])
  new_df[ ,aux + 2] = a
}
colnames(new_df)[seq_along(positions) + 2] = sample_id

new_df$mz<-round(new_df$mz,5)
new_df$rt <- round(new_df$rt, 1)
new_df$featureid<-as.character(paste(new_df$mz,new_df$rt,sep = "-"))


# replace those with values rel intensity <1000 to 0 because they are noise. 

#x <- new_df[1:100, 10:15]
# apply(x,2, function(y) sum(y==0))

# for(j in 10:ncol(new_df)){
#   new_df[new_df[, j]<1000, j] <- 0
# }

#merge the individual intensities with the summary data

# df_plot <- merge(new_df, df5_neg , by = "featureid")

# to plot the limited variables
df_plot<- merge(new_df, lim_df5_neg, by = "featureid")

# head(new_df$featureid)
# head(df4_neg$featureid)
# plots describing features detected in positive  mode



# select_df <- new_df%>%
#   filter(featureid %in% df5_pos$featureid) #select those feature ids that are in the summarized and selected for MB correction N= 113 so only unique matches

# select_df <- merge(new_df, df5_pos, by  = "featureid") # to include all matches N= 182



dtsc <- select(dtsc, c("dtsc_id", "FF_id", "visit", "new_visit", "newFF_id"))
dtsc <- unique(dtsc)
df_plot_long <-  gather(df_plot, "sample_id", "intensity", BD01346:DIWB2.B )
df_plot_long$type<-substr(df_plot_long$sample_id,1,2)
df_plot <- merge(df_plot_long, dtsc, by.x= "sample_id", by.y = "dtsc_id", all.x = TRUE)

# df_plot <- merge(df_plot, annotation_pos, by = "featureid")

# df_plot <- df_plot%>%
#   mutate(new_visit_blk = ifelse(is.na(new_visit) & type =="MB", "mb", new_visit),
#          group= factor(new_visit_blk))


df_plot <- df_plot%>%
  mutate(new_visit_blk = ifelse(newFF_id == "F197"|newFF_id== "F223"|newFF_id== "F231"|newFF_id== "F232", "FB", new_visit),
         new_visit_blk = ifelse(is.na(new_visit) & type =="DI", "DIWB", new_visit_blk),
         group= factor(new_visit_blk), 
         class = ifelse(PFAS==1, "pfas", ifelse(FRs ==1, "fr", (ifelse(nitroPAH_bin==1, "nitroPAH", ifelse(ERactive_bin==1, "er_act", ifelse(E2Up_bin==1, "e2up", ifelse(P4Up_bin==1, "p4up", ifelse(MC==1, "mc", NA)))))))))

class_name <- df_plot%>% 
  group_by(featureid)%>%
  summarize(class = ifelse(PFAS==1, "pfas", ifelse(FRs ==1, "fr", (ifelse(nitroPAH_bin==1, "nitroPAH", ifelse(ERactive_bin==1, "er_act", ifelse(E2Up_bin==1, "e2up", ifelse(P4Up_bin==1, "p4up", ifelse(MC==1, "mc", NA)))))))))

# df_plot_long$newFF_id[is.na(df_plot_long$newFF_id)] <- "FBlank"

#add name to feature id. 

# make new COl with metabolite



# Feature candidates for confirmation based on if mammary carcinogen, Nitro-compounds, not drugs. 
#also meets FB filter and is significant between visits 1 and 2 or visits 1 and 3


count <- df_plot%>%
  group_by(featureid)%>%
  filter(of_interest=="yes")%>%
  count(featureid)


df_plot%>%
  group_by(featureid)%>%
  filter(of_interest=="yes")%>%
  # filter(mz.x < 200)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Selected features of interest (MC, Nitro compound, not a UCSF drug)", 
       caption = "(neg-ionization)") 

# ggsave(filename = "Boxplot_selected_features_neg.png")


# 1. All freatures detected in at least 3 participants in visit 1, has significant difference between V1 vs V2 or V3, and meets FB filter 
#* Meets FB filter, Min DF and has statistically significant difference between visit 1 vs Visit 2 OR visit 1 vs Visit 3. 
#* Wilcoxon rank-sum test, p-value<=0.10.1
#* n = 172


count <- df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  count(featureid)

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x < 150)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization) mz<150") 

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=150, mz.x<200)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=150, mz<200") 

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=200, mz.x<225)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=200, mz<225")

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=225, mz.x<250)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=225, mz<250")

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=250, mz.x<280)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=250, mz<280") 

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=280, mz.x<300)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=280, mz<300") 

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=300, mz.x<350)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=300, mz<350") 


df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=350, mz.x<400)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization); mz>=350, mz<400") 


df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,
         PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  filter(mz.x >=400)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter, detected in at least 3 samples of visit 1 and \n has different intensity between visit 1 and 3 or visit 1 and 2", 
       caption = "wilcoxon p-value <= 0.1 (neg-ionization), mz>=400") 

# 2. Features that are any mammary carcinogen AND detected in at least 3 participants in visit 1, has significant difference between V1 vs V2 or V3, and meets FB filter 
#* Is MC (using a sum of all variables associated with carcinogenicity)
#* Includes those with a FB<50% of Visit 1, detected in at least 3 of visit 1. MC_sum >=1 and is significant difference between visit 1 and Visit 2 or visit 3
#* n = 33



count_mc <- df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,  MC_sum>=1)%>%
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  count(featureid, name_original, DTXSID)

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50,  DFv1_min==1,  
         MC_sum>=1, PvalV1_V2 <=0.1 | PvalV1_v3<0.1 )%>%
  ggplot(aes(x = DTXSID, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = max(intensity) , group = featureid, label = DTXSID))+
  labs(title = "Mammary Carcinogen- meets FB filter, min-DF, & significant difference V1,V2 orV3", 
       caption = "(neg-ionization)") 

# 3 Features that are metabolites AND detected in at least 3 participants in visit 1, has significant difference between V1 vs V2 or V3, and meets FB filter 
#* is metabolite (is_metab==1)
#* Includes those with a FB<50% of Visit 1, detected in at least 3 of visit 1, is_metab==1, and is significant difference between visit 1 and Visit 2 or visit 3. 
#* n = 18


count_met <- df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,  is_metab==1)%>%
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  count(featureid, DTXSID)

df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  # filter(mz.x<240)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = DTXSID, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, meets FB filter and has dif V1, V2 or V3", 
       caption = "(neg-ionization)") 


# 4 Eposure predicted and meets FB filter, min detection frequency and significant p-value v1 v v2 or v3. 
#* n= 82 

count_exposurePred <- df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min ==1,  exposurepred_bin==1)%>%
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  count(featureid)

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  filter(mz<150)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  # geom_text(aes(x= featureid, y = 10, label = DTXSID, group = featureid), size= 2, angle= 90)+
  
  labs(title = "Exposure predicted", 
       caption = "(neg-ionization), mz<150") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  filter(mz>=150, mz<200)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted", 
       caption = "(neg-ionization), mz>=150, mz<200") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  filter(mz>=200, mz<250)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted", 
       caption = "(neg-ionization), mz>=200, mz<250")

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  filter(mz>=250, mz<300)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted", 
       caption = "(neg-ionization), mz>=250, mz<300") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(PvalV1_V2 <=0.1 | PvalV1_v3<=0.1)%>%
  filter(mz>=300)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted", 
       caption = "(neg-ionization), mz>=300") 


# Plots that meet one of the additional criteria (phase 2), meets FB filter, min detection in visit 1 (n= 3)
## Mammary Carcinogen (sum)

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50, DFv1_min==1,  
         MC_sum>=1 )%>%
  filter(mz.x<150 )%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter and is a mammary carcinogen", 
       caption = "(neg-ionization) mz<150") 

df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50,  DFv1_min==1,  
         MC_sum>=1)%>%
  filter(mz.x >=150, mz.x<200)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter and is a mammary carcinogen", 
       caption = "(neg-ionization)mz>=150, mz<200") 



df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50,  DFv1_min==1,  
         MC_sum>=1)%>%
  filter(mz.x >=200, mz.x<250)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter and is a mammary carcinogen", 
       caption = "(neg-ionization)mz>=200, mz<250") 


df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50,  DFv1_min==1,  
         MC_sum>=1)%>%
  filter(mz.x >=250, mz.x<300)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter and is a mammary carcinogen", 
       caption = "(neg-ionization) mz>=250, mz<300") 


df_plot%>%
  group_by(featureid)%>%
  filter(percent_FBvV1<=50,  DFv1_min==1,  
         MC_sum>=1)%>%
  filter(mz.x >=300)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  # geom_text(aes(x= featureid, y = 1000000, label = class))+
  labs(title = "Features that meet FB filter and is a mammary carcinogen", 
       caption = "(neg-ionization) mz>=300") 


## PFAS

df_plot %>% 
  filter(percent_FBvV1<=50 & PFAS==1)%>%  
  # filter(PvalV1_V2 <=0.1 | PvalV1_v3<0.1)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "PFAS and FB<50", 
       caption = "(neg-ionization)") 

## Flame retardants-- None plotted
df_plot %>% 
  filter(FRs==1)%>%  
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Flame Retardant chemicals", 
       caption = "BUT does not meet FB filter (neg-ionization)") 


# nitroPAH-- does not meet field blank filter
df_plot %>% 
  filter(nitroPAH_bin==1)%>%  
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "NItroPAH", 
       caption = "BUT does not meet FB filter (neg-ionization)") 

## Metabolite; fits FB Filter, and has a significant difference between V1, V2 or V3
df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%
  filter(mz.x<240)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz<240", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%  
  filter(mz.x>=240 & mz.x<300)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz>=300 & mz<400", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%  
  filter(mz.x>=300 & mz.x<343)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz>=300 & mz<343", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%  
  filter(mz.x>=343 & mz.x<400)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz>=343 & mz<400", 
       caption = "(neg-ionization)") 
df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%  
  filter(mz.x>=400 & mz.x<500)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz>=400 & mz<500", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & is_metab==1)%>%  
  filter(mz.x>=500)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is metabolite, mz>=500", 
       caption = "(neg-ionization)") 

## exposure predicted-
df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x<125)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz <125", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=125 & mz.x<150)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title ="Exposure predicted and mz>=125 & mz<150", 
       caption = "(neg-ionization)") 


df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=150 & mz.x<175)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=150 & mz>175", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=175 & mz.x<200)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=175 & mz<200", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=200& mz.x<250)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=200 & mz<250", 
       caption = "(neg-ionization)") 


df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=250 & mz.x<300)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=250 & mz<300", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=300 & mz.x<400 )%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=300 & mz<400", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & exposurepred_bin==1)%>%  
  filter(mz.x>=400)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  labs(title = "Exposure predicted and mz>=400", 
       caption = "(neg-ionization)") 


## ERactive
df_plot %>% 
  filter(percent_FBvV1<=50 & ERactive_bin==1)%>%  
  filter(mz.x<250)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  # theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+
  labs(title = "Is ERactive", 
       caption = "(neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & ERactive_bin==1)%>%  
  filter(mz.x>=250)%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  # theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+
  labs(title = "Is ERactive", 
       caption = "(neg-ionization)") 


## E2Up
df_plot %>% 
  filter(percent_FBvV1<=50 & E2Up_bin==1)%>%  
  
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  # facet_wrap(.~name_original, scales = "free" )+
  # theme(axis.title.x=element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank())+
  # ylim(NA, 400000)+
  labs(title = "Is E2Up bin", 
       caption = "(neg-ionization)") 


## P4Up
df_plot %>% 
  filter(percent_FBvV1<=50 & P4Up_bin==1)%>%  
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~name_original, scales = "free" )+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(trans= "log10")+
  labs(title = "Is P4Upbin", 
       caption = "(neg-ionization)") 

## MC (just MC_bin)
df_plot %>% 
  filter(percent_FBvV1<=50 & MC==1)%>%  
  filter(mz.x<280)%>%
  
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~featureid, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  # theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "Mammary Carcinogens" , 
       caption = "Is MC (neg-ionization)") 

df_plot %>% 
  filter(percent_FBvV1<=50 & MC==1)%>%  
  filter(mz.x>=280)%>%
  # filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = featureid, y = intensity,  color = new_visit_blk))+
  geom_boxplot()+
  # facet_wrap(.~featureid, scales = "free" )+
  scale_y_continuous(trans= "log10")+
  # theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "Mammary Carcinogens" , 
       caption = "Is MC (neg-ionization)") 

dfnomatch_1v2<- df_nomatch_pos %>%
  filter(conc_p_V2 <=0.05)%>%
  mutate(absMED = abs(med_diff1v2), 
         absDF = abs(df1v2_difference),
         compared = "1v2")%>%
  arrange(desc(absDF))
selected_absmed_1v2 <- med_nomatch_1v2[1:12, ]

#select those with largest absolute difference in MEDIAN between visit 1 and visit 3

med_nomatch_1v3<- df_nomatch_pos %>%
  filter(conc_p<=0.05)%>%
  mutate(absMED = abs(med_diff1v3), 
         absDF = abs(df1v3_difference), 
         compared = "1v3")%>%
  arrange(desc(absMED))

selected_absmed_1v3 <- med_nomatch_1v3[1:12, ]


#join the selected median lists. 
test2 <- rbind(selected_absmed_1v2, selected_absmed_1v3)
test2$featureid[duplicated(test2$featureid)]
test2 <- test2%>%distinct(test2$featureid, .keep_all = TRUE)
test2 <- test2%>%
  mutate(feat_dup = round(mz, 0))

nomatch_select_df <- new_df%>%
  filter(featureid %in% df_nomatch_pos$featureid) #select those feature ids that are in the summarized and selected for MB correction N= 113 so only unique matches

# select_df <- merge(new_df, df5_pos, by  = "featureid") # to include all matches N= 182

#merge again with dtsc. 
dtsc <- dtsc %>%
  mutate(time = ifelse(new_visit =="Visit 1", 1, ifelse(new_visit =="Visit 2", 2, ifelse(new_visit =="Visit 3", 3, NA))))

notmatch_plot_long <-  gather(nomatch_select_df, "sample_id", "intensity", BD01346:MB2..2. )


notmatch_plot_long$type<-substr(notmatch_plot_long$sample_id,1,2)

# df_plot_long$newFF_id[is.na(df_plot_long$newFF_id)] <- "FBlank"


# long_df2 <- long_df%>%
#   mutate(new_visit_blk = ifelse(is.na(new_visit) & type =="MB", "mb", new_visit),
#          group= factor(new_visit_blk))
# table(long_df2$new_visit_bl)


notmatch_plot <- merge(notmatch_plot_long, dtsc, by.x= "sample_id", by.y = "dtsc_id")


sigp <- df_nomatch_pos$featureid[df_nomatch_pos$conc_p<0.05]
sigp2 <- df_nomatch_pos$featureid[df_nomatch_pos$conc_p_V2<0.05]


# df_plot %>% 
#   # filter(featureid %in% cv_keep)%>% 
#   filter(featureid %in% gm_diff_12| featureid %in% df_diff)%>%
#   filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
#   ggplot(aes(x = visit, y = intensity, group= featureid, color = featureid))+
#   geom_line()+
#   facet_wrap(.~FF_id, scales = "free" )+
#   theme(legend.position = "none")+
#   scale_y_log10()+
#   labs(title = "Relative intensity of features different at time 1 vs time 2 for each study ID n=18", 
#        caption = "filtered by: GMv1 >GMv3 or DFv1 > DFv3") 



# ggsave("n20_POSFFfeatures_GMorDF_Diff.png", height = 7, width = 14 )

notmatch_plot %>% 
  # filter(featureid %in% selectedmed)%>%
  filter(featureid %in% sigp | featureid %in% sigp2)%>%
  filter(featureid %in% test2$featureid)%>%
  filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = new_visit, y = intensity, group= featureid, color = featureid))+
  geom_line()+
  facet_wrap(.~newFF_id, scales = "free" )+
  theme(legend.position = "none")+
  scale_y_log10()+
  labs(title = "Relative intensity of features with a significant difference at T1 vs T2 or T1 vs T3  for each study ID n= 20", 
       caption = "p value < 0.05; (neg-ionization)") 

# ggsave("n20_POSFFfeatures_sig_difference.png", height = 7, width = 14 )

#Plot by feature rather than id: 
notmatch_plot %>% 
  filter(featureid %in% sigp | featureid %in% sigp2)%>%
  filter(featureid %in% test2$featureid)%>%
  filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = new_visit, y = intensity, group= newFF_id, color = newFF_id))+
  geom_line()+
  facet_wrap(.~featureid, scales = "free" )+
  theme(legend.position = "none")+
  scale_y_log10()+
  labs(title = "Relative intensity of non-matched selected features", 
       caption = "p-value <=0.05. includes 21 features with largest absolute median difference between visit 1 and visit 2 or 3") 


notmatch_plot %>% 
  filter(featureid %in% sigp | featureid %in% sigp2)%>%
  filter(featureid %in% test2$featureid)%>%
  filter(FF_id != "F197" & FF_id != "F223" & FF_id != "F231"& FF_id != "F232")%>%
  ggplot(aes(x = new_visit, y = intensity))+
  geom_boxplot()+
  facet_wrap(.~featureid, scales = "free" )+
  theme(legend.position = "none")+
  scale_y_log10()+
  labs(title = "Relative intensity of non-matched selected features", 
       caption = "p-value <=0.05. includes 21 features with largest absolute median difference between visit 1 and visit 2 or 3") 

