####Processing non-targeted analysis data for nurses and OW#######
library(IPO)
library(xcms)
library(xMSannotator)
devtools::install_github("yufree/enviGCMS")
library(enviGCMS)
library(dplyr)
library(reshape2)
library(xMSanalyzer)
library(ggplot2)
library(mice)
library(stringr)
library(plyr)


##Find the best parameters for xcms - 
setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results")
peakpickingParameters<-getDefaultXcmsSetStartingParams("centWave")
path<-list.files("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/ESI_neg/QC/",full.names = T, recursive = T)
peakpickingParameters$ppm<-10
resultPeakpicking<-
  optimizeXcmsSet(files = path,
                  params = peakpickingParameters,
                  subdir = NULL)
optimizedXcmsSetObject<-resultPeakpicking$best_settings$xset
retcorGroupParameters<-getDefaultRetGroupStartingParams()
resultRetcorGroup<-
  optimizeRetGroup(xset = optimizedXcmsSetObject,
                   params = retcorGroupParameters,
                   subdir = NULL)
writeRScript(resultPeakpicking$best_settings$parameters,
             resultRetcorGroup$best_settings)
sessionInfo()

###get the data
getopqedata<-function(path,
                      index=F,
                      xsmethod="centWave",
                      peakwidth=c(12.8,54.5),
                      ppm=10,
                      noise=0,
                      snthresh=10,
                      mzdiff=0.00395,
                      prefilter=c(3,100),
                      mzCenterFun="wMean",
                      integrate=1,
                      fitgauss=FALSE,
                      verbose.columns=FALSE,
                      BPPARAM=BiocParallel::SnowParam(workers = 4),
                      rmethod="obiwarp",
                      plottype="none",
                      distFunc="cor_opt",
                      profStep=1,
                      center=4,
                      response=1,
                      gapInit=0.2,
                      gapExtend=2.4,
                      factorDiag=2,
                      factorGap=1,
                      localAlignment=0,
                      gmethod="density",
                      bw=22,
                      mzwid=0.015,
                      minfrac=0.1,
                      minsamp=1,
                      gmax=50,
                      ...){
  cdffiles<-list.files(path, recursive = TRUE, full.names = TRUE)
  if (index) {
    cdffiles<-cdffiles[index]
  }
  xset<-xcms::xcmsSet(
    cdffiles,
    method=xsmethod,
    snthresh=snthresh,
    mzdiff=mzdiff,
    BPPARAM = BPPARAM,
    peakwidth=peakwidth,
    ppm=ppm,
    noise=noise,
    prefilter=prefilter,
    mzCenterFun=mzCenterFun,
    integrate=integrate,
    fitgauss=fitgauss,
    verbose.columns=verbose.columns,
    ...
  )
  if (index & length(index)==1) {
    xset3<-xset
  } else {
    xset<-xcms::group(
      xset,
      method=gmethod,
      bw=bw,
      mzwid=mzwid,
      minfrac=minfrac,
      minsamp=minsamp,
      max=gmax
    )
    xset2<-xcms::retcor(
      xset,
      method=rmethod,
      plottype=plottype,
      distFunc=distFunc,
      profStep=profStep,
      center=center,
      response=response,
      gapInit=gapInit,
      gapExtend=gapExtend,
      factorDiag=factorDiag,
      factorGap=factorGap,
      localAlignment=localAlignment
    )
    
    #you need group the peaks again for this corrected data
    xset2<-xcms::group(
      xset2,
      method=gmethod,
      bw=bw,
      mzwid=mzwid,
      minfrac=minfrac,
      minsamp=minsamp,
      max=gmax
    )
    xset3<-xcms::fillPeaks(xset2, BPPARAM=BPPARAM)
  }
  return(xset3)
}

### get the data
path<-"D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/ESI_neg"
nurses_ow<-getopqedata(path)
xset3<-peakTable(nurses_ow)


write.csv(xset3, file="Nurses_OW_GSS__neg_ESI_preprocessed.csv")
df<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/Nurses_OW_GSS__neg_ESI_preprocessed.csv", header = T)


### Annotate ions with WWBC MS database (annotate log-transformed data not batch corrected. Batch effects will be evaluated later using batch # as covariate) 
##ESI pos
df<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/Nurses_OW_GSS__pos_ESI_preprocessed.csv", header = T)
df<-df[,c(2,5,18:289)]
colnames(df)[2]<-"time"
df<-unique(df)
df$QC2.2.r2_B5<-NULL

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
queryadductlist=c("M+H","M+NH4","M+Na","M+H-H2O","M+H-2H2O") #ESI pos 
queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H") #ESI neg
adduct_weights<-cbind.data.frame(Adduct=c("M+H","M-H"),Weight=c(5,5))

outloc<-"D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_annotation_03.21/ESI_pos"

WWBC_MS_library<- read.csv("C:/Dropbox (Silent Spring)/8. External Collaboration/building the wwbc database/DB_prepared/WWBC_MS_database_3.27.21_prepforanot.csv", sep="\t")
WWBC_MS_library<-WWBC_MS_library[,c("ID","name_cleaned","formula_cleaned","Molweight_cleaned")]
names(WWBC_MS_library)[names(WWBC_MS_library)=="name_cleaned"]<-"Name"
names(WWBC_MS_library)[names(WWBC_MS_library)=="formula_cleaned"]<-"Formula"
names(WWBC_MS_library)[names(WWBC_MS_library)=="Molweight_cleaned"]<-"MonoisotopicMass"

#library(tidyverse)
#WWBC_MS_library$MonoisotopicMass<-as.character(WWBC_MS_library$MonoisotopicMass)
#WWBC_MS_library$Formula<-as.character(WWBC_MS_library$Formula)
#WWBC_MS_library1<-WWBC_MS_library %>%
transform(MonoisotopicMass=strsplit(MonoisotopicMass,","),
          Formula=strsplit(Formula,",")) %>%
  unnest(MonoisotopicMass, Formula)
#WWBC_MS_library1<-WWBC_MS_library1[!duplicated(WWBC_MS_library1$DTXSID),]
#WWBC_MS_library$MonoisotopicMass<-as.numeric(WWBC_MS_library$MonoisotopicMass)
#WWBC_DB<-filter(WWBC_MS_library1, MonoisotopicMass>=80 & MonoisotopicMass<=1200)

WWBC_MS_library$MonoisotopicMass<-as.character(WWBC_MS_library$MonoisotopicMass)
WWBC_MS_library$MonoisotopicMass<-as.numeric(WWBC_MS_library$MonoisotopicMass)
WWBC_MS_library$ID<-as.character(WWBC_MS_library$ID)
WWBC_MS_library$Name<-as.character(WWBC_MS_library$Name)
WWBC_MS_library$Formula<-as.character(WWBC_MS_library$Formula)
WWBC_MS_library<-WWBC_MS_library[!is.na(WWBC_MS_library$MonoisotopicMass),]


#customDB<-as.data.frame(WWBC_MS_library)
#customDB$DTXSID<-as.character(customDB$DTXSID)
#customDB$Name<-as.character(customDB$Name)
#customDB$Formula<-as.character(customDB$Formula)

annotres<-multilevelannotation(df, max.mz.diff = max.mz.diff, max.rt.diff = max.rt.diff, 
                               cormethod = "pearson", 
                               num_nodes = num_nodes, queryadductlist = queryadductlist, mode = "pos",
                               outloc=outloc, db_name = db_name, adduct_weights = NA, num_sets = num_sets,
                               allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, customIDs = NA,
                               missing.value = NA, deepsplit = 2, 
                               networktype = "unsigned",
                               minclustsize = 10, module.merge.dissimilarity = 0.2, filter.by = c("M+H"),
                               redundancy_check = TRUE,
                               min_ions_perchem = 1, biofluid.location = NA, origin = NA,
                               status = status, boostIDs = NA, max_isp = 5,
                               MplusH.abundance.ratio.check = FALSE, 
                               customDB = WWBC_MS_library, HMDBselect = "union", mass_defect_window = 0.01,
                               mass_defect_mode = "pos",  pathwaycheckmode = "pm")






###Calculate detection frequency for each feature for each group including field blanks
df[df == 0] <- NA
df$count_n<-apply(df[,c(23:44,47:68,71:114,117:148)],1, function(z) sum(!is.na(z)))
df$count_n<-as.numeric(df$count_n)
df$DF_n<-(df$count_n/120)*100

df$count_ow<-apply(df[,c(169:236,241:252)],1, function(z) sum(!is.na(z)))
df$DF_ow<-(df$count_ow/80)*100

df$count_blank<-apply(df[,c(45,46,69,70,115,116,237:240,253,254)],1, function(z) sum(!is.na(z)))
df$DF_blank<-(df$count_blank/12)*100

df$DF_diff<-(df$DF_n-df$DF_ow)

####Reshape dataset (df) to have sample file names and intensity in two columns
df[is.na(df)] <- 0
df$time<-round(df$time,1)
df$mz<-round(df$mz,5)
dsn<-melt(df,id.vars = c("mz","time","count_n","DF_n","count_ow","DF_ow","count_blank","DF_blank","DF_diff"),
          variable.name = "filename",
          value.name = "intensity")

###Create id, group and replicate columns
dsn$id<-substr(dsn$filename,2,4)
dsn$id<-as.numeric(dsn$id)
dsn$replicate<-substr(dsn$filename,7,7)
dsn$group<-substr(dsn$filename,1,1)
dsn$group<-as.factor(dsn$group)
dsn$group<-as.character(dsn$group)
dsn$id<-as.character(dsn$id)

dsn$group[dsn$id %in% c("147","150","163","165","194","196","197",
                        "223","231","232")]<-"FB"
dsn$group<-as.factor(dsn$group)

###Keep only OW, nurses and field blanks samples
dsn<-filter(dsn, group=="N" | group=="W" | group=="FB")


##Calculate mean, SD and CV peak area for nurses and OW and MB
library(EnvStats)
cv<-function(x) 100*(sd(x)/mean(x))
dsn$featureid<-as.character(paste(dsn$mz,dsn$time,sep = "-"))
dsn_N<-filter(dsn, group=="N")

dsn_N_summary<-dsn_N %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity, na.rm = TRUE),gMean=geoMean(intensity, na.rm=TRUE))

#dsn_N_summary_t<-dsn_N %>%
group_by(featureid,id) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))


dsn_N<-merge(dsn_N,dsn_N_summary,by="featureid")
names(dsn_N)[names(dsn_N)=="mean"]<-"mean_N"
names(dsn_N)[names(dsn_N)=="sd"]<-"sd_N"
names(dsn_N)[names(dsn_N)=="cv"]<-"cv_N"
names(dsn_N)[names(dsn_N)=="median"]<-"median_N"
names(dsn_N)[names(dsn_N)=="gMean"]<-"gMean_N"

df_N<-dsn_N[,c(1:3,5,7,9,10,16:20)]
df_N<-df_N[!duplicated(df_N$featureid),]

dsn_OW<-filter(dsn, group=="W")

dsn_OW_summary<-dsn_OW %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity, na.rm = TRUE),gMean=geoMean(intensity,na.rm = TRUE))
dsn_OW<-merge(dsn_OW,dsn_OW_summary,by="featureid")
names(dsn_OW)[names(dsn_OW)=="mean"]<-"mean_OW"
names(dsn_OW)[names(dsn_OW)=="sd"]<-"sd_OW"
names(dsn_OW)[names(dsn_OW)=="cv"]<-"cv_OW"
names(dsn_OW)[names(dsn_OW)=="median"]<-"median_OW"
names(dsn_OW)[names(dsn_OW)=="gMean"]<-"gMean_OW"
df_OW<-dsn_OW[,c(1:3,16:20)]
df_OW<-df_OW[!duplicated(df_OW$featureid),]

dsn_FB<-filter(dsn, group=="FB")

dsn_FB_summary<-dsn_FB %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity, na.rm = TRUE),gMean=geoMean(intensity,na.rm = TRUE))
dsn_FB<-merge(dsn_FB,dsn_FB_summary,by="featureid")
names(dsn_FB)[names(dsn_FB)=="mean"]<-"mean_FB"
names(dsn_FB)[names(dsn_FB)=="sd"]<-"sd_FB"
names(dsn_FB)[names(dsn_FB)=="cv"]<-"cv_FB"
names(dsn_FB)[names(dsn_FB)=="median"]<-"median_FB"
names(dsn_FB)[names(dsn_FB)=="gMean"]<-"gMean_FB"
df_FB<-dsn_FB[,c(1:3,16:20)]
df_FB<-df_FB[!duplicated(df_FB$featureid),]

df1<-merge(df_N,df_OW,by=c("featureid","mz","time"))
df1<-merge(df1,df_FB,by=c("featureid","mz","time"))

####Merge with annotation
setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos")

annotation<- read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_annotation_03.21/ESI_pos/Stage5.csv")
annotation<-annotation[,c("chemical_ID","Confidence","score","mz","time","MatchCategory","theoretical.mz",
                          "delta_ppm","Name","Formula","MonoisotopicMass","Adduct")]
annotation$mz<-round(annotation$mz,5)
annotation<-filter(annotation, Confidence>0)

df1<-merge(df1,annotation, by=c("mz","time"), all.x=T)
df1$Adduct<-as.character(df1$Adduct)

df1[["Adduct"]][is.na(df1[["Adduct"]])]<-"No match"
df1<-filter(df1, Adduct=="M+H" | Adduct=="No match")

##Check number of unique features and unique matches
df1 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ###1491 features match and 20619 no match

df1 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(chemical_ID)) ###1573 unique matches


###Check for percent intensity in FB compared to mean N and OW
df1$mean_N_OW<-(df1$mean_N+df1$mean_OW)/2
df1$percent_FB<-(df1$mean_FB/df1$mean_N_OW)*100

###Merge with source and tox data
#WWBC_MS_library<- read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_list_annotation/WWBC_MS_database_7.9.20.csv", header = T)
WWBC_MS_library<-read.csv("C:/Dropbox (Silent Spring)/8. External Collaboration/building the wwbc database/DB_prepared_4.07.21/WWBC_MS_database_4.7.21_prepForAnotation.csv", sep="\t")

WWBC_MS_library<-WWBC_MS_library[,c(1:3,5,8:11,13:38)]
WWBC_MS_library$carcino<-ifelse(WWBC_MS_library$MC==1 | WWBC_MS_library$MGDev==1,1,0)

WWBC_MS_library %>%
  group_by(Drug_UCSF_PXYS) %>%
  summarise(count=n_distinct(ID)) ##306 unique suspects
WWBC_MS_library %>%
  group_by(Drug_most.comon.and.haz) %>%
  summarise(count=n_distinct(ID)) ##17 unique suspects

names(df1)[names(df1)=="chemical_ID"]<-"ID"
names(df1)[names(df1)=="Name"]<-"name_cleaned"
df2<-merge(df1,WWBC_MS_library,by="ID",all.x=T)

df2 %>%
  group_by(PFAS)%>%
  summarise(count=n_distinct(ID)) ##19 unique PFASs

df2 %>%
  group_by(FRs)%>%
  summarise(count=n_distinct(ID)) ##7 unique FRs

df2 %>%
  group_by(Drug_UCSF_PXYS)%>%
  summarise(count=n_distinct(ID)) ##66 unique UCSF drugs

df2 %>%
  group_by(Drug_most.comon.and.haz)%>%
  summarise(count=n_distinct(ID)) ##3 unique UCSF drugs

df2 %>%
  group_by(Disinfectant)%>%
  summarise(count=n_distinct(ID)) ##12 unique UCSF drugs

df2 %>%
  group_by(phthalate)%>%
  summarise(count=n_distinct(ID)) ##21 unique phthalates

df2$metabolite<-ifelse(df2$source=="metabolite",1,0)

df2 %>%
  group_by(metabolite)%>%
  summarise(count=n_distinct(ID)) ##107 unique phase II metabolites



######keep only all features  
df2$DF_high_N<-as.factor(ifelse(df2$DF_diff>=10,1,0))

df2 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ###1491 features match and 20619 no match

df2 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(ID)) ###1573 unique matches

#####Test for differences in intensity by group controlling for batch number
dsn<-filter(dsn, group=="N" | group=="W")

##Merge metadata files with batch number
batch1_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/batch1.txt",header = TRUE)
batch2_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/batch2.txt",header = TRUE)
batch3_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/batch3.txt",header = TRUE)
batch4_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/batch4.txt",header = TRUE)
batch5_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/batch5.txt",header = TRUE)
QC_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/QC.txt",header = TRUE)
MB_pos<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_pos/MB.txt",header = TRUE)


metadata_pos<-rbind(batch1_pos,batch2_pos,batch3_pos,batch4_pos,batch5_pos,QC_pos,MB_pos)
metadata_pos$sample_id<-str_replace(metadata_pos$filename, ".d|.mzXML", "") ##remove .d and .mzXML from filenames

metadata_pos$sample_id<-gsub("-",".",metadata_pos$sample_id,fixed=TRUE) ##Replace - by . in sample_id to match sample_id in dsn

colnames(dsn)[10]<-"sample_id"
dsn1<-merge(dsn,metadata_pos,by="sample_id")
dsn1$batch<-as.factor(dsn1$batch)

###Quantify associations between intensity~group+batch using lm for each feature
dsn1$log2_intensity<-log2(dsn1$intensity+1)
dsn1<-filter(dsn1, group=="N" | group=="W")
dsn1$nurses<-as.factor(ifelse(dsn1$group=="N",1,0))

###Perform glm
library(broom)
library(tidyr)
library(purrr)
glm.model_pos<-dsn1 %>%
  nest_by(featureid) %>%
  mutate(model=list(glm(log2_intensity~nurses+batch, data= data))) %>%
  ungroup() %>%
  transmute(featureid, modelcoef=map(model, tidy)) %>%
  unnest(modelcoef)


glm_models_nurses<-filter(glm.model_pos, term=="nurses1" & p.value<=0.1)

####Alternative without controlling for batch effect
wilcox<-dsn %>%
  group_by(featureid) %>%
  summarise(p=wilcox.test(intensity~group)$p.value)

wilcox.test(dsn$intensity~dsn$group)

######Merge with df2
df2<-merge(df2,glm_models_nurses,by="featureid",all.x = TRUE)
df2$higher_intensity_N<-as.factor(ifelse(df2$estimate>0,1,0))
df2$higher_intensity_N[is.na(df2$higher_intensity_N)]<-0

df2$estimate[is.na(df2$estimate)] <- 0
df2$metabolite<-as.factor(ifelse(df2$source=="metabolite",1,0))
df2$metabolite[is.na(df2$metabolite)] <- 0

###Add tox info for parent compounds of suspected metabolites
##New df with only tentatively identified metabolites
df3<-filter(df2,metabolite==1)
df3<-df3[,c(1:34,36,69:75)]

#Extract the DTXSID of the parent compound from the metabolite information 
df3<-df3 %>%
  mutate(parent_id=str_extract(DB_name, ".*(?=-)"))

##Remove DB_name and ID, and rename DTXSID
WWBC_MS_library_m<-WWBC_MS_library
WWBC_MS_library_m<-WWBC_MS_library_m[,c(2,4:35)]
colnames(WWBC_MS_library_m)[1]<-"parent_id"
df3<-merge(df3,WWBC_MS_library_m,by="parent_id",all.x = TRUE)

##Reformat df3 to match df2 structure
df3<-df3[,c(2:35,1,36,44:75,37:43)]

##Remove all metabolites from df2
df4<-filter(df2,metabolite==0)
##Rename DTXSID as parent_id
colnames(df4)[35]<-"parent_id"
##Rbind df3 and df4
df5<-rbind(df3,df4)


df2  %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(featureid)) #370 features matches 4,434 no matches

df2 %>%
  group_by(Adduct,higher_intensity_N)%>% #588 unique matches
  summarise(count=n_distinct(ID)) 


df2 %>%
  group_by(Drug_UCSF_PXYS,higher_intensity_N)%>% #25 unique drugs
  summarise(count=n_distinct(ID)) 

df2 %>%
  group_by(PFAS,higher_intensity_N)%>%
  summarise(count=n_distinct(ID))

df2 %>%
  group_by(FRs,higher_intensity_N)%>%
  summarise(count=n_distinct(DTXSID))

df2  %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(featureid)) 

df2 %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(DTXSID)) 

df2 %>%
  group_by(Drug,DF_high_N)%>%
  summarise(count=n_distinct(DTXSID)) 

df2 %>%
  group_by(PFAS,DF_high_N)%>%
  summarise(count=n_distinct(DTXSID))

df2 %>%
  group_by(FRs,DF_high_N)%>%
  summarise(count=n_distinct(DTXSID))



###Check for DF diff w/o blank correction
df1$DF_high_N<-as.factor(ifelse(df1$DF_diff>=10,1,0))

df1 %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(featureid)) ###175 features match and 2806 no match

df1 %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(DTXSID)) ###197 unique matches

###Check for intensity differences w/o blank correction
df2<-merge(df2,glm_models_nurses,by="featureid",all.x = TRUE)
df2$higher_intensity_N<-as.factor(ifelse(df2$estimate>0,1,0))
df2$higher_intensity_N[is.na(df2$higher_intensity_N)]<-0

df1  %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(featureid)) 

df1 %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(DTXSID)) 

#names(df1)[names(df1)=="chemical_ID"]<-"DTXSID"
#names(df1)[names(df1)=="Name"]<-"Compound"
df2<-merge(df1,WWBC_MS_library,by="DTXSID",all.x=T)

df2 %>%
  group_by(Drug)%>%
  summarise(count=n_distinct(featureid)) 

####Plot CV nurses and OW
df2_cv<-df2[,c("featureid","mz","time","cv_N","cv_OW","Adduct")]
df2_cv<-melt(df2_cv, id.vars = c("featureid", "mz","time","Adduct"))

df2_cv[df2_cv == "NaN"] <- NA

mu<-df2_cv %>%
  group_by(Adduct,variable) %>%
  summarise(mean=mean(value,na.rm=T))

ggplot(df2_cv, aes(x=value, color=variable))+
  geom_histogram(fill="white", position = "dodge", binwidth = 10)+
  geom_vline(data=mu, aes(xintercept=mean, color=variable),
             linetype="dashed")+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45))+
  facet_grid(vars(Adduct),scales = "free")+
  scale_x_continuous(breaks=seq(0,1200,40))

###Flag most hazardous drugs
hazard_drug<-data.frame(list(Compound=c("OXCARBAZEPINE","CLONAZEPAM","METRONIDAZOLE",
                                        "PREDNISONE","HALOPERIDOL","MYCOPHENOLIC ACID",
                                        "WARFARIN","TACROLIMUS","FUROSEMIDE","SPIRONOLACTONE",
                                        "TESTOSTERONE","TOPIRAMATE","VALGANCICLOVIR",
                                        "FOSPHENYTOIN","ZONISAMIDE"),drug_hazard=1))

df2$Compound<-toupper(df2$Compound)##Convert compound name from lower to upper case

##Merge hazard drug list with df2
df2<-merge(df2,hazard_drug,by="Compound",all.x = TRUE)


###Select top 10 unknown features with the highest DF diff (positive value)
no_match_flag<-filter(df5,Adduct=="No match")
no_match_flag<-arrange(no_match_flag,desc(DF_diff))%>%
  mutate(rank_DF=1:nrow(no_match_flag)) 

###Select top 10 unknown features with the lowest p values
no_match_flag<-no_match_flag %>%
  mutate(rank_pval=dense_rank(p.value))%>%
  arrange(rank_pval)

no_match_pval5<-filter(no_match_flag,p.value<=0.05)##5755 features p <= 0.05
#no_match_pval1<-filter(no_match_flag,p.value<=0.01)##873 features p <= 0.01

no_match_pval5$diff_N_OW<-no_match_pval5$median_N-no_match_pval5$median_OW ##Rank top 10 based on diff in median
no_match_pval5<-arrange(no_match_pval5,desc(diff_N_OW))%>%
  mutate(rank_pval_median=1:nrow(no_match_pval5))

no_match_pval5$ratio_N_OW<-no_match_pval5$median_N/no_match_pval5$median_OW ##Rank top 10 based on ratio
no_match_pval5<-arrange(no_match_pval5,desc(ratio_N_OW))%>%
  mutate(rank_pval_ratio=1:nrow(no_match_pval5))

write.csv(no_match_pval5,file="all_no_match_pos_FB.06.21.csv",row.names = FALSE)

no_match_MSMS_DF<-filter(no_match_pval5,rank_DF<=20)
write.csv(no_match_MSMS_DF,file="top10_no_match_DF_pos_FB.06.21.csv",row.names = FALSE)

no_match_MSMS_intensity<-filter(no_match_pval5,rank_pval_median<=20)
write.csv(no_match_MSMS_intensity,file="top10_no_match_diff_intensity_pos_FB.06.21.csv",row.names = FALSE)


match_flag<-filter(df5, Adduct=="M+H")

match_flag %>%
  summarise(n_distinct(featureid, na.rm = T)) ###1491 features


df2_hist<-df2[!duplicated(df2$featureid),]
ggplot(subset(df2_hist, DF_high_N %in% 1), aes(x=DF_diff)) + 
  geom_histogram(binwidth=1, color="black", fill="white")+
  theme_bw()+
  ggtitle("DF > 10% in nurses")

ggplot(subset(df2_hist, stat_diff %in% 1), aes(x=p)) + 
  geom_histogram(binwidth=0.005, color="black", fill="white")+
  theme_bw()+
  ggtitle("Significant differences in intensity of features between nurses and office workers \n (Wilcoxon test)")

match_flag$ESI<-"pos"
setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos")

write.csv(match_flag,file="List_matched_no_filter_pos_FB_06.21.csv")
#write.csv(no_match_MSMS,file="no_match_MSMS_pos.csv")

###Perform MWAS of occupational group to identify discriminating features
#library(MWASTools)
#mwas<-dsn[,c("filename","featureid","intensity","group")]
#mwas<-dcast(mwas,filename+group~featureid, value.var = "intensity") ###Not working fix this
#mwas$sample_type<-0
#mwas$group_mwas<-ifelse(mwas$Group=="FF",1,0)
#mwas[,3:8321]<-log10(mwas[3:8321])

#####Calculate mean intensity over 2 replicates for each id (Functions not working - Need to be fixed)
df3$mean<- df3 %>%
  group_by(mz,time,id) %>%
  summarize(Mean=mean(intensity,na.rm=TRUE))

##ESI neg
df_neg<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/Nurses_OW_GSS__neg_ESI_preprocessed.csv", header = T)
df_neg<-df_neg[,c(2,5,14:283)]
colnames(df_neg)[2]<-"time"
df_neg<-unique(df_neg)


outloc<-"D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_annotation_03.21/ESI_neg"

max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.7 #correlation threshold between adducts/isotopes
max_isp=5
mass_defect_window=0.01

num_nodes<-2   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption

db_name="Custom" #other options: KEGG, LipidMaps, T3DB, HMDB
status=NA 
num_sets<-3000 #number of sets into which the total number of database entries should be split into;

mode<-"neg" #ionization mode
queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H") #ESI neg
adduct_weights<-cbind.data.frame(Adduct=c("M+H","M-H"),Weight=c(5,5))

customDB<-as.data.frame(WWBC_DB)
customDB$DTXSID<-as.character(customDB$DTXSID)
customDB$Name<-as.character(customDB$Name)
customDB$Formula<-as.character(customDB$Formula)

annotres<-multilevelannotation(df_neg, max.mz.diff = max.mz.diff, max.rt.diff = max.rt.diff, 
                               cormethod = "pearson", 
                               num_nodes = num_nodes, queryadductlist = queryadductlist, mode = mode,
                               outloc=outloc, db_name = db_name, adduct_weights = NA, num_sets = 3000,
                               allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, customIDs = NA,
                               missing.value = NA, deepsplit = 2, 
                               networktype = "unsigned",
                               minclustsize = 10, module.merge.dissimilarity = 0.2, filter.by = c("M-H"),
                               redundancy_check = TRUE,
                               min_ions_perchem = 1, biofluid.location = NA, origin = NA,
                               status = status, boostIDs = NA, max_isp = 5,
                               MplusH.abundance.ratio.check = FALSE, 
                               customDB = WWBC_MS_library, HMDBselect = "union", mass_defect_window = 0.01,
                               mass_defect_mode = "pos",  pathwaycheckmode = "pm")

###Calculate detection frequency for each feature for each group
df_neg[df_neg == 0] <- NA
df_neg$count_n<-apply(df_neg[,c(23:44,47:68,71:114,117:148)],1, function(z) sum(!is.na(z)))
df_neg$count_n<-as.numeric(df_neg$count_n)
df_neg$DF_n<-(df_neg$count_n/120)*100

df_neg$count_ow<-apply(df_neg[,c(169:236,241:252)],1, function(z) sum(!is.na(z)))
df_neg$DF_ow<-(df_neg$count_ow/80)*100

df_neg$count_blank<-apply(df_neg[,c(45,46,69,70,115,116,237:240,253,254)],1, function(z) sum(!is.na(z)))
df_neg$DF_blank<-(df_neg$count_blank/12)*100

df_neg$DF_diff<-(df_neg$DF_n-df_neg$DF_ow)

####Reshape dataset (df_neg) to have sample file names and intensity in two columns
df_neg[is.na(df_neg)] <- 0
df_neg$time<-round(df_neg$time,1)
df_neg$mz<-round(df_neg$mz,5)
dsn_neg<-melt(df_neg,id.vars = c("mz","time","count_n","DF_n","count_ow","DF_ow","count_blank","DF_blank","DF_diff"),
              variable.name = "filename",
              value.name = "intensity")

###Create id, group and replicate columns
dsn_neg$id<-substr(dsn_neg$filename,2,4)
dsn_neg$replicate<-substr(dsn_neg$filename,7,7)
dsn_neg$group<-substr(dsn_neg$filename,1,1)

dsn_neg$group[dsn_neg$id %in% c("147","150","163","165","194","196","197",
                                "223","231","232")]<-"FB"
dsn_neg$group<-as.factor(dsn_neg$group)

###Keep OW, nurses and MB samples
dsn_neg<-filter(dsn_neg, group=="N" | group=="W" | group=="FB")


##Calculate mean, SD and CV peak area for nurses and OW and MB
cv<-function(x) 100*(sd(x)/mean(x))
dsn_neg$featureid<-as.character(paste(dsn_neg$mz,dsn_neg$time,sep = "-"))
dsn_neg_N<-filter(dsn_neg, group=="N")

dsn_neg_N_summary<-dsn_neg_N %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity,na.rm = TRUE),gMean=geoMean(intensity,na.rm = TRUE))

#dsn_N_summary_t<-dsn_N %>%
group_by(featureid,id) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))


dsn_neg_N<-merge(dsn_neg_N,dsn_neg_N_summary,by="featureid")
names(dsn_neg_N)[names(dsn_neg_N)=="mean"]<-"mean_N"
names(dsn_neg_N)[names(dsn_neg_N)=="sd"]<-"sd_N"
names(dsn_neg_N)[names(dsn_neg_N)=="cv"]<-"cv_N"
names(dsn_neg_N)[names(dsn_neg_N)=="median"]<-"median_N"
names(dsn_neg_N)[names(dsn_neg_N)=="gMean"]<-"gMean_N"
df_neg_N<-dsn_neg_N[,c(1:3,5,7,9,10,16:20)]
df_neg_N<-df_neg_N[!duplicated(dsn_neg_N$featureid),]

dsn_neg_OW<-filter(dsn_neg, group=="W")

dsn_neg_OW_summary<-dsn_neg_OW %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity,na.rm=TRUE),gMean=geoMean(intensity,na.rm = TRUE))
dsn_neg_OW<-merge(dsn_neg_OW,dsn_neg_OW_summary,by="featureid")
names(dsn_neg_OW)[names(dsn_neg_OW)=="mean"]<-"mean_OW"
names(dsn_neg_OW)[names(dsn_neg_OW)=="sd"]<-"sd_OW"
names(dsn_neg_OW)[names(dsn_neg_OW)=="cv"]<-"cv_OW"
names(dsn_neg_OW)[names(dsn_neg_OW)=="median"]<-"median_OW"
names(dsn_neg_OW)[names(dsn_neg_OW)=="gMean"]<-"gMean_OW"
df_neg_OW<-dsn_neg_OW[,c(1:3,16:20)]
df_neg_OW<-df_neg_OW[!duplicated(df_neg_OW$featureid),]

dsn_neg_FB<-filter(dsn_neg, group=="FB")

dsn_neg_FB_summary<-dsn_neg_FB %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity),median=median(intensity,na.rm = TRUE),gMean=geoMean(intensity,na.rm = TRUE))
dsn_neg_FB<-merge(dsn_neg_FB,dsn_neg_FB_summary,by="featureid")
names(dsn_neg_FB)[names(dsn_neg_FB)=="mean"]<-"mean_FB"
names(dsn_neg_FB)[names(dsn_neg_FB)=="sd"]<-"sd_FB"
names(dsn_neg_FB)[names(dsn_neg_FB)=="cv"]<-"cv_FB"
names(dsn_neg_FB)[names(dsn_neg_FB)=="median"]<-"median_FB"
names(dsn_neg_FB)[names(dsn_neg_FB)=="gMean"]<-"gMean_FB"
df_neg_FB<-dsn_neg_FB[,c(1:3,16:20)]
df_neg_FB<-df_neg_FB[!duplicated(df_neg_FB$featureid),]

df1_neg<-merge(df_neg_N,df_neg_OW,by=c("featureid","mz","time"))
df1_neg<-merge(df1_neg,df_neg_FB,by=c("featureid","mz","time"))

###Merge df1_neg with annotation output (stage 5) and detection frequency
setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_neg")

annotation_neg<- read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_annotation_03.21/ESI_neg/Stage5.csv")
annotation_neg<-annotation_neg[,c("chemical_ID","Confidence","score","mz","time","MatchCategory","theoretical.mz",
                                  "delta_ppm","Name","Formula","MonoisotopicMass","Adduct")]
annotation_neg$mz<-round(annotation_neg$mz,5)
annotation_neg<-filter(annotation_neg, Confidence>0)

df1_neg<-merge(df1_neg,annotation_neg, by=c("mz","time"), all.x=T)

df1_neg$Adduct<-as.character(df1_neg$Adduct)

df1_neg[["Adduct"]][is.na(df1_neg[["Adduct"]])]<-"No match"

df1_neg<-filter(df1_neg, Adduct=="M-H" | Adduct=="No match")

##Check number of unique features and unique matches
df1_neg %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ###1983 features match and 11012 no match

df1_neg %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(chemical_ID)) ###2928 unique matches


###Check for percent intensity in MB compared to mean N and OW
df1_neg$mean_N_OW<-(df1_neg$mean_N+df1_neg$mean_OW)/2
df1_neg$percent_FB<-(df1_neg$mean_FB/df1_neg$mean_N_OW)*100

###Merge with source and tox data
#WWBC_MS_library<- read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/WWBC_new_list_annotation/WWBC_MS_database_6.24.20.csv", header = T)
#WWBC_MS_library<-WWBC_MS_library[,c(2:3,9:30)]
WWBC_MS_library<-read.csv("C:/Dropbox (Silent Spring)/8. External Collaboration/building the wwbc database/DB_prepared_4.07.21/WWBC_MS_database_4.7.21_prepForAnotation.csv", sep="\t")
names(df1_neg)[names(df1_neg)=="chemical_ID"]<-"ID"
#names(df1_neg)[names(df1_neg)=="Name"]<-"Compound"
df2_neg<-merge(df1_neg,WWBC_MS_library,by="ID",all.x=T)

df2_neg %>%
  group_by(PFAS)%>%
  summarise(count=n_distinct(ID)) ###29 PFASs

df2_neg %>%
  group_by(FRs)%>%
  summarise(count=n_distinct(ID)) ###2 FRs

df2_neg %>%
  group_by(Drug_UCSF_PXYS)%>%
  summarise(count=n_distinct(ID)) ###40 UCSF drugs

df2_neg %>%
  group_by(Drug_most.comon.and.haz)%>%
  summarise(count=n_distinct(ID)) ###2 UCSF drugs


###blank correction 10% (10x)
df_percentMB_10<-filter(df2_neg,percent_MB<=10) ###10% blank correction 

df_percentMB_10 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ##103 unique features match and 2456 no match

df_percentMB_10 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(DTXSID)) ##179 unique matches

df_percentMB_10 %>%
  group_by(PFAS)%>%
  summarise(count=n_distinct(DTXSID)) ##2 unique PFASs

df_percentMB_10 %>%
  group_by(FRs)%>%
  summarise(count=n_distinct(DTXSID)) ##0 unique FRs

df_percentMB_10 %>%
  group_by(Drug)%>%
  summarise(count=n_distinct(DTXSID)) ##11 unique UCSF drugs

####Blank correction 20% (5x)
df_percentMB_20<-filter(df2_neg,percent_MB<=20) ###20% blank correction 

df_percentMB_20 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ##151 unique M-H and 3249 no match

df_percentMB_20 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(DTXSID)) ##246 unique matches

df_percentMB_20 %>%
  group_by(PFAS)%>%
  summarise(count=n_distinct(DTXSID)) ##9 unique PFASs

df_percentMB_20 %>%
  group_by(FRs)%>%
  summarise(count=n_distinct(DTXSID)) ##0 unique FRs

df_percentMB_20 %>%
  group_by(Drug)%>%
  summarise(count=n_distinct(DTXSID)) ##13 unique UCSF drugs

###blank correction 33% (3x)
df_percentMB_33<-filter(df2_neg,percent_MB<=33) ###33% blank correction 

df_percentMB_33 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(featureid)) ##234 unique M+H and 4196 no match

df_percentMB_33 %>%
  group_by(Adduct)%>%
  summarise(count=n_distinct(DTXSID)) ##401 unique matches

df_percentMB_33 %>%
  group_by(PFAS)%>%
  summarise(count=n_distinct(DTXSID)) ##13 unique PFASs

df_percentMB_33 %>%
  group_by(FRs)%>%
  summarise(count=n_distinct(DTXSID)) ##0 unique FRs

df_percentMB_33 %>%
  group_by(Drug)%>%
  summarise(count=n_distinct(DTXSID)) ##16 unique UCSF drugs


######Keep only features that pass blank correction 5x
#df2_neg<-df_percentMB_20
df2_neg$DF_high_N<-as.factor(ifelse(df2_neg$DF_diff>=10,1,0))

df2_neg %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(featureid)) ###75 features match and 1221 no match

df2_neg %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(ID)) ###140 unique matches

df2_neg %>%
  group_by(Drug_UCSF_PXYS,DF_high_N)%>%
  summarise(count=n_distinct(ID)) ###5 unique matches

#####Test for differences in intensity
dsn_neg<-filter(dsn_neg, group=="N" | group=="W")

##Merge metadata files with batch number
batch1_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/batch1_neg.txt",header = TRUE)
batch2_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/batch2_neg.txt",header = TRUE)
batch3_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/batch3_neg.txt",header = TRUE)
batch4_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/batch4_neg.txt",header = TRUE)
batch5_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/batch5_neg.txt",header = TRUE)
QC_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/QC_neg.txt",header = TRUE)
MB_neg<-read.table("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Data/Metadata/ESI_neg/MB_neg.txt",header = TRUE)


metadata_neg<-rbind(batch1_neg,batch2_neg,batch3_neg,batch4_neg,batch5_neg,QC_neg,MB_neg)
metadata_neg$sample_id<-str_replace(metadata_neg$filename, ".d|.mzXML", "") ##remove .d and .mzXML from filenames

metadata_neg$sample_id<-gsub("-",".",metadata_neg$sample_id,fixed=TRUE) ##Replace - by . in sample_id to match sample_id in dsn

colnames(dsn_neg)[10]<-"sample_id"
dsn1_neg<-merge(dsn_neg,metadata_neg,by="sample_id")
dsn1_neg$batch<-as.factor(dsn1_neg$batch)

###Quantify associations between intensity~group+batch using lm for each feature
dsn1_neg$log2_intensity<-log2(dsn1_neg$intensity+1)
dsn1_neg<-filter(dsn1_neg, group=="N" | group=="W")
dsn1_neg$nurses<-as.factor(ifelse(dsn1_neg$group=="N",1,0))

###Perform glm
library(broom)
library(tidyr)
library(purrr)
glm.model_neg<-dsn1_neg %>%
  nest_by(featureid) %>%
  mutate(model=list(glm(log2_intensity~nurses+batch, data= data))) %>%
  ungroup() %>%
  transmute(featureid, modelcoef=map(model, tidy)) %>%
  unnest(modelcoef)

glm_models_nurses_neg<-filter(glm.model_neg, term=="nurses1" & p.value<=0.1)

##Alternative to glm without controlling for batch number
wilcox_neg<-dsn_neg %>%
  group_by(featureid) %>%
  summarise(p=wilcox.test(intensity~group)$p.value)

df2_neg<-merge(df2_neg,glm_models_nurses_neg,by="featureid",all.x = TRUE)
df2_neg$higher_intensity_N<-as.factor(ifelse(df2_neg$estimate>0,1,0))
df2_neg$higher_intensity_N[is.na(df2_neg$higher_intensity_N)]<-0

df2_neg$metabolite<-as.factor(ifelse(df2_neg$source=="metabolite",1,0))
df2_neg$metabolite[is.na(df2_neg$metabolite)] <- 0

###Add tox info for parent compounds of suspected metabolites
##New df with only tentatively identified metabolites
df3_neg<-filter(df2_neg,metabolite==1)
df3_neg$Name<-NULL
df3_neg<-df3_neg[,c(1:33,68,35,69:75)]

#Extract the DTXSID of the parent compound from the metabolite information 
df3_neg<-df3_neg %>%
  mutate(parent_id=str_extract(DB_name, ".*(?=-)"))

##Remove DB_name and ID, and rename DTXSID
WWBC_MS_library_m<-WWBC_MS_library
WWBC_MS_library_m<-WWBC_MS_library_m[,c(2,4:35)]
colnames(WWBC_MS_library_m)[1]<-"parent_id"
df3_neg<-merge(df3_neg,WWBC_MS_library_m,by="parent_id",all.x = TRUE)

##Reformat df3 to match df2 structure
df3_neg<-df3_neg[,c(2:35,1,36,44:75,37:43)]

##Remove all metabolites from df2
df4_neg<-filter(df2_neg,metabolite==0)
df4_neg$Name<-NULL
df4_neg<-df4_neg[,c(1:33,68,34:67,69:75)]

##Rename DTXSID as parent_id
colnames(df4_neg)[35]<-"parent_id"
##Rbind df3_neg and df4_neg
df5_neg<-rbind(df3_neg,df4_neg)


df2_neg  %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(featureid)) 

df2_neg %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(ID)) 


df2_neg %>%
  group_by(Drug_UCSF_PXYS,higher_intensity_N)%>%
  summarise(count=n_distinct(ID)) 


###Check for DF diff w/o blank correction
df1_neg$DF_high_N<-as.factor(ifelse(df1_neg$DF_diff>=10,1,0))

df1_neg %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(featureid)) ###66 features match and 1237 no match

df1_neg %>%
  group_by(Adduct,DF_high_N)%>%
  summarise(count=n_distinct(chemical_ID)) ###121 unique matches

###Check for intensity differences w/o blank correction
df1_neg<-merge(df1_neg,glm_models_nurses_neg,by="featureid",all.x = TRUE)
df1_neg$higher_intensity_N<-as.factor(ifelse(df1_neg$estimate>0,1,0))
df1_neg$higher_intensity_N[is.na(df1_neg$higher_intensity_N)]<-0

df1_neg %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(featureid)) 

df1_neg %>%
  group_by(Adduct,higher_intensity_N)%>%
  summarise(count=n_distinct(chemical_ID)) 

names(df1_neg)[names(df1_neg)=="chemical_ID"]<-"DTXSID"
names(df1_neg)[names(df1_neg)=="Name"]<-"Compound"
df2_neg<-merge(df1_neg,WWBC_MS_library,by="DTXSID",all.x=T)

df2_neg %>%
  group_by(Drug)%>%
  summarise(count=n_distinct(DTXSID)) 



df2_neg_cv<-df2_neg[,c("featureid","mz","time","cv_N","cv_OW","Adduct")]
df2_neg_cv<-melt(df2_neg_cv, id.vars = c("featureid", "mz","time","Adduct"))

df2_neg_cv[df2_neg_cv == "NaN"] <- NA

mu_neg<-df2_neg_cv %>%
  group_by(Adduct,variable) %>%
  summarise(mean=mean(value,na.rm=T))

ggplot(df2_neg_cv, aes(x=value, color=variable))+
  geom_histogram(fill="white", position = "dodge", binwidth = 10)+
  geom_vline(data=mu_neg, aes(xintercept=mean, color=variable),
             linetype="dashed")+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45))+
  facet_grid(vars(Adduct),scales = "free")+
  scale_x_continuous(breaks = seq(0,1000,40))



###Select top 10 unknown features with the highest DF diff (positive value)
no_match_flag_neg<-filter(df5_neg,Adduct=="No match")
no_match_flag_neg<-arrange(no_match_flag_neg,desc(DF_diff))%>%
  mutate(rank_DF=1:nrow(no_match_flag_neg))

###Select top 10 unknown features with the lowest p values
no_match_flag_neg<-no_match_flag_neg %>%
  mutate(rank_pval=dense_rank(p.value))%>%
  arrange(rank_pval)

no_match_pval5_neg<-filter(no_match_flag_neg,p.value<=0.05)##886 features p <= 0.05
#no_match_pval1_neg<-filter(no_match_flag_neg,p.value<=0.01)##567 features p <= 0.01

no_match_pval5_neg$diff_N_OW<-no_match_pval5_neg$median_N-no_match_pval5_neg$median_OW
no_match_pval5_neg<-arrange(no_match_pval5_neg,desc(diff_N_OW))%>%
  mutate(rank_pval_median=1:nrow(no_match_pval5_neg))

no_match_pval5_neg$ratio_N_OW<-no_match_pval5_neg$median_N/no_match_pval5_neg$median_OW ##Rank top 10 based on ratio
no_match_pval5_neg<-arrange(no_match_pval5_neg,desc(ratio_N_OW))%>%
  mutate(rank_pval_ratio=1:nrow(no_match_pval5_neg))

write.csv(no_match_pval5_neg,file="all_no_match_neg_06.21.csv",row.names = FALSE)

no_match_MSMS_DF_neg<-filter(no_match_pval5_neg,rank_DF<=20)
write.csv(no_match_MSMS_DF_neg,file="top10_no_match_DF_neg_06.21.csv",row.names = FALSE)

no_match_MSMS_intensity_neg<-filter(no_match_pval5_neg,rank_pval_median<=20)
write.csv(no_match_MSMS_intensity_neg,file="top10_no_match_diff_intensity_neg_06.21.csv",row.names = FALSE)


match_flag_neg<-filter(df5_neg, Adduct=="M-H")

match_flag_neg %>%
  summarise(n_distinct(featureid, na.rm = T)) ###1036 unique features

match_flag_neg %>%
  group_by(MS_MS) %>%
  summarise(count=n_distinct(featureid)) ###23 unique features meet criteria and 78 not

match_flag_neg %>%
  group_by(MS_MS) %>%
  summarise(count=n_distinct(DTXSID))###37 matches meet criteria and 147 not


df2_neg_hist<-df2_neg[!duplicated(df2_neg$featureid),]
ggplot(subset(df2_neg_hist, DF_high_N %in% 1), aes(x=DF_diff)) + 
  geom_histogram(binwidth=1, color="black", fill="white")+
  theme_bw()+
  ggtitle("DF > 10% in nurses")

ggplot(subset(df2_neg_hist, stat_diff %in% 1), aes(x=p)) + 
  geom_histogram(binwidth=0.005, color="black", fill="white")+
  theme_bw()+
  ggtitle("Significant differences in intensity of features between nurses and office workers \n (Wilcoxon test)")

match_flag_neg$ESI<-"neg"

write.csv(match_flag_neg,file="List_matched_no_filter_neg_FB_06.21.csv")
write.csv(no_match_MSMS_neg,file="no_matches_MSMS_neg_05.21.csv")

###Look for chemicals detected in both ESI pos and ESI neg
df7<-rbind(match_flag,match_flag_neg)
df7$ESI<-as.factor(df7$ESI)

both_mode<-as.data.frame(table(df7$DTXSID,df7$ESI))
both_mode<-dcast(both_mode,Var1~Var2, value.var = "Freq")
names(both_mode)[names(both_mode)=="Var1"]<-"DTXSID"
both_mode$pos_neg<-as.factor(ifelse(both_mode$pos>0 & both_mode$neg>0,1,0))
df7<-merge(df7,both_mode,by="DTXSID")

df7_unique_DTXSID<-df7[!duplicated(df7$DTXSID),] ###363 unique chemical matches
summary(df7_unique_DTXSID$pos_neg) ###15 unique chemical matches found in both ESI pos and ESI neg

write.csv(df7,file="List_500features_MSMS_pos_neg.csv")

##########################################################################################
##FF data
ff_pos<-read.csv("C:/Users/VincentBessonneau/Downloads/list_featuresFFpos.csv", header=T)


##########################################################################################
##Prepare data for GGM analysis to identify associations between exogenous and endogenous molecules
##Add HMDB Annotation 
annotation_hmdb_pos<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/HMDB_annotation/ESI_pos/Stage5.csv", header = TRUE)
annotation_hmdb_pos<-annotation_hmdb_pos[,c("chemical_ID","Confidence","score","mz","time","MatchCategory","theoretical.mz",
                                            "delta_ppm","Name","Formula","MonoisotopicMass","Adduct")]
annotation_hmdb_pos$mz<-round(annotation_hmdb_pos$mz,5)
annotation_hmdb_pos<-filter(annotation_hmdb_pos, Confidence>0)
annotation_pos<-rbind(annotation_hmdb_pos,annotation)

df3<-merge(df2,annotation_pos,by=c("mz","time"),all.x = T)

names(df3)[names(df3)=="chemical_ID"]<-"DTXSID"
names(df3)[names(df3)=="Name"]<-"Compound"
df4<-merge(df3,WWBC_MS_library,by="DTXSID",all.x=T)

df4$Adduct<-as.character(df4$Adduct)
df4[["Adduct"]][is.na(df4[["Adduct"]])]<-"No match"

###Keep only annotated features
df5<-filter(df4, Adduct=="M+H")
df5<-filter(df5, mean_N>=1000)
df5$MS_MS<-as.factor(ifelse(df5$DF_high_N==1 | df5$higher_intensity_N==1 | df5$source=="Drug",1,0))

###Exclude features with DF<=50 in both nurses and ow samples
df6<-filter(df5,DF_n>=50 & DF_ow>=50)

##Remove duplicated featureid
df6<-df6[!duplicated(df6$featureid),]
df6<-df6[,2:4]

###Merge with samples info
df7<-merge(dsn1,df6,by=c("featureid","mz","time"))

#Compute Partial Correlations and Select Relevant Edges
#data matrix of metabolites
library(GeneNet)

df7<-df7[,c("featureid","sample_id","log2_intensity","batch")]

##Reshape df7 to have sample_id as rows and featureid with intensity as columns
df8<-dcast(df7,sample_id+batch~featureid, value.var = "log2_intensity")

data.x<-as.matrix(df8[,3:221])

#Perform GGM
pcor.dyn<-ggm.estimate.pcor(data.x,method="dynamic")

#estimating optimal shrinkage
met.edges<-network.test.edges(pcor.dyn,direct=TRUE,plot = FALSE)

#Extract network containing edges with prob >0.9 (i.e. FDR<0.1)
net<-extract.network(met.edges, cutoff.ggm = 0.9)
net$node1<-as.character(net$node1)
net$node2<-as.character(net$node2)

node.labels<-as.data.frame(colnames(data.x))
node.labels$node1<-as.character(1:nrow(node.labels))
node.labels$node2<-as.character(1:nrow(node.labels))
colnames(node.labels)[1]<-"featureid"

df9<-merge(net,node.labels,by="node1")
df9<-merge(df9,node.labels,by.x = "node2.x",by.y = "node2")

write.csv(df9,file="ggm_significant_edges_pos.csv")

##Prepare table for annotation on Cytoscape
node_annotation<-df5[,c(4,1,27,32,36)]

node_annotation[-1]<-apply(node_annotation[-1],2,as.character)
node_cytoscape<-aggregate(node_annotation[-1],by=list(node_annotation$featureid),c)

node_cytoscape$DTXSID<-unlist(lapply(node_cytoscape$DTXSID,paste,collapse=", "))
node_cytoscape$Compound<-unlist(lapply(node_cytoscape$Compound,paste,collapse=", "))
node_cytoscape$Confidence<-unlist(lapply(node_cytoscape$Confidence,paste,collapse=", "))
node_cytoscape$source<-unlist(lapply(node_cytoscape$source,paste,collapse=", "))

colnames(node_cytoscape)[1]<-"featureid"

node_feature<-df5[,c(4,2:3,5:6,8:9,12,20,26,60)]
node_feature<-node_feature[!duplicated(node_feature$featureid),]
node_feature<-merge(node_feature,node.labels,by="featureid")

node_cytoscape<-merge(node_cytoscape,node_feature,by="featureid")
write.csv(node_cytoscape,file="node_info_ggm_cytoscape.csv", row.names = FALSE)


#####Merge FF and nurses list features and de-duplicate
##ESI pos
FF_pos<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/list_matched_POSfeaturesFF29082020.csv", header = TRUE)

FF_pos %>%
  group_by(Adduct) %>%
  summarise(count=n_distinct(featureid)) ###113 unique features meet criteria 

FF_pos %>%
  group_by(Adduct) %>%
  summarise(count=n_distinct(DTXSID)) ###168 unique matches

FF_pos<-FF_pos[,c("featureid","mz","rt","DTXSID","Compound","Formula")]
names(FF_pos)[names(FF_pos)=="rt"]<-"time"
FF_pos$population<-"firefighters"


n_pos<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos/List_500features_MSMS_pos.csv", header=TRUE)
n_pos%>%
  group_by(Adduct)%>%
  summarise(n_distinct(DTXSID))

n_pos<-n_pos[,c("featureid","mz","time","DTXSID","Compound","Formula")]
n_pos$population<-"nurses"

wwbc_pos<-rbind(n_pos,FF_pos)
wwbc_pos$DTXSID<-as.character(wwbc_pos$DTXSID)


##Venn diagram
library(ggVennDiagram)
x<-list(
  Nurses=wwbc_pos%>%
    filter(population=="nurses") %>%
    select(DTXSID)%>%
    unlist(),
  Firefighters=wwbc_pos%>%
    filter(population=="firefighters")%>%
    select(DTXSID)%>%
    unlist()
)

ggVennDiagram(x)

wwbc_pos<-wwbc_pos[!duplicated(wwbc_pos$DTXSID),]

setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos")
no_match_n_pos_int<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos/top10_no_match_diff_intensity_pos.csv", header=TRUE)
no_match_n_pos_int$DTXSID<-"unknown"
no_match_n_pos_int$Compound<-"unknown"
no_match_n_pos_int$Formula<-"unknown"
no_match_n_pos_int<-no_match_n_pos_int[,c("featureid","mz","time","DTXSID","Compound","Formula")]
no_match_n_pos_int$population<-"nurses"

no_match_n_pos_DF<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_pos/top10_no_match_DF_pos.csv", header=TRUE)
no_match_n_pos_DF$DTXSID<-"unknown"
no_match_n_pos_DF$Compound<-"unknown"
no_match_n_pos_DF$Formula<-"unknown"
no_match_n_pos_DF<-no_match_n_pos_DF[,c("featureid","mz","time","DTXSID","Compound","Formula")]
no_match_n_pos_DF$population<-"nurses"

no_match_FF_pos_DF<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/2020_09POS_nomatch_absDF_select.csv", header=TRUE)
no_match_FF_pos_DF<-no_match_FF_pos_DF[,c("featureid","mz","rt","DTXSID","Compound","Formula")]
names(no_match_FF_pos_DF)[names(no_match_FF_pos_DF)=="rt"]<-"time"
no_match_FF_pos_DF$population<-"firefighters"
no_match_FF_pos_DF$DTXSID<-"unknown"
no_match_FF_pos_DF$Compound<-"unknown"
no_match_FF_pos_DF$Formula<-"unknown"

no_match_FF_pos_int<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/POS_nomatch_abs_median_difference_selected.csv", header=TRUE)
no_match_FF_pos_int<-no_match_FF_pos_int[,c("featureid","mz","rt","DTXSID","Compound","Formula")]
names(no_match_FF_pos_int)[names(no_match_FF_pos_int)=="rt"]<-"time"
no_match_FF_pos_int$population<-"firefighters"
no_match_FF_pos_int$DTXSID<-"unknown"
no_match_FF_pos_int$Compound<-"unknown"
no_match_FF_pos_int$Formula<-"unknown"

wwbc_pos<-rbind(no_match_FF_pos_DF,no_match_n_pos_int,wwbc_pos,no_match_FF_pos_int,no_match_n_pos_DF)

write.csv(wwbc_pos,file="list_DTXSID_fragmentation_nurses_FF_pos.csv", row.names = FALSE)

wwbc_pos<-wwbc_pos[!duplicated(wwbc_pos$featureid),]
write.csv(wwbc_pos,file="list_mz_fragmentation_nurses_FF_pos.csv", row.names = FALSE)


##ESI neg
FF_neg<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/list_matched_NEGfeaturesFF28082020.csv", header = TRUE)
FF_neg<-FF_neg[,c("featureid","mz.x","rt.x","DTXSID","Compound","Formula")]
names(FF_neg)[names(FF_neg)=="mz.x"]<-"mz"
names(FF_neg)[names(FF_neg)=="rt.x"]<-"time"
FF_neg$population<-"firefighters"

n_neg<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_neg/List_500features_MSMS_neg.csv", header=TRUE)
n_neg%>%
  group_by(Adduct)%>%
  summarise(n_distinct(DTXSID))

n_neg<-n_neg[,c("featureid","mz","time","DTXSID","Compound","Formula")]
n_neg$population<-"nurses"

wwbc_neg<-rbind(n_neg,FF_neg)
wwbc_neg$DTXSID<-as.character(wwbc_neg$DTXSID)

##Venn diagram
library(ggVennDiagram)
x<-list(
  Nurses=wwbc_neg%>%
    filter(population=="nurses") %>%
    select(DTXSID)%>%
    unlist(),
  Firefighters=wwbc_neg%>%
    filter(population=="firefighters")%>%
    select(DTXSID)%>%
    unlist()
)

ggVennDiagram(x)

wwbc_neg<-wwbc_neg[!duplicated(wwbc_neg$DTXSID),]

no_match_n_neg_int<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_neg/top10_no_match_diff_intensity_neg.csv", header=TRUE)
no_match_n_neg_int$DTXSID<-"unknown"
no_match_n_neg_int$Compound<-"unknown"
no_match_n_neg_int$Formula<-"unknown"
no_match_n_neg_int<-no_match_n_neg_int[,c("featureid","mz","time","DTXSID","Compound","Formula")]
no_match_n_neg_int$population<-"nurses"

no_match_n_neg_DF<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_neg/top10_no_match_DF_neg.csv", header=TRUE)
no_match_n_neg_DF$DTXSID<-"unknown"
no_match_n_neg_DF$Compound<-"unknown"
no_match_n_neg_DF$Formula<-"unknown"
no_match_n_neg_DF<-no_match_n_neg_DF[,c("featureid","mz","time","DTXSID","Compound","Formula")]
no_match_n_neg_DF$population<-"nurses"

no_match_FF_neg_DF<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/2020_09_16NEGabs_DF_difference_selected.csv", header=TRUE)
no_match_FF_neg_DF<-no_match_FF_neg_DF[,c("featureid","mz.x","rt.x","DTXSID","Compound","Formula")]
names(no_match_FF_neg_DF)[names(no_match_FF_neg_DF)=="mz.x"]<-"mz"
names(no_match_FF_neg_DF)[names(no_match_FF_neg_DF)=="rt.x"]<-"time"
no_match_FF_neg_DF$population<-"firefighters"
no_match_FF_neg_DF$DTXSID<-"unknown"
no_match_FF_neg_DF$Compound<-"unknown"
no_match_FF_neg_DF$Formula<-"unknown"

no_match_FF_neg_int<-read.csv("D:/MYBOOK/Silent Spring Institute/WWBC/FF/Data/NEG_nomatch_abs_median_difference_selected.csv", header=TRUE)
no_match_FF_neg_int<-no_match_FF_neg_int[,c("featureid","mz.x","rt.x","DTXSID","Compound","Formula")]
names(no_match_FF_neg_int)[names(no_match_FF_neg_int)=="mz.x"]<-"mz"
names(no_match_FF_neg_int)[names(no_match_FF_neg_int)=="rt.x"]<-"time"
no_match_FF_neg_int$population<-"firefighters"
no_match_FF_neg_int$DTXSID<-"unknown"
no_match_FF_neg_int$Compound<-"unknown"
no_match_FF_neg_int$Formula<-"unknown"

wwbc_neg<-rbind(no_match_FF_neg_DF,no_match_n_neg_int,wwbc_neg,no_match_FF_neg_int,no_match_n_neg_DF)


setwd("D:/MYBOOK/Silent Spring Institute/WWBC/Nurses/GSS_processing/Results/MSMS_selection/ESI_neg")
write.csv(wwbc_neg,file="list_DTXSID_fragmentation_nurses_FF_neg.csv", row.names = FALSE)

wwbc_neg<-wwbc_neg[!duplicated(wwbc_neg$featureid),]
write.csv(wwbc_neg,file="list_mz_fragmentation_nurses_FF_neg.csv", row.names = FALSE)

