install.packages("remotes")


install.packages()


BiocManager::install(c('MASS', 'rgl', 'mzR' , 'splines', 'doParallel', 'foreach', 'iterators', 'snow','ROCR', 'ROCS','e1071','randomForest','gbm'))
BiocManager::install("GO.db")
BiocManager::install("xcms")
install.packages(c('BiocManager','data.table','digest'))
remotes::install_github('omegahat/XMLSchema')
remotes::install_github('cran/SSOAP')
BiocManager::install(c("SSOAP","KEGGREST","pcaMethods","Rdisop","GO.db","matrixStats","WGCNA"))
devtools::install_github("yufree/xMSannotator")
remotes::install_github("jaspershen/MSannotator")
BiocManager::install("sva")
remotes::install_github("yufree/xMSanalyzer", dependencies = TRUE)

devtools::install_github("yufree/enviGCMS")

#####Processing non-targeted analysis data for FF#######
library(IPO)
library(xcms)
# library(xMSannotator)
library(enviGCMS)
library(dplyr)
library(reshape2)


library(ggplot2)
library(mice)


##Find the best parameters for xcms using the R-package IPO
#setwd("E:/MYBOOK/Silent Spring Institute/WWBC/Nurses/ESI_pos/Data preprocessing")
setwd("c:/Users/jesstro/Box")
getwd()

peakpickingParameters<-getDefaultXcmsSetStartingParams("centWave")
path<-list.files("C:/Users/jesstro/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_neg/preprocessing_neg/Neg_ionization_rawdata/QC_neg/",full.names = T, recursive = T)##Folder that contain QC raw data (.mzXML)
peakpickingParameters$ppm<-10
resultPeakpicking<-
  optimizeXcmsSet(files = path,
                  params = peakpickingParameters,
                  subdir = NULL)
optimizedXcmsSetObject<-resultPeakpicking$best_settings$xset
retcorGroupParameters<-getDefaultRetGroupStartingParams()


#-----------this command was to save what I had already optimized and had to close the session.-------------

save(optimizedXcmsSetObject,peakpickingParameters, resultPeakpicking, retcorGroupParameters, file = "SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_neg/preprocessing_neg/20210423neg_optimizedpeak_picking.RData" )



#----------Start here next time you open R----
load("SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ neg_optimizedpeak_picking.RData")

resultRetcorGroup<-
  optimizeRetGroup(xset = optimizedXcmsSetObject,
                   params = retcorGroupParameters,
                   subdir = NULL, nSlaves = 2)

# save(resultRetcorGroup, file = "neg_resultretcorgroup.RData" )

writeRScript(resultPeakpicking$best_settings$parameters,
             resultRetcorGroup$best_settings)
sessionInfo()

###get the data
getopqedata<-function(path,
                      index=F,
                      xsmethod="centWave",
                      peakwidth=c(13.6,69.5),
                      ppm=10,
                      noise=0,
                      snthresh=10,
                      mzdiff=9.99999999999994e-05,
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
                      center=1,
                      response=1,
                      gapInit=0.928,
                      gapExtend=1.968,
                      factorDiag=2,
                      factorGap=1,
                      localAlignment=0,
                      gmethod="density",
                      bw=0.8799999999,
                      mzwid=0.0265,
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

### In case the function "getopqedata" doesn't work
#xset1<-group(xset, bw=0.8799999999,mzwid=0.0265,minfrac=0.1,minsamp=1)
#xset2<-retcor(xset1, method="obiwarp",plottype="none",distFunc="cor_opt",profStep=1,center=1,response=1,gapInit=0.928,gapExtend=1.968,factorDiag=2,factorGap=1,localAlignment=0)
#xset3<-group(xset2, bw=0.8799999999,mzwid=0.0265,minfrac=0.1,minsamp=1,max=50)
#xset4<-fillPeaks(xset3)
#FF<-peakTable(xset4)


### get the data

path<-"c:/Users/jesstro/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_neg/preprocessing_neg/Neg_ionization_rawdata/SerumQCBK_samples_neg/" ###Folder that contain all raw data files, including FF, QCs and blanks (.mzXML)
FF<-getopqedata(path)
xset3<-peakTable(FF)

write.csv(xset3, "20210423neg_xset.csv")

getwd(
)
# xset3 <- read.csv("Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/Code/neg_xset15072020.csv", header =  TRUE)


### Code Ends Here. The following steps are replicated in file 3_FF_NEG_NTA...

####Import FF id coding DTSC
dtsc<-read.csv("C:/Users/jesstro/Desktop/FF_NTA_processing/FF_id_visit.csv", header = T)

###Extract FF id and visit
dtsc$dtsc_id<-substr(dtsc$id,1,7)
dtsc$dtsc_run <- substr(dtsc$id, 30,33)
dtsc$FF_id<-substr(dtsc$id,9,12)
dtsc$visit<-substr(dtsc$id,14,20)

dtsc$dtsc_fid <- as.character(paste(dtsc$dtsc_id,dtsc$dtsc_run,sep = "."))

#two participants had visit #4,5 & 6 will give a new number later

dtsc <- dtsc%>%
  mutate(new_visit= ifelse(FF_id == "F227" & visit == "Visit 3", "Visit 1", 
                           ifelse(FF_id =="F227" & visit == "Visit 4", "Visit 2", 
                                  ifelse(FF_id =="F227" &visit == "Visit 5", "Visit 3", 
                                         ifelse( FF_id =="F216" & visit =="Visit 4", "Visit 1", 
                                                 ifelse(FF_id =="F216" & visit == "Visit 5", "Visit 2", 
                                                        ifelse(FF_id =="F216" & visit == "Visit 6", "Visit 3", visit)))))), 
         newFF_id = ifelse(FF_id == "F227" & visit == "Visit 3", "F502", 
                           ifelse(FF_id == "F227" & visit == "Visit 4", "F502", 
                                  ifelse(FF_id == "F227" & visit == "Visit 5", "F502",
                                         ifelse(FF_id == "F216" & visit == "Visit 4", "F501",
                                                ifelse(FF_id == "F216" & visit == "Visit 5", "F501",
                                                       ifelse(FF_id == "F216" & visit == "Visit 6", "F501", FF_id)))))))




# table(dtsc$newFF_id)


# #remove blank sample IDs
# dtsc <- dtsc[-(dtsc$newFF_id=="F197") , ]
# dtsc <- dtsc[-(dtsc$newFF_id=="F223") , ]
# dtsc <- dtsc[-(dtsc$newFF_id=="F231") , ]



###Check for missing data and exclude if more than 50% in visit 1. 
pMiss<-function(x){sum(is.na(x))/length(x)*100}
dataA<-xset3[,c(2,5,10:141)]
dataA[dataA==0]<-NA
names(dataA)[names(dataA)=="rt"]<-"time"
dataA$mz<-round(dataA$mz,5)
dataA$time<-round(dataA$time,2)
dataA$featureid<-as.character(paste(dataA$mz,dataA$time,sep = "-"))


dataB<-dataA[,-135]
row.names(dataB)<-dataA[,135]
miss_samples<-data.frame(apply(dataB,2,pMiss))
miss_features<-data.frame(apply(dataB,1,pMiss))


##subset by Visit
# colnames(new_df)

v1_df <- dtsc %>%
  filter(new_visit =="Visit 1")

v3_df <- dtsc %>%
  filter(new_visit =="Visit 3")
v1 <- v1_df$dtsc_id

id_vars <- dataA[,c(1,2)]

visit1_df_1 <- dataA[ , v1_df$dtsc_fid]
visit1_df <- cbind(id_vars,visit1_df_1)
visit1_df$feature <- as.character(paste(visit1_df$mz,visit1_df$time,sep = "-"))
dataB_visit1<-visit1_df[,-51]
row.names(dataB_visit1)<-visit1_df[,51]

visit3_df_1 <- dataA[, v3_df$dtsc_fid]
visit3_df <- cbind(id_vars,visit3_df_1)
visit3_df$feature <- as.character(paste(visit3_df$mz,visit3_df$time,sep = "-"))
dataB_visit3<-visit3_df[,-41]
row.names(dataB_visit3)<-visit3_df[,41]

miss_features_visit1<-data.frame(apply(dataB_visit1,1,pMiss))
miss_features_visit1<-tibble::rownames_to_column(miss_features_visit1, "featureid")
miss_features_visit3<-data.frame(apply(dataB_visit3,1,pMiss))
miss_features_visit3<-tibble::rownames_to_column(miss_features_visit3, "featureid")

miss_features_class<-merge(miss_features_visit1,miss_features_visit3,by="featureid",all=TRUE)
miss_features_class$exclude<-ifelse(miss_features_class$apply.dataB_visit1..1..pMiss.>=50 ,1,0)
miss_features_class$exclude<-as.factor(miss_features_class$exclude)
summary(miss_features_class$exclude)

dataA<-merge(dataA,miss_features_class,by="featureid")
data_clean<-filter(dataA, exclude==0)


###Calculate detection frequency for each feature for each group
data_clean$count<-apply(data_clean[,4:127],1, function(z) sum(!is.na(z))) #All. 
data_clean$DF<-(data_clean$count/127)*100


colnames(data_clean)
data_clean$count_n<-apply(data_clean[,v1_df$dtsc_fid],1, function(z) sum(!is.na(z))) #visit 1
data_clean$DF_v1<-(data_clean$count_n/48)*100

data_clean$count_ow<-apply(data_clean[,v3_df$dtsc_fid],1, function(z) sum(!is.na(z))) #visit 3
data_clean$DF_v3<-(data_clean$count_ow/38)*100

see<-data_clean[,c("featureid","DF","DF_v1","DF_v3","exclude")]

colnames(data_clean)
##Replace missing values (0 by 10) and log2 transfromation 
data<-data_clean[,c(4:127)]
min<-data.frame(apply(data,1,min, na.rm=TRUE)) ###find min intensity for each features.

data_clean[is.na(data_clean)]<-63
data_clean[,4:127]<-log(data_clean[4:127],2)
check<-data_clean[,c(1,4:127)] ##check normality of distribution 
tcheck<-setNames(data.frame(t(check[,-1])), check[,1])


summary(tcheck$`1004.67289-1042.2`)
hist(tcheck$`1005.67594-1042.15`)


write.csv(data_clean, "Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/Code/ff_gss_neg_esi_preprocessed.csv")
###Look for batch effects using PCA plots
library(FactoMineR)
tcheck<-tibble::rownames_to_column(tcheck, "rowid")
tcheck$rowid<-lapply(tcheck$rowid, function(x) {
  gsub("\\.", "\\-", x)
})
tcheck$rowid<-as.character(tcheck$rowid)
# metadata$rowid<-lapply(metadata$rowid, function(x) {
#   gsub("\\.", "\\-", x)
# })
# metadata$rowid<-as.character(metadata$rowid)

pca<-tcheck
pca$batch<-as.factor(pca$batch)
pca$group<-as.factor(pca$group)
pca1<-pca[,-c(1:3)]
row.names(pca1)<-pca[,1]
res.pca<-PCA(pca1,quali.sup = 1)
plot(res.pca, cex=0.7)
plot(res.pca, cex=0.7,axes=3:4)

aa<-cbind.data.frame(pca1[,1], res.pca$ind$coord[,1:2])
bb<-coord.ellipse(aa, bary = TRUE)
plot.PCA(res.pca, axes=1:2, ellipse = bb, graph.type = "ggplot", ggoptions = list(size=2))
aa<-cbind.data.frame(pca1[,1], res.pca$ind$coord[,3:4])
bb<-coord.ellipse(aa, bary = TRUE)
plot.PCA(res.pca, axes=3:4, ellipse = bb, graph.type = "ggplot", ggoptions = list(size=2))

###PCA using ggplot2
pca1$group<-NULL
pca1_pca<-prcomp(pca1)
plot(pca1_pca$x[,1], pca1_pca$x[,2])
pca1_pca_out<-as.data.frame(pca1_pca$x)
pca1_pca_out<-tibble::rownames_to_column(pca1_pca_out, "rowid")
# pca1_pca_out<-merge(metadata,pca1_pca_out,by="rowid")
# pca1_pca_out$batch<-as.factor(pca1_pca_out$batch)

p<-ggplot(pca1_pca_out,aes(x=PC1,y=PC2))

p+geom_point()

percentage<-round(pca1_pca$sdev / sum(pca1_pca$sdev) * 100, 2)
percentage<- paste(colnames(pca1_pca_out[,5:251]), "(",paste(as.character(percentage), "%", ")", sep=""))

p<-ggplot(pca1_pca_out,aes(x=PC1,y=PC2, color=as.factor(group), shape=batch))
p<-p+geom_point()+theme()+xlab(percentage[1])+ylab(percentage[2])
p<-p+stat_ellipse()
p

p<-ggplot(pca1_pca_out,aes(x=PC3,y=PC4, color=as.factor(group), shape=batch))
p<-p+geom_point()+theme()+xlab(percentage[3])+ylab(percentage[4])
p<-p+stat_ellipse()
p

###--- simple PCA on untransformed data. 

library('FactoMineR')
colanmaes()

obj_wide <- PCA(xset3)
data_tall <- t(xset3[,c(10:140)])
x <- log(data_tall)


for(i in 1:ncol(x)){
  x[x[,i]<0, i] <- 0
}


ob_tall <- PCA(x)

#PCA shows that there is no batch effect bewetween run 6 and run 7. 

############################################################# 
### Annotate ions with WWBC MS database (annotate log-transformed data not batch corrected. Batch effects will be evaluated later using batch # as covariate) 
##ESI pos

library(xMSannotator)
library(dplyr)
getwd()

###Upload the xcms-preprocessed files 
setwd("c:/Users/jesstro/")

dataC <- read.csv("Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/neg_xset03072020.csv")

df<-dataC[,c(2,5,1:141)] ###Keep only columns: mz, rt and each samples
Ccolnames(dataC)
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

mode<-"neg" #ionization mode
#queryadductlist=c("M+H","M+NH4","M+Na","M+H-H2O","M+H-2H2O") #for positive ionization 
queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H") #for negative ionization
adduct_weights<-cbind.data.frame(Adduct=c("M+H","M-H"),Weight=c(5,5))

outloc<-"c:/Users/jesstro/Desktop/neg_annotator_07032020/" ##Folder where to save annotation results

###Upload MS library
# dsstox <- read.csv("c:/Users/jesstro/Box/SHE lab/WWBC/Raw data [de-ID]/WWBC NTA Results & Database/CompTox_17March2019_SelectMetaData.csv", header = T)
# dsstox<-dsstox[,c(3,1,6,7)] ###Keep only columns: DTXSID, chemical name, formula and monoisotopic mass  #### Try with only 5000 vars first. see if error goes away. 
# colnames(dsstox)[1]<-"DTXSID"
# colnames(dsstox)[2]<-"Name"
# colnames(dsstox)[3]<-"Formula"
# colnames(dsstox)[4]<-"MonoisotopicMass"
# dsstox$MonoisotopicMass<-as.numeric(as.character(dsstox$MonoisotopicMass)) ##Should be numeric
# dsstox<-filter(dsstox, MonoisotopicMass>=80 & MonoisotopicMass<=1200) ##Search chemical with mass ranging from 80 to 1,200 DA
# 
# colnames(dsstox)
# customDB<-dsstox

wwbcDB <- read.csv("c:/Users/jesstro/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/Chemicals_of_Interest_DB/WWBC_MS_database_6.24.20.csv")
wwbcDB <- wwbcDB[, c(2,3,5,8)] ###Keep only columns: DTXSID, chemical name, formula and monoisotopic mass
colnames(wwbcDB)
colnames(wwbcDB)[2] <- "Name"
colnames(wwbcDB)[4] <- "MonoisotopicMass"
#remove chemical formulas after the comma since that was giving me an error. 
library(stringr)
wwbcDB$Formula <- str_extract(wwbcDB$Formula,"\\w+")
wwbcDB$MonoisotopicMass <- str_extract(wwbcDB$MonoisotopicMass,"\\w+")


wwbcDB$MonoisotopicMass<-as.numeric(as.character(wwbcDB$MonoisotopicMass)) ##Should be numeric

is.numeric(wwbcDB$MonoisotopicMass)

customDB<-wwbcDB


#if ran previous with error delete the files in the folder. 
# neg_annotaordirectory <- "c:/Users/jesstro/Desktop/FF_non-targeted/neg_annotator/"
# removethese <- dir(neg_annotaordirectory)
# for(i in removethese) {
#   file.remove(paste0(neg_annotaordirectory, i))
# }

##Annotation - Make sure you are setting the correct ionization mode (pos or neg)
annotres<-multilevelannotation(df, max.mz.diff = max.mz.diff, max.rt.diff = max.rt.diff, 
                               cormethod = "pearson", 
                               num_nodes = num_nodes, queryadductlist = queryadductlist, mode = "neg",
                               outloc=outloc, db_name = db_name, adduct_weights = NA, num_sets = num_sets,
                               allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, customIDs = NA,
                               missing.value = NA, deepsplit = 2, 
                               networktype = "unsigned",
                               minclustsize = 10, module.merge.dissimilarity = 0.2, filter.by = c("M-H"),
                               redundancy_check = TRUE,
                               min_ions_perchem = 1, biofluid.location = NA, origin = NA,
                               status = status, boostIDs = NA, max_isp = 5,
                               MplusH.abundance.ratio.check = FALSE, 
                               customDB = customDB, HMDBselect = "union", mass_defect_window = 0.01,
                               mass_defect_mode = "neg",  pathwaycheckmode = "pm")



###Calculate detection frequency for each feature for each group
df[df == 0] <- NA
df$count_n<-apply(df[,23:148],1, function(z) sum(!is.na(z)))
df$count_n<-as.numeric(df$count_n)
df$DF_n<-(df$count_n/126)*100

df$count_ow<-apply(df[,169:254],1, function(z) sum(!is.na(z)))
df$DF_ow<-(df$count_ow/86)*100

df$count_blank<-apply(df[,13:22],1, function(z) sum(!is.na(z)))
df$DF_blank<-(df$count_blank/10)*100

df$DF_diff<-(df$DF_n-df$DF_ow)

####Reshape dataset (df) to have sample file names and intensity in two columns
df[is.na(df)] <- 0
df$time<-round(df$time,1)
df$mz<-round(df$mz,5)
dsn<-melt(df,id.vars = c("mz","time","count_n","DF_n","count_ow","DF_ow","count_blank","DF_blank","DF_diff"))
names(dsn)[names(dsn)=="variable"]<-"filename"
names(dsn)[names(dsn)=="value"]<-"intensity"

###Create id, group and replicate columns
dsn$id<-substr(dsn$filename,2,4)
dsn$id<-as.numeric(dsn$id)
dsn$replicate<-substr(dsn$filename,7,7)
dsn$group<-substr(dsn$filename,1,1)
dsn$group<-as.factor(dsn$group)

###Keep only OW and nurses samples
dsn<-filter(dsn, group=="N" | group=="W" | group=="M")


##Calculate mean, SD and CV peak area for nurses and OW and MB
cv<-function(x) 100*(sd(x)/mean(x))
dsn$featureid<-as.character(paste(dsn$mz,dsn$time,sep = "-"))
dsn_N<-filter(dsn, group=="N")

dsn_N_summary<-dsn_N %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))

dsn_N_summary_t<-dsn_N %>%
  group_by(featureid,id) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))


dsn_N<-merge(dsn_N,dsn_N_summary,by="featureid")
names(dsn_N)[names(dsn_N)=="mean"]<-"mean_N"
names(dsn_N)[names(dsn_N)=="sd"]<-"sd_N"
names(dsn_N)[names(dsn_N)=="cv"]<-"cv_N"
df_N<-dsn_N[,c(1:3,5,7,9,10,16:18)]
df_N<-df_N[!duplicated(df_N$featureid),]

dsn_OW<-filter(dsn, group=="W")

dsn_OW_summary<-dsn_OW %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))
dsn_OW<-merge(dsn_OW,dsn_OW_summary,by="featureid")
names(dsn_OW)[names(dsn_OW)=="mean"]<-"mean_OW"
names(dsn_OW)[names(dsn_OW)=="sd"]<-"sd_OW"
names(dsn_OW)[names(dsn_OW)=="cv"]<-"cv_OW"
df_OW<-dsn_OW[,c(1:3,16:18)]
df_OW<-df_OW[!duplicated(df_OW$featureid),]

dsn_MB<-filter(dsn, group=="M")

dsn_MB_summary<-dsn_MB %>%
  group_by(featureid) %>%
  summarise(mean=mean(intensity), sd=sd(intensity), cv=cv(intensity))
dsn_MB<-merge(dsn_MB,dsn_MB_summary,by="featureid")
names(dsn_MB)[names(dsn_MB)=="mean"]<-"mean_MB"
names(dsn_MB)[names(dsn_MB)=="sd"]<-"sd_MB"
names(dsn_MB)[names(dsn_MB)=="cv"]<-"cv_MB"
df_MB<-dsn_MB[,c(1:3,16:18)]
df_MB<-df_MB[!duplicated(df_MB$featureid),]

df1<-merge(df_N,df_OW,by=c("featureid","mz","time"))
df1<-merge(df1,df_MB,by=c("featureid","mz","time"))

###Check for percent intensity in MB compared to mean N and OW
df1$mean_N_OW<-(df1$mean_N+df1$mean_OW)/2
df1$percent_MB<-(df1$mean_MB/df1$mean_N_OW)*100

##Check for percent intensity in MB (median) compared to median N and OW
median_N_ow<-dsn %>%
  group_by(featureid) %>%
  summarise(Median=median(intensity))

median_MB<-dsn_MB %>%
  group_by(featureid) %>%
  summarise(Median_MB=median(intensity))

df1<-merge(df1,median_N_ow,by="featureid")
df1<-merge(df1,median_MB,by="featureid")
df1$percent_MB<-(df1$Median_MB/df1$Median.x)*100

df1_percentMB<-filter(df1,percent_MB<300)

ggplot(df1_percentMB, aes(x=percent_MB))+
  geom_histogram(binwidth = 5)+
  scale_x_continuous(breaks = seq(0, 300, by = 10))

df_percentMB_10<-filter(df1,percent_MB<=10)


######Exclude feature with cv<30%
df2<-filter(df_percentMB_10, cv_N>30 & cv_OW>30)
df2$DF_high_N<-as.factor(ifelse(df2$DF_diff>=10,1,0))


#####Test for differences in intensity
dsn<-filter(dsn, group=="N" | group=="W")
wilcox<-dsn %>%
  group_by(featureid) %>%
  summarise(p=wilcox.test(intensity~group)$p.value)

wilcox.test(dsn$intensity~dsn$group)

df2<-merge(df2,wilcox,by="featureid")
df2$stat_diff<-as.factor(ifelse(df2$p<=0.05,1,0))

###Merge df2 with annotation output (stage 5) and detection frequency
