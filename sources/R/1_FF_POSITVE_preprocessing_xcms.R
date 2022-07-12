# install.packages("remotes")
# install.packages(c('BiocManager','data.table','digest'))
# 
# BiocManager::install(c('MASS', 'rgl', 'mzR' , 'splines', 'doParallel', 'foreach', 'iterators', 'snow','ROCR', 'ROCS','e1071','randomForest','gbm'))
# BiocManager::install("GO.db")
# BiocManager::install("xcms")
# BiocManager::install("IPO")
# BiocManager::install("multtest")
# 
# remotes::install_github('omegahat/XMLSchema')
# remotes::install_github('cran/SSOAP')
# BiocManager::install(c("SSOAP","KEGGREST","pcaMethods","Rdisop","GO.db","matrixStats","WGCNA"))
# devtools::install_github("yufree/xMSannotator")
# remotes::install_github("jaspershen/MSannotator")
# BiocManager::install("sva")
# remotes::install_github("yufree/xMSanalyzer", dependencies = TRUE)
# 
# devtools::install_github("yufree/enviGCMS")

#####Processing non-targeted analysis data for FF#######
library(IPO)
library(xcms)
library(enviGCMS)
library(dplyr)
library(reshape2)
#library(xMSanalyzer)
#library(xMSannotator)

library(ggplot2)
library(mice)


getwd()
##Find the best parameters for xcms using the R-package IPO
#setwd("E:/MYBOOK/Silent Spring Institute/WWBC/Nurses/ESI_pos/Data preprocessing")
# setwd("c:/Users/jesstro/Desktop/FF_NTA_processing/")
# setwd("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/")


peakpickingParameters<-getDefaultXcmsSetStartingParams("centWave")
path<-list.files("SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processingPOS/Pos_ionization/QC-pos/",full.names = T, recursive = T)##Folder that contain QC raw data (.mzXML)
peakpickingParameters$ppm<-10
resultPeakpicking<-
  optimizeXcmsSet(files = path,
                  params = peakpickingParameters,
                  subdir = NULL)
optimizedXcmsSetObject<-resultPeakpicking$best_settings$xset
retcorGroupParameters<-getDefaultRetGroupStartingParams()


#-----------this command was to save what I had already optimized and had to close the session.-------------

save(optimizedXcmsSetObject,peakpickingParameters, resultPeakpicking, retcorGroupParameters, file = "SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processingPOS/20210325POS_optimizedpeak_picking.RData" )


#----------Start here next time you open R----

load("20210325POS_optimizedpeak_picking.RData")

resultRetcorGroup<-
  optimizeRetGroup(xset = optimizedXcmsSetObject,
                   params = retcorGroupParameters,
                   subdir = NULL)

save(resultRetcorGroup, file = "POS_resultRetcorGroup.RData" )
writeRScript(resultPeakpicking$best_settings$parameters,
             resultRetcorGroup$best_settings)
sessionInfo()

###get the data
getopqdata<-function(path,
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
xset1<-group(xset, bw=0.8799999999,mzwid=0.0265,minfrac=0.1,minsamp=1)
xset2<-retcor(xset1, method="obiwarp",plottype="none",distFunc="cor_opt",profStep=1,center=1,response=1,gapInit=0.928,gapExtend=1.968,factorDiag=2,factorGap=1,localAlignment=0)
xset3<-group(xset2, bw=0.8799999999,mzwid=0.0265,minfrac=0.1,minsamp=1,max=50)
xset4<-fillPeaks(xset3)
FF<-peakTable(xset4)





####Import FF id coding DTSC
# dtsc<-read.csv("C:/Users/jesstro/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/FF_id_visit.csv", header = T)
dtsc<-read.csv("Code/data_files_csv/FF_id_visit.csv", header = T)

###Extract FF id and visit
dtsc$dtsc_id<-substr(dtsc$id,1,7)
dtsc$FF_id<-substr(dtsc$id,9,12)
dtsc$visit<-substr(dtsc$id,14,20)

### get the data
# path<-"c:/Users/jesstro/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processing/Pos_ionization/Serum_samples_pos/" ###Folder that contain all raw data files, including FF, QCs and blanks (.mzXML)
path<-"SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processingPOS/Pos_ionization/Serum_samples_pos/" ###Folder that contain all raw data files, including FF, QCs and blanks (.mzXML)

FF<-getopqdata(path)
xset3<-peakTable(FF)

write.csv(xset3, "SHE lab/WWBC/Data analysis/Non-targeted analysis/FF_non-targeted/ESI_pos/pre-processingPOS/20210423POS_xset.csv")

###---PCA
install.packages("FactoMineR")
library('FactoMineR')


obj_wide <- PCA(xset3)
data_tall <- t(xset3[,c(9:148)])
x <- log(data_tall)


for(i in 1:ncol(x)){
  x[x[,i]<0, i] <- 0
}


ob_tall <- PCA(x)



library(tidyr)
new_df <- gather(xset3, "sample_id", "rel_conc", BD01346.r001:BD01407.r002 )
head(new_df)
new_df$dtsc_id <- substr(new_df$sample_id, 1, 7)
new_df$run <- substr(new_df$sample_id, 9,12)
new_df <- merge(new_df, dtsc, 'dtsc_id', all = TRUE)
head(new_df)
# write.csv(new_df, "pos_long_df.csv")


