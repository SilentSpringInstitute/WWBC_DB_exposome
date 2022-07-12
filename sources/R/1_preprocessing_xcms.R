#!/usr/bin/env Rscript
library(IPO)
library(xcms)
library(enviGCMS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(mice)


###get the data
getopqedata<-function(p_samples,
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
  cdffiles<-list.files(p_samples, recursive = TRUE, full.names = TRUE)
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



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pr_row_data_QC = args[1]
pr_row_data_sample = args[2]
pr_out = args[3]

pr_row_data_QC = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/FF_raw_data/ESI_neg/QC_neg/"
pr_row_data_sample = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/FF_raw_data/ESI_neg/SerumQCBK_samples_neg/"
pr_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/XCMS/FF_ESI_neg/"


####
# check if optimization is already done
# if not do the peak picking

p_env = paste(pr_out, "optimizedpeak_picking.RData", sep ="")


if(file.exists(p_env)){
  load(p_env)
}else{
  peakpickingParameters<-getDefaultXcmsSetStartingParams("centWave")
  p_QC_files<-list.files(pr_row_data_QC,full.names = T, recursive = T)##Folder that contain QC raw data (.mzXML)
  peakpickingParameters$ppm<-10
  
  resultPeakpicking<-
    optimizeXcmsSet(files = p_QC_files,
                    params = peakpickingParameters,
                    subdir = NULL)
  
  optimizedXcmsSetObject<-resultPeakpicking$best_settings$xset
  retcorGroupParameters<-getDefaultRetGroupStartingParams()
  

  
  # save in a RData
  save(optimizedXcmsSetObject,peakpickingParameters, resultPeakpicking, retcorGroupParameters, file = p_env)
}


########
# Optimize set -> process mzML

# also long process -> Add a save environment here too

p_env_resultRetcorGroup = paste(pr_out, "resultRetGroup.RData")

if(file.exists(p_env_resultRetcorGroup)){
  load(p_env_resultRetcorGroup)
}else{
  resultRetcorGroup<-
    optimizeRetGroup(xset = optimizedXcmsSetObject,
                     params = retcorGroupParameters,
                     subdir = NULL, nSlaves = 2)
  
  save(resultRetcorGroup, file = p_env_resultRetcorGroup)
  
}


### get the data
FF<-getopqedata(pr_row_data_sample)
xset3<-peakTable(FF)
write.csv(xset3, paste(pr_out, "xset.csv", sep = ""))




