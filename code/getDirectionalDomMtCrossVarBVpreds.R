# Wrapper function for runMtCrossVarPredsA.
# For each rep-fold (==unique set of marker effects), predict the relevant cross variances.
# This version is for a directional dominance model.
# The only difference from getMtCrossVarPreds are:
# 1. that the inbreeding effect for each trait is extract from the BGLR output,
# divided by N snps and added to the vector of SNP effects
# 2. the allele substitution effects are computed as: a+d(q-p)
# runMtCrossVarPredsA() is run and the
# output predicted variances should be the predVar(BV) for each family
getDirectionalDomMtCrossVarBVpreds<-function(outprefix,outpath,predType,nIter,burnIn,thin,
                                              CrossesToPredict,recombFreqMat,haploMat,doseMat,sampleIDs,ncores,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  # estimate genome-wide effect of homozygosity
  homEffects<-mtbrrFit$mtbrrFit$ETA$GmeanD$beta
  # mean dominance effect = -1*genome-wide homEffects divided by number of SNPs
  meanDeffects<--1*homEffects/length(snpIDs)
  colnames(meanDeffects)<-traits
  meanDeffects %<>% array_branch(.,2)
  recombFreqMat<-recombFreqMat[snpIDs,snpIDs]
  haploMat<-haploMat[,snpIDs]
  # Load posterior effects estimates
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),"_ETA_Gd_beta.bin")))
  DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)

  # Add mean dominance effect to dom effects vector
  DomEffectList<-map2(DomEffectList,meanDeffects,~.x+.y)

  # Allele substitution effects = a+d(q-p)
  p<-getAF(doseMat[sampleIDs,snpIDs])
  q<-1-p
  AlleleSubEffectList<-map2(AddEffectList,DomEffectList,~.x+.y*(q-p))
  rm(doseMat) ; gc()

  predvars<-runMtCrossVarPredsA(outprefix=outprefix,outpath=outpath,predType="PMV",
                                CrossesToPredict=CrossesToPredict,
                                AlleleSubEffectList=AlleleSubEffectList,
                                haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores)

  runtime<-proc.time()[3]-start
  print(paste0("Predicted the vars and covars for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossVars=predvars,runtime=runtime))
}



# # debug - getDirectionalDomMtCrossVarPreds -------------------
# outprefix<-parentfolds$outprefix[[1]]
# predType<-"PMV"
# CrossesToPredict<-parentfolds$CrossesToPredict[[1]][1:10,]
# doseMat<-snps
# sampleIDs<-parentfolds$sampleIDs[[1]]
