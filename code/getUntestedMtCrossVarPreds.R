getUntestedMtCrossVarPreds<-function(inprefix,outpath,outprefix,predType,Model,nIter,burnIn,thin,
                                     CrossesToPredict,recombFreqMat,haploMat,ncores,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(inprefix,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  recombFreqMat<-recombFreqMat[snpIDs,snpIDs]
  haploMat<-haploMat[,snpIDs]
  # Load posterior effects estimates
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  if(Model=="A"){
    predvars<-runMtCrossVarPredsA(outprefix=outprefix,outpath=outpath,predType=predType,
                                  CrossesToPredict=CrossesToPredict,
                                  AlleleSubEffectList=AddEffectList,
                                  haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores) }
  if(Model=="AD"){
    DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Gd_beta.bin")))
    DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
    predvars<-runMtCrossVarPredsAD(outprefix=outprefix,outpath=outpath,predType=predType,
                                   CrossesToPredict=CrossesToPredict,
                                   AddEffectList=AddEffectList,DomEffectList=DomEffectList,
                                   haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores) }
  runtime<-proc.time()[3]-start
  print(paste0("Predicted the vars and covars for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossVars=predvars,runtime=runtime))
}

getDirDomUntestedMtCrossVarTGVpreds<-function(inprefix,outpath,outprefix,predType,nIter,burnIn,thin,
                                              CrossesToPredict,recombFreqMat,haploMat,ncores,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(inprefix,".rds")))
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
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Gd_beta.bin")))
  DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)

  # Add mean dominance effect to dom effects vector
  DomEffectList<-map2(DomEffectList,meanDeffects,~.x+.y)

  predvars<-runMtCrossVarPredsAD(outprefix=outprefix,outpath=outpath,predType=predType,
                                 CrossesToPredict=CrossesToPredict,
                                 AddEffectList=AddEffectList,DomEffectList=DomEffectList,
                                 haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores)

  runtime<-proc.time()[3]-start
  print(paste0("Predicted the vars and covars for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossVars=predvars,runtime=runtime))
}

getDirDomUntestedMtCrossVarBVpreds<-function(inprefix,outpath,outprefix,predType,nIter,burnIn,thin,
                                             CrossesToPredict,recombFreqMat,haploMat,doseMat,sampleIDs,ncores,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(inprefix,".rds")))
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
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(inprefix,"_ETA_Gd_beta.bin")))
  DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)

  # Add mean dominance effect to dom effects vector
  DomEffectList<-map2(DomEffectList,meanDeffects,~.x+.y)

  # Allele substitution effects = a+d(q-p)
  p<-getAF(doseMat[sampleIDs,snpIDs])
  q<-1-p
  AlleleSubEffectList<-map2(AddEffectList,DomEffectList,~.x+.y*(q-p))
  rm(doseMat,AddEffectList,DomEffectList) ; gc()

  predvars<-runMtCrossVarPredsA(outprefix=outprefix,outpath=outpath,predType=predType,
                                CrossesToPredict=CrossesToPredict,
                                AlleleSubEffectList=AlleleSubEffectList,
                                haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores)

  runtime<-proc.time()[3]-start
  print(paste0("Predicted the vars and covars for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossVars=predvars,runtime=runtime))
}
