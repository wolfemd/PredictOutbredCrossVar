# Wrapper function for runMtCrossVarPredsA / runMtCrossVarPredsAD.
# For each rep-fold-Model (==unique set of marker effects), predict the relevant cross variances.
getMtCrossVarPreds<-function(outprefix,outpath,predType,Model,nIter,burnIn,thin,
                             CrossesToPredict,recombFreqMat,haploMat,ncores,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  recombFreqMat<-recombFreqMat[snpIDs,snpIDs]
  haploMat<-haploMat[,snpIDs]
  # Load posterior effects estimates
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  if(Model=="A"){
    predvars<-runMtCrossVarPredsA(outprefix=outprefix,outpath=outpath,predType="PMV",
                                  CrossesToPredict=CrossesToPredict,
                                  AlleleSubEffectList=AddEffectList,
                                  haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores) }
  if(Model=="AD"){
    DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(gsub("_ReDoSelfs","",outprefix),"_ETA_Gd_beta.bin")))
    DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
    predvars<-runMtCrossVarPredsAD(outprefix=outprefix,outpath=outpath,predType="PMV",
                                   CrossesToPredict=CrossesToPredict,
                                   AddEffectList=AddEffectList,DomEffectList=DomEffectList,
                                   haploMat=haploMat,recombFreqMat=recombFreqMat,ncores=ncores) }
  runtime<-proc.time()[3]-start
  print(paste0("Predicted the vars and covars for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossVars=predvars,runtime=runtime))
}



# # debug - getMtCrossVarPreds -------------------
# outprefix<-parentfolds$outprefix[[1]]
# predType<-"PMV"
# Model<-parentfolds$Model[[1]]
# CrossesToPredict<-parentfolds$CrossesToPredict[[1]][1:30,]
# #
# # # debug - runMtCrossVarPredsA -------------------
# # # debug - predCrossVarsA -------------------
# Trait1<-varcovars$Trait1[[1]]
# Trait2<-varcovars$Trait2[[1]]
# varcovars<-varcovars[1:2,]
# # debug - predOneCrossVarA -------------------
# sireID<-CrossesToPredict$sireID[[1]]
# damID<-CrossesToPredict$damID[[1]]
