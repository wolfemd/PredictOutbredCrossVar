# getVarComps function -----------
## Wrapper function for getMultiTraitPMVs_A and getMultiTraitPMVs_AD
## For a given Model / data chunk, load stored posterior marker effects
## Compute vars/covars
getVarComps<-function(blups,outName,Model,snps,nIter,burnIn,thin,...){

  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outName,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)

  # Compute genoVarCovarMat
  snps<-snps[rownames(snps) %in% blups$germplasmName,snpIDs]
  Madd<-centerDosage(snps)
  genoVarCovarMat<-genoVarCovarMatFunc(Madd)
  rm(Madd, snps); gc()

  # Load posterior effects estimates
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(outName,"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  if(Model=="AD"){
    DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(outName,"_ETA_Gd_beta.bin")))
    DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  }
  if(Model=="A"){ pmv<-getMultiTraitPMVs_A(AddEffectList=AddEffectList, genoVarCovarMat=genoVarCovarMat) }
  if(Model=="AD"){ pmv<-getMultiTraitPMVs_AD(AddEffectList=AddEffectList,DomEffectList=DomEffectList,genoVarCovarMat=genoVarCovarMat) }
  runtime<-proc.time()[3]-start
  print(paste0("Computed PMVs for",outName," in ",runtime/60," mins"))
  return(list(pmv=pmv,runtime=runtime))
}
