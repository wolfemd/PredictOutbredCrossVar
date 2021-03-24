# getDirectionalDomVarComps function
## Wrapper function for getMultiTraitPMVs_A and getMultiTraitPMVs_AD
## For a given Model / data chunk, load stored posterior marker effects
## Compute vars/covars
getDirectionalDomVarComps<-function(blups,outName,Model,snps,nIter,burnIn,thin,...){

  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outName,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  # estimate genome-wide effect of homozygosity
  homEffects<-mtbrrFit$mtbrrFit$ETA$GmeanD$beta
  # mean dominance effect = -1*genome-wide homEffects divided by number of SNPs
  meanDeffects<--1*homEffects/length(snpIDs)
  colnames(meanDeffects)<-traits
  meanDeffects %<>% array_branch(.,2)

  # Compute genoVarCovarMat
  snps<-snps[rownames(snps) %in% blups$germplasmName,snpIDs]
  Madd<-centerDosage(snps)
  genoVarCovarMat<-genoVarCovarMatFunc(Madd)

  # Load posterior effects estimates
  AddEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(outName,"_ETA_Ga_beta.bin")))
  AddEffectList<-effectsArray2list(AddEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  DomEffectList<-readBinMatMultitrait(filename = here::here("output/mtMarkerEffects",paste0(outName,"_ETA_Gd_beta.bin")))
  DomEffectList<-effectsArray2list(DomEffectList,snpIDs,traits,nIter=30000, burnIn=5000,thin=5)
  # Add mean dominance effect to dom effects vector
  DomEffectList<-map2(DomEffectList,meanDeffects,~.x+.y)

  if(Model=="DirDomA"){
    # Allele substitution effects = a+d(q-p)
    p<-getAF(snps) # frequencies in the training data
    q<-1-p
    AlleleSubEffectList<-map2(AddEffectList,DomEffectList,~.x+.y*(q-p))
    rm(Madd, snps,AddEffectList,DomEffectList); gc()
    pmv<-getMultiTraitPMVs_A(AddEffectList=AlleleSubEffectList, genoVarCovarMat=genoVarCovarMat)
  }
  if(Model=="DirDomAD"){
    rm(Madd, snps); gc()
    pmv<-getMultiTraitPMVs_AD(AddEffectList=AddEffectList,DomEffectList=DomEffectList,genoVarCovarMat=genoVarCovarMat)
  }
  runtime<-proc.time()[3]-start
  print(paste0("Computed PMVs for",outName," in ",runtime/60," mins"))
  return(list(pmv=pmv,runtime=runtime))
}
