# fitMtBRR function -------------
## Wrapper function for BGLR::Multitrait()
## For a given set of training blups+snps for either model "A" or "AD"
fitmtBRR <- function(blups, snps, Model, outPath, outName, nIter=30000, burnIn=5000,thin=5,...){
  require(tidyverse); require(magrittr); require(BGLR); require(predCrossVar)
  starttime<-proc.time()[3]
  blupsMat<-blups %>%
    filter(germplasmName %in% rownames(snps)) %>%
    column_to_rownames(var = "germplasmName") %>% as.matrix
  nTraits<-ncol(blupsMat)

  snps<-snps[rownames(snps) %in% rownames(blupsMat),] %>%
    maf_filter(.,0.01) %>%
    remove_invariant(.)
  snpIDs<-colnames(snps) # keep the SNP IDs since each dataset separately MAF filtered
  Madd<-centerDosage(snps)
  eta<-list(Ga=list(X = Madd, model = "BRR",
                    saveEffects=TRUE, Cov=list(df0 = nTraits)))

  if(Model=="AD"){
    Mdom<-dose2domDev(snps)
    eta[["Gd"]]<-list(X = Mdom, model = "BRR",
                      saveEffects=TRUE, Cov=list(df0 = nTraits))
  }
  rm(snps); gc()
  mtbrrFit<-Multitrait(y=blupsMat,
                       ETA=eta,
                       nIter=nIter,burnIn=burnIn,thin=thin,
                       saveAt=here::here(outPath,paste0(outName,"_")))
  runtime<-proc.time()[3]-starttime
  mtbrrFit<-list(mtbrrFit=mtbrrFit,runtime=runtime,snpIDs=snpIDs)
  saveRDS(mtbrrFit,file=here::here(outPath,paste0(outName,".rds")))
}
