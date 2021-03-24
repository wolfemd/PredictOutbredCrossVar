# fitDirectionalDomMtBRR function -------------
## Wrapper function for BGLR::Multitrait()
## For a given set of training blups+snps, fit a directional dominance model as in Xiang et al. 2016.
## using "biologically" partitioned additive and dominance effects + a mean effect for overall proportion homozygous (inbreeding).
fitDirectionalDomMtBRR <- function(blups, snps, outPath, outName, nIter=30000, burnIn=5000,thin=5,...){
  require(tidyverse); require(magrittr); require(BGLR); require(predCrossVar)
  starttime<-proc.time()[3]

  blups %<>%
    filter(germplasmName %in% rownames(snps))
  blupsMat<-blups %>%
    column_to_rownames(var = "germplasmName") %>% as.matrix
  nTraits<-ncol(blupsMat)

  snps<-snps[rownames(snps) %in% rownames(blupsMat),] %>%
    maf_filter(.,0.01) %>%
    remove_invariant(.)
  snpIDs<-colnames(snps) # keep the SNP IDs since each dataset separately MAF filtered

  blups %<>%
    dplyr::mutate(germplasmName=factor(germplasmName,levels=rownames(snps)))
  Xmat<-model.matrix(~germplasmName-1,data=blups)
  # Additive effects
  Madd<-centerDosage(snps)
  eta<-list(Ga=list(X = Xmat%*%Madd, model = "BRR",
                    saveEffects=TRUE, Cov=list(df0 = nTraits)))
  # Dominance effects
  Mdom<-dose2domDevGenotypic(snps)
  eta[["Gd"]]<-list(X = Xmat%*%Mdom, model = "BRR",
                    saveEffects=TRUE, Cov=list(df0 = nTraits))
  # Genome-wide mean dominance effect
  f<-getPropHom(snps)
  eta[["GmeanD"]]<-list(X = Xmat%*%f, model = "FIXED",
                        saveEffects=TRUE, Cov=list(df0 = nTraits))
  rm(snps); gc()
  mtbrrFit<-Multitrait(y=blupsMat,
                       ETA=eta,
                       nIter=nIter,burnIn=burnIn,thin=thin,
                       saveAt=here::here(outPath,paste0(outName,"_")))

  runtime<-proc.time()[3]-starttime
  mtbrrFit<-list(mtbrrFit=mtbrrFit,runtime=runtime,snpIDs=snpIDs)
  saveRDS(mtbrrFit,file=here::here(outPath,paste0(outName,".rds")))
}

# Testing -------------
# blups<-parentfolds$blups[[1]]
# snps<-snps
# outPath<-"output/mtMarkerEffects"
# outName<-"test_test"
# nIter=20; burnIn=5; thin=5
