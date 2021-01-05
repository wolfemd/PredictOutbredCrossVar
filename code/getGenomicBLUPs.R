# getGenomicBLUPs function -----------------
## Wrapper function for either Model=="A" or "AD"
## For a given set of multi-trait posterior mean marker effects
## and a given set of SNPs, load the effects and compute GEBV and GETGV
getGenomicBLUPs<-function(outName,Model,sampleIDs,snps,...){
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outName,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  snps<-snps[,snpIDs]

  # Load posterior effects estimates
  postMeanAddEffects<-mtbrrFit$mtbrrFit$ETA$Ga$beta
  colnames(postMeanAddEffects)<-traits;
  rownames(postMeanAddEffects)<-snpIDs
  postMeanAddEffects<-purrr::array_branch(postMeanAddEffects,2)

  ## Center predictors on overall pop freqs
  ### This probably isn't necessary...?
  Madd<-centerDosage(snps)
  if(Model=="A"){
    # gblups = GEBV
    gblups<-tibble(Trait=names(postMeanAddEffects),
                   GEBV=postMeanAddEffects) %>%
      mutate(GEBV=map(GEBV,
                      ~Madd%*%. %>%
                        as.data.frame %>%
                        rownames_to_column(var = "germplasmName") %>%
                        as_tibble),
             predOf="BV") %>%
      unnest(GEBV) %>%
      pivot_wider(values_from = "V1", names_from = "Trait") %>%
      select(germplasmName,predOf,all_of(traits))
  }
  if(Model=="AD"){
    postMeanDomEffects<-mtbrrFit$mtbrrFit$ETA$Gd$beta
    colnames(postMeanDomEffects)<-traits;
    rownames(postMeanDomEffects)<-snpIDs
    postMeanDomEffects<-purrr::array_branch(postMeanDomEffects,2)
    Mdom<-dose2domDev(snps)
    # gblups = GETGV
    gblups<-tibble(Trait=names(postMeanAddEffects),
                   Ga=postMeanAddEffects,
                   Gd=postMeanDomEffects) %>%
      mutate(Ga=map(Ga, # GEBV
                    ~Madd%*%. %>%
                      as.data.frame %>%
                      rownames_to_column(var = "germplasmName") %>%
                      as_tibble),
             Gd=map(Gd, # GEDD
                    ~Mdom%*%. %>%
                      as.data.frame %>%
                      rownames_to_column(var = "germplasmName") %>%
                      as_tibble),
             Gtot=map2(Ga,Gd, # total genomic values
                       ~left_join(.x %>% rename(Ga=V1),
                                  .y %>% rename(Gd=V1)))) %>%
      select(Trait,Gtot) %>%
      unnest(Gtot) %>%
      mutate(Gtot=Ga+Gd,
             predOf="TGV") %>%
      # Format to (wide-form) matrix of GETGV
      ## make it easy to compute selection indices downstream
      select(-Ga,-Gd) %>% # keep only TGVs
      pivot_wider(values_from = "Gtot", names_from = "Trait") %>%
      select(germplasmName,predOf,all_of(traits))
  }
  return(gblups)
}
