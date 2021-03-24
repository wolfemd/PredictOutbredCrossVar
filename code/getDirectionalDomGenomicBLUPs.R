# getDirectionalDomGenomicBLUPs function
## Similar to "getGenomicBLUPs.R"
## For a given set of multi-trait posterior mean marker effects
## and a given set of SNPs, load the effects and compute GEBV and GETGV
# 1. inbreeding effect for each trait is extracted from the BGLR output,
### divided by N snps and added to the vector of SNP effects
# 2. Allele substitution effects are computed as a+d(q-p) and used to predict GEBV
# 3. GETGV = sum(X_a*a + X_d*d)

getDirectionalDomGenomicBLUPs<-function(outName,sampleIDs,snps,...){
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outName,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  snps<-snps[,snpIDs]

  # estimate genome-wide effect of homozygosity
  homEffects<-mtbrrFit$mtbrrFit$ETA$GmeanD$beta
  # mean dominance effect = -1*genome-wide homEffects divided by number of SNPs
  meanDeffects<--1*homEffects/length(snpIDs)
  colnames(meanDeffects)<-traits
  meanDeffects %<>% array_branch(.,2)

  # Load posterior effects estimates
  postMeanAddEffects<-mtbrrFit$mtbrrFit$ETA$Ga$beta
  colnames(postMeanAddEffects)<-traits;
  rownames(postMeanAddEffects)<-snpIDs
  postMeanAddEffects<-purrr::array_branch(postMeanAddEffects,2)
  postMeanDomEffects<-mtbrrFit$mtbrrFit$ETA$Gd$beta
  colnames(postMeanDomEffects)<-traits;
  rownames(postMeanDomEffects)<-snpIDs
  postMeanDomEffects<-purrr::array_branch(postMeanDomEffects,2)
  # Add mean dominance effect to dom effects vector
  postMeanDomEffects<-map2(postMeanDomEffects,meanDeffects,~.x+.y)

  # Allele substitution effects = a+d(q-p)
  p<-getAF(snps[sampleIDs,]) # frequencies in the training data
  q<-1-p
  postMeanAlleleSubEffects<-map2(postMeanAddEffects,postMeanDomEffects,~.x+.y*(q-p))

  ## Center predictors on overall pop freqs
  ### This probably isn't necessary...
  Madd<-centerDosage(snps)
  Mdom<-dose2domDevGenotypic(snps)

  # GEBV
  gebv<-tibble(Trait=names(postMeanAlleleSubEffects),
               GEBV=postMeanAlleleSubEffects) %>%
    mutate(GEBV=map(GEBV,
                    ~Madd%*%. %>%
                      as.data.frame %>%
                      rownames_to_column(var = "germplasmName") %>%
                      as_tibble),
           predOf="BV") %>%
    unnest(GEBV) %>%
    pivot_wider(values_from = "V1", names_from = "Trait") %>%
    select(germplasmName,predOf,all_of(traits))
  # GETGV
  getgv<-tibble(Trait=names(postMeanAddEffects),
                Ga=postMeanAddEffects,
                Gd=postMeanDomEffects) %>%
    mutate(Ga=map(Ga, # additive genomic values
                    ~Madd%*%. %>%
                      as.data.frame %>%
                      rownames_to_column(var = "germplasmName") %>%
                      as_tibble),
           Gd=map(Gd, # dominance genomic values
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
  gblups<-bind_rows(gebv,getgv)
  return(gblups)
}

# Test ----------
#outName<-parentfolds$outName[[1]]
#sampleIDs<-parentfolds$sampleIDs[[1]]
