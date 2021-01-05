# getDirectionalDomMtCrossMeanPreds function
# Wrapper function: for each rep-fold-Model (==unique set of marker effects), predict the relevant cross means.
# This version is for a directional dominance model.
# Differences from getMtCrossMeanPreds:
# 1. inbreeding effect for each trait is extracted from the BGLR output,
### divided by N snps and added to the vector of SNP effects
# 2. Predicted cross mean GEBV:
### Compute allele sub effects as a+d(q-p) and multiply by allelic dosages of parents.
### Cross mean GEBV = 0.5*(GEBV_P1 + GEBV+P2)
# 3. Predicted cross mean GETGV: G = sum( 𝑎(𝑝 − 𝑞 − 𝑦) + 𝑑[2𝑝𝑞 + 𝑦(𝑝 − 𝑞)] )
### a and d being the additive and dominance effects
### p and q being the allele frequencies of one parent
### y is the difference of freq. between the two parents

getDirectionalDomMtCrossMeanPreds<-function(outprefix,CrossesToPredict,doseMat,sampleIDs,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outprefix,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)

  # estimate genome-wide effect of homozygosity
  homEffects<-mtbrrFit$mtbrrFit$ETA$GmeanD$beta
  # mean dominance effect = -1*genome-wide homEffects divided by number of SNPs
  meanDeffects<--1*homEffects/length(snpIDs)
  colnames(meanDeffects)<-traits
  meanDeffects %<>% array_branch(.,2)

  doseMat<-doseMat[,snpIDs]

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
  p<-getAF(doseMat[sampleIDs,])
  q<-1-p
  postMeanAlleleSubEffects<-map2(postMeanAddEffects,postMeanDomEffects,~.x+.y*(q-p))
  P <- matrix(rep(p,nrow(doseMat)),byrow=T,ncol=ncol(doseMat))
  Z <- doseMat-2*P

  predmeanBVs<-predCrossMeanBVs(CrossesToPredict=CrossesToPredict,
                                postMeanAlleleSubEffects=postMeanAlleleSubEffects,
                                doseMat=Z)

  predmeanTGVs<-predCrossMeanTGVs(CrossesToPredict=CrossesToPredict,
                                  postMeanAddEffects=postMeanAddEffects,
                                  postMeanDomEffects=postMeanDomEffects,
                                  doseMat=doseMat)
  predmeans<-left_join(predmeanBVs,predmeanTGVs)
  runtime<-proc.time()[3]-start
  print(paste0("Predicted the means for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossMeans=predmeans,runtime=runtime))
}

# Testing ---------
# outprefix<-parentfolds$outprefix[[2]]
# CrossesToPredict<-parentfolds$CrossesToPredict[[2]]#[1:30,]
# doseMat<-snps
# sampleIDs<-parentfolds$sampleIDs[[2]]
# rm(outprefix,CrossesToPredict,doseMat,sampleIDs)

# sireID<-CrossesToPredict$sireID[[21]]
# damID<-CrossesToPredict$damID[[21]]
