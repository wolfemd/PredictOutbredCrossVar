# getMtCrossMeanPreds function
# Wrapper function: for each rep-fold-Model (==unique set of marker effects), predict the relevant cross means.
getMtCrossMeanPreds<-function(outprefix,Model,CrossesToPredict,doseMat,sampleIDs,...){
  start<-proc.time()[3]
  mtbrrFit<-readRDS(here::here("output/mtMarkerEffects",paste0(outprefix,".rds")))
  snpIDs<-mtbrrFit$snpIDs
  traits<-colnames(mtbrrFit$mtbrrFit$yHat)
  doseMat<-doseMat[,snpIDs]

  # Load posterior effects estimates
  postMeanAddEffects<-mtbrrFit$mtbrrFit$ETA$Ga$beta
  colnames(postMeanAddEffects)<-traits;
  rownames(postMeanAddEffects)<-snpIDs
  postMeanAddEffects<-purrr::array_branch(postMeanAddEffects,2)
  if(Model=="A"){
    p<-getAF(doseMat[sampleIDs,])
    P <- matrix(rep(p,nrow(doseMat)),byrow=T,ncol=ncol(doseMat))
    Z <- doseMat-2*P
    predmeans<-predCrossMeanBVs(CrossesToPredict=CrossesToPredict,
                                  postMeanAlleleSubEffects=postMeanAddEffects,
                                  doseMat=Z)
  }
  if(Model=="AD"){
    postMeanDomEffects<-mtbrrFit$mtbrrFit$ETA$Gd$beta
    colnames(postMeanDomEffects)<-traits;
    rownames(postMeanDomEffects)<-snpIDs
    postMeanDomEffects<-purrr::array_branch(postMeanDomEffects,2)
    predmeans<-predCrossMeanTGVs(CrossesToPredict=CrossesToPredict,
                                    postMeanAddEffects=postMeanAddEffects,
                                    postMeanDomEffects=postMeanDomEffects,
                                    doseMat=doseMat)
    }
  runtime<-proc.time()[3]-start
  print(paste0("Predicted the means for all crosses: ",outprefix,". Compute Time: ",runtime/60," mins"))
  return(list(predictedCrossMeans=predmeans,runtime=runtime))
}

# Testing ---------
# outprefix<-parentfolds$outprefix[[2]]
# Model<-parentfolds$Model[[2]]
# CrossesToPredict<-parentfolds$CrossesToPredict[[2]]
# doseMat<-snps
# sireID<-CrossesToPredict$sireID[[21]]
# damID<-CrossesToPredict$damID[[21]]
