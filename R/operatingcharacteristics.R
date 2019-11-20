#operatingcharacteristics calculates the probability of recommending the treatment in the positive group and combined group after the adaptive enrichment study
#Arguments:
#pvalues - output from simulateadaptiveenrichment
#futilitypvalue - threshold for stopping for futility if the p-value is above this
#finalcriticalvalue - threshold on p-value at end of trial: if p-value is below this the hypothesis is rejected
#nperarmperstage - same as in simulateadaptiveenrichment, only used for getting expected sample size
#proportion.positive - same as in simualteadaptiveenrichment, only used for getting expected sample size.


#Returns a list containing the probability of recommending treatment in overall group and in positive subgroup together with expected sample size used.
operatingcharacteristics=function(pvalues,futilitypvalue,finalcriticalpvalue,nperarmperstage,proportion.positive)
{
  #set up vector to record whether null hypothesis for positive group was rejected
  rejecth0.positive=rep(0,dim(pvalues)[1])

  #similar, for combined group
  rejecth0.combined=rep(0,dim(pvalues)[1])

  #reject combined group hypothesis if: 1) interim combined pvalue is below futilitythreshold; 2) final combined pvalue is below finalcriticalpvalue

  rejecth0.combined=ifelse(pvalues[,3]<=futilitypvalue & pvalues[,6]<=finalcriticalpvalue,1,0)

  #reject positive group hypothesis if: 1a) interim combined pvalues is above futilitypvalue AND interim positive pvalue is below futilitypvalue 2a) final positive pvalue is below finalcriticalpvalue
  #OR 1b) interim combined pvalue is below futilitypvalue AND 2b) final combined pvalue is above finalcriticalpvalue AND 3b) final postive pvalue is below finalcriticalpvalue
  rejecth0.positive=ifelse(pvalues[,3]>futilitypvalue & pvalues[,1]<=futilitypvalue & pvalues[,4]<=finalcriticalpvalue,1,0)+
    ifelse(pvalues[,3]<=futilitypvalue & pvalues[,6]>finalcriticalpvalue & pvalues[,4]<=finalcriticalpvalue,1,0)



  #approximate sample size used in each case. This is 2*nperarmperstage if trial stops for futility, 2*nperarmperstage+2*proportion.positive*nperarmperstage if trial continues in positive subgroup only, and 4*nperarmperstage otherwise

  samplesize=rep(4*nperarmperstage,dim(pvalues)[1])
  samplesize=replace(samplesize,pvalues[,1]>futilitypvalue & pvalues[,3]>futilitypvalue,2*nperarmperstage)
  samplesize=replace(samplesize,pvalues[,1]<=futilitypvalue & pvalues[,3]>futilitypvalue,2*nperarmperstage+2*proportion.positive*nperarmperstage)

  return(list(Power.combined=mean(rejecth0.combined),Power.positive.only=mean(rejecth0.positive),Expected_sample_size=mean(samplesize)))
}