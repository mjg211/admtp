#two-stage adaptive enrichment


#example:
#set.seed(1)
#pvalues=simulateadaptiveenrichment(nperarmperstage=50,proportion.positive=0.4,mean.positive=0.5,mean.negative=0.5,sigma=1,niterations=1000)
#print(operatingcharacteristics(pvalues=pvalues,futilitypvalue=0.5,finalcriticalpvalue=0.05,nperarmperstage=50,proportion.positive=0.4))
#$Power.combined
#[1] 0.93

#$Power.positive.only
#[1] 0.005

#$Expected_sample_size
#[1] 197.52




#simulateadaptiveenrichment simulates a user-specified number of realisations of p-values from a two-arm adaptive enrichment study. P-values for testing the positive, negative and combined groups are recorded for both stages
#Arguments are:
#nperarmperstage - the number of patients to be recruited per arm in each stage. The total sample size is therefore 4*nperarmperstage
#proportion.positive - the probability of each recruited individual being in the postive subgroup
#mean.positive - the mean difference between arms in the positive subgroup
#mean.negative - the mean difference between arms in the negative subgroup
#sigma - the standard deviation of the outcome
#niterations - the number of simulation replicates to perform
#the function returns a matrix of dimension niterations x 6. The first column is the p-values for first stage positive subgroup test;
#the second column is p-values for the first stage negative subgroup; third column is p-values for the first stage combined set of patients
#columns 4-6 are the p-values calculated at the end of the second stage for positive, negative and combined respectively.

simulateadaptiveenrichment=function(nperarmperstage=50,proportion.positive=0.4,mean.positive=0.75,mean.negative=0.5,sigma=1,niterations=1000)
{
  pvalues=matrix(0,niterations,6)
  for(iteration in 1:niterations)
  {
    #simulate the subgroup status and allocation of 4nperarmperstage individuals assuming 1:1 allocation and a probability of being in the positive subgroup of proportion.positive

    subgroup=rbinom(4*nperarmperstage,1,proportion.positive)
    arm=rbinom(4*nperarmperstage,1,0.5)
    #split patients into two stages
    stage=c(rep(1,2*nperarmperstage),rep(2,2*nperarmperstage))
    #simulate outcome data with mean difference being mean.positive in the postive subgroup and mean.negative in the negative subgroup

    outcome=rnorm(4*nperarmperstage,mean.positive*subgroup*arm+mean.negative*(1-subgroup)*arm,sigma)

    #calculate pvalues from test statistics for interim analysis for each group

    test.positive.interim=t.test(outcome[arm==1 & subgroup==1 & stage==1],outcome[arm==0 & subgroup==1 & stage==1])$p.value
    test.negative.interim=t.test(outcome[arm==1 & subgroup==0 & stage==1],outcome[arm==0 & subgroup==0 & stage==1])$p.value
    test.all.interim=t.test(outcome[arm==1 & stage==1],outcome[arm==0 & stage==1])$p.value

    #calculate pvalues from test statistics for final analysis for each group:

    test.positive=t.test(outcome[arm==1 & subgroup==1],outcome[arm==0 & subgroup==1])$p.value
    test.negative=t.test(outcome[arm==1 & subgroup==0],outcome[arm==0 & subgroup==0])$p.value
    test.all=t.test(outcome[arm==1],outcome[arm==0])$p.value

    pvalues[iteration,]=c(test.positive.interim,test.negative.interim,test.all.interim,test.positive,test.negative,test.all)
  }

  return(pvalues)
}