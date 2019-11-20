# IncSimes function is a secondary function which computes a truncated version of the weighted
# incomplete Simes p-value for an intersection hypothesis
incsimes<-function(p,w,gamma)
  # P, Vector of raw p-values
  # W, Vector of hypothesis weights
  # GAMMA, Truncation parameter
{
  # Number of null hypotheses included in the intersection
  k<-length(w[w!=0])
  if (k>1)
  {
    temp<-matrix(0,3,k)
    # P-values
    temp[1,]<-p[w!=0]
    # Weights
    temp[2,]<-w[w!=0]
    # Normalized weights
    temp[3,]<-w[w!=0]/sum(w[w!=0])
    # Sort by p-values
    sorted<-temp[,order(temp[1,])]
    modw<-w[w!=0]
    modw[1]<-0
    modw[2:k]<-sorted[3,1:k-1]
    numer<-sorted[1,]
    denom<-gamma*sorted[3,]/(1-cumsum(modw))+(1-gamma)*sorted[2,]
    incsimes<-min(numer/denom)
  }
  if (k==1)
  {
    numer<-p[w!=0]
    denom<-gamma+(1-gamma)*w[w!=0]
    incsimes<-numer/denom
  }
  if (k==0) incsimes<-1
  return(incsimes)
}
# End of incsimes
