# Simes function is a secondary function which computes a truncated version of the
# weighted Simes p-value for an intersection hypothesis
simes<-function(p,w,gamma)
  # P, Vector of raw p-values
  # W, Vector of hypothesis weights
  # GAMMA, Truncation parameter
{
  # Number of null hypotheses included in the intersection
  k<-length(w[w!=0])
  if (k>1)
  {
    temp<-matrix(0,3,k)# TODO: Add comment
    # P-values
    temp[1,]<-p[w!=0]
    # Weights
    temp[2,]<-w[w!=0]
    # Normalized weights
    temp[3,]<-w[w!=0]/sum(w[w!=0])
    # Sort by p-values
    sorted<-temp[,order(temp[1,])]
    numer<-sorted[1,]
    denom<-gamma*cumsum(sorted[3,])+(1-gamma)*sorted[2,]
    simes<-min(numer/denom)
  }
  if (k==1)
  {
    numer<-p[w!=0]
    denom<-gamma+(1-gamma)*w[w!=0]
    simes<-numer/denom
  }
  if (k==0) simes<-1
  return(simes)
}
# End of simes