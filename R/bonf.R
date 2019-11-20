# Bonf function is a secondary function which computes the weighted Bonferroni
# p-value for an intersection hypothesis

# End of bonf
bonf<-function(p,w)
  # P, Vector of raw p-values
  # W, Vector of hypothesis weights
{
  # Number of null hypotheses included in the intersection
  k<-length(w[w!=0])
  if (k>0) bonf<-min(p[w!=0]/w[w!=0])
  if (k==0) bonf<-1
  return(bonf)
}