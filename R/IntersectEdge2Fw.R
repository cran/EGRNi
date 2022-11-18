#'To compute F score for significant edges from individual methods
#'@description IntsctEdg2Fw computes the Fw score using Fisher's weighted method for the significant edges obtained in k numbers of individual methods. The probability values are combined to compute the Fw score which follows chi-square distribution. The significant edges (>fdr) are selected using intersection.
#'@export
#'@param s Matrix containing edge scores obtained from k numbers of methods with gene pairs in row and edge scores in column
#'@param w Matrix containing weight for all individual methods
#'@param k Numbers of independent methods considered for computing edge scores
#'@param fdr Cut-off for selecting significant edges
#'@return Fw matrix containing F score for significant gene pairs
#'@details IntsctEdg2Fw function generates mixture distribution based on edge score for each method given column wise in s matrix. The probability value for each pair of gene obtained from mixture distribution are combined using Fisher's weighted method. The combined score Fw follows chi-square distribution with 2k degrees of freedom.
#'@author Chiranjib Sarkar(cschiranjib9@gmail.com)
#'@references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A. (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#'@examples
#'#load "EGRNi" library and Edge score data
#'library(EGRNi)
#'data(Edgescore)
#'data(weight)
#'IntsctEdg2Fw(Edgescore,weight, k=4, fdr=0.1)
#'@importFrom fdrtool fdrtool
#'@importFrom gdata upperTriangle
IntsctEdg2Fw<-function(s,w,k,fdr){
  requireNamespace("gdata")
  requireNamespace("fdrtool")
  p<-as.data.frame(matrix(0,nrow = nrow(s),ncol = k))
  FDR<-as.data.frame(matrix(0,nrow = nrow(s),ncol = k))
  for (i in 1:k) {
    mxdist<-fdrtool(s[,i+2],plot=FALSE)  #### fdr calculation
    p[,i]<-mxdist$pval
    FDR[,i]<-mxdist$lfdr
  }
  p<-cbind(s[,1:2],p)
  for (j in 1:k) {
    p=subset(p,FDR[,j]<fdr)
  }

  Fw=0
  for (l in 1:k) {
    F1=-2*(w[,l]*(log(p[,l+2])))
    Fw=Fw+F1
  }

  return(Fw)
}
