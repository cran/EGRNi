#'To compute F score from edge scores
#'@description Edg2Fw computes the F-score from edge score using Fisher's weighted method. One pair of genes with k numbers of edge scores obtained from k numbers of independent method are combined using the probability value. The weight matrix w contains the weight for k number of methods.
#'@export
#'@param s Matrix containing edge scores obtained from k numbers of methods with gene pairs in row and edge scores in column
#'@param w Matrix containing weight for all individual methods
#'@param k Numbers of independent methods considered for computing edge scores
#'@return Fw_sum matrix containing F score for significant gene pairs
#'@details Edg2Fw function generates mixture distribution based on edge score for each method given column wise in s matrix. The probability value for each pair of gene obtained from mixture distribution are combined using Fisher's weighted method. The combined score Fw follows chi-square distribution with 2k degrees of freedom.
#'@author Chiranjib Sarkar(cschiranjib9@gmail.com)
#'@references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A. (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#'@examples
#'#load "EGRNi" library and Edge score data
#'library(EGRNi)
#'data(Edgescore)
#'data(weight)
#'Edg2Fw(Edgescore, weight, 4)
#'@importFrom fdrtool fdrtool
#'@importFrom gdata upperTriangle
Edg2Fw<-function(s,w,k){
  requireNamespace("gdata")
  requireNamespace("fdrtool")
  p<-as.data.frame(matrix(0,nrow = nrow(s),ncol = k))

  for (i in 3:ncol(s)) {
    mxdist<-fdrtool(s[,i],plot=FALSE)  #### fdr calculation
    p[,i]<-mxdist$pval
  }

  Fw=as.data.frame(matrix(0,nrow = nrow(p),ncol = k))
  for (i in 1:k) {
    Fw[,i]=-2*(w[,i]*(log(p[,i])))
  }
  Fw_sum=rowMeans(Fw)
  Fw_sum<-cbind(s[,1:2],Fw_sum)
  return(Fw_sum)
}
