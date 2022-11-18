#'To compute F score using probability value
#'@description F_score computes the combined edge score from multiple probability values of edges obtained from independent methods. The probability score follows uniform distribution [0,1].The F score follows chi-square distribution with 2k degrees of freedom.
#'@export
#'@param p Matrix containing probability values of edges for each method column wise having gene pairs in row
#'@param w Matrix containing weight for all individual methods
#'@param k Numbers of independent methods considered for computing edge scores
#'@return Fw_sum matrix containing F score for significant gene pairs
#'@details F_score function generates mixture distribution based on probability value for each method given column wise in p matrix. The probability value for each pair of gene are combined using Fisher's weighted method. The combined score Fw follows chi-square distribution with 2k degrees of freedom. The F_score provides the network file as output containing the interacting pair of genes in row with final score.
#'@author Chiranjib Sarkar(cschiranjib9@gmail.com)
#'@references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A. (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#'@examples
#'#load "EGRNi" library and probability value data
#'library(EGRNi)
#'data(pvalue)
#'data(weight)
#'F_score(pvalue, weight, 4)
F_score<-function(p,w,k){
  Fw=0

  for (i in 1:k) {
    F1=-2*(w[,i]*(log(p[,i+2])))
    Fw=Fw+F1
  }

  Fw_sum<-cbind(p[,1:2],Fw)
  return(Fw_sum)
}
