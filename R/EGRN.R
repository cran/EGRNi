#'Ensemble Gene Regulatory Network Inference
#'@description EGRN computes F-score from probability values obtained individual method for each edge. The F-score follows ch-square distribution with 2k degrees of freedom, where k is the number of individual methods consider for ensemble study.The EGRN combines the outcomes obtained from four methods i.e. correlation, principal component regression, partial least square regression and ridge regression.The function EGRN has been implemented using Fisher's weighted method.
#'@export
#'@param x Matrix containing gene expression data with genes in row and samples in column
#'@param n Number of Bootstrap samples to obtain estimate of mean connectivity score and mean square error
#'@param w Matrix containing weight for all individual methods
#'@return Fw_sum matrix containing F score for significant gene pairs
#'@details The function works step-by-step as follows:The input gene expression data is considered for withdrawing n number of bootstrap samples to obtain the estimate of pairwise connectivity score for all possible pairs of genes in the dataset. The t-test statistic is calculated for each pair of genes and performed probability value and false discovery rate calculation from mixture distribution. The p-values for each edge are further used for computing F-score using fisher's weighted method. The fisher's weighted method provides the F-score which follows chi-square distribution with degrees of freedom twice the number of individual methods considered for ensemble study. The EGRN provides the network file as output containing the interacting pair of genes in row with final score.
#'@author Chiranjib Sarkar(cschiranjib9@gmail.com)
#'@references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A. (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#'@examples
#'#load "EGRNi" library and gene expression data
#'library(EGRNi)
#'data(gene_exp)
#'data(weight)
#'EGRN(gene_exp[1:50,], 2, weight)
#'@importFrom fdrtool fdrtool
#'@importFrom stats setNames
#'@importFrom gdata upperTriangle
EGRN <- function(x,n,w) {

  requireNamespace("gdata")
  requireNamespace("fdrtool")
  requireNamespace("stats")
  #### Data set
  y<-t(x)
  #### Connectivity score based on Bootstrap samples

  r<-nrow(y) #### Sample size
  c<-ncol(y) #### No. of genes
  E<-r*c     #### Elemnets in the matrix
  s<-n     #### No. of samples

  listscorr<-vector("list",s)
  listsPC<-vector("list",s)
  listsPLS<-vector("list",s)
  listsridge<-vector("list",s)
  for(i in 1:s){

    Z<-matrix(sample(1:nrow(y),E,replace = T),ncol=c,byrow=T)
    y1<-as.matrix(y[Z[,1],])
    listscorr[[i]]<-CRN(y1)
    listsPC[[i]]<-PCN(y1)
    listsPLS[[i]]<-PLSN(y1)
    listsridge[[i]]<-RidgN(y1)
  }

  ##### MSE calculation
  corsum<-0
  PCsum<-0
  PLSsum<-0
  ridgesum<-0
  for(j in 1:s){
    r1<-abs(listscorr[[j]])
    corsum<-r1+corsum
    r2<-abs(listsPC[[j]])
    PCsum<-r2+PCsum
    r3<-abs(listsPLS[[j]])
    PLSsum<-r3+PLSsum
    r4<-abs(listsridge[[j]])
    ridgesum<-r4+ridgesum
  }
  mean_corr<-corsum/s
  mean_PC<-PCsum/s
  mean_PLS<-PLSsum/s
  mean_ridge<-ridgesum/s
  dvsq_cor<-0
  dvsq_PC<-0
  dvsq_PLS<-0
  dvsq_ridge<-0
  for(j in 1:s){
    r1<-abs(listscorr[[j]])
    r11<-((r1-mean_corr)^2)
    dvsq_cor<-(dvsq_cor+r11)

    r2<-abs(listsPC[[j]])
    r22<-((r2-mean_PC)^2)
    dvsq_PC<-(dvsq_PC+r22)

    r3<-abs(listsPLS[[j]])
    r33<-((r3-mean_PLS)^2)
    dvsq_PLS<-(dvsq_PLS+r33)

    r4<-abs(listsridge[[j]])
    r44<-((r4-mean_ridge)^2)
    dvsq_ridge<-(dvsq_ridge+r44)
  }
  MSE_cor<-sqrt(dvsq_cor/(s-1))
  MSE_PC<-sqrt(dvsq_PC/(s-1))
  MSE_PLS<-sqrt(dvsq_PLS/(s-1))
  MSE_ridge<-sqrt(dvsq_ridge/(s-1))

  ### Converting Matrix into Datafrme with NAME
  X<-setNames(as.data.frame.table(mean_corr),c("rownames","colnames","values")) ##### setNames arrange by column
  X<-subset(X,!duplicated(X[,3])) #### Removing Duplicate values
  X<-as.matrix(X[-1,])

  ### t-test statistic and mixture distribution
  cor_mean<-upperTriangle(mean_corr,diag = FALSE,byrow = TRUE) #### Mean vector
  cor_MSE<-upperTriangle(MSE_cor,diag = FALSE,byrow = TRUE)   #### MSE vector
  cor_t<-as.vector(cor_mean/cor_MSE)                                #### t-test statistic
  mxdist_cor<-fdrtool(abs(cor_t), statistic = "normal",plot=FALSE)  #### fdr calculation

  PC_mean<-upperTriangle(mean_PC,diag = FALSE,byrow = TRUE) #### Mean vector
  PC_MSE<-upperTriangle(MSE_PC,diag = FALSE,byrow = TRUE)   #### MSE vector
  PC_t<-as.vector(PC_mean/PC_MSE)                                #### t-test statistic
  mxdist_PC<-fdrtool(PC_mean, statistic = "normal",plot=FALSE)  #### fdr calculation

  PLS_mean<-upperTriangle(mean_PLS,diag = FALSE,byrow = TRUE) #### Mean vector
  PLS_MSE<-upperTriangle(MSE_PLS,diag = FALSE,byrow = TRUE)   #### MSE vector
  PLS_t<-as.vector(PLS_mean/PLS_MSE)                                #### t-test statistic
  mxdist_PLS<-fdrtool(PLS_mean, statistic = "normal",plot=FALSE)  #### fdr calculation

  ridge_mean<-upperTriangle(mean_ridge,diag = FALSE,byrow = TRUE) #### Mean vector
  ridge_MSE<-upperTriangle(MSE_ridge,diag = FALSE,byrow = TRUE)   #### MSE vector
  ridge_t<-as.vector(ridge_mean/ridge_MSE)                                #### t-test statistic
  mxdist_ridge<-fdrtool(abs(ridge_mean), statistic = "normal",plot=FALSE)  #### fdr calculation

  t_fdr_pval<-cbind(X[,1:2],cor_t,mxdist_cor$lfdr,mxdist_cor$pval,PC_t,mxdist_PC$lfdr,mxdist_PC$pval,PLS_t,mxdist_PLS$lfdr,mxdist_PLS$pval, ridge_t,mxdist_ridge$lfdr,mxdist_ridge$pval)

  colnames(t_fdr_pval)<-c("Gene_name1","Gene_name2","t(cor)","FDR(cor)","P(cor)","t(PC)","FDR(PC)","P(PC)","t(PLS)","FDR(PLS)","P(PLS)","t(ridge)","FDR(ridge)","P(ridge)")

  ######### Fisher's weighted method
  p<-cbind(mxdist_cor$pval,mxdist_PC$pval,mxdist_PLS$pval,mxdist_ridge$pval)
  Fw=matrix(0,nrow = nrow(t_fdr_pval),ncol = 4)
  for (i in 1:4) {
    Fw[,i]=-2*(w[,i]*(log(p[,i])))
  }
  Fw_sum=rowMeans(Fw)
  return(Fw_sum)

}
