#' Correlation based network
#'
#' @param x microarray dataset with genes in columns and samples in rows.
#'
#' @return s matrix containing connectivity scores
#' @export
#'
#' @examples
#'library(EGRNi)
#'data(gene_exp)
#'t_geneexp<-t(gene_exp)
#'CRN(t_geneexp)
#'@importFrom stats cor
CRN = function(x){
requireNamespace("stats")
m<-ncol(x)
s<-cor(x)
colnames(s)<-paste0("Gene",1:m)
rownames(s)<-paste0("Gene",1:m)
return(s)
}

#' Principal component regression based network
#'
#' @param x microarray dataset with genes in columns and samples in rows.
#'
#' @return s matrix containing connectivity scores
#' @export
#' @author Chiranjib Sarkar(cschiranjib9@gmail.com)
#' @references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#' @examples
#'library(EGRNi)
#'data(gene_exp)
#'t_geneexp<-t(gene_exp)
#'PCN(t_geneexp)
#'@importFrom stats cov
#'@importFrom stats lm
PCN = function(x){
requireNamespace("stats")
m<-ncol(x)
s<-as.data.frame(matrix(0,nrow = m,ncol = m))
colnames(s)<-paste0("Gene",1:m)
for (i in 1:m) {
  y<-scale(x[,i],center = TRUE,scale = FALSE)
  k<-scale(x[,-i],center = TRUE,scale = TRUE)
  A<-cov(k)
  EV<-eigen(A,only.values = FALSE)
  Evec<-EV$vectors
  p<-10
  z<-k%*%EV$vectors[,1:p]
  model<-lm(y~z)
  b<-model$coefficients[-1]
  s[i,-i]<-abs(EV$vectors[,1:p]%*%b)
}
s <- 0.5*(s+t(s))
return(s)
}

#' Partial least square based network
#'
#' @param x microarray dataset with genes in columns and samples in rows.
#'
#' @return s matrix containing connectivity scores
#' @export
#'
#' @examples
#'library(EGRNi)
#'data(gene_exp)
#'t_geneexp<-t(gene_exp)
#'PLSN(t_geneexp)
#'@importFrom stats lm
PLSN = function(x){
  requireNamespace("stats")
s<-as.data.frame(matrix(0,nrow = ncol(x),ncol = ncol(x)))
n<-ncol(x)
ncom<-10
for(i in 1:n){
  y<-scale(x[,i],center = TRUE,scale = TRUE)
  k<-scale(x[,-i],center = TRUE,scale = TRUE)

  c <- matrix(0, n-1, ncom)
  TT <- matrix(0, nrow(x), ncom)
  tX<- k
  for(j in 1:ncom){
    c[,j] <- t(tX)%*%y
    TT[,j] <- tX%*%c[,j]
    tX <- tX-TT[,j]%*%solve(t(TT[,j])%*%TT[,j])%*%t(TT[,j])%*%tX
  }
  model<-lm(y~TT)
  b<-model$coefficients[-1]
  s[i,-i]<-c%*%b
}
s <- 0.5*(s+t(s))
return(s)
}

#' Ridge regression based network
#'
#' @param x microarray dataset with genes in columns and samples in rows.
#'
#' @return s matrix containing connectivity scores
#' @export
#'
#' @examples
#'library(EGRNi)
#'data(gene_exp)
#'t_geneexp<-t(gene_exp)
#'RidgN(t_geneexp)
#'@importFrom MASS lm.ridge
RidgN = function(x){
  requireNamespace("MASS")
  s<-as.data.frame(matrix(0,nrow = ncol(x),ncol = ncol(x)))
  n<-ncol(x)
  for (i in 1:n) {
    y<-scale(x[,i],center = TRUE,scale = TRUE)
    k<-scale(x[,-i],center = TRUE,scale = TRUE)
    model<-MASS::lm.ridge(y~k,lambda = 0.0001)
    b<-model$coef
    s[i,-i]<-b
  }
s <- 0.5*(s+t(s))
return(s)
}
