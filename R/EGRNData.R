#'Gene expression data for Ensemble Gene Regulatory Network Inference
#'@param path path to file name
#'
#'@return a \code{tibble}
#'@export
#'@importFrom readr read_csv
#'@author Chiranjib Sarkar(cschiranjib9@gmail.com)
#'@references Sarkar, C., Parsad, R., Mishra, D.C. and Rai, A (2020). An ensemble approach for gene regulatory network study in rice blast. Journal of Crop and Weed , 16 , 1-8.
#@examples
#csv=system.file("inst","inst/gene_exp.csv",package = "EGRNi")
#gene_exp(csv)

gene_exp<-function(path){
  readr::read_csv(path)
}
