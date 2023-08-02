
#'@title CalDEscore
#'@description Function "CalDEscore" uses gene expression to calculate differential expression level.
#'@param exp A gene expression profile of interest (rows are genes, columns are samples).
#'@param Label A character vector consist of "0" and "1" which represent sample class in the gene expression profile. "0" means normal sample and "1" means disease sample.
#'@return A matrix with one column of zscore.
#'@importFrom stats pt
#'@importFrom stats qnorm
#'@usage CalDEscore(exp, Label)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-Gettest("GEP")
#' label<-Gettest("label")
#' # Run the function
#' DEscore<-CalDEscore(GEP,label)

CalDEscore<-function(exp, Label){
  expdata<-as.matrix(exp)
  test<-matrix(nrow=nrow(expdata),ncol=1)
  rownames(test)<-rownames(expdata)
  colnames(test)<-c("Zscore")
  logfc.vector<-apply(expdata, 1, function(x){
    ind1 <- which(Label == 1)
    ind2 <- which(Label == 0)
    m <- length(ind1 <- which(Label == 1))
    n <- length(ind2 <- which(Label == 0))
    expdata1 <- x[ind1]
    expdata2 <- x[ind2]
    rmean1 <- mean(expdata1)
    rmean2 <- mean(expdata2)
    ss1 <- sum((expdata1 - rmean1)^2)
    ss2 <- sum((expdata2 - rmean2)^2)
    Tscore <- (m + n - 2)^0.5*(rmean1 - rmean2)/((1/m + 1/n)*(ss1 + ss2))^0.5
    Pvalue <- pt(abs(Tscore),lower.tail=FALSE,df=m+n-2)
    Zvalue <- qnorm(Pvalue, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    return(Zvalue)
  })
  test[,1]<-logfc.vector
  test <- as.matrix(test[!is.infinite(test[,1]),])
  test <- as.matrix(test[!is.na(test[,1]),])
  test <- as.matrix(test[!is.nan(test[,1]),])

  colnames(test)<-c("Zscore")
  return(test)
}
