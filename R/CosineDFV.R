#' calculate disease similarity by using feature vectors
#' 
#' given two (lists of) disease names, this function will calculate cosine 
#' similarity between these diseases' feature vectors.
#' 
#' @param D1 a vector consists of disease ids/names
#' @param D2 another vector consists of disease ids/names
#' @param d2f data.frame, contains term co-occurrences between features and diseases
#' @param dcol integer, disease column number in d2f
#' @param fcol integer, feature column number in d2f
#' @param ccol integer, co-occurrences column number in d2f
#' @return a matrix of disease disease similarity which rownames and colnames 
#' are the disease names
#' @useDynLib dSimer
#' @export
#' @author Zhihui Fei, Peng Ni, Min Li
#' @references Zhou X Z, Menche J, Barabasi A L, et al. Human symptoms-disease 
#' network[J]. Nature communications, 2014, 5.
#' 
#' Van Driel M A, Bruggeman J, Vriend G, et al. A text-mining analysis of the human 
#' phenome[J]. European journal of human genetics, 2006, 14(5): 535-542.
#' @examples
#' ### this is a disease-symptom-cooccurrence sample, if you want to use 
#' ### the complete data, please use "data(d2s_hsdn)" command
#' data(d2s_hsdn_sample)
#' ds <- sample(unique(d2s_hsdn_sample[,2]), 10)
#' simmat <- CosineDFV(ds, ds, d2s_hsdn_sample)
CosineDFV<-function(D1, D2, d2f, dcol=2, fcol=1, ccol=3){
  message("preparing before calculating..")
  features_names <- unique(d2f[, fcol])
  diseases_names <- unique(d2f[, dcol])
  
  D1<-intersect(D1,diseases_names)
  D2<-intersect(D2,diseases_names)
  stopifnot(length(D1)>0 & length(D2)>0)
  
  dlen<-length(diseases_names)
  flen<-length(features_names)
  #======================================================================
  #1.obtain the term frequency matrix from symptoms_occurrence document
  tf <- matrix(data = 0, nrow=flen, ncol =dlen, 
               dimnames = list(features_names,diseases_names) )
  tf <- settf(tf, d2f, dcol = dcol, fcol = fcol, ccol = ccol)
  #======================================================================
  #2.obtain the inverse-document frequency vector from  matrix tf
  idf <- vector(mode = "numeric",length =length(features_names) )
  names(idf) <- features_names
  idf<-setidf(idf,tf)
  #======================================================================
  #3.calculate the dj which is described by a vector of symptoms,then 
  #save in a matrix
  tfidf <- idf*tf
  #======================================================
  #4.calculate the cosine value of each pair of diseases, and save in 
  #a square matrix
  dsim <- setdsim(D1,D2,tfidf)
  return(dsim)
}

settf <- function(tf, d2f, dcol, fcol, ccol){
  d2fnr<-nrow(d2f)
  for(i in 1:d2fnr){
    feature_term <- as.character(d2f[i,fcol])
    disease_term <- as.character(d2f[i,dcol])
    occurrence <-as.numeric(d2f[i,ccol])
    tf[feature_term,disease_term] <- occurrence
  }
  return(tf)
}

setidf <- function(idf, tf){
  N=ncol(tf)
  for(i in 1:nrow(tf)){
    idf[i] <- log(N/length(which(tf[i,]>0)))
  }
  return(idf)
}

setdsim<-function(D1, D2, tfidf){
  tfidfr<-tfidf[,D1]
  if(identical(D1, D2)){
    vmod <- apply(tfidfr,2,modular)
    dsim<-t(tfidfr) %*% tfidfr
    N<-ncol(dsim)
    message("calculating symptom similarity of each disease pair.. this may take a while..")
    for(i in 1:N){
      for(j in i:N){
        dsim[i,j]<-dsim[i,j]/(vmod[i]*vmod[j])
        dsim[j,i]<-dsim[i,j]
      }
    }
    message("calculation completed..")
    return(dsim)
  }else{
    tfidfc<-tfidf[,D2]
    vmodr <- apply(tfidfr,2,modular)
    vmodc <- apply(tfidfc,2,modular)
    dsim<-t(tfidfr) %*% tfidfc
    message("calculating symptom similarity of each disease pair.. this may take a while..")
    dsimnr<-nrow(dsim)
    dsimnc<-ncol(dsim)
    for(i in 1:dsimnr){
      for(j in 1:dsimnc){
        dsim[i,j]<-dsim[i,j]/(vmodr[i]*vmodc[j])
      }
    }
    message("calculation completed..")
    return(dsim)
  }
}

modular <- function(x){
  return(sqrt(sum(x*x)))
}
