

sim_mat<-function(d2g){
  #d2g<-d2g[!duplicated(names(d2g))]
  totalGeneNum<-length(unique(unlist(d2g)))
  Genenumlist<-lapply(d2g,function(x) length(unique(x)))
  names(Genenumlist)<-names(d2g)
  
  return(BOG_simmat_cpp(d2g,totalGeneNum,Genenumlist))
}

#' calculate disease similarity by BOG
#' 
#' given two vectors of diseases and a list of disease-gene associations,
#' this function will calculate disease similarity by method BOG.
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations
#' @return a matrix of disease disease simialrity 
#' which rownames is D1 and colnames is D2
#' @importFrom Rcpp sourceCpp
#' @useDynLib dSimer
#' @export
#' @author Peng Ni, Min Li
#' @references Mathur S, Dinakarpandian D. Automated ontological gene annotation for computing disease similarity[J]. 
#'  AMIA Summits on Translational Science Proceedings, 2010, 2010: 12
#' @seealso \code{\link{Normalize}}
#' @examples
#' data(d2g_separation) #get disease-gene associations
#' ds<-sample(names(d2g_separation),5)
#' sim<-BOG(ds,ds,d2g_separation)
#' Normalize(sim) #normalize BOG sim scores
BOG<-function(D1,D2,d2g){
  d2g<-d2g[!duplicated(names(d2g))]
  
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("calculating Information Content of diseases..")
  IC<-InformationContent(d2g)
  
  message("calculating similarity matrix of diseases.. this may take a while..")
  simmat<-sim_mat(d2g)
  
  message("normalize the similarity matrix.. this may also take a while..")
  maxsim<-apply(simmat,1,max)
  names(maxsim)<-row.names(simmat)
  
  snames<-as.list(c(0:(nrow(simmat)-1)))
  names(snames)<-row.names(simmat)
  
  sim<-BOG_normat_cpp(D1,D2,simmat,snames,maxsim,IC);
  message("done..")
  return(sim)
}
