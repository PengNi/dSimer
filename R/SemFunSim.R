
#' calculate disease similarity by FunSim
#'
#' given two vectors of diseases, a list of disease-gene associations ,
#' and a list of gene-gene log-likelihood score from HumanNet, 
#' this function will calculate disease similarity by method FunSim 
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations, while gene ids should be entrez id.
#' @param LLSnList a list of gene-gene log-likelihood score from HumanNet
#' @return a matrix of disease disease simialrity
#' which rownames is D1 and colnames is D2
#' @export
#' @author Peng Ni, Min Li
#' @references Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring 
#' disease similarity by integrating semantic and gene functional association[J]. 
#' PloS one, 2014, 9(6): e99415.
#' @seealso \code{\link{LLSn2List}}
#' @examples
#' ## in this method, we must use disease-gene associations 
#' ## which genes are represented by entrez ids because of
#' ## HumanNet
#' data(d2g_fundo_entrezid)
#' data(HumanNet_sample)
#' ## we specified 5 DOIDs to match Human_sample
#' ds<-c("DOID:8176","DOID:2394","DOID:3744","DOID:8466","DOID:5679")
#' llsnlist<-LLSn2List(HumanNet_sample)
#' FunSim(ds,ds,d2g_fundo_entrezid,llsnlist)
FunSim<-function(D1,D2,d2g,LLSnList){
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  lend1<-length(D1)
  lend2<-length(D2)
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("calculating simialrity matrix of diseases.. this may take a lot of time..")
  if(identical(D1,D2)){
    result.matrix<-matrix(nrow = length(D1),ncol = length(D1),dimnames = list(D1,D1))
    for(i in 1:length(D1)){
      for(j in i:length(D1)){
        sim<-FunSim_onepair(D1[i],D1[j],d2g,LLSnList)
        result.matrix[i,j]<-sim
        result.matrix[j,i]<-sim
      }
    }
    message("done..")
    return(result.matrix)
  }else{
    result.matrix<-matrix(mapply(FunSim_onepair,
                         rep(D1,lend2),
                         rep(D2,each=lend1),
                         MoreArgs = list(d2g=d2g,LLSnList=LLSnList)),
                  dimnames = list(D1,D2),ncol=lend2)
    message("done..")
    return(result.matrix)
  }
}


FunSim_onepair<-function(d1,d2,d2g,LLSnList){
  
  return((sum(vapply(d2g[[d1]],
                     funsim,numeric(1),
                     d2g[[d2]],
                     LLSnList))+
            sum(vapply(d2g[[d2]],
                       funsim,numeric(1),
                       d2g[[d1]],
                       LLSnList)))/
           (length(d2g[[d1]])+length(d2g[[d2]])))
}


funsim<-function(g1,G2,LLSnList){
  if(g1 %in% G2) return(1)
  
  simvalues<-.subset2(LLSnList,g1)[intersect(G2,names(LLSnList[[g1]]))]
  
  if(length(simvalues)==0) return(0)
  return(max(simvalues))
}



