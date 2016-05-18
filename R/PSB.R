#' calculate disease similarity by PSB
#'
#' given two vectors of diseases, a list of disease-GO term associations and a 
#' list of GO term-gene associations, this function will calculate disease 
#' similarity by method PSB 
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2go a list of disease-go associations
#' @param go2g a list of go-gene associations
#' @importFrom Rcpp sourceCpp
#' @useDynLib dSimer
#' @return a matrix of disease disease simialrity 
#' which rownames is D1 and colnames is D2
#' @export
#' @author Peng Ni, Min Li
#' @references Mathur S, Dinakarpandian D. Finding disease similarity based on implicit semantic similarity[J]. 
#' Journal of biomedical informatics, 2012, 45(2): 363-371.
#' @seealso \code{\link{get_GOterm2GeneAssos}}, \code{\link{HypergeometricTest}}, \code{\link{Normalize}}
#' @examples

#' ## these are samples of GO-gene associations and disease-GO associations
#' data(go2g_sample)
#' data(d2go_sample)
#' 
#' ##### the entire associations can be obtained by follows:
#' ## go2g<-get_GOterm2GeneAssos(GOONTOLOGY = "BP") #get go-gene associations
#' ## ## in this method, we must use disease-gene associations
#' ## ## which genes are represented by gene symbols
#' ## data(d2g_fundo_symbol)
#' ## d2go<-HypergeometricTest(d2g = d2g_fundo_symbol,go2g = go2g)
#' ##### ###################################################
#' 
#' ds<-names(d2go_sample)
#' sim<-PSB(ds,ds,d2go_sample,go2g_sample)
#' Normalize(sim)
PSB<-function(D1,D2,d2go,go2g){
  d2go<-d2go[!duplicated(names(d2go))]
  
  D1<-intersect(D1,names(d2go))
  D2<-intersect(D2,names(d2go))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  go<-unique(unlist(d2go[union(D1,D2)]))
  
  if(length(go)==0){
    sim<-matrix(data = rep(0,length(D1)*length(D2)),
                ncol = length(D2),
                dimnames = list(D1,D2))
    return(sim)
  }
  
  message("calculating similarity matrix of go terms.. this may need a lot of time..")
  gosim<-termsim(go,go2g)
  
  message("calculating normalize factor of each go term..")
  norFactor<-NormalizingFactor(go,go2g,x2y_conv2_y2x(d2go))
  gosim<-t(t(gosim)*norFactor)
  
  message("calculating simialrity matrix of diseases..")
  len_d1<-length(D1)
  if(identical(D1,D2)){
    matlist<-lapply(D1,psb_line,D1,d2go=d2go,gosim=gosim)
    sim<-matrix(nrow = len_d1,ncol = len_d1,dimnames = list(D1,D1))
    for(i in 1:len_d1){
      sim[i,i:len_d1]<-matlist[[i]]
      sim[i:len_d1,i]<-matlist[[i]]
    }
    message("done..")
    return(sim)
  }else{
    len_d2<-length(D2)
    sim<-matrix(mapply(psb_onepair,
                       rep(D1,len_d2),
                       rep(D2,each=len_d1),
                       MoreArgs = list(d2go,gosim)),
                dimnames=list(D1,D2),
                ncol=len_d2)
    message("done..")
    return(sim)
  }
}

psb_line<-function(d,ds,d2go,gosim){
  idloc=match(d,ds)
  n<-length(ds)
  return(unlist(lapply(ds[idloc:n],psb_onepair,d2=d,d2go,gosim)))
}

psb_onepair<-function(d1,d2,d2go,gosim){
  go1<-d2go[[d1]]
  go2<-d2go[[d2]]
  if(length(go1)==0||length(go2)==0){
    return(0)
  }
  return((sum(apply(as.matrix(gosim[go1,go2]),1,max))/length(go1)+sum(apply(as.matrix(gosim[go2,go1]),1,max))/length(go2))/2)
}


NormalizingFactor<-function(go,go2g,go2d){
  ## adding this sentence or not?
  go2g<-go2g[names(go2g) %in% names(go2d)]
  
  ICgo<-InformationContent(T2G = go2g)
  ICdis<-InformationContent(T2G = go2d)
  MaxICgo<-max(unlist(ICgo))
  MaxICdis<-max(unlist(ICdis))
  
  nf<-sapply(go,NormalizingFactor_single,ICgo,ICdis,MaxICgo,MaxICdis,simplify = TRUE)
  names(nf)<-go
  return(nf)
}

NormalizingFactor_single<-function(termName,ICgo,ICdis,maxICgo,maxICdis){
  return(ICgo[[termName]]*ICdis[[termName]]/maxICgo/maxICdis)
}




termsim<-function(ts,t2g){
  len<-length(ts)
  
  stopifnot(len!=0)
  
  IC<-InformationContent(t2g)
  
  matlist<-lapply(ts,jaccardindex_line,xs=ts,x2y=t2g)
  sim_jc<-matrix(nrow = len,ncol = len,dimnames = list(ts,ts))
  for(i in 1:len){
    sim_jc[i,i:len]<-matlist[[i]]
    sim_jc[i:len,i]<-matlist[[i]]
  }
  
  return(PSB_termsim_cpp(sim_jc,IC))
}






