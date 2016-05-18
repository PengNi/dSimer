#' Sun's annotation measure of disease similarity calculating
#' 
#' given two vectors of diseases and a list of disease-gene associations,
#' this function will calculate disease similarity by method Sun_annotation
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations
#' @return a matrix of disease disease simialrity
#' which rownames is D1 and colnames is D2
#' @export
#' @author Peng Ni, Min Li
#' @references Sun K, Goncalves JP, Larminie C. Predicting disease associations 
#' via biological network analysis[J]. BMC bioinformatics, 2014, 15(1): 304.
#' @examples
#' data(d2g_separation)
#' ds<-sample(names(d2g_separation),5)
#' Sun_annotation(ds,ds,d2g_separation)
Sun_annotation<-function(D1,D2,d2g){
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("calculating similarity of each disease pair..")
  if(identical(D1,D2)){
    scorelist<-lapply(D1,jaccardindex_line,xs=D1,d2g)
    
    result.matrix<-matrix(nrow = length(D1),
                          ncol = length(D1),
                          dimnames = list(D1,D1))
    n<-length(D1)
    for(i in 1:n){
      result.matrix[i,i:n]<-scorelist[[i]]
      result.matrix[i:n,i]<-scorelist[[i]]
    }
    message("done..")
    return(result.matrix)
  }else{
    
    result.matrix<-matrix(mapply(jaccardindex,
                         rep(D1,length(D2)),
                         rep(D2,each=length(D1)),
                         MoreArgs=list(x2y=d2g)),
                  dimnames=list(D1,D2),
                  ncol=length(D2))
    message("done..")
    return(result.matrix)
  }
}

#' Sun's function measure of disease similarity calculating
#' 
#' given two vectors of diseases and a list of disease-go term associations,
#' this function will calculate disease similarity by method Sun_function
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2go a list of disease-go term associations
#' @return a matrix of disease disease simialrity
#' which rownames is D1 and colnames is D2
#' @export
#' @author Peng Ni, Min Li
#' @references Sun K, Goncalves JP, Larminie C. Predicting disease associations 
#' via biological network analysis[J]. BMC bioinformatics, 2014, 15(1): 304.
#' @seealso \code{\link{get_GOterm2GeneAssos}}, \code{\link{HypergeometricTest}}
#' @examples
#' ## get a sample of disease-GO associations
#' data(d2go_sample)
#' 
#' ##### the entire disease-GO associations can be obtained by follows:
#' ## go2g<-get_GOterm2GeneAssos(GOONTOLOGY = "BP") #get go-gene associations
#' ## ## in this method, we must use disease-gene associations
#' ## ## which genes are represented by gene symbols
#' ## data(d2g_fundo_symbol)
#' ## d2go<-HypergeometricTest(d2g = d2g_fundo_symbol,go2g = go2g)
#' ##### ###################################################
#' 
#' ds<-names(d2go_sample)
#' Sun_function(ds,ds,d2go_sample)
Sun_function<-function(D1,D2,d2go){
  D1<-intersect(D1,names(d2go))
  D2<-intersect(D2,names(d2go))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("calculating similarity of each disease pair..")
  if(identical(D1,D2)){
    scorelist<-lapply(D1,jaccardindex_line,xs=D1,d2go)
    
    result.matrix<-matrix(nrow = length(D1),
                          ncol = length(D1),
                          dimnames = list(D1,D1))
    n<-length(D1)
    for(i in 1:n){
      result.matrix[i,i:n]<-scorelist[[i]]
      result.matrix[i:n,i]<-scorelist[[i]]
    }

    message("done..")
    return(result.matrix)
  }else{
    result.matrix<-matrix(mapply(jaccardindex,
                         rep(D1,length(D2)),
                         rep(D2,each=length(D1)),
                         MoreArgs=list(x2y=d2go)),
                  dimnames=list(D1,D2),
                  ncol=length(D2))
    
    message("done..")
    return(result.matrix)
  }
}

#' Sun's topology measure of disease similarity calculating
#' 
#' given two vectors of diseases, a list of disease-gene associations,
#' a matrix of genes' graphlet signature in a PPI network and a weight vector of
#' 73 orbits in graphlet theory, this function will calculate disease similarity by 
#' method Sun_function
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations
#' @param graphlet_sig_mat matrix of graphlet signature of nodes in a ppi network calculated by orca, see examples below.
#' @param weight a vector which elements are weight factors to each orbit in graphlet theory 
#' @return a disease disease similarity matrix 
#' @useDynLib dSimer
#' @importFrom Rcpp sourceCpp
#' @export
#' @author Peng Ni, Min Li
#' @references Sun K, Goncalves JP, Larminie C. Predicting disease associations 
#' via biological network analysis[J]. BMC bioinformatics, 2014, 15(1): 304.
#' @examples
#' data(d2g_fundo_symbol)
#' data(graphlet_sig_hprd) #get graphlet signatures of genes in HPRD PPI network
#' data(weight)
#' ds<-sample(names(d2g_fundo_symbol),5)
#' Sun_topology(ds,ds,d2g_fundo_symbol,graphlet_sig_hprd,weight)
Sun_topology<-function(D1,D2,d2g,graphlet_sig_mat,weight){
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("calculating similarity of each disease pair..")
  result.matrix<-Sun_topology_cpp(D1,D2,d2g,graphlet_sig_mat,weight)
  message("done..")
  return(result.matrix)
}




#' set weight factor
#' 
#' set weight factor of 73-orbits in graphlet theory
#' @param orbit_dependency_count a vector which each element are the dependency count of each orbit
#' @return a vector which contains weight factors to each orbit
#' @export
#' @author Peng Ni
#' @references Milenkovic T, Przulj N. Uncovering biological network function via 
#' graphlet degree signatures[J]. Cancer informatics, 2008, 6: 257.
#' @examples
#' data(orbit_dependency_count)
#' setWeight(orbit_dependency_count)
setWeight<-function(orbit_dependency_count){
  return(sapply(orbit_dependency_count,function(x){1-log(x,73)},simplify = TRUE))
}







