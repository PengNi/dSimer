#' calculate disease similarity by ICod
#' 
#' given two vectors of diseases, a list of disease-gene associations 
#' and a PPI neteork, this function will calculate disease similarity by method ICod 
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations
#' @param graph an igraph graph object of PPI network
#' @param A a parameter used in ICod to calculate transformed distance of node pair, default 0.9
#' @param b a parameter used in ICod to calculate transformed distance of node pair, default 1
#' @param C a parameter used in ICod to calculate disease similarity, default 0
#' @return a matrix of disease disease simialrity
#' which rownames is D1 and colnames is D2
#' @useDynLib dSimer
#' @importFrom Rcpp sourceCpp
#' @importFrom igraph V
#' @importFrom igraph "V<-"
#' @importFrom igraph neighborhood
#' @importFrom igraph shortest.paths
#' @export
#' @author Peng Ni, Min Li
#' @references Paik H, Heo HS, Ban H, et al. Unraveling human protein interaction 
#' networks underlying co-occurrences of diseases and pathological conditions[J]. 
#' Journal of translational medicine, 2014, 12(1): 99.
#' @examples
#' data(d2g_fundo_symbol)
#' data(PPI_HPRD)
#' 
#' graph_hprd<-graph.data.frame(PPI_HPRD,directed=FALSE) #get a igraph object based on HPRD data
#' ds<-sample(names(d2g_fundo_symbol),5)
#' ICod(ds,ds,d2g_fundo_symbol,graph_hprd)
ICod<-function(D1,D2,d2g,graph,A=0.9,b=1,C=0){
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  message("finding first neighbours of the disease-related genes..")
  newD<-union(D1,D2)
  newD2g<-list()
  
  newGenes<-unique(unlist(d2g[newD]))
  newGenes<-intersect(newGenes,V(graph)$name)
  neighbours<-neighborhood(graph,1,nodes = newGenes,mode="all")
  names(neighbours)<-newGenes
  
  for(i in 1:length(neighbours)){
    neighbours[[i]]<-neighbours[[i]]$name
  }
  
  for(i in 1:length(newD)){
    tmp_v<-intersect(d2g[[newD[i]]],newGenes)
    allneighbours<-Reduce(union,neighbours[tmp_v])
    newD2g[[newD[i]]]<-union(tmp_v,allneighbours)
  }
  
  message("calculating shortest path of each gene pair.. this may take a while..")
  newGenes<-unique(unlist(newD2g))
  #print(length(newGenes))
  SDistance<-shortest.paths(graph,newGenes,newGenes,mode = c("all"))
  
  message("calculating transformed distance of each gene pair.. this may take a while..")
  TDistance<-TDistanceMat_cpp(SDistanceMat = SDistance,A = A,b = b)
  
  message("calculating similarity of each disease pair..")
  mat<-matrix(nrow = length(D1),ncol=length(D2),dimnames = list(D1,D2))
  if(identical(D1,D2)){
    for(i in 1:nrow(mat)){
      for(j in i:ncol(mat)){
        val<-ICod_onepair(D1[i],D2[j],newD2g,SDistance,TDistance,C)
        mat[i,j]<-val
        mat[j,i]<-val
      }
    }
    message("done..")
    return(mat)
  }else{
    mat<-matrix(mapply(ICod_onepair,
                         rep(D1,length(D2)),
                         rep(D2,each=length(D1)),
                         MoreArgs = list(d2g=newD2g,
                                         SDistance=SDistance,
                                         TDistance=TDistance,
                                         C=C)),
                  dimnames = list(D1,D2),
                  ncol=length(D2))
    message("done..")
    return(mat)
  }
}


ICod_onepair<-function(d1,d2,d2g,SDistance,TDistance,C=0){
  if(identical(d1,d2)){
    return(1)
  }else{
    if(length(d2g[[d1]])==0||length(d2g[[d2]])==0){
      return(0)
    }
    
    return(ICod_onepair_cpp(SDistance[d2g[[d1]],d2g[[d2]]],TDistance[d2g[[d1]],d2g[[d2]]],C))
  }
}



