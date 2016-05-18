#graph should be a connected graph or shoud be the max connected subgraph

#'calculating network-based separation of disease pairs 
#'
#' given two vectors of diseases, a list of disease-gene associations
#' and a PPI network, this function will calculate network-based separation
#' by method Separation.
#' @param D1 a vector consists disease ids
#' @param D2 another vector consists disease ids
#' @param d2g a list of disease-gene associations
#' @param graph an igraph graph object of PPI network
#' @return a matrix of disease disease network-based separation
#' which rownames is D1 and colnames is D2
#' @importFrom igraph components
#' @importFrom igraph groups
#' @importFrom igraph shortest.paths
#' @export
#' @author Peng Ni, Min Li
#' @references Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease 
#' relationships through the incomplete interactome[J]. Science, 2015, 
#' 347(6224): 1257601.
#' @seealso \code{\link{Separation2Similarity}}
#' @examples
#' data(d2g_separation)
#' data(interactome)
#' 
#' graph_interactome<-graph.data.frame(interactome,directed=FALSE)
#' ds<-sample(names(d2g_separation),5)
#' sep<-Separation(ds,ds,d2g_separation,graph_interactome)
#' sim<-Separation2Similarity(sep)
#' sim
Separation<-function(D1,D2,d2g,graph){
  D1<-intersect(D1,names(d2g))
  D2<-intersect(D2,names(d2g))
  
  stopifnot(length(D1)>0 & length(D2)>0)
  
  # first filter D1 and D2's gene set, use the max connected component of graph 
  message("filtering disease-gene associations..")
  connectedcom<-groups(components(graph))
  con_com_size<-sapply(connectedcom,function(x){return(length(x))},simplify=TRUE)
  maxcc_loci<-match(max(con_com_size),con_com_size)
  maxcc_nodes<-connectedcom[maxcc_loci][[1]]
  
  newdnames<-union(D1,D2)
  newd2g<-lapply(newdnames,
                 function(d,d2g,graph_nodes){return(intersect(d2g[[d]],graph_nodes))},
                 d2g=d2g,graph_nodes=maxcc_nodes)
  names(newd2g)<-newdnames

  message("calculating separation distance between diseases.. this may need a lot of time.. be patient..")
  single_set_distance<-sapply(newd2g,cal_single_set_distance,graph,simplify = TRUE)
  names(single_set_distance)<-newdnames
  
  if(identical(D1,D2)){
    len_d1<-length(D1)
    matlist<-lapply(D1,cal_disease_oneline_separation,D1,
                    d2g=newd2g,graph=graph,single_set_distance=single_set_distance)
    sim<-matrix(nrow = len_d1,ncol = len_d1,dimnames = list(D1,D1))
    for(i in 1:len_d1){
      sim[i,i:len_d1]<-matlist[[i]]
      sim[i:len_d1,i]<-matlist[[i]]
    }
    message("done..")
    return(sim)
  }
  else{
    sim<-matrix(mapply(cal_disease_pair_separation,
                       rep(D1,length(D2)),
                       rep(D2,each=length(D1)),
                       MoreArgs = list(d2g=newd2g,
                                       graph=graph,
                                       single_set_distance=single_set_distance)),
                dimnames=list(D1,D2),
                ncol=length(D2))
    message("done..")
    return(sim)
  }

}


cal_disease_oneline_separation<-function(d,ds,d2g,graph,single_set_distance){
  idloc=match(d,ds)[1]
  n<-length(ds)
  
  #print(idloc)
  
  return(unlist(lapply(ds[idloc:n],cal_disease_pair_separation,d2=d,
                       d2g=d2g,graph=graph,single_set_distance=single_set_distance)))
}

# cal separation between two diseases
cal_disease_pair_separation<-function(d1,d2,d2g,graph,single_set_distance){
  Dab<-cal_set_pair_distance(d2g[[d1]],d2g[[d2]],graph)
  #Daa<-cal_single_set_distance(d2g[[d1]],graph)
  #Dbb<-cal_single_set_distance(d2g[[d2]],graph)
  Daa<-single_set_distance[d1]
  Dbb<-single_set_distance[d2]
  
  return(Dab-(Daa+Dbb)/2)
}

# cal set pair distance 
cal_set_pair_distance<-function(gset1,gset2,graph){
  # genes in gset must be all in graph
  STpath<-shortest.paths(graph,gset1,gset2,mode="all")
  all_min_length_row<-apply(STpath,1,min)
  all_min_length_col<-apply(STpath,2,min)
  return((sum(all_min_length_row)+
            sum(all_min_length_col))/
           (length(all_min_length_row)+
              length(all_min_length_col)))
}

# cal single set distance
cal_single_set_distance<-function(gset,graph){
  # genes in gset must be all in graph
  STpath<-shortest.paths(graph,gset,gset,mode = "all")
  diag(STpath)<-Inf
  all_min_length<-apply(STpath,1,min)
  return(mean(all_min_length))
}

#' a method which convert separation to similarity
#' 
#' convert a separation matrix to a similarity matrix
#' @param data a numeric/integer matrix calculated by method Separation
#' @return a similarity matrix
#' @export
#' @author Peng Ni
#' @seealso \code{\link{Separation}}
#' @examples
#' a<-matrix(c(-4:4),3,3)
#' Separation2Similarity(a)
Separation2Similarity<-function(data){
  stopifnot(class(data)=="matrix")
  #diag(data)<-0
  data<-Normalize(data)
  #sdiag(data)<-0
  return(1-data)
}





