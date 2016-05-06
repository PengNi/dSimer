#' calculating Jaccard Index
#' 
#' calculate Jaccard Index of two terms by using their annotated genes
#' @param x1  a disease id
#' @param x2  another disease id
#' @param x2y a list of disease-gene associations which consists x1 and x2
#' @return numeric value of a jaccard index of x1 and x2
#' @export
#' @author Peng Ni
#' @examples
#' ## this function is not just for disease-gene associations
#' data(d2go_sample)
#' d1<-names(d2go_sample)[1]
#' d2<-names(d2go_sample)[2]
#' jaccardindex(d1,d2,d2go_sample)
jaccardindex<-function(x1,x2,x2y){
  if(identical(x1,x2))return(1)
  x2y<-as.list(x2y)
  i=length(intersect(x2y[[x1]],x2y[[x2]]))
  u=length(union(x2y[[x1]],x2y[[x2]]))
  if(u==0)return(0)
  return(i/u)
}


jaccardindex_line<-function(x,xs,x2y){
  idloc=match(x,xs)
  n<-length(xs)
  return(unlist(lapply(xs[idloc:n],jaccardindex,x2=x,x2y=x2y)))
}


informationContent<-function(x,T2G,totalGeneNum){
  p=as.numeric(length(T2G[[x]])/as.numeric(totalGeneNum))
  ICx=-log2(p)
  return(ICx)
}

#' calculating information content
#' 
#' calculate information content of all term ids in a term list
#' @param T2G a list of Term-Gene associations which names are term ids
#' @return a list of IC values of inputted term ids
#' @export
#' @author Peng Ni
#' @examples
#' data(d2g_fundo_symbol)
#' InformationContent(d2g_fundo_symbol[1:5])
InformationContent<-function(T2G){
  totalGeneNum<-length(unique(unlist(T2G)))
  IC<-lapply(names(T2G),informationContent,T2G,totalGeneNum)
  names(IC)<-names(T2G)
  return(IC)
}

#' Hypergeometric test and multiple testing
#' 
#' given disease-gene associations and go-gene associations, return disease-go associations by using hypergeometric test and fdr mulitiple testing
#' @param d2g a list of disease-gene associations
#' @param go2g a list of GOterm-gene associations
#' @param method multiple testing method, the same as parameter in method p.adjust
#' @param cutoff multiple testing cut off value
#' @return a list of disease-GO term associations
##' @importFrom stats phyper
##' @importFrom stats p.adjust
#' @export
#' @author Peng Ni
#' @seealso \code{\link{PSB}}, \code{\link{Sun_function}}
#' @examples
#' ## see examples in function PSB or Sun_function
#' data(d2go_sample)
#' data(go2g_sample)
#' data(d2g_fundo_symbol)
#' HypergeometricTest(d2g_fundo_symbol[names(d2go_sample)],go2g_sample)
HypergeometricTest<-function(d2g,go2g,method="BH",cutoff=0.05){
  allgene_assogo_num<-length(unique(unlist(go2g)))
  
  g2go<-x2y_conv2_y2x(go2g)
  
  d2go<-lapply(d2g,FindOneDiseaseAssoGOterms,
               go2g,g2go,allgene_assogo_num,
               method,cutoff)
  names(d2go)<-names(d2g)
  return(d2go)
}

FindOneDiseaseAssoGOterms<-function(d,go2g,g2go,allgene_assogo_num,method="BH",cutoff=0.05){
  go_candidate<-unique(unlist(g2go[d]))
  #print(d)
  #print(go_candidate)
  if(length(go_candidate)==0)return(character(0))
  dpvalue<-HypergeometricTest_onedisease(d,go2g[go_candidate],allgene_assogo_num)
  return(MultipleTesting(dpvalue,method,cutoff))
}


HypergeometricTest_onedisease<-function(d,go2g,allgene_assogo_num){
  dlen=length(d)
  #print(length(go2g))
  dpvalue<-unlist(lapply(go2g,HypergeometricTest_onedisease_onego,d=d,dlen=dlen,allgene_assogo_num=allgene_assogo_num))
  names(dpvalue)<-names(go2g)
  return(dpvalue)
}


HypergeometricTest_onedisease_onego<-function(go,d,dlen,allgene_assogo_num){
  intersectnum<-length(intersect(d,go))
  pvalue<-1.0-phyper(intersectnum-1,length(go),allgene_assogo_num-length(go),dlen)
  return(pvalue)
}

MultipleTesting<-function(dpvalue,method="BH",cutoff=0.05){
  dqvalue<-p.adjust(dpvalue,method = method)
  return(names(dqvalue[which(dqvalue<=cutoff)]))
}

#' normalize data
#' 
#' normalize a vector or a matrix based on the formula from SemFunSim
#' @param data a numeric/integer vector or matrix
#' @return normalized vector or matrix
#' @export
#' @author Peng Ni
#' @references Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring 
#' disease similarity by integrating semantic and gene functional association[J]. 
#' PloS one, 2014, 9(6): e99415.
#' @examples
#' sim<-matrix(1:9,3,3)
#' Normalize(sim)
Normalize <-function(data){
  stopifnot(class(data)=="numeric"|class(data)=="integer"|class(data)=="matrix")
  minValue<-min(data)
  delta<-max(data)-minValue
  return((data-minValue)/delta)
}






