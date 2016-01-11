#' convert x-y associations
#' 
#' concert x-y associations (e.g. disease-gene associations) from data.frame 
#' to list
#' @param x2ydf data.frame of x-y associations
#' @param xcol col of x in x2ydf
#' @param ycol col of y in x2ydf
#' @return a list of x-y associations
#' @export
#' @author Peng Ni
#' @examples
#' options(stringsAsFactors = FALSE)
#' 
#' d2g_fundo_sample<-read.table(text = "DOID:5218    IL6
#' DOID:8649  EGFR
#' DOID:8649	PTGS2
#' DOID:8649	VHL
#' DOID:8649	ERBB2
#' DOID:8649	PDCD1
#' DOID:8649	KLRC1
#' DOID:5214	MPZ
#' DOID:5214	EGR2
#' DOID:5210	AMH")
#' 
#' d2g_fundo_list<-x2y_df2list(d2g_fundo_sample)
x2y_df2list<-function(x2ydf,xcol=1,ycol=2){
  colnames(x2ydf)[c(xcol,ycol)]<-c('V1','V2')
  
  class(x2ydf[,xcol])<-"character"
  class(x2ydf[,ycol])<-"character"
  
  xname<-unique(as.data.frame(x2ydf)[,xcol])
  x2y_list<-list(length=length(xname))
  for(i in 1:length(xname)){
    x2y_list[[i]]<-unique(c(x2ydf[x2ydf$V1==xname[i],]$V2))
  }
  names(x2y_list)<-xname
  return(x2y_list)
}

#' convert x2ylist to y2xlist
#' 
#' convert list of x-y associations to list of y-x associations
#' @param x2ylist a list which the names are xs and the elements are ys of each x
#' @return a list of y2x
#' @export
#' @author Peng Ni
#' @examples
#' data(go2g_sample)
#' g2go_sample<-x2y_conv2_y2x(go2g_sample[1:100])
x2y_conv2_y2x<-function(x2ylist){
  xnames<-names(x2ylist)
  ynames<-unique(unlist(x2ylist))
  y2xlist<-list()
  for(n in xnames){
   gofn<-x2ylist[[n]]
   for(g in gofn){
     y2xlist[[g]]<-union(y2xlist[[g]],n)
   }
  }
  return(y2xlist)
}

## AnnotationDbi using org.Hs.eg.db, GO-gene association remove IEAs
## only BP
## return a list
#' get GO-gene associations
#' 
#' get GO-gene associations from GO.db and org.Hs.eg.db
#' @param GOONTOLOGY "BP" or "MF" or "CC
#' @param rm.IEAs logical value, remove GO terms with evidence "IEA" or not
#' @param rm.termlessthan3genes logical value, remove terms whose number of annotated genes are less than 3 or not
#' @return a list which names are GO term IDs and elements are gene symbols annotated with GO terms
#' @importMethodsFrom AnnotationDbi as.list
#' @importFrom GO.db GOBPOFFSPRING
#' @importFrom GO.db GOMFOFFSPRING
#' @importFrom GO.db GOCCOFFSPRING
#' @importMethodsFrom AnnotationDbi keys
#' @importMethodsFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @export
#' @author Peng Ni
#' @references Mathur S, Dinakarpandian D. Finding disease similarity based on implicit semantic similarity[J]. 
#' Journal of biomedical informatics, 2012, 45(2): 363-371.
#' @seealso \code{\link{PSB}}, \code{\link{Sun_function}}
#' @examples
#' ## see examples in function PSB or Sun_function
get_GOterm2GeneAssos<-function(GOONTOLOGY=c("BP","MF","CC"),
                               rm.IEAs=TRUE,
                               rm.termlessthan3genes=TRUE){
  entKeys <- keys(org.Hs.eg.db, keytype="ENTREZID")
  cols<-c("ENTREZID","SYMBOL","ONTOLOGY","GO","EVIDENCE")
  asso<-select(org.Hs.eg.db, keys=entKeys, columns=cols, keytype="ENTREZID")
  
  GOONTOLOGY<-match.arg(GOONTOLOGY)
  asso<-asso[which(asso$ONTOLOGY==GOONTOLOGY),]
  
  if(rm.IEAs==TRUE){
    asso<-asso[which(asso$EVIDENCE!="IEA"),]
  }
  
  bpnames<-unique(asso$GO)
  switch(GOONTOLOGY,BP={
    go_offsp<-as.list(GOBPOFFSPRING)
  },MF={
    go_offsp<-as.list(GOMFOFFSPRING)
  },CC={
    go_offsp<-as.list(GOCCOFFSPRING)
  })
  
  go_offsp<-go_offsp[which(names(go_offsp) %in% bpnames)]
  
  go2g<-x2y_df2list(asso[,c(4,2)])
  
  go2g_full<-list()
  
  for(i in 1:length(go2g)){
    go2g_full[[i]]<-go2g[[i]]
    if(length(go_offsp[[names(go2g)[i]]])==1&&is.na(go_offsp[[names(go2g)[i]]][1])){
      
    }else{
      offsp=go_offsp[[names(go2g)[i]]]
      tmp<-Reduce(union,go2g[offsp])
      go2g_full[[i]]<-unique(union(tmp,go2g[[i]]))
    }
    go2g_full[[i]]<-go2g_full[[i]][!(is.na(go2g_full[[i]]))]
  }
  names(go2g_full)<-names(go2g)
  if(rm.termlessthan3genes==TRUE){
    go2g_full<-filterGO2Glists(go2g_full)
  }
  return(go2g_full)
}

filterGO2Glists<-function(go2g){
  p_i=1
  resultlist<-list()
  resultname<-list()
  for(i in 1:length(go2g)){
    if(length(go2g[[i]])>2){
      resultlist[[p_i]]<-go2g[[i]]
      resultname[[p_i]]<-names(go2g)[i]
      p_i=p_i+1
    }
  }
  names(resultlist)<-unlist(resultname)
  return(resultlist)
}


predict_data_size = function(numeric_size, number_type = "numeric") {
  if(number_type == "integer") {
    byte_per_number = 4
  } else if(number_type == "numeric") {
    byte_per_number = 8
  } else {
    stop(sprintf("Unknown number_type: %s", number_type))
  }
  estimate_size_in_bytes = (numeric_size * byte_per_number)
  class(estimate_size_in_bytes) = "object_size"
  print(estimate_size_in_bytes, units = "auto")
}


LLSn2List_c<-function(LLSn){
  class(LLSn[,1])<-"character"
  class(LLSn[,2])<-"character"
  
  reslist<-vector("list",length = nrow(LLSn))
  for(i in 1:nrow(LLSn)){
    reslist[[i]]<-LLSn[i,3]
  }
  names(reslist)<-apply(LLSn[,c(1,2)],1,function(v){paste(v[1],v[2],sep = "|")})
  return(reslist)
}

#' convert data.frame of HumanNet log-likelihood Score to list
#' 
#' convert HumanNet normalized log-likelihood score from data.frame to list,
#'  which will be used in FunSim method
#' @param LLSn data.frame of gene-gene normalized log-likelihood score in HumanNet
#' @return  a list of normalized log-likelihood score
#' @export
#' @author Peng Ni
#' @references Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring 
#' disease similarity by integrating semantic and gene functional association[J]. 
#' PloS one, 2014, 9(6): e99415.
#' @seealso \code{\link{FunSim}}
#' @examples
#' ## see examples in function FunSim
#' data(HumanNet_sample)
#' llsnlist<-LLSn2List(HumanNet_sample[1:100,])
#' llsnlist
LLSn2List<-function(LLSn){
  class(LLSn[,1])<-"character"
  class(LLSn[,2])<-"character"
  
  names<-unique(union(LLSn[,1],LLSn[,2]))
  
  LLSnSMat<-vector("list",length = length(names))
  names(LLSnSMat)<-names
  
  for(i in 1:length(LLSnSMat)){
    LLSnSMat[[i]]<-vector("numeric")
  }
  for(n in names){
    subsetLLSn<-rbind(LLSn[LLSn[,1]==n,],setNames(LLSn[LLSn[,2]==n,][,c(2,1,3)],names(LLSn)))
    LLSnSMat[[n]]<-subsetLLSn[,3]
    names(LLSnSMat[[n]])<-subsetLLSn[,2]
  }
  return(LLSnSMat)
}
















