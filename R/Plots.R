#' similarity matrix heatmap plotting
#' 
#' plot heatmap of a disease similarity matrix
#' @param simmat a similarity matrix
#' @param xlab xlab
#' @param ylab ylab
#' @param color.low color of low value
#' @param color.high color of high value
#' @param labs logical, add text label or not
#' @param digits round digit numbers
#' @param labs.size lable size
#' @param font.size font size
#' @return a ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 scale_fill_gradient
##' @importFrom ggplot2 scale_x_discrete
##' @importFrom ggplot2 scale_y_discrete
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 %+replace%
##' @importFrom reshape2 melt
#' @export
#' @author Peng Ni
#' @references Yu G, Wang L G, Yan G R, et al. DOSE: an R/Bioconductor package for 
#' disease ontology semantic and enrichment analysis[J]. 
#' Bioinformatics, 2015, 31(4): 608-609.
#' @examples
#' 
#' data(d2g_separation)
#' data(interactome)
#' 
#' graph_interactome<-graph.data.frame(interactome,directed=FALSE)

#' ds<-c("myocardial ischemia","myocardial infarction","coronary artery disease",
#'  "cerebrovascular disorders","arthritis, rheumatoid","diabetes mellitus, type 1",
#'  "autoimmune diseases of the nervous system","demyelinating autoimmune diseases, cns",
#'  "respiratory hypersensitivity","asthma","retinitis pigmentosa",
#'  "retinal degeneration","macular degeneration")
#'  
#' sep<-Separation(ds,ds,d2g_separation,graph_interactome)
#' sim<-Separation2Similarity(sep)
#' plot_heatmap(sim)
plot_heatmap<-function(simmat,xlab="", ylab="", 
                       color.low="white", color.high="red", labs=TRUE, 
                       digits=2, labs.size=3, font.size=14){
  simplot(sim = simmat,xlab, ylab, color.low, color.high, labs, 
          digits, labs.size, font.size, readable=FALSE )
}

#' plot a network based on a symmetric disease similarity matrix
#'
#' plot a network/graph of a symmetric disease similarity matrix, note that a unsymmetric
#' matrix can't be visualized into a network by this method.
#' @param simmat a symmetric similarity matrix
#' @param cutoff a cutoff value, only disease pairs have similarity scores 
#' no less than cutoff will be visualized in the network
#' @param vertex.label.font font size
#' @param vertex.label.dist font dist
#' @param vertex.label.color font text color
#' @param vertex.label.cex cex of vertex label
#' @param vertex.shape vertex shape
#' @param vertex.color vertex color
#' @param vertex.size vertex size
#' @param vertexsizefold a fold of vertex size
#' @param edge.color edge color
#' @param layout layout
#' @return a igraph plot object
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph graph.adjacency
#' @importFrom igraph simplify
#' @importFrom igraph induced.subgraph
#' @importFrom igraph degree
#' @export
#' @author Peng Ni
#' @examples
#' 
#' data(d2g_separation)
#' data(interactome)
#' 
#' graph_interactome<-graph.data.frame(interactome,directed=FALSE)

#' ds<-c("myocardial ischemia","myocardial infarction","coronary artery disease",
#'  "cerebrovascular disorders","arthritis, rheumatoid","diabetes mellitus, type 1",
#'  "autoimmune diseases of the nervous system","demyelinating autoimmune diseases, cns",
#'  "respiratory hypersensitivity","asthma","retinitis pigmentosa",
#'  "retinal degeneration","macular degeneration")
#'  
#' sep<-Separation(ds,ds,d2g_separation,graph_interactome)
#' sim<-Separation2Similarity(sep)
#' plot_net(sim,cutoff=0.2)
plot_net<-function(simmat,
                   cutoff=1,
                   vertex.label.font=2,
                   vertex.label.dist=0.5,
                   vertex.label.color='black',
                   vertex.label.cex=0.8,
                   vertex.shape="circle",
                   vertex.size=1,
                   vertexsizefold=5,
                   vertex.color="paleturquoise",
                   edge.color="red",
                   layout=layout.fruchterman.reingold){

  
  simmat<-as.matrix(simmat)
  simmat<-simmat[order(row.names(simmat)),order(colnames(simmat)),drop=FALSE]
  
  stopifnot(isSymmetric(simmat))
  
  diag(simmat)<-0
  
  #cutoff
  stopifnot(cutoff<=max(simmat) & cutoff>0)
  for(i in 1:nrow(simmat)){
    for(j in 1:ncol(simmat)){
      if(simmat[i,j]>=cutoff){
        simmat[i,j]=1
      }else{
        simmat[i,j]=0
      }
    }
  }
  
  g<-graph.adjacency(simmat,mode = "undirected",weighted = NULL,diag = FALSE)
  g<-simplify(g)
  dgee<-degree(g)
  g<-induced.subgraph(g,names(dgee[dgee>0]))
  
  plot.igraph(g,
              vertex.size=vertexsizefold*degree(g),
              vertex.label=V(g)$name,
              vertex.label.font=vertex.label.font,
              vertex.label.dist=vertex.label.dist,
              vertex.label.color=vertex.label.color,
              vertex.label.cex=vertex.label.cex,
              vertex.shape=vertex.shape,
              vertex.color=vertex.color,
              edge.color=edge.color,
              layout=layout)
}

## simplot author Yu Guangchuang
simplot <- function(sim, xlab="", ylab="", color.low="white", color.high="red", labs=TRUE, digits=2, labs.size=3, font.size=14, readable=FALSE ) {
  sim.df <- as.data.frame(sim)
  if(readable == TRUE) {
    ##rownames(sim.df) <- TERM2NAME.DO(rownames(sim.df))
    ##colnames(sim.df) <- TERM2NAME.DO(colnames(sim.df))
  }
  rn <- row.names(sim.df)
  
  sim.df <- cbind(ID=rownames(sim.df), sim.df)
  sim.df <- melt(sim.df)
  
  sim.df[,1] <- factor(sim.df[,1], levels=rev(rn))
  if (labs == TRUE) {
    ## lbs <- c(apply(round(sim, digits), 2, as.character))
    sim.df$label <- as.character(round(sim.df$value, digits))
  }
  variable <- ID <- value <- label <- NULL ## to satisfy codetools
  if (labs == TRUE)
    p <- ggplot(sim.df, aes(variable, ID, fill=value, label=label))
  else
    p <- ggplot(sim.df, aes(variable, ID, fill=value))
  
  p <- p + geom_tile(color="black")+
    scale_fill_gradient(low=color.low, high=color.high) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))+
    theme(axis.ticks=element_blank())
  if (labs == TRUE)
    p <- p+geom_text(size=labs.size)
  p <- p+theme_dose(font.size)
  p <- p + theme(axis.text.x=element_text(hjust=0, angle=-90)) +
    theme(axis.text.y=element_text(hjust=0))
  p <- p+theme(legend.title=element_blank())
  ##geom_point(aes(size=value))
  p <- p+xlab(xlab)+ylab(ylab)
  
  if (readable == TRUE) {
    p <- p + theme(axis.text.y = element_text(hjust=1))
  }
  p <- p + theme(axis.text.x = element_text(vjust=0.5))
  return(p)
}

##  ggplot theme of DOSE
##
##  theme_dose
##  param font.size font size
##  author Yu Guangchuang
theme_dose <- function(font.size=14) {
  theme_bw() %+replace%
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.title.x = element_text(colour="black",
                                      size = font.size),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title.y = element_text(colour="black",
                                      size = font.size, angle=90)
    )
}








