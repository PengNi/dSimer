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
#' @author Peng Ni, Min Li
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
#' @param vertex.label.font label text font
#' @param vertex.label.dist label text dist
#' @param vertex.label.color label text color
#' @param vertex.label.cex label text cex
#' @param vertex.shape vertex shape
#' @param vertex.color vertex color
#' @param vertex.size vertex size
#' @param edge.color edge color
#' @param layout layout
#' @return an igraph plot object
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph graph.adjacency
#' @importFrom igraph simplify
#' @importFrom igraph induced.subgraph
#' @importFrom igraph degree
#' @export
#' @author Peng Ni, Min Li
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
                   vertex.label.color="black",
                   vertex.label.cex=0.8,
                   vertex.shape="circle",
                   vertex.color="paleturquoise",
                   vertex.size=20,
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
              vertex.label=V(g)$name,
              vertex.label.font=vertex.label.font,
              vertex.label.dist=vertex.label.dist,
              vertex.label.color=vertex.label.color,
              vertex.label.cex=vertex.label.cex,
              vertex.shape=vertex.shape,
              vertex.color=vertex.color,
              vertex.size=vertex.size,
              edge.color=edge.color,
              layout=layout)
}

## simplot author Yu Guangchuang from DOSE
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

#' plot topological relationship of two gene sets
#' 
#' plot topological relationship of two gene sets (which are associated 
#' with two diseases respectively).
#' @param geneset1 a character vector contains gene ids
#' @param geneset2 another character vector contains gene ids
#' @param graph an igraph graph object which represents a gene network
#' @param vertexcolor a character vector contains 3 colors for vertexs
#' @param vertex.shape vertex shape
#' @param vertex.size vertex size
#' @param vertex.label.font label text font
#' @param vertex.label.dist label text dist
#' @param vertex.label.color label text color
#' @param vertex.label.cex label text cex
#' @param edge.color edge color
#' @param layout layout
#' @return an igraph plot object
#' @importFrom igraph V
#' @importFrom igraph "V<-"
#' @importFrom igraph induced.subgraph
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout.auto
#' @importFrom graphics legend
#' @export
#' @author Peng Ni, Min Li
#' @examples
#' data("PPI_HPRD")
#' g<-graph.data.frame(PPI_HPRD,directed = FALSE) #get an igraph graph 
#' 
#' data(d2g_fundo_symbol)
#' a<-d2g_fundo_symbol[["DOID:8242"]] # get gene set a
#' b<-d2g_fundo_symbol[["DOID:4914"]] # get gene set b
#' 
#' plot_topo(a,b,g)
plot_topo<-function(geneset1, geneset2, graph,
                    vertexcolor=c("tomato","orange","lightsteelblue"), 
                    vertex.shape="circle",
                    vertex.size=14,
                    vertex.label.font=1,
                    vertex.label.dist=0,
                    vertex.label.color='black',
                    vertex.label.cex=0.5,
                    edge.color="black",
                    layout=layout.auto){
  geneset1<-intersect(geneset1, V(graph)$name)
  geneset2<-intersect(geneset2, V(graph)$name)
  stopifnot(length(geneset1)>0 & length(geneset2)>0)
  
  g<-induced.subgraph(graph,union(geneset1, geneset2))
  shareg<-intersect(geneset1, geneset2)
  diffg1<-setdiff(geneset1, geneset2)
  diffg2<-setdiff(geneset2, geneset1)
  gname<-c(shareg, diffg1, diffg2)
  V(g)[gname]$type<-c(rep(1, length(shareg)), 
                      rep(2, length(diffg1)), 
                      rep(3, length(diffg2)))
  vertex.color<-vertexcolor[V(g)$type]
  plot.igraph(g, 
              vertex.color=vertex.color, 
              vertex.shape=vertex.shape, 
              vertex.size=vertex.size, 
              vertex.label.font=vertex.label.font, 
              vertex.label.dist=vertex.label.dist, 
              vertex.label.color=vertex.label.color, 
              vertex.label.cex=vertex.label.cex, 
              edge.color=edge.color, 
              layout=layout, 
              vertex.shape=vertex.shape,
              layout=layout)
  legend(x=-1.1, y=-1.0, 
         legend = c("shared genes","geneset 1", "geneset2"), 
         pch = 21, pt.bg=vertexcolor, pt.cex=2, cex=.8,bty="n", ncol=3)
}

#' plot disease-gene (or GO term etc.) associations as a bipartite graph
#' 
#' plot a bipartite graph which visualizes associations between diseases and genes 
#' (or GO terms etc.)
#' @param xylist a named list object which names are diseases and each element of the
#' list is a gene set with respect to each disease.
#' @param vertex.size vertex size
#' @param vertex.shape1 shape for one kind of vertex
#' @param vertex.shape2 shape for another kind of vertex
#' @param vertex.color1 color for one kind of vertex
#' @param vertex.color2 color for another kind of vertex
#' @param vertex.label.font label text font
#' @param vertex.label.dist label text dist
#' @param vertex.label.color label text color
#' @param vertex.label.cex label text cex
#' @param edge.color edge color
#' @param layout layout
#' @return an igraph plot object
#' @importFrom igraph V
#' @importFrom igraph "V<-"
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout.kamada.kawai
#' @importFrom igraph graph_from_incidence_matrix
#' @export
#' @author Peng Ni, Min Li
#' @examples
#' data(d2g_fundo_symbol)
#' d2g_sample<-sample(d2g_fundo_symbol, 3)
#' plot_bipartite(d2g_sample)
plot_bipartite<-function(xylist,
                         vertex.size=12,
                         vertex.shape1="circle",
                         vertex.shape2="square",
                         vertex.color1="darkseagreen", 
                         vertex.color2="turquoise1", 
                         vertex.label.font=2,
                         vertex.label.dist=0,
                         vertex.label.color='black',
                         vertex.label.cex=0.8,
                         edge.color="black",
                         layout=layout.kamada.kawai){
  llen<-length(xylist)
  eles<-unique(unlist(xylist))
  elelen<-length(eles)
  mat<-matrix(data = 0, nrow = llen, 
              ncol = elelen, 
              dimnames = list(names(xylist), eles))
  for(i in 1:llen){
    ele<-xylist[[i]]
    lname<-names(xylist[i])
    for(e in ele){
      mat[lname, e]=1
    }
  }
  g<-graph_from_incidence_matrix(mat, weighted=TRUE)
  vertex_shape<-c(vertex.shape1, vertex.shape2)[1+(V(g)$type==FALSE)]
  vertex_color<-c(vertex.color1, vertex.color2)[1+(V(g)$type==FALSE)]
  plot.igraph(g, 
              vertex.label.font=vertex.label.font, 
              vertex.label.dist=vertex.label.dist, 
              vertex.label.color=vertex.label.color, 
              vertex.label.cex=vertex.label.cex, 
              edge.color=edge.color, 
              layout=layout, 
              vertex.shape=vertex_shape, 
              vertex.color=vertex_color, 
              vertex.size=vertex.size)
}
