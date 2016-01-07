#' Integration of Disease Similarity Methods
#' 
#' dSimer is an R package which provides eight function-based methods 
#' for disease similarity calculation. These eight methods take advantage of 
#' diverse biological data which calculate disease similarity from different 
#' perspectives. The disease similarity matrix obtained from these eight methods 
#' can also be visualized by dSimer.
#'
#' \tabular{ll}{ Package: \tab dSimer\cr Type: \tab Package\cr Version: \tab
#' 0.99.0\cr Date: \tab 12-10-2015\cr biocViews:\tab Bioinformatics,
#' Similarity\cr Depends: \tab R (>= 3.2.0)\cr Imports: \tab Rcpp (>= 0.11.3), 
#' igraph (>= 1.0.1), ggplot2, reshape2, org.Hs.eg.db, GO.db, AnnotationDbi\cr 
#' Suggests: \tab knitr, rmarkdown, BiocStyle\cr LinkingTo: \tab Rcpp\cr License: \tab
#' GPL (>= 2)\cr }
#'
#' @name dSimer-package
#' @aliases dSimer-package dSimer
#' @docType package
#' @author Peng Ni, Min Li
#' @keywords package
NULL

#' disease-gene associations
#'
#' a list of disease-gene associations from FunDO which diseases are DOIDs and
#' genes are entrez ids
#'
#' @name d2g_fundo_entrezid
#' @aliases d2g_fundo_entrezid
#' @docType data
#' @keywords dataset
#' @references Osborne J D, Flatow J, Holko M, et al. Annotating the human genome 
#' with Disease Ontology[J]. BMC genomics, 2009, 10(Suppl 1): S6.
NULL

#' disease-gene associations
#' 
#' a list of disease-gene associations from FunDO which diseases are DOIDs and
#' genes are gene symbols
#' 
#' @name d2g_fundo_symbol
#' @aliases d2g_fundo_symbol
#' @docType data
#' @keywords dataset
#' @references Osborne J D, Flatow J, Holko M, et al. Annotating the human genome 
#' with Disease Ontology[J]. BMC genomics, 2009, 10(Suppl 1): S6.
NULL

#' disease-gene associations
#' 
#' a list of disease-gene associations from the reference paper (see below) 
#' which diseases are disease names and genes are entrez ids
#' @name d2g_separation
#' @aliases d2g_separation
#' @docType data
#' @keywords dataset
#' @references Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease 
#' relationships through the incomplete interactome[J]. Science, 2015, 
#' 347(6224): 1257601.
NULL


#' a sample of GO term-gene associations
#' 
#' a sample list of GO term-gene associations, the entire data of GO term-gene
#' can be obtained by function get_GOterm2GeneAssos
#' @name go2g_sample
#' @aliases go2g_sample
#' @docType data
#' @keywords dataset
#' @seealso \code{\link{get_GOterm2GeneAssos}}
NULL

#' a sample of disease-GO term associations
#' 
#' a sample list of disease-GO term associations, the entire data of disease-
#' GO term associations can be obtained by function HypergeometricTest
#' @name d2go_sample
#' @aliases d2go_sample
#' @docType data
#' @keywords dataset
#' @seealso \code{\link{HypergeometricTest}}
NULL

#' a sample of HumanNet likelihood score data, will be used in method FunSim
#' 
#' data.frame, where each row indicates a pair of genes and 
#' their normalized likelihood score in HumanNet, will be used in
#' method FunSim after being converted to list by method LLSn2List
#' @name HumanNet_sample
#' @aliases HumanNet_sample
#' @docType data
#' @keywords dataset
#' @references Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring 
#' disease similarity by integrating semantic and gene functional association[J]. 
#' PloS one, 2014, 9(6): e99415.
#' @seealso \code{\link{FunSim}}, \code{\link{LLSn2List}}
NULL

#' PPI data from HPRD
#' 
#' a data.frame, fetched from HPRD which genes are gene symbols
#' @name PPI_HPRD
#' @aliases PPI_HPRD
#' @docType data
#' @keywords dataset
NULL

#' interactome data
#' 
#' a data.frame, contains interactome data in reference below
#' @name interactome
#' @aliases interactome
#' @docType data
#' @keywords dataset
#' @references Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease 
#' relationships through the incomplete interactome[J]. Science, 2015, 
#' 347(6224): 1257601.
NULL

#' graphlet signature
#' 
#' matrix, graphlet signatures of nodes in HPRD PPI network calculated by 
#' ORCA tool, will be using in method Sun_topology
#' @name graphlet_sig_hprd
#' @aliases graphlet_sig_hprd
#' @docType data
#' @keywords dataset
#' @references Hocevar T, Demsar J. A combinatorial approach to graphlet 
#' counting[J]. Bioinformatics, 2014, 30(4): 559-565.
#' @seealso \code{\link{Sun_topology}}
NULL

#' orbit dependency count
#' 
#' a 73-dim vector, indicating 73 orbits' dependency count in graphlet theory, 
#' used to calculate weight factor in method setWeight
#' @name orbit_dependency_count
#' @aliases orbit_dependency_count
#' @docType data
#' @keywords dataset
#' @references Milenkovic T, Przulj N. Uncovering biological network function via 
#' graphlet degree signatures[J]. Cancer informatics, 2008, 6: 257.
#' @seealso \code{\link{setWeight}}
NULL

#' weight factor 
#' 
#' a 73-dim vector, indicating 73 orbits' weight factor, will be used in
#' method Sun_topology
#' @name weight
#' @aliases weight
#' @docType data
#' @keywords dataset
#' @references Sun K, Goncalves JP, Larminie C. Predicting disease associations 
#' via biological network analysis[J]. BMC bioinformatics, 2014, 15(1): 304.
#' @seealso \code{\link{setWeight}}, \code{\link{Sun_topology}}
NULL



