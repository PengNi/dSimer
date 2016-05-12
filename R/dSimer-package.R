#' Integration of Disease Similarity Methods
#' 
#' dSimer is an R package which provides eight function-based methods 
#' for disease similarity calculation. These eight methods take advantage of 
#' diverse biological data which calculate disease similarity from different 
#' perspectives. The disease similarity matrix obtained from these eight methods 
#' can also be visualized by dSimer.
#'
#' \tabular{ll}{ Package: \tab dSimer\cr Type: \tab Package\cr Version: \tab
#' 0.99.2\cr Date: \tab 12-10-2015\cr biocViews:\tab Bioinformatics,
#' Similarity\cr Depends: \tab R (>= 3.3.0), igraph (>= 1.0.1)\cr Imports: \tab stats, Rcpp (>= 0.11.3), 
#' ggplot2, reshape2, GO.db, AnnotationDbi, org.Hs.eg.db\cr 
#' Suggests: \tab knitr, rmarkdown, BiocStyle\cr LinkingTo: \tab Rcpp\cr License: \tab
#' GPL (>= 2)\cr }
#'
#' @name dSimer-package
#' @aliases dSimer-package dSimer
#' @docType package
#' @author Peng Ni, Min Li
#' @keywords package
NULL

#' d2g_fundo_entrezid
#'
#' a list of disease-gene associations from FunDO.
#' @name d2g_fundo_entrezid
#' @aliases d2g_fundo_entrezid
#' @docType data
#' @keywords dataset
#' @return d2g_fundo_entrezid is a named list of length 1855 which stored disease-gene 
#' associations from FunDO. The names are the DOIDs (DOIDs are ids of terms 
#' in Disease Ontology, e.g. "DOID:4" ) and list elements are vectors of Entrez gene IDs.
#' @examples
#' data(d2g_fundo_entrezid)
#' @references Osborne J D, Flatow J, Holko M, et al. Annotating the human genome 
#' with Disease Ontology[J]. BMC genomics, 2009, 10(Suppl 1): S6.
NULL

#' d2g_fundo_symbol
#' 
#' a list of disease-gene associations from FunDO.
#' @name d2g_fundo_symbol
#' @aliases d2g_fundo_symbol
#' @docType data
#' @keywords dataset
#' @return d2g_fundo_symbol is a named list of length 1855 which stored disease-gene 
#' associations from FunDO. The names are the DOIDs (DOIDs are ids of terms 
#' in Disease Ontology, e.g. "DOID:4" ) and list elements are vectors of gene symbols.
#' @examples
#' data(d2g_fundo_symbol)
#' @references Osborne J D, Flatow J, Holko M, et al. Annotating the human genome 
#' with Disease Ontology[J]. BMC genomics, 2009, 10(Suppl 1): S6.
NULL

#' d2g_separation
#' 
#' a list of disease-gene associations from the reference paper (see below).
#' @name d2g_separation
#' @aliases d2g_separation
#' @docType data
#' @keywords dataset
#' @return d2g_separation is a named list of length 299 which stored disease-gene
#' associations from the reference paper (see below). The names are diseases and 
#' list elements are vectors of gene entrez ids.
#' @examples 
#' data(d2g_separation)
#' @references Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease 
#' relationships through the incomplete interactome[J]. Science, 2015, 
#' 347(6224): 1257601.
NULL


#' go2g_sample
#' 
#' a sample list of GO term-gene associations.
#' @name go2g_sample
#' @aliases go2g_sample
#' @docType data
#' @keywords dataset
#' @return go2g_sample is a named list of length 465. The names are GO term ids (GOIDs)
#' and list elements are vectors of gene symbols. 
#' The entire data of GO term-gene assos can be obtained by function get_GOterm2GeneAssos.
#' @examples 
#' data(go2g_sample)
#' @seealso \code{\link{get_GOterm2GeneAssos}}
NULL

#' d2go_sample
#' 
#' a sample list of disease-GO term associations.
#' @name d2go_sample
#' @aliases d2go_sample
#' @docType data
#' @keywords dataset
#' @return d2go_sample is a named list of length 3. The names are are the DOIDs 
#' (DOIDs are ids of terms in Disease Ontology, e.g. "DOID:4" ) and list 
#' elements are vectors of GO term ids.
#' The entire data of disease-GO term associations can be obtained by function HypergeometricTest.
#' @examples 
#' data(d2go_sample)
#' @seealso \code{\link{HypergeometricTest}}
NULL

#' HumanNet_sample
#' 
#' a sample of HumanNet likelihood score data which will be used in method FunSim.
#' @name HumanNet_sample
#' @aliases HumanNet_sample
#' @docType data
#' @keywords dataset
#' @return HumanNet_sample is a data.frame has 22708 rows and 3 columns. Each row 
#' indicates a pair of genes and their normalized likelihood score in HumanNet. 
#' HumanNet_sample will be used in method FunSim after being converted to list 
#' by method LLSn2List.
#' The entire data of HumanNet can be downloaded from 
#' the website http://www.functionalnet.org/humannet/ .
#' @examples 
#' data(HumanNet_sample)
#' @references Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring 
#' disease similarity by integrating semantic and gene functional association[J]. 
#' PloS one, 2014, 9(6): e99415.
#' @seealso \code{\link{FunSim}}, \code{\link{LLSn2List}}
NULL

#' PPI_HPRD
#' 
#' PPI data from HPRD
#' @name PPI_HPRD
#' @aliases PPI_HPRD
#' @docType data
#' @keywords dataset
#' @return PPI_HPRD is a data.frame of 36867 rows and 2 columns. Each rows indicates
#' an interaction of two gene symbols. It was fetched from HPRD.
#' @examples 
#' data(PPI_HPRD)
#' @references Prasad T S K, Goel R, Kandasamy K, et al. Human protein reference 
#' database-2009 update[J]. Nucleic acids research, 2009, 37(suppl 1): D767-D772.
NULL

#' interactome
#' 
#' interactome data
#' @name interactome
#' @aliases interactome
#' @docType data
#' @keywords dataset
#' @return interactome is a data.frame of 141296 rows and 2 columns. Each row indicates
#' an interaction of two gene entrez ids.
#' It was obtained from the reference below.
#' @examples 
#' data(interactome)
#' @references Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease 
#' relationships through the incomplete interactome[J]. Science, 2015, 
#' 347(6224): 1257601.
NULL

#' graphlet_sig_hprd
#' 
#' graphlet signature of nodes in HPRD PPI network.
#' @name graphlet_sig_hprd
#' @aliases graphlet_sig_hprd
#' @docType data
#' @keywords dataset
#' @return #' graphlet_sig_hprd is a matrix of 9270 rows and 73 rows. The rownames
#' of graphlet_sig_hprd are gene symbols of nodes from HPRD. Each row indicates a graphlet
#' signature of one node.
#' Graphlet signatures of nodes in HPRD PPI network were calculated 
#' by ORCA tool, will be used in method Sun_topology.
#' @examples 
#' data(graphlet_sig_hprd)
#' @references Hocevar T, Demsar J. A combinatorial approach to graphlet 
#' counting[J]. Bioinformatics, 2014, 30(4): 559-565.
#' @seealso \code{\link{Sun_topology}}
NULL

#' orbit_dependency_count
#' 
#' orbit dependency count
#' @name orbit_dependency_count
#' @aliases orbit_dependency_count
#' @docType data
#' @keywords dataset
#' @return orbit_dependency_count is a 73-dim vector, indicating 73 orbits' 
#' dependency count in graphlet theory, used to calculate weight factor in 
#' method setWeight.
#' @examples 
#' data(orbit_dependency_count)
#' @references Milenkovic T, Przulj N. Uncovering biological network function via 
#' graphlet degree signatures[J]. Cancer informatics, 2008, 6: 257.
#' @seealso \code{\link{setWeight}}
NULL

#' weight
#' 
#' weight factor 
#' @name weight
#' @aliases weight
#' @docType data
#' @keywords dataset
#' @return weight is a 73-dim vector, indicating 73 orbits' weight factor, will be used in
#' method Sun_topology.
#' @examples 
#' data(weight)
#' @references Sun K, Goncalves JP, Larminie C. Predicting disease associations 
#' via biological network analysis[J]. BMC bioinformatics, 2014, 15(1): 304.
#' @seealso \code{\link{setWeight}}, \code{\link{Sun_topology}}
NULL



