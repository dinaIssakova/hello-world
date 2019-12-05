#' getSeq
#'
#' @description Get the corresponding species' gene sequence.
#'
#' @param tree a phylogenetic tree
#' @param phydat phydat object containing corresponding sequences
#' @param spe name of species
#' @param speNum optional: node (for ancestral reconstruction of internal nodes w/o species names)
#'
#' @return sequence
#'
#' @examples
#' getSeq(tree, primates, "Human")
#' @export

getSeq <- function(tree, phydat, spe, speNum){

  species <- tree$tip.label
  if (missing(speNum)){
    speNum = which(species == spe)
  }
  anc.acctran <- phangorn::ancestral.pars(tree, phydat, "ACCTRAN")
  # Get matrix corresponding to the sequence
  geneMat = convertToAA(anc.acctran[[speNum]])

  # Convert into AA sequence
  geneSeq = toString(matrixtoAAString(geneMat))
  return(geneSeq)
}


#' getAlignment
#'
#' @description Get MSA ClustalOmega alignment of the sequences of spe1 and spe2.
#'
#' @param tree phylogenetic tree
#' @param phydat phydat object containing corresponding sequences
#' @param spe1 name of species 1
#' @param spe2 name of species 2
#'
#' @return alignment
getAlignment <- function(tree, phydat, spe1, spe2){

  geneSeq1 <- getSeq(tree, phydat, spe1)
  geneSeq2 <- getSeq(tree, phydat, spe2)
  ancSeq <- getSeq(tree, phydat, "anc", getMostRecentCommonAncestor(tree, spe1, spe2))

  mySeq <- c(geneSeq1, geneSeq2, ancSeq)
  phyloSet1 <- Biostrings::AAStringSet(mySeq)
  msaAnc <-  msa::msaClustalOmega(phyloSet1)

  return(Biostrings::unmasked(msaAnc))
}
