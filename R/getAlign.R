#' Title
#'
#' @param tree
#' @param phydat
#' @param spe
#' @param speNum
#'
#' @return

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


#' Title
#'
#' @param tree
#' @param phydat
#' @param spe1
#' @param spe2
#'
#' @return
getAlignment <- function(tree, phydat, spe1, spe2){

  geneSeq1 <- getSeq(tree, phydat, spe1)
  geneSeq2 <- getSeq(tree, phydat, spe2)
  ancSeq <- getSeq(tree, phydat, "anc", getMostRecentCommonAncestor(tree, spe1, spe2))

  mySeq <- c(geneSeq1, geneSeq2, ancSeq)
  phyloSet1 <- Biostrings::AAStringSet(mySeq)
  msaAnc <-  msa::msaClustalOmega(phyloSet1)

  return(Biostrings::unmasked(msaAnc))
}
