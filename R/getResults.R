
#' Find convergent species
#'
#' @description For a given species, find others in the tree that satisfy the conditions of convergent evolution at a certain position with a less than 0.05 chance of this having occured by chance.
#'
#' @import Biostrings
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe The species to compare others in the tree to
#' @param pos The position at which to compare the genes
#' @param type Type of analysis: 'abs' for basic model or 'score' for by convergence score model
#' @param t threshold
#'
#' @return A vector of species names which satisfy the conditions listed above.
#' @examples
#' \dontrun{
#' getConvergent(tree, primates, "Human", 1, 5)
#' }
#' @export
getConvergent <- function(tree, phydat, spe, pos, type=c("abs, score"), t=NA){
  species <- tree$tip.label
  convSpe <- c(spe)
  for (s in species){
      #library(Biostrings)
      #data(BLOSUM62)

      #Is this species convergent with spe?
      cond = areCondSatisfied(tree, phydat, s, spe, pos, type, threshold=t, BLOSUM62)

    if (s != spe && cond){
      # If it is (and it isn't spe)
      # Find the probability that this is by chance
      p = probOfSiteConfig(tree, primates, s, spe, pos)
      print(sprintf("Species %s is potentially convergent with p %e.", s, p))
      if (p < 0.05){
        convSpe = c(convSpe, s)
      }
    }
  }

  return (convSpe)
}

#' Get gene length
#'
#' Get the length of two genes. (Genes must be of equal length.)
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @return The length of each gene.
#'
getm <- function(tree, phydat, spe1, spe2){

  # Get length of gene sequences being compared.

  geneSeq1 <- getSeq(tree, phydat, spe1)
  geneSeq2 <- getSeq(tree, phydat, spe2)

  if (width(geneSeq1) != width(geneSeq2)){
    print("Unequal gene lengths")
    stop()
  } else {
    return (width(geneSeq1))
  }
}

#' Get number of convergent sites
#'
#' Get the number of potentially evolved convergent sites and print the probability that this occured by chance.
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param m The length of each gene (Up to what is desired to be evaluated). Default is entire gene
#' @param t threshold
#' @param type Type of analysis: 'abs' for basic model or 'score' for by convergence score model
#'
#' @return The number of potentially convergent sites
#' @examples
#' \dontrun{
#' convSiteData(smallTree, primates, "Human", "Chimp", 5)
#' }
#' @export
convSiteData <- function(tree, phydat, spe1, spe2, t, type=c("abs", "score"), m=getm(tree, phydat, spe1, spe2)){

  numSites = 0
  for (i in 1:m){
    # Get number of convergent sites
    if (areCondSatisfied(tree, phydat, spe1, spe2, i, type=type, threshold=t)) {
      numSites = numSites + 1
      }
  }
  print(sprintf("%i potentially convergent sites with a %e probability of this occuring by chance.", numSites, probOfNSitesByChance(tree, spe1, spe2, m, n=numSites, t=t, type=type)))

  return (numSites)

}


