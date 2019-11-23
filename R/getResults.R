
#' Find convergent species
#'
#' For a given species, find others in the tree that satisfy the conditions of convergent evolution at a certain position with a less than 0.05 chance of this having occured by chance.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe The species to compare others in the tree to
#' @param pos The position at which to compare the genes
#'
#' @return A vector of species names which satisfy the conditions listed above.
#' @examples
#' \dontrun{
#' getConvergent(tree, primates, "Human", 1)
#' }
#' @export
getConvergent <- function(tree, phydat, spe, pos){
  species <- tree$tip.label
  convSpe <- c(spe)
  for (s in species){
      cond = areCondSatisfied(tree, phydat, s, spe, pos)


    if (s != spe && cond){

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
#' getm(tree, primates, "Mouse", "Bovine")
getm <- function(tree, phydat, spe1, spe2){

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- ancestral.pars(tree, phydat, "ACCTRAN")

  geneSeq1 = convertToAA(anc.acctran[[speNum1]])
  geneSeq2 = convertToAA(anc.acctran[[speNum2]])

  if (nrow(geneSeq1) != nrow(geneSeq2)){
    print("Unequal gene lengths")
    stop()
  } else {
    return (nrow(geneSeq1))
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
#'
#' @return The number of potentially convergent sites
#' @examples
#' \dontrun{
#' convSiteData(smallTree, primates, "Human", "Chimp", 5)
#' }
#' @export
convSiteData <- function(tree, phydat, spe1, spe2, m=getm(tree, phydat, spe1, spe2)){

  numSites = 0
  for (i in 1:m){
    if (areCondSatisfied(tree, phydat, spe1, spe2, i)) {
      numSites = numSites + 1
      }
  }
  print(sprintf("%i potentially convergent sites with a %e probability of this occuring by chance.", numSites, probOfNSitesByChance(tree, spe1, spe2, m, n=numSites)))

  return (numSites)

}


