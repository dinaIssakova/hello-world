
#' Calculate likihood of chance condition satisfaction
#'
#' Calculate the probability that an amino acid site will satisfy the conditions of convergent evolution by chance.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat.
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param pos The position at which to evaluate if conditions are satisfied
#' @param p The fraction of the amino acid at the target position in the evaluated genome. (Default is 1/20)
#'
#' @return Scaled likelihood that the amino acids satisfying the conditions of convergent evolution is by chance.
#' @examples
#' probOfSiteConfig(tree, primates, "Human", "Chimp", 1)
#' @export
probOfSiteConfig <- function(tree, phydat, spe1, spe2, pos, p=(1/20)){
  if (spe1 == spe2){
    stop("Species cannot be the same.")
  }

  ancNode = getMostRecentCommonAncestor(tree, spe1, spe2)

  b1 = getBranchLength(tree, spe1, ancNode)
  b2 = getBranchLength(tree, spe2, ancNode)

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- ancestral.pars(tree, phydat, "ACCTRAN")

  geneSeq1 = convertToAA(anc.acctran[[speNum1]])
  geneSeq2 = convertToAA(anc.acctran[[speNum2]])
  ancSeq = convertToAA(anc.acctran[[getMostRecentCommonAncestor(tree, spe1, spe2)]])

  gene1Val <- colnames(geneSeq1)[which(geneSeq1[pos,] == 1)]
  gene2Val <- colnames(geneSeq2)[which(geneSeq2[pos,] == 1)]
  ancVal <- colnames(ancSeq)[which(ancSeq[pos,]==1)]

  prob1 <- probOfChange(pam, gene1Val, ancVal, b1)
  prob2 <- probOfChange(pam, gene2Val, ancVal, b2)

  totalProb = p * prob1 * prob2

  return (totalProb)
}

#' Calculate probability of N chance 'convergent' sites
#'
#' Calculate the probability that n sites between the two species will satisfy the conditions of convergent evolution by chance.
#'
#' @importFrom ape is.binary
#'
#' @param tree A phylogenetic tree
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param m The length of the genes. (Must of equal length)
#' @param n The number of potential convergent sites
#' @param p The fraction of the amino acid at the target position in the evaluated genome.
#'  (Default is 1/20)
#' @return Scaled likelihood that the amino acids satisfying the conditions of convergent evolution is by chance.
#' @examples
#' \dontrun{
#' probOfNSitesByChance(tree, "Human", "Chimp", 20, 2)
#' }
#' @export
probOfNSitesByChance <- function(tree, spe1, spe2, m, n, p=(1/20)){

  if (n == 0){
    return (0)
  }

  avg=0
  for (j in 1:m){
    avg = avg + probOfSiteConfig(tree, primates, spe1, spe2, j)
  }

  f = avg/m

  prob = 0
  for (i in 1:n-1){

    prob = prob + ( factorial(m) / factorial(i)*factorial(m-i) ) * f**i * (1-f)**(m-i)
  }

  return (1 - prob)
}
