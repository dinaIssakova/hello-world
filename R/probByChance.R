
getf <- function(tree, phydat, spe1, spe2, p=(1/20), type=c("abs", "score"), threshold=NA){

  alph <- Biostrings::AA_ALPHABET[1:26]
  ancNode = getMostRecentCommonAncestor(tree, spe1, spe2)
  b1 = getBranchLength(tree, spe1, ancNode)
  b2 = getBranchLength(tree, spe2, ancNode)

  total = 0
  count = 0

  if (type=="abs"){

    for (aa in alph){

      x = aa
      y = aa

      for (anc in alph[which(alph != aa)]){

        total = total + probOfSiteConfig(tree, phydat, spe1, spe2, anc=anc, x = x, y=y)
      }
    }

  } else if (type=="score"){
    for (x in alph){
      for (y in alph){
        for (anc in alph){
          count = count + 1
          data(BLOSUM62)

          score = BLOSUM62[which(colnames(BLOSUM62) == x), which(rownames(BLOSUM62 == y))]
          -  BLOSUM62[which(colnames(BLOSUM62) == x), which(rownames(BLOSUM62 == anc))]

          if (is.na(threshold)){
            print("Please supply threshold.")
            stop()
          }
          total = total + (score > threshold)

        }

      }
    }
    return (total/count)
  }

  return (total)
}

#' probOfSiteConfig
#'
#' Calculate the probability of the given configuration
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat.
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param pos The position at which to evaluate if conditions are satisfied
#' @param p The fraction of the amino acid at the target position in the evaluated genome. (Default is 1/20)
#'
#' @return Probability of given configuration
#' @examples
#' probOfSiteConfig(tree, primates, "Human", "Chimp", pos=1)
#' @export
probOfSiteConfig <- function(tree, phydat, spe1, spe2, pos, p=(1/20), anc, x, y){

  if (spe1 == spe2){
    stop("Species cannot be the same.")
  }



  if(missing(x)){
    ancM <- getAlignment(tree, phydat, spe1, spe2)
    x = toString(ancM[1][[1]][pos])
  }
  if(missing(y)){
    ancM <- getAlignment(tree, phydat, spe1, spe2)
    y = toString(ancM[2][[1]][pos])
  }

  if (missing(anc)){
    ancM <- getAlignment(tree, phydat, spe1, spe2)
    anc = toString(ancM[3][[1]][pos])
  }

  ancNode <- getMostRecentCommonAncestor(tree, spe1, spe2)

  b1 <- getBranchLength(tree, spe1, ancNode)
  b2 <- getBranchLength(tree, spe2, ancNode)

  prob1 <- probOfChange(pam, x, anc, b1)
  prob2 <- probOfChange(pam, y, anc, b2)

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
probOfNSitesByChance <- function(tree, phydat, spe1, spe2, m, n, p=(1/20), t=NA, type=c("abs", "score")){

  if (n == 0){
    return (0)
  }

  if (type == "abs"){

    f = getf(tree, phydat, spe1, spe2, p, type)
    prob = 0
    for (i in 1:n-1){

      prob = prob + ( factorial(m) / factorial(i)*factorial(m-i) ) * f**i * (1-f)**(m-i)
    }

    return (1 - prob)

  } else if (type == "score"){

    f = getf(tree, phydat, spe1, spe2, p, type="score", threshold=t)

    prob = (m - n)*(1 - f) + f*n

    return (prob)

  }
  print("Please select type: either abs or score.")
  stop()
}
