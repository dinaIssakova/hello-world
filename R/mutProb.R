

#' Calculate mutation probability over 1 PAM
#'
#' Calculate the probability of one amino acid mutating into another over a distance of 1 PAM.
#'
#' @param pam The PAM matrix (given in ./data)
#' @param AA1 The origin amino acid
#' @param AA2 The target amino acid
#'
#' @return A probability of mutation.
#' probOfChange(pam, "L", "K")
probOfChange1PAM <- function(pam, AA1, AA2){

  i <- which(rownames(pam) == AA1)
  j <- which(rownames(pam) == AA2)

  probability = pam[i, j]/10000
  return (probability)

}

#' Calculate mutation probability
#'
#' Calculate the probability of one amino acid mutating into another over a given branch length.
#'
#' @import foreach
#'
#' @param pam The PAM matrix given in ./data.
#' @param AA1 The origin amino acid
#' @param AA2 The target amino acid
#' @param d The branch length between the two sequences.
#'
#' @return A scaled probability of mutation.
#' @examples
#' probOfChange(pam, "L", "K", 3)
#' @export
probOfChange <- function(pam, AA1, AA2, d){
  AAlist <- Biostrings::AA_ALPHABET[-c(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)]
  if (!(AA1%in%AAlist) | !(AA2%in%AAlist)){
    return(0)
  }

  prob = 0
  if (d == 1){
    prob = probOfChange1PAM(pam, AA1, AA2)
    return(prob)
  }
  foreach::foreach (j=1:length(AAlist)) %do% {
    aa = AAlist[[j]]
    prob = prob + (probOfChange1PAM(pam, AA1, aa) * probOfChange(pam, aa, AA2, d-1))
  }
  return(prob)
}
