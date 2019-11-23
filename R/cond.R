
#' Convert DNA base format
#'
#' Map DNA bases to a number vector with 1=A, 2=C, 3=G, 4=T.
#'
#' @param numVec The vector to be mapped.
#'
#' @return A vector of equivalent length with each value mapped to the appropriate base.
#'

mapLetters <- function(numVec){

  retVec = c()

  for (num in 1:length(numVec)){

    if (numVec[num] == 1){
      retVec <- c(retVec, 'A')
    } else if (numVec[num] == 2){
      retVec <- c(retVec, 'C')
    } else if (numVec[num] == 3){
      retVec <- c(retVec, 'G')
    } else if (numVec[num] == 4){
      retVec <- c(retVec, 'T')
    } else{
      stop("Nucleotide must be A,C,G or T.")
    }
  }
  return (retVec)
}

#' Convert genomic to AA matrix
#'
#' Convert a binary matrix representing the genomic sequence to one representing the amino acid sequence.
#'
#' @param acctranData A matrix of binary data as returned by phangorn function acctran().
#'
#' @return A matrix of amino acid data.
#'
#' @import Biostrings

convertToAA <- function(acctranData){

  numrow = nrow(acctranData)

  while (numrow %% 3 != 0){
    numrow = numrow - 1
  }

  AAmatrix <- matrix(data=0, numrow, ncol=30)
  colnames(AAmatrix) = Biostrings::AA_ALPHABET

  DNArow = 0
  AArow = 1
  while (AArow <= numrow && DNArow+3 <= nrow(acctranData)){
    codon = acctranData[DNArow:(DNArow+3),]
    firstP = mapLetters(as.vector(which(codon[1,]!=0)))
    secondP = mapLetters(as.vector(which(codon[2,]!=0)))
    thirdP = mapLetters(as.vector(which(codon[3,]!=0)))

    codonList <- expand.grid(firstP, secondP, thirdP)
    AAlist <- c()
    for (i in nrow(codonList)){
      seq = paste(c(as.character(codonList[i,1]), as.character(codonList[i,2]), as.character(codonList[i,3])), collapse='')
      AAlist <- c(AAlist, toString(Biostrings::translate(DNAString(seq))))
    }
    for (item in codonList){
      AAmatrix[AArow,][which(colnames(AAmatrix)%in%AAlist)] = 1
    }
    AArow = AArow + 1
    DNArow = DNArow + 3
  }
  return (AAmatrix)
}

matrixtoAAString <- function(m){

  seqString = c()

  for (i in 1:nrow(m)){

    seqString <- paste(seqString, colnames(m)[which(m[i,] == 1)])

  }

  ret = Biostrings::AAString(seqString)

  return(ret)
}


#' areCondSatisfied
#'
#' @description  Determine whether conditions are satisfied for two species to be convergently evolved at a given position.
#'
#' @import msa
#' @import Biostrings
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param pos The position at which to evaluate if conditions are satisfied
#' (AA position with ref to species 1; others are aligned to species 1)
#' @param simMatrix similarity matrix to quantify similarity between amino acids (default is BLOSUM62)
#' @param threshold score threshold above which a position is considered convergent
#'
#' @return TRUE if conditions are satisfied; FALSE if they are not.
#'
#' @examples
#' \dontrun{
#' data(BLOSUM62)
#' areCondSatisfied(smallTree, primates, "Human", "Chimp", 1, 10, BLOSUM62)
#' }
#'
#' @import phangorn
#'
#' @export
areCondSatisfied <- function(tree, phydat, spe1, spe2, pos, simMatrix=BLOSUM62, threshold=1){

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- phangorn::ancestral.pars(tree, phydat, "ACCTRAN")

  geneMat1 = convertToAA(anc.acctran[[speNum1]])
  geneMat2 = convertToAA(anc.acctran[[speNum2]])

  if (geneMat1[pos,which(colnames(geneMat1) == '*')] | geneMat2[pos,which(colnames(geneMat2) == '*')]){
    print(sprintf("Cannot evaluate convergence for location of stop codon at position %i. Returning FALSE.", pos))
    return (FALSE)
  }

  geneSeq1 = matrixtoAAString(convertToAA(anc.acctran[[speNum1]]))
  geneSeq2 = matrixtoAAString(convertToAA(anc.acctran[[speNum2]]))

  phyloSet <- Biostrings::AAStringSet(list(geneSeq1, geneSeq2))

  phyloSet@ranges@NAMES <- c(spe1, spe2)

  # Align both sequences together.

  msaM <-  msa::msaMuscle(phyloSet, order = "aligned")

   # fetch the BLOSUM62 package from the Biostrings package
  msaMScores <- msa::msaConservationScore(msaM, substitutionMatrix = simMatrix)

  simscore = msaMScores[pos]

  ancSeq = matrixtoAAString(convertToAA(anc.acctran[[getMostRecentCommonAncestor(tree, spe1, spe2)]]))

  phyloSet1 <- Biostrings::AAStringSet(list(geneSeq1, ancSeq))

  phyloSet1@ranges@NAMES <- c(spe1, "ancestral")

  msaAnc <-  msa::msaMuscle(phyloSet1, order = "aligned")

  msaAncScores <- msa::msaConservationScore(msaAnc, substitutionMatrix = simMatrix)

  ancScore <- -msaAncScores[pos]

  totalscore = ancScore + simscore

  cond = totalscore >= threshold

  return (cond)
}

