#' matrixtoAAString
#'
#' @import Biostrings
#'
#' @description Convert matrix generated during workflow of areCondSatisfied to AA sequence.
#'
#' @param m matrix to be converted
#'
#' @return AAString object
matrixtoAAString <- function(m){

  seqString = c()

  for (i in 1:nrow(m)){

    seqString <- paste(seqString, colnames(m)[which(m[i,] == 1)], sep='')

  }

  #ret = Biostrings::AAString(seqString)

  return(seqString)
}

#' Convert DNA base format
#'
#' @description Map DNA bases to a number vector with 1=A, 2=C, 3=G, 4=T.
#'
#' @param numVec The vector to be mapped.
#'
#' @return A vector of equivalent length with each value mapped to the appropriate base.
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

#' convertToAA
#'
#' @import Biostrings
#'
#' @description Convert a binary matrix representing the genomic sequence to one representing the amino acid sequence.
#'
#' @param acctranData A matrix of binary data as returned by phangorn function acctran().
#'
#' @return A matrix of amino acid data.
#'
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

      AAlist <- c(AAlist, toString(Biostrings::translate(Biostrings::DNAString(x=seq))))
    }
    for (item in codonList){
      AAmatrix[AArow,][which(colnames(AAmatrix)%in%AAlist)] = 1
    }
    AArow = AArow + 1
    DNArow = DNArow + 3
  }
  return (AAmatrix)
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
#' areCondSatisfied(smallTree, primates, "Human", "Chimp", pos=1, type="abs", 1, BLOSUM62)
#' areCondSatisfied(tree, primates, "Human", "Chimp", pos=2, type="score", threshold=0, BLOSUM62)
#' areCondSatisfied(tree, primates, "Human", "Chimp", pos=3, type="score", threshold=0, BLOSUM62)
#' }
#'
#' @import phangorn
#'
#' @export
areCondSatisfied <- function(tree, phydat, spe1, spe2, pos, type=c("abs", "score"), threshold, simMatrix=BLOSUM62){


  msaAnc <- getAlignment(tree, phydat, spe1, spe2)
  x = toString(msaAnc[1][[1]][pos])
  y = toString(msaAnc[2][[1]][pos])
  anc = toString(msaAnc[3][[1]][pos])


  if (type=="score"){
    # Get the score
    # Match score: how similar are the two AAs based on the simMatrix

    mScore <- simMatrix[which(rownames(simMatrix) == x), which(colnames(simMatrix) == y)]
    # Difference score: how different is the reference species from the ancestor
    dScore <- simMatrix[which(rownames(simMatrix) == x), which(colnames(simMatrix) == anc)]
    score = mScore - dScore

    #Condition is satisfied if the score is higher or equal to the threshold.
    cond <- score >= threshold

  } else if (type=="abs"){

    #Condition is satisfied if the two species' AA's are identical and different from the ancestral state
    # Based on:
    # Zhang, J. and Kumar, S. (1997) Detection of Convergent and Parallel Evolution at the Amino Acid Sequence
    # Level. Mol. Biol. Evol. 14(5):527-536.

    cond = ((x == y) && (x != anc) && (y != anc))

  }

  return (cond)
}

