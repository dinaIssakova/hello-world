#' Get length of branch in tree.
#'
#' Find the (rounded) length of a branch from a given species to a given internal node.
#'
#' @param tree A phylogenetic tree
#' @param spe The name of the species
#' @param nodeNum The number of the internal node desired.
#' @return The length of the branch between nodeNum and spe as given by tree.
#' @examples
#' getBranchLength(tree, "Mouse", 15)
#' @export
getBranchLength<- function(tree, spe, nodeNum){
  leafNum <- which(tree$tip.label == spe)
  if (tree$edge.length[1]> 10){
    tree$edge.length = tree$edge.length/10
  }
  if(leafNum == nodeNum){
    return(0)
  }

  rowNum = which(tree$edge[ ,2] == leafNum)
  ancNum = as.numeric((tree$edge[rowNum, ])[1])

  len = 0
  len = len + tree$edge.length[rowNum]

  while(!is.na(ancNum) && ancNum != nodeNum){
    leafNum = ancNum
    rowNum = which(tree$edge[ ,2] == leafNum)
    ancNum = as.numeric((tree$edge[rowNum, ])[1])
    len = len + tree$edge.length[rowNum]
  }
  return(round(len))
}

#' Get most recent common ancestor
#'
#' Get the most recent common ancestor between two nodes as given by the given tree.
#'
#' @param tree A phylogenetic tree
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @return The number of the internal node that is the most recent point of divergence.
#' @examples
#' getMostRecentCommonAncestor(tree, "Mouse", "Bovine")
#' @export
getMostRecentCommonAncestor<- function(tree, spe1, spe2){

  if(!is.recursive(tree)){
    print(tree)
    stop()
  }

  leafNum1 <- which(tree$tip.label == spe1)
  leafNum2 <- which(tree$tip.label == spe2)

  if (spe1 == spe2){
    return (leafNum1)
  }

  visited <- c()
  row1 = tree$edge[which(tree$edge[,2] == leafNum1),]
  anc1 = as.numeric(row1[1])

  while (!is.na(anc1)){

    visited <- c(visited, anc1)
    row1 = tree$edge[which(tree$edge[,2] == anc1),]
    anc1 = as.numeric(row1[1])
  }

  row2 = tree$edge[which(tree$edge[,2] == leafNum2),]
  anc2 = as.numeric(row2[1])
  nodeNum = leafNum2

  while (!is.na(anc2)){
    if (anc2%in%visited){
      return (anc2)
    }
    row2 = tree$edge[which(tree$edge[,2] == anc2),]
    nodeNum = anc2
    anc2 = as.numeric(row2[1])
  }
  if (is.na(anc2)){
    return (nodeNum)
  } else {
    stop("Tree parsed incorrectly. Please check that input is a complete binary tree.")
  }
}

#' Convert DNA base format
#'
#' Map DNA bases to a number vector with 1=A, 2=C, 3=G, 4=T.
#'
#' @param numVec The vector to be mapped.
#' @return A vector of equivalent length with each value mapped to the appropriate base.
#' @examples
#' mapLetters(c(1,2,3,4))
#' @export
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

#' Determine if convergently evolved
#'
#' Determine whether conditions are satisfied for two species to be convergently evolved at a given position.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param pos The position at which to evaluate if conditions are satisfied
#' @param anc Optional: the value of the ancestral amino acid at position pos (if not supplied, will be predicted using ancestral reconstruction)
#' @return TRUE if conditions are satisfied; FALSE if they are not.
#' @examples
#' areCondSatisfied(tree, primates, "Mouse", "Bovine", 17)
#' areCondSatisfied(tree, primates, "Mouse", "Bovine", 1, "L")
#' @export
areCondSatisfied <- function(tree, phydat, spe1, spe2, pos, anc){

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- ancestral.pars(tree, phydat, "ACCTRAN")

  geneSeq1 = convertToAA(anc.acctran[[speNum1]])
  geneSeq2 = convertToAA(anc.acctran[[speNum2]])

  gene1Val <- which(geneSeq1[pos,] == 1)
  gene2Val <- which(geneSeq2[pos,] == 1)

  if (length(gene1Val) == 0 | length(gene2Val) == 0){
    print(sprintf("No AA at location %d.", pos))
    return (FALSE)
  }
  if (missing(anc)){
    ancSeq = convertToAA(anc.acctran[[getMostRecentCommonAncestor(tree, spe1, spe2)]])
    ancVal = ancSeq[pos,]
    cond2 = (ancSeq[pos,][gene1Val] == 0)
    cond3 = (ancSeq[pos,][gene2Val] == 0)
  } else {
    ancVal = anc
    cond2 = (anc != gene1Val)
    cond3 = (anc != gene2Val)
  }

  cond1 = (gene1Val==gene2Val)

  return (cond1 && cond2 && cond3)
}

#' Convert genomic to AA matrix
#'
#' Convert a binary matrix representing the genomic sequence to one representing the amino acid sequence.
#'
#' @param acctranData A matrix of binary data as returned by phangorn function acctran().
#' @return A matrix of amino acid data.

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

#' Calculate mutation probability over 1 PAM
#'
#' Calculate the probability of one amino acid mutating into another over a distance of 1 PAM.
#'
#' @param pam The PAM matrix (given in ./data)
#' @param AA1 The origin amino acid
#' @param AA2 The target amino acid
#' @return A probability of mutation.
#' probOfChange(pam, "L", "K")
probOfChange1PAM <- function(pam, AA1, AA2){

  i <- which(rownames(pam) == AA1)
  j <- which(rownames(pam) == AA2)
  return (pam[i, j]/10000)

}

#' Calculate mutation probability
#'
#' Calculate the probability of one amino acid mutating into another over a given branch length.
#'
#' @param pam The PAM matrix given in ./data.
#' @param AA1 The origin amino acid
#' @param AA2 The target amino acid
#' @param d The branch length between the two sequences.
#' @return A scaled probability of mutation.
# @examples
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
  for (aa in AAlist){
    prob = prob + (probOfChange1PAM(pam, AA1, aa) * probOfChange(pam, aa, AA2, d-1))
  }
  return(prob)
}

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
#' @return Scaled likelihood that the amino acids satisfying the conditions of convergent evolution is by chance.
# @examples
#' probOfSiteConfig(tree, primates, "Human", "Chimp", 1)
#' @export
probOfSiteConfig <- function(tree, phydat, spe1, spe2, pos, p=(1/20)){

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

  return (p * prob1 * prob2)
}

#' Calculate likelihood of N chance 'convergent' sites
#'
#' Calculate the probability that n sites between the two species will satisfy the conditions of convergent evolution by chance.
#'
#' @param tree A phylogenetic tree
#' @param spe1 The name of species 1
#' @param spe2 The name of species 2
#' @param m The length of the genes. (Must of equal length)
#' @param n The number of potential convergent sites
#' @param p The fraction of the amino acid at the target position in the evaluated genome. (Default is 1/20)
#' @return Scaled likelihood that the amino acids satisfying the conditions of convergent evolution is by chance.
# @examples
#' probOfNSitesByChance(tree, "Human", "Chimp", 20, 2)
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

#' Find convergent species
#'
#' For a given species, find others in the tree that satisfy the conditions of convergent evolution at a certain position with a less than 0.05 chance of this having occured by chance.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe The species to compare others in the tree to
#' @param pos The position at which to compare the genes
#' @param anc Optional: the value of the ancestral amino acid at position pos (if not supplied, will be predicted using ancestral reconstruction)
#' @return A vector of species names which satisfy the conditions listed above.
# @examples
#' getConvergent(tree, primates, "Human", 1, "L")
#' @export
getConvergent <- function(tree, phydat, spe, pos, anc){
  species <- tree$tip.label
  convSpe <- c(spe)
  for (s in species){
    if (missing(anc)){
      cond = areCondSatisfied(tree, phydat, s, spe, pos)
    } else {
      cond = areCondSatisfied(tree, phydat, s, spe, pos, anc)
    }

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
#' @return The number of potentially convergent sites
#' @examples
#' convSiteData(tree, primates, "Mouse", "Bovine", 20)
#' @export
convSiteData <- function(tree, phydat, spe1, spe2, m=getm(tree, phydat, spe1, spe2)){

  numSites = 0
  for (i in 1:m){
    if (areCondSatisfied(tree, phydat, spe1, spe2, i)) {numSites = numSites + 1}
  }
  print(sprintf("%d potentially convergent sites with a %d probability of this occuring by chance.", numSites, probOfNSitesByChance(tree, spe1, spe2, m, n=numSites)))

  return (numSites)

}


