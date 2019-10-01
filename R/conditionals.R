#library(msa)
library(phangorn)

getBranchLength<- function(tree, spe, nodeNum){
  leafNum <- which(tree$tip.label == spe)
  rowNum = which(tree$edge[ ,2] == leafNum)
  ancNum = as.numeric((tree$edge[rowNum, ])[1])

  len = 0
  len = len + tree$edge.length[rowNum]

  while(!is.na(ancNum) && ancNum != nodeNum){
    leafNum = ancNum
    rowNum = which(subtree$edge[ ,2] == leafNum)
    ancNum = as.numeric((tree$edge[rowNum, ])[1])
    len = len + tree$edge.length[rowNum]
  }
  return(len)
}

getMostRecentCommonAncestor<- function(tree, spe1, spe2, areLeaves=TRUE){

  if(areLeaves){
    leafNum1 <- which(tree$tip.label == spe1)
    leafNum2 <- which(tree$tip.label == spe2)
  } else {
    leafNum1 <- as.numeric(spe1)
    leafNum2 <- as.numeric(spe2)
  }

  row1 = tree$edge[which(tree$edge[,2] == leafNum1),]
  row2 = tree$edge[which(tree$edge[,2] == leafNum2),]
  anc1 = as.numeric(row1[1])
  anc2 = as.numeric(row2[1])

  if (is.na(anc1) | is.na(anc2)){
    return (leafNum1)
  }else if (anc1==anc2){
    return(anc1)
  }
  else {
    return (max(getMostRecentCommonAncestor(tree, leafNum1, anc2, areLeaves = FALSE), getMostRecentCommonAncestor(tree, leafNum2, anc1, areLeaves = FALSE)))
  }
}

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

areCondSatisfied <- function(tree, phydat, spe1, spe2, pos){

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- ancestral.pars(tree, phydat, "ACCTRAN")

  geneSeq1 = convertToAA(anc.acctran[[speNum1]])
  geneSeq2 = convertToAA(anc.acctran[[speNum2]])
  ancSeq = convertToAA(anc.acctran[[getMostRecentCommonAncestor(tree, spe1, spe2)]])

  # For now : without MSA.

  gene1Val <- which(geneSeq1[pos,] == 1)
  gene2Val <- which(geneSeq2[pos,] == 1)

  # Position is convergently evolved if:
  #the amino acids at the descendant nodes are identical with each other (x1 = x3)
  cond1 = (gene1Val==gene2Val)

  # And different from their respective ancestral amino acids
  cond2 = (ancSeq[pos,][gene1Val] == 0)
  cond3 = (ancSeq[pos,][gene2Val] == 0)

  return (cond1 && cond2 && cond3)
}

convertToAA <- function(acctranData){

  numrow = nrow(acctranData)

  while (numrow %% 3 != 0){
    numrow = numrow - 1
  }

  AAmatrix <- matrix(data=0, numrow, ncol=30)
  colnames(AAmatrix) = Biostrings::AA_ALPHABET

  DNArow = 0
  AArow = 1
  while (AArow < numrow && DNArow+3 < nrow(acctranData)){
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

probOfChange <- function(pam, AA1, AA2, b){

  i <- which(rownames(pam) == AA1)
  j <- which(rownames(pam) == AA2)
  return (pam[i+1, j+1] * b)

}

# Future - refine these numbers by considering 'inside' nodes.
probOfSiteConfig <- function(p=(1/20), tree, spe1, spe2, pos){

  ## Probability that both just mutated to get there.
  ancNode = getMostRecentCommonAncestor(tree, spe1, spe2)

  b1 = getBranchLength(tree, spe1, ancNode)
  b2 = getBranchLength(tree, spe2, ancNode)

  # TODO: Probably remove out with helper.

  species <- tree$tip.label

  speNum1 = which(species == spe1)
  speNum2 = which(species == spe2)

  anc.acctran <- ancestral.pars(tree, phydat, "ACCTRAN")

  geneSeq1 = convertToAA(anc.acctran[[speNum1]])
  geneSeq2 = convertToAA(anc.acctran[[speNum2]])
  ancSeq = convertToAA(anc.acctran[[getMostRecentCommonAncestor(tree, spe1, spe2)]])

  gene1Val <- which(geneSeq1[pos,] == 1)
  gene2Val <- which(geneSeq2[pos,] == 1)
  ancVal <- which(ancSeq[pos,]==1)

  prob1 <- probOfChange(pam, gene1Val, ancVal, b1)
  prob2 <- probOfChange(pam, gene2Val, ancVal, b2)

  return (p * prob1 * prob2)
}

probOfNSitesByChance <- function(p=(1/20), tree, spe1, spe2, m, n){
  # probability of observing n sites or more convergent-change sites by chance
  # Future - consider them independently. Currently - assume uniform model (on average)
  avg=0
  for (j in 1:m){
    avg = avg + probOfSiteConfig(p, tree, spe1, spe2, j)
  }

  f = avg/m

  prob = 0
  for (i in 0:n-1){

    prob = prob + ( factorial(m) / factorial(i)*factorial(m-i) ) * f**i * (1-f)**(m-i)
  }

  return (1 - prob)
}

getConvergent <- function(tree, phyDat, spe, pos){
  species <- tree$tip.label
  convSpe <- c()
  for (s in species){
    if (s != spe && areCondSatisfied(tree, phyDat, s, spe, pos)){

      p = probOfSiteConfig(tree, s, spe, pos)
      print(sprintf("Species %s is potentially convergent with p %d.", spe, p))
    }

    if (p < 0.05){
      convSpe = c(convSpe, s)
    }

  }
  return (convSpe)
}

convSiteData <- function(tree, phydat, spe1, spe2, m){

  numSites = 0
  for (i in 1:m){
    if (areCondSatisfied(tree, phydat, spe1, spe2, i)) {numSites = numSites + 1}
  }
  print(sprintf("%d potentially convergent sites with a %d probability of this occuring by chance.", numSites, probOfNSitesByChance(tree, spe1, spe2, m, numSites)))

  return (numSites)

}

