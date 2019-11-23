#' Get length of branch in tree.
#'
#' Find the (rounded) length of a branch from a given species to a given internal node.
#'
#' @param tree A phylogenetic tree
#' @param spe The name of the species
#' @param nodeNum The number of the internal node desired.
#'
#' @return The length of the branch between nodeNum and spe as given by tree.
#'
#' @examples
#' getBranchLength(tree, "Mouse", 15)
#' @export
#' @importFrom ape is.binary
getBranchLength<- function(tree, spe, nodeNum){
  if (!ape::is.binary(tree)){
    stop("Tree must be binary!")
  }
  if (!is.numeric(nodeNum)){
    stop("nodeNum must be numeric!")
  }
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
#'
#' @return The number of the internal node that is the most recent point of divergence.
#'
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
    return(invisible(NULL))
    stop("Tree parsed incorrectly. Please check that input is a complete binary tree.")
  }
}
