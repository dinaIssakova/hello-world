#' Plot the convergent network associated with the given gene at the given position on top of the phylogenetic tree.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe The name of the species to compare to
#' @param pos The position of interest
#' @return None
# @examples
#' rgenesconvergedPlot(tree, primates, "Mouse", 17)
#' @export
rgenesconvergedPlot <- function(tree, phydat, spe, pos) {

  convNodes <- getConvergent(tree, phydat, spe, pos)
  convNodes <- c(convNodes, spe)

  groupInfo <- split(tree$tip.label, tree$tip.label%in%convNodes)
  tree <- groupOTU(tree, groupInfo)
  ggtree(tree, aes(color=group)) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

}
