#' rgenesconvergedPlot
#'
#' Plot tree associated with the given gene, accentuating convergent lineages.
#'
#' @param tree A phylogenetic tree
#' @param phydat An object of class phydat
#' @param spe The name of the species to compare to
#' @param pos The position of interest
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' rgenesconvergedPlot(tree, primates, "Mouse", 17)
#' }
#'
#' @importFrom ggtree ggtree
#' @importFrom ggtree groupOTU
#' @importFrom ggtree geom_tiplab
#' @importFrom ggtree geom_text2
#' @import ggplot2
#' @export
rgenesconvergedPlot <- function(tree, phydat, spe, pos) {

  convNodes <- getConvergent(tree, phydat, spe, pos)
  convNodes <- c(convNodes, spe)

  groupInfo <- split(tree$tip.label, tree$tip.label%in%convNodes)
  tree <- ggtree::groupOTU(tree, groupInfo)

  # Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
  # ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates
  # and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628
  ggtree::ggtree(tree, aes(color=group)) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

  return(invisible(NULL))
}
