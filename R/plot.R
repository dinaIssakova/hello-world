### Graphical output.
library(ggtree)

rgenesconvergedPlot <- function(tree, phyDat, spe, pos) {

  convNodes <- getConvergent(tree, phyDat, spe, pos)
  convNodes <- c(convNodes, spe)

  groupInfo <- split(tree$tip.label, tree$tip.label%in%convNodes)
  tree <- groupOTU(tree, groupInfo)
  ggtree(tree, aes(color=group)) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

}
