#' Standard PAM matrix.
#'
#'A matrix for computing probabilities of common AA mutations.
#'
#' @format A matrix with 20 rows and 20 columns.
#'
#' @source Dayhoff, M.O. et al. (1978) A model of evolutionary change in proteins. Atlast of Protein Sequence and Structure. p. 345-352
"pam"

#' Phylogenetic tree with 14 tips and 12 internal nodes.
#' Tip labels:
#'  Mouse, Bovine, Lemur, Tarsier, Squir Monk, Jpn Macaq, ...
#' Unrooted; includes branch lengths.
#'
#' @format Contains leaf names, edge matrix, edge lengths, and number of internal nodes.
#'
#' @source phangorn
"tree"

#' Phylogenetic tree with 2 tips and 1 internal nodes. Subtree to tree.
#' Tip labels:
#'  Human, Chimp
#' Rooted; includes branch lengths.
#'
#' @format Contains leaf names, edge matrix, edge lengths, and number of internal nodes.
#'
"smallTree"

#' Object of class phyDat containing sequences for species corresponding to tree.
#'
#' 14 sequences with 232 character and 217 different site patterns.
#' The states are a c g t
#'
#' @format phyDat
#'
"primates"
