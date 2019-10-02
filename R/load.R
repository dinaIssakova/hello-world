library(devtools)

pam <- read.table("PAM.txt")


if(!require("ape")){
  install.packages("ape")
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (! require("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
if (! require("phangorn")) {
  BiocManager::install("phangorn")
}
if(!require("seqLogo")){
  BiocManager::install("seqLogo")
}
if(!require("ggtree")){
  BiocManager::install("ggtree")
}



