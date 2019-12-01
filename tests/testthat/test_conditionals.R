library(Biostrings)
data(BLOSUM62)
context("getBranchLength")
test_that("Distance is Correct", {
  expect_equal(getBranchLength(tree, "Mouse", 15), 7)
  expect_equal(getBranchLength(tree, "Mouse", 1), 0)
})

context("getMostRecentCommonAncestor")
test_that("Ancestor is Correct", {
  expect_equal(getMostRecentCommonAncestor(tree, "Mouse", "Human"), 15)
  expect_equal(getMostRecentCommonAncestor(tree, "Mouse", "Mouse"), 1)
  expect_equal(getMostRecentCommonAncestor(tree, "BarbMacaq", "Gorilla"), 18)
})

context("mapLetters")
test_that("translates correctly", {
  expect_equal(mapLetters(c(1,2,3,4)), c("A", "C", "G", "T"))
  expect_error(mapLetters(c(1,2,6)))
})

context("areCondSatisfied")
test_that("True is returned correctly", {
  expect_true(areCondSatisfied(tree, primates, "Human", "Chimp", pos=2, type="score", threshold=0, BLOSUM62))
})


context("convertToAA")
test_that("matrices converted correctly", {
  anc.acctran <- ancestral.pars(tree, primates, "ACCTRAN")
  expect_equal(convertToAA(anc.acctran[[1]])[,"T"][1], 1)
  expect_equal(convertToAA(anc.acctran[[1]])[,"T"][2], 0)
})

context("probOfChange1PAM")

test_that("Math works", {
  expect_equal(probOfChange1PAM(pam, "K", "L"), 1e-04)
  expect_equal(probOfChange1PAM(pam, "A", "A"), 0.9867)
})

context("probOfChange")
test_that("Math works", {
  expect_equal(probOfChange(pam, "L", "K", 3), 0.0006039205)
  expect_equal(probOfChange(pam, "K", "K", 3), 0.9780055, tolerance = 0.01)
  expect_equal(probOfChange(pam, "K", "K", 1), 0.9926)
})

context("probOfChange")
test_that("Functions are equivalent", {
  expect_equal(probOfChange(pam, "A", "K", 1), probOfChange1PAM(pam, "A", "K"))
})

context("probOfSiteConfig")
test_that("Function works", {
  expect_equal(probOfSiteConfig(tree, primates, "Human", "Chimp", pos=1), 0.04686169, tolerance = 0.01)
  expect_error(probOfSiteConfig(tree, primates, "Human", "Human", pos=1))
})


context("probOfNSitesByChance")
test_that("Basic case works", {
  expect_equal(probOfNSitesByChance(tree, primates, spe1="Human", spe2="Chimp", m=20, n=2, type="abs"), -5.905249e+36, tolerance = 0.01)
})



context("probOfNSitesByChance")
test_that("zero gives 0", {
  expect_equal(probOfNSitesByChance(tree, primates, spe1="Human", spe2="Chimp", m=20, n=0, type="abs"), 0)
})


context("getConvergent")
test_that("function works on small dataset", {
  expect_equal(getConvergent(smallTree, primates, "Human", 1, type="abs"), c("Human"), ignore.order=TRUE)
})


context("getm")
test_that("function returns correct gene length", {
  expect_equal(getm(tree, primates, "Mouse", "Bovine"), 216)
})

context("convSitedata")
test_that("function correctly returns", {
  expect_equal(convSiteData(smallTree, primates, "Human", "Chimp", t=5, type="abs", m=5), 0)
})



