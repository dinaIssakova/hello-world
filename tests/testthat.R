library(testthat)
library(rgenesconverged)

context("getBranchLength")
test_that("Distance is Correct", {
  expect_equal(getBranchLength(tree, "Mouse", 15), 67.7, tolerance = 0.1)
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
test_that("False is parsed correctly", {
  expect_false(areCondSatisfied(tree, primates, "Mouse", "Bovine", 17))
})



















#test_check("rgenesconverged")

