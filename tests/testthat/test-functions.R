library(testthat)
library(mrhevo)

test_that("pnorm.extreme returns string in scientific notation", {
  result <- pnorm.extreme(10)
  expect_type(result, "character")
  expect_match(result, "E")
})

test_that("format.scinot.pvalue formats correctly", {
  result <- format.scinot.pvalue("5E-10")
  expect_type(result, "character")
  expect_match(result, "times")
})

test_that("format.z.aspvalue returns formatted p-values", {
  result <- format.z.aspvalue(2)
  expect_type(result, "character")
})

test_that("get_coeffratios calculates ratios correctly", {
  coeffs.dt <- data.table::data.table(
    alpha_hat = c(1, 2, 3),
    gamma_hat = c(0.5, 1, 1.5),
    se.alpha_hat = c(0.1, 0.1, 0.1),
    se.gamma_hat = c(0.1, 0.1, 0.1)
  )
  result <- get_coeffratios(coeffs.dt)
  expect_true("theta_IV" %in% colnames(result))
  expect_true("se.theta_IV" %in% colnames(result))
  expect_true("inv.var" %in% colnames(result))
})

test_that("get_estimatorsMR returns valid estimates", {
  coeffs.dt <- data.table::data.table(
    alpha_hat = c(1, 2, 3),
    gamma_hat = c(0.5, 1, 1.5),
    se.alpha_hat = c(0.1, 0.1, 0.1),
    se.gamma_hat = c(0.1, 0.1, 0.1),
    theta_IV = c(0.5, 0.5, 0.5),
    se.theta_IV = c(0.1, 0.05, 0.03),
    inv.var = c(100, 400, 1111)
  )
  result <- get_estimatorsMR(coeffs.dt)
  expect_true("Estimate" %in% colnames(result))
  expect_true("SE" %in% colnames(result))
  expect_true("z" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$Estimator, "Inverse variance weighted")
})

test_that("set.tau0 calculates correctly", {
  result <- set.tau0(fraction_pleio = 0.5, nu_global = 1, J = 10, info = 1)
  expect_type(result, "double")
  expect_gt(result, 0)
})
