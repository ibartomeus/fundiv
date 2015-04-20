#very silly to check if I get the concept

context("tests on inputs")

test_that("tests for length_to_mass variable",{
  x <- c(0,1,"cosa")
  taxa <- "Bees"
  
  expect_that(length_to_mass(x, taxa),condition=throws_error())
  })

context("test on outputs")

test_that("test values are numeric and non zero",{
  x <- c(2,1,NA)
  taxa <- "Bees"
  
  expect_that(length_to_mass(x, taxa),is_a("numeric"))
  expect_that(all(length_to_mass(x, taxa) > 0),is_true())
})

#print("All tests OK")