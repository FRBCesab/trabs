test_that('simulate_assemblage input tests works OK', {
  expect_error(simulate_assemblage(S = "abc"))
  expect_error(simulate_assemblage(S = c(1, 35)))
  expect_error(simulate_assemblage(S = NA))
  expect_error(simulate_assemblage(S = -5))
  expect_error(simulate_assemblage(S = 10, N = 5))
  expect_error(simulate_assemblage(trait_mean = NA))
  expect_error(simulate_assemblage(trait_sd = "abc"))
})
