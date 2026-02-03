# tests/testthat/test-seq_scale.R

test_that("seq_scale rejects invalid inputs", {
  # n_points must be > 1
  expect_error(
    seq_scale(n_points = 1),
    "`n_points` must be a number higher than 1"
  )
  expect_error(seq_scale(n_points = -5))
  expect_error(seq_scale(n_points = "a"))

  # max_indiv > n_points
  expect_error(
    seq_scale(max_indiv = 5, n_points = 10),
    "`max_indiv` must be a number higher than `n_points`"
  )
  expect_error(seq_scale(max_indiv = Inf))

  # min_indiv positive and < max_indiv
  expect_error(seq_scale(min_indiv = 0))
  expect_error(seq_scale(min_indiv = -1))
  expect_error(seq_scale(min_indiv = 100, max_indiv = 50))

  expect_error(seq_scale(step_scale = "banana"))
})

test_that("seq_scale returns increasing unique integers", {
  x <- seq_scale(1, 1000, 10)
  expect_type(x, "integer")
  expect_true(all(diff(x) >= 0)) # monotonic non-decreasing
  expect_equal(length(x), length(unique(x))) # no duplicates
})

test_that("seq_scale includes the expected range endpoints", {
  x <- seq_scale(1, 1000, 10)
  expect_true(min(x) == 1)
  expect_true(max(x) == 1000)
})
