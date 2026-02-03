## Make steps
#' Returns a uniformly distributed steps along scale
#'
#' @description
#' Make a monotonically increasing integer vector of steps in \[min_indiv, max_indiv\],
#' with denser spacing near small m when mode = "log" or "sqrt".
#'
#' @param min_indiv minimum number of
#' @param max_indiv maximum number of
#' @param n_points number of points estimated in the curves
#' @param step_scale scale of the steps to calculate the curves (not clear)
#'
#' @returns A `data.frame` with estimated trait-abundance relationship
#'
#' @export
#'
#' @examples
#' plot(seq_scale(1, 1000, 10, "log"), type="l")
#' lines(seq_scale(1, 1000, 10, "sqrt"), col="red")
#' lines(seq_scale(1, 1000, 10, "linear"), col="blue")
#'
seq_scale <- function(
  min_indiv = 1,
  max_indiv = 1000,
  n_points = 10,
  step_scale = c("log", "sqrt", "linear")
) {
  # Checking the inputs ---------------------------------------------
  stopifnot(
    "`n_points` must be a number higher than 1." = {
      is.numeric(n_points) & length(n_points) == 1 & n_points > 1
    }
  )
  stopifnot(
    "`max_indiv` must be a number higher than `n_points`." = {
      is.numeric(max_indiv) &
        length(max_indiv) == 1 &
        max_indiv > n_points &
        is.finite(max_indiv)
    }
  )
  stopifnot(
    "`min_indiv` must be a positive number lower than `max_indiv`." = {
      is.numeric(min_indiv) &
        length(min_indiv) == 1 &
        min_indiv > 0 &
        min_indiv < max_indiv
    }
  )
  step_scale <- match.arg(
    step_scale,
    c("log", "sqrt", "linear"),
    several.ok = FALSE
  )

  # Creating the steps ---------------------------------------------
  vals <- switch(
    step_scale,
    "log" = {
      # Robust near 1 using log1p/expm1; maps 0..log(m_max-1+1) to 0..(m_max-1)
      expm1(seq(log1p(min_indiv), log1p(max_indiv), length.out = n_points))
    },
    "sqrt" = {
      # Heavier near low m via concave power transform
      seq(sqrt(min_indiv), sqrt(max_indiv), length.out = n_points)**2
    },
    "linear" = {
      # default: linear
      seq(min_indiv, max_indiv, length.out = n_points)
    }
  )

  # transform to integer
  steps <- as.integer(round(vals))
  # and remove the duplicates
  steps <- sort(unique(steps))

  # not needed anymore
  # steps <- steps[steps >= 1 & steps <= max_indiv]
  # ensure endpoints present
  # if (1L %notin% steps) {
  #   steps <- c(1L, steps)
  # }
  # if (ensure_last && max_indiv %notin% steps) {
  #   steps <- c(steps, max_indiv)
  # }

  return(steps)
}

#' Get the first element of a string
#'
#' @param x a vector of character
#' @param split separator. By default, '_'
#'
#' @returns the first element separated by '_'
#'
#' @keywords internal
#' @noRd
first_str <- function(x, split = "_") {
  strsplit(x, split = "_") |>
    sapply(function(y) y[[1]])
}

#' Get the second element of a string
#'
#' @param x a vector of character
#' @param split separator. By default, '_'
#'
#' @returns the second element separated by '_'
#'
#' @keywords internal
#' @noRd
second_str <- function(x, split = "_") {
  strsplit(x, split = "_") |>
    sapply(function(y) ifelse(length(y) > 1, y[[2]], NA))
}

#' Get the last element of a string
#'
#' @param x a vector of character
#' @param split separator. By default, '_'
#'
#' @returns the last element separated by '_'
#'
#' @keywords internal
#' @noRd
last_str <- function(x, split = "_") {
  strsplit(x, split = "_") |>
    sapply(function(y) y[[length(y)]])
}


#' Get number of unique elements
#'
#' @param x string vector
#'
#' @returns a number of unique elements
#'
#' @keywords internal
#' @noRd
nodup <- function(x) {
  return(sum(!duplicated(x)))
}


# Scale traits
#' Scale variables with option to keep mean and set objective sd
#'
#' @description
#' scaling function
#'
#' @param min_indiv minimum number of
#' @param max_indiv maximum number of
#' @param n_points number of points estimated in the curves
#' @param step_scale scale of the steps to calculate the curves (not clear)
#' @param ensure_last scale of the steps to calculate the curves (not clear)
#'
#' @returns A `data.frame` with estimated trait-abundance relationship
#'
#' @keywords internal
#' @noRd
scale_traits <- function(Tr, sd_ref = 1, keep_mean = TRUE) {
  mu <- mean(Tr)
  s <- stats::sd(Tr)
  if (!is.finite(s) || s <= 0) {
    return(Tr)
  }
  Tz <- (Tr - mu) / s
  if (keep_mean) mu + sd_ref * Tz else sd_ref * Tz
}
