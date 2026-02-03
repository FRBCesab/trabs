# Can be improved with:
# - different null models
# - different distribution

## Simple simulator to validate pipeline against known cov(T,A)
#' Clean and homogenized city names
#'
#' @description
#' Create random communities with:
#' - normally distributed trait values
#' - negatively distributed individuals among species
#' - trait abundance coupling
#'
#' @param S number of species
#' @param N number of individuals
#' @param trait_mean mean value of the trait
#' @param trait_sd standard deviation of the trait
#' @param sad_k dispersion of the neagitve binomial, closer to 0 means more unequal SAD
#' @param beta_cov standard deviation of the trait
#'
#' @returns A `data.frame` with community attributes
#'
#' @export
#'
#' @examples
#' simulate_assemblage(S=20, N=1000)
#'
simulate_assemblage <- function(
  S = 20,
  N = 1000,
  trait_mean = 1,
  trait_sd = 1,
  sad_k = 1, # NB shape, smaller -> more uneven SAD
  beta_cov = 0.8 # strength of trait–abundance covariance
) {
  # Checking the inputs ------------
  stopifnot(
    "`S` must be a single positive number." = {
      is.numeric(S) & length(S) == 1 & S > 0
    }
  )
  stopifnot(
    "`N` must be a number strictly higher than `S`." = {
      is.numeric(N) & length(N) == 1 & N > S
    }
  )
  stopifnot(
    "`trait_mean` must be a single number." = {
      is.numeric(trait_mean) & length(trait_mean) == 1
    }
  )
  stopifnot(
    "`trait_sd` must be a single positive number." = {
      is.numeric(trait_sd) & length(trait_sd) == 1 & trait_sd >= 0
    }
  )
  stopifnot(
    "`sad_k` must be a single positive number." = {
      is.numeric(sad_k) & length(sad_k) == 1 & sad_k > 0
    }
  )
  stopifnot(
    "`beta_cov` must be a single number." = {
      is.numeric(beta_cov) & length(beta_cov) == 1
    }
  )

  # Simulation -----------
  # Species traits
  Tsp <- stats::rnorm(S, trait_mean, trait_sd)

  # Base SAD via NegBin (overdispersed, nonneg)
  base_mu <- N / S
  A0 <- stats::rnbinom(S, mu = base_mu, size = sad_k)

  # Trait–abundance coupling as weights (positive)
  w <- (A0 + 1e-6) * exp(beta_cov * scale(Tsp))
  w <- pmax(w, 0) # just in case sad_k is odd or scale(Tsp) is extreme
  if (!is.finite(sum(w)) || sum(w) <= 0) {
    w <- rep(1, S)
  }

  # Final counts: integer, ≥0, sum exactly to N
  A <- as.vector(stats::rmultinom(1, size = N, prob = w / sum(w)))

  out <- data.frame(Sp = paste0("sp", seq_len(S)), Ab = A, Tr = Tsp)
  return(new_trasm(out))
}
