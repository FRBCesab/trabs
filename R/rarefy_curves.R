#' Calculate the rarefaction curves
#'
#' @description
#' from assemblages in long format
#' calculate the rarefaction curves based on different null models
#' ...
#' to be completed
#'
#' @param asm a `trasm` data object
#' @param n_perm number of permutations per null model
#' @param n_points number of estimates per curve (needed?)
#' @param max_indiv maximum individuals sampled for the rarefaction curves
#' @param step_scale scale of the steps to build the rarefaction curves
#' @param probs quantile of the permutation probabilities
#' @param indices trait indices to be calculated per null model
#' @param null_models null models to be calculated
#'
#' @returns A `data.frame` with community attributes
#'
#' @export
#'
#' @examples
#' asm <- simulate_assemblage(S = 20, N = 1000)
#' cur <- get_rarefy_curves(asm, n_perm = 100, n_points = 100)
#' gplot_curves(cur, indices = c("S", "CWM"), null_models = c("emp", "flat"))
#'
get_rarefy_curves <- function(
  asm,
  n_perm = 999,
  n_points = 100,
  max_indiv = 10000,
  step_scale = c("log", "sqrt", "linear"),
  probs = c(0.1, 0.5, 0.9),
  indices = c("CWM", "sumT", "rangeT", "S"),
  null_models = c("emp", "null_Tr", "null_SAD", "null_all", "flat")
) {
  # Checking the inputs ------------
  stopifnot("`asm` must be a `trasm` object." = {
    "trasm" %in% class(asm)
  })
  stopifnot(
    "`n_perm` must be a single positive number." = {
      is.numeric(n_perm) & length(n_perm) == 1 & n_perm > 0
    }
  )
  stopifnot(
    "`n_points` must be a single number higher than 1." = {
      is.numeric(n_points) & length(n_points) == 1 & n_points > 1
    }
  )
  stopifnot(
    "`max_indiv` must be a single positive number." = {
      is.numeric(max_indiv) & length(max_indiv) == 1 & max_indiv > 0
    }
  )

  # create steps
  step_scale <- match.arg(
    step_scale,
    c("log", "sqrt", "linear"),
    several.ok = FALSE
  )

  # Create steps
  # Number of individuals to use
  m_use <- min(sum(asm$Ab), max_indiv)
  steps <- seq_scale(2L, m_use, n_points, step_scale = step_scale)

  # characteristics of the dataset
  tot_Ab <- sum(asm$Ab)
  N <- nrow(asm)
  # random draw of individuals based on abundance
  if (any(c("emp", "null_Tr") %in% null_models)) {
    p_emp <- asm$Ab / tot_Ab
    i_emp <- sapply(1:n_perm, function(x) {
      sample(1:N, m_use, replace = TRUE, prob = p_emp)
    })
  }
  # random traits
  if (any(c("null_Tr", "null_all", "flat") %in% null_models)) {
    perm_tr <- sapply(1:n_perm, function(x) {
      sample(asm$Tr, N, replace = FALSE)
    })
  }
  # random abundance
  if (any(c("null_SAD", "null_all") %in% null_models)) {
    i_nullsad <- sapply(1:n_perm, function(x) {
      perm_ab <- sample(asm$Ab, N, replace = FALSE)
      p_nullsad <- perm_ab / tot_Ab
      sample(1:N, m_use, replace = TRUE, prob = p_nullsad)
    })
  }
  # flat abundance
  if ("flat" %in% null_models) {
    i_flat <- sapply(1:n_perm, function(x) {
      sample(1:N, m_use, replace = TRUE)
    })
  }

  out <- list()
  if ("emp" %in% null_models) {
    out[["emp"]] <- rarefy_curves(
      asm$Sp,
      asm$Tr,
      i_emp,
      steps = steps,
      scn = "emp",
      probs = probs
    )
  }
  if ("null_Tr" %in% null_models) {
    out[["null_Tr"]] <- rarefy_curves(
      asm$Sp,
      perm_tr,
      i_emp,
      steps = steps,
      scn = "null_Tr",
      probs = probs
    )
  }
  if ("null_SAD" %in% null_models) {
    out[["null_SAD"]] <- rarefy_curves(
      asm$Sp,
      asm$Tr,
      i_nullsad,
      steps = steps,
      scn = "null_SAD",
      probs = probs
    )
  }
  if ("null_all" %in% null_models) {
    out[["null_all"]] <- rarefy_curves(
      asm$Sp,
      perm_tr,
      i_nullsad,
      steps = steps,
      scn = "null_all",
      probs = probs
    )
  }
  if ("flat" %in% null_models) {
    out[["flat"]] <- rarefy_curves(
      asm$Sp,
      perm_tr,
      i_flat,
      steps = steps,
      scn = "flat",
      probs = probs
    )
  }
  return(do.call(rbind, out))
}

#' Calculate the rarefaction curves
#'
#' @description
#' from assemblages permutation matrices of species or traits
#'
#' @param species a vector of species id
#' @param trait a matrix with species in row and permutation in column
#' @param samples a premutation matrix with row id and permutation in column
#' @param steps steps to be kept for the rarefaction curves
#' @param probs quantile of the permutation probabilities
#' @param indices trait indices to be calculated per null model
#' @param scn name of the null model
#'
#' @returns A `data.frame` with rarefaction curves
#'
#' @keywords internal
#' @noRd
rarefy_curves <- function(
  species,
  trait,
  samples,
  steps,
  probs = c(0.1, 0.5, 0.9),
  indices = c("CWM", "sumT", "rangeT", "S"),
  scn = "TEST"
) {
  trait <- as.matrix(trait)
  samples <- as.matrix(samples)
  # calculate indices across permutations
  perm_ind <- list()
  # transform with lapply to make it faster?
  for (i in 1:ncol(samples)) {
    sample_i <- samples[, i]
    j <- ifelse(ncol(trait) == 1, 1, i)
    perm_ind[[i]] <- get_single_curve(
      species[sample_i],
      trait[sample_i, j],
      steps = steps,
      indices = indices
    )
  }
  # transform to array
  perm_ar <- abind::abind(perm_ind, along = 3)

  # compute stats
  perm_quant <- apply(perm_ar[, -1, ], c(1, 2), stats::quantile, probs = probs)

  dimnames(perm_quant)[[1]] <- paste0("p", round(probs * 100))

  # format as original script
  # perm_quant <- aperm(perm_quant, c(2, 3, 1))
  # as.matrix(stats::ftable(perm_quant, row.vars = 1))

  # format for easy plotting
  out <- as.matrix(stats::ftable(perm_quant, row.vars = c(3, 2)))
  # transform to matrix
  out <- data.frame(
    N = steps,
    scn = scn,
    ind = first_str(row.names(out)),
    out,
    row.names = NULL
  )
  return(out)
}

#' Calculate the rarefaction curves
#'
#' @description
#' from a single permutation
#'
#' @param sample_Sp a vector of species id
#' @param sample_Tr a matrix with species in row and permutation in column
#' @param steps steps to be kept for the rarefaction curves
#' @param indices trait indices to be calculated per null model
#'
#' @returns A `data.frame` with rarefaction curves
#'
#' @keywords internal
#' @noRd
get_single_curve <- function(
  sample_Sp,
  sample_Tr,
  steps = 1,
  indices = c("CWM", "sumT", "rangeT", "S")
) {
  n <- seq_along(sample_Tr)
  cumT <- cumsum(sample_Tr)
  run_min <- cummin(sample_Tr)
  run_max <- cummax(sample_Tr)
  cumS <- cumsum(!duplicated(sample_Sp))
  out <- data.frame(
    "n" = n,
    "CWM" = cumT / n,
    "sumT" = cumT,
    "rangeT" = run_max - run_min,
    "S" = cumS
  )
  return(out[steps, c("n", indices)])
}
