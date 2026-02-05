#' Partition the effects from the rarefaction curves
#'
#' @description
#' partition the effect between covTA, SAD, traits, SAD_flat
#' ...
#' to be completed
#'
#' @param curves a `trasmout` data object
#' @param conservative_ci calculate the confidence intervals
#' @param add_percent calculate percentages
#' @param eps threshold value for calculation percentages
#' @param output format of the output, either "wide", "long", or "effect"
#'
#' @returns A `data.frame` with partition effects in long or wide format
#'
#' @export
#' @examples
#' asm <- simulate_assemblage(S = 20, N = 1000)
#' cur <- get_rarefy_curves(asm, n_perm = 100, n_points = 100)
#' partition(cur)
#'
partition <- function(
  curves,
  conservative_ci = TRUE,
  add_percent = TRUE,
  eps = 1e-8,
  output = c("wide", "long", "effect")
) {
  output <- match.arg(output, several.ok = FALSE)

  # transformation to be moved in rarefy()
  ind_p <- grep("^p[1-9]", names(curves))
  p_order <- gsub("p", "", names(curves)[ind_p]) |> as.numeric() |> order()
  names(curves)[ind_p[p_order]] <- c("lo", "med", "hi")

  stopifnot(
    "`curves` must contain columns `N`, `ind`, `scn`, `med`, `lo`, `hi`" = {
      all(c("N", "ind", "scn", "med", "lo", "hi") %in% names(curves))
    }
  )

  keepC <- names(curves)[!names(curves) %in% c("scn", "lo", "med", "hi")]
  # to be checked : unique(curves$scn)

  # check if duplicated rows before reshape
  wcur <- stats::reshape(
    curves,
    idvar = c("N", "ind"),
    timevar = "scn",
    direction = "wide",
    v.names = c("med", "lo", "hi")
  )
  # equivalent to tidyr pivot_wider()
  # wcur <- tidyr::pivot_wider(
  #   curves,
  #   names_from = c(scn),
  #   values_from = c(med, lo, hi)
  # )

  out <- data.frame(
    wcur[, c(keepC, grep("med", names(wcur), value = TRUE))],
    # simple effects
    "effect_covTA_med" = wcur$med.emp - wcur$med.null_Tr, # isolates T–A coupling
    "effect_SAD_med" = wcur$med.null_Tr - wcur$med.null_all, # restores SAD only (T random)
    "effect_traits_med" = wcur$med.null_SAD - wcur$med.null_all, # restores trait marginals only (A random)
    "effect_SADflat_med" = wcur$med.null_Tr - wcur$med.flat
  )

  if (conservative_ci) {
    # confidence interval
    ci <- data.frame(
      "effect_covTA_lo" = wcur$lo.emp - wcur$hi.null_Tr,
      "effect_covTA_hi" = wcur$hi.emp - wcur$lo.null_Tr,
      "effect_SAD_lo" = wcur$lo.null_Tr - wcur$hi.null_all,
      "effect_SAD_hi" = wcur$hi.null_Tr - wcur$lo.null_all,
      "effect_traits_lo" = wcur$lo.null_SAD - wcur$hi.null_all,
      "effect_traits_hi" = wcur$hi.null_SAD - wcur$lo.null_all,
      "effect_SADflat_lo" = wcur$lo.null_Tr - wcur$hi.flat,
      "effect_SADflat_hi" = wcur$hi.null_Tr - wcur$lo.flat
    )
    out <- cbind(out, ci)
  }

  if (add_percent & output != "effect") {
    is_emp_positive <- is.finite(wcur$med.emp) & abs(wcur$med.emp) > eps
    wcur$denom_emp = ifelse(is_emp_positive, wcur$med.emp, NA)

    # Additive path (EMP - NULL_ALL)
    delta_all = wcur$med.emp - wcur$med.null_all
    is_delta_positive <- is.finite(delta_all) & abs(delta_all) > eps
    wcur$denom_delta = ifelse(is_delta_positive, delta_all, NA)

    pct <- data.frame(
      "pct_covTA_of_EMP" = 100 * out$effect_covTA_med / wcur$denom_emp,
      "pct_SAD_of_EMP" = 100 * out$effect_SAD_med / wcur$denom_emp,
      "pct_traits_of_EMP" = 100 * out$effect_traits_med / wcur$denom_emp,
      "pct_covTA_of_EMP_abs" = 100 * abs(out$effect_covTA_med / wcur$denom_emp),
      "pct_SAD_of_EMP_abs" = 100 * abs(out$effect_SAD_med / wcur$denom_emp),
      "pct_traits_of_EMP_abs" = 100 *
        abs(out$effect_traits_med / wcur$denom_emp),
      "pct_covTA_viaSAD" = 100 * out$effect_covTA_med / wcur$denom_delta,
      "pct_SAD_only" = 100 * out$effect_SAD_med / wcur$denom_delta,
      "pct_TRAITS_only" = 100 * out$effect_traits_med / wcur$denom_delta,
      "pct_covTA_viaSAD_abs" = 100 *
        abs(out$effect_covTA_med / wcur$denom_delta),
      "pct_SAD_only_abs" = 100 * abs(out$effect_SAD_med / wcur$denom_delta),
      "pct_TRAITS_only_abs" = 100 *
        abs(out$effect_traits_med / wcur$denom_delta)
    )
    out <- cbind(out, pct)
  }
  if (output == "wide") {
    return(out)
  } else {
    var <- names(out)[!names(out) %in% keepC]
    long <- stats::reshape(
      out,
      varying = var,
      v.names = "value",
      timevar = "partition",
      times = var,
      direction = "long",
      idvar = c("N", "ind")
    )
    if (output == "long") {
      return(long)
    } else {
      # select only the effect
      effect <- long[grep("^effect_", long$partition), ]
      # build new variables
      effect$part <- second_str(effect$partition)
      effect$level <- last_str(effect$partition)
      # reshape to wide
      w_effect <- stats::reshape(
        effect[names(effect) != "partition"],
        idvar = c("N", "ind", "part"),
        timevar = "level",
        direction = "wide",
        v.names = "value"
      )
      return(w_effect)
    }
  }
}
