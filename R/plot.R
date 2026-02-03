#' Plot the rarefaction curves
#'
#' @description
#' plot the metrics and the different
#' calculate the rarefaction curves based on different null models
#' ...
#' to be completed
#'
#' @param curves a `data.frame` object output of rarefy_curves
#' @param indices number of permutations per null model
#' @param null_models number of estimates per curve (needed?)
#' @param title maximum individuals sampled for the rarefaction curves
#' @param show_ci logical argument to show confidence interval
#'
#' @returns A `ggplot2` object
#'
#' @export
#'
#' @examples
#' asm <- simulate_assemblage(S = 20, N = 1000)
#' cur <- get_rarefy_curves(asm, n_perm = 100, n_points = 100)
#' gplot_curves(cur)
#'
gplot_curves <- function(
  curves,
  indices = c("CWM", "sumT", "rangeT", "S"),
  null_models = c("emp", "null_Tr", "null_SAD", "null_all", "flat"),
  title = "Accumulation curves",
  show_ci = TRUE
) {
  indices <- match.arg(indices, several.ok = TRUE)
  null_models <- match.arg(null_models, several.ok = TRUE)

  curves <- curves[curves$ind %in% indices, ]
  if (length(indices) > 1) {
    curves$ind <- factor(curves$ind, levels = indices, ordered = TRUE)
  }

  curves <- curves[curves$scn %in% null_models, ]
  if (length(null_models) > 1) {
    curves$scn <- factor(curves$scn, levels = null_models, ordered = TRUE)
  }

  ind_p <- grep("^p[1-9]", names(curves))
  p_order <- gsub("p", "", names(curves)[ind_p]) |> as.numeric() |> order()
  names(curves)[ind_p[p_order]] <- c("lo", "med", "hi")

  p <- ggplot2::ggplot(curves, ggplot2::aes(x = N, group = scn))

  if (isTRUE(show_ci)) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lo, ymax = hi, fill = scn),
        alpha = 0.10,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(y = med, colour = scn),
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~ind, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Individuals sampled (m)",
      y = NULL,
      colour = NULL,
      fill = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12)
  return(p)
}

#' Plot the effect curves
#'
#' @description
#' plot the effect curves
#' ...
#' to be completed
#'
#' @param curves a `data.frame` object output of rarefy_curves
#' @param indices number of permutations per null model
#' @param effects number of estimates per curve (needed?)
#' @param show_ci logical argument to show confidence interval
#' @param palette color palette
#'
#' @returns A `ggplot2` object
#'
#' @export
#'
#' @examples
#' asm <- simulate_assemblage(S = 20, N = 1000)
#' cur <- get_rarefy_curves(asm, n_perm = 100, n_points = 100)
#' gplot_effects(cur)
#'
gplot_effects <- function(
  curves,
  indices = c("CWM", "sumT", "rangeT", "S"),
  effects = c("covTA", "SAD", "traits"),
  show_ci = TRUE,
  palette = "Set1"
) {
  indices <- match.arg(indices, several.ok = TRUE)
  effects <- match.arg(effects, several.ok = TRUE)

  # select indices
  curves <- curves[curves$ind %in% indices, ]
  if (length(indices) > 1) {
    curves$ind <- factor(curves$ind, levels = indices, ordered = TRUE)
  }
  # partition effects
  curves_part <- partition(
    curves,
    conservative_ci = show_ci,
    long_format = TRUE
  )
  # select only the effect
  effect <- curves_part[grep("^effect_", curves_part$partition), ]
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

  w_effect <- w_effect[w_effect$part %in% effects, ]
  if (length(effects) > 1) {
    w_effect$part <- factor(w_effect$part, levels = effects, ordered = TRUE)
  }

  # Build plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60")

  if (show_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        data = w_effect,
        ggplot2::aes(x = N, ymin = value.lo, ymax = value.hi, fill = part),
        alpha = 0.1
      )
  }

  p <- p +
    ggplot2::geom_line(
      data = w_effect,
      ggplot2::aes(x = N, y = value.med, colour = part),
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_brewer(palette = palette, drop = FALSE) +
    ggplot2::scale_fill_brewer(palette = palette, drop = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Individuals sampled (m)",
      y = "Effect (EMP - null)",
      colour = NULL,
      fill = NULL,
      title = glue::glue(
        "Effects vs m - {paste(indices, collapse = ', ')}"
      ),
      subtitle = if (show_ci) {
        "Ribbon = 10-90% quantiles (conservative bounds)"
      } else {
        NULL
      }
    ) +
    ggplot2::guides(fill = "none") + # hide duplicate legend
    ggplot2::theme_minimal(base_size = 12)

  return(p)
}
