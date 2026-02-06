# The objectif of this script is to test whether the duplicates change the trait-abundance relationship
# Conclusion:
#   issue with summary_trasm() meanT, sdT and beta_cov should be weighted version
#   issue with null_SAD and S accumulation curves that change with duplicates
devtools::load_all()
library(ggplot2)

# Normal run with no duplicates ----------------
asm1 <- simulate_assemblage(
  S = 200,
  N = 10000,
)

system.time({
  cur1 <- get_rarefy_curves(
    asm1,
    n_perm = 999,
    n_points = 1000
  )
}) # 9 sec

gplot_curves(cur1)
gplot_effects(cur1)

# With duplicates ------------------------------
asm2 <- asm1[rep(1:nrow(asm1), 2), ]
asm2$Ab <- asm2$Ab / 2
summary_trasm(asm2) == summary_trasm(asm1)
nrow(asm2) == nrow(asm1) * 2
# not sdT et beta_hat !! should be weighted indivators?

system.time({
  cur2 <- get_rarefy_curves(
    asm2,
    n_perm = 999,
    n_points = 1000
  )
}) # 9 sec
# same time to compute the curves
gplot_curves(cur2)
# in S, null all different from empirical !
gplot_effects(cur2)

# With quintuplates ------------------------------
asm5 <- asm1[rep(1:nrow(asm1), 5), ]
asm5$Ab <- asm5$Ab / 5
summary_trasm(asm5) == summary_trasm(asm1)
# not sdT et beta_hat !
nrow(asm5) == nrow(asm1) * 5


system.time({
  cur5 <- get_rarefy_curves(
    asm5,
    n_perm = 999,
    n_points = 1000
  )
}) # 9 sec
# same time to compute the curves
gplot_curves(cur5)
# in S, null all different from empirical !
gplot_effects(cur5)


# Comparison metrics ----------------------------
# focus on species richness
met = "rangeT"

# merge into a single data.frame
cur1$simu = 1
cur2$simu = 2
cur5$simu = 5
scurves <- rbind(
  cur1[cur1$ind == met, ],
  cur2[cur2$ind == met, ],
  cur5[cur5$ind == met, ]
)

# ggplot with facet
p <- ggplot2::ggplot(scurves, ggplot2::aes(x = N, group = scn))
p <- p +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = p10, ymax = p90, fill = scn),
    alpha = 0.10,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = p50, colour = scn),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~simu, nrow = 1) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Individuals sampled (m)",
    y = met,
    colour = NULL,
    fill = NULL,
    title = "Comparison with duplicated data"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggsave(
  p,
  filename = here::here("analysis", "figures", "simulations_duplicate_S.png"),
  width = 8,
  height = 6
)
