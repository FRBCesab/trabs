## ----setup--------------------------------------------------------------------
devtools::load_all()
library(ggplot2)

# Optional global S and N for all runs
S_global <- 200
N_global <- 10000
trait_mean_global <- 6
N_perm <- 999
N_points <- 100
N_rep <- 20
metrics <- c("S", "rangeT", "CWM", "sumT")

# Simulation parameters (edit to taste)
grid <- expand.grid(
  trait_sd = c(0.1, 1, 2),
  sad_k = c(0.5, 1, 2), # smaller -> more uneven SAD, needed?
  beta_cov = c(-1.0, 0.0, 1.0),
  rep = 1:N_rep # replicates per combination, needed?
)

# run in xxmins

## ----simulations ------------------------------------------------------------
library(futurize)
plan(multisession)

curves <- lapply(1:nrow(grid), function(i) {
  asm <- simulate_assemblage(
    S = S_global,
    N = N_global,
    trait_mean = trait_mean_global,
    trait_sd = grid$trait_sd[i],
    sad_k = grid$sad_k[i],
    beta_cov = grid$beta_cov[i]
  )
  out <- data.frame(
    grid[i, ],
    summary_trasm(asm),
    get_rarefy_curves(
      asm,
      n_perm = N_perm,
      n_points = N_points,
      indices = metrics
    ),
    row.names = NULL
  )
  return(out)
}) |>
  futurize::futurize(future.seed = TRUE)

# merge the curves into a large data.frame
curves_df <- do.call(rbind, curves)

# calculate the effect size per simulation
part_list <- lapply(curves, partition, output = "effect")
part_df <- do.call(rbind, part_list)

# save output
simu <- list(
  "curves" = curves_df,
  "effects" = part_df
)
# save the simulation
saveRDS(simu, "analysis/data/simulations_curves_effects.rds")

# load the previously run simulations
# simu <- readRDS("analysis/data/simulations_curves_effects.rds")

## Plot curves -----------------------------------------------------------------
for (i in metrics) {
  sub_curves <- simu$curves[simu$curves$ind == i & simu$curves$sad_k == 1, ]

  p <- ggplot(
    sub_curves,
    aes(x = N, y = p50, group = interaction(rep, scn, trait_sd))
  ) +
    # geom_errorbar(aes(ymin = p10, ymax = p90, y = p50, colour = as.factor(sad_k))) +
    # geom_ribbon(aes(ymin = p10, ymax = p90,y = p50, fill = as.factor(scn)), alpha = 0.1) +
    geom_line(aes(y = p50, col = as.factor(scn)), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    labs(x = "Individuals sampled (m)", y = i) +
    facet_grid(
      sad_k ~ beta_cov,
      labeller = labeller(.rows = label_both, .cols = label_both)
    )
  nameout <- paste0("simulations_curves_", i, ".png")
  ggsave(
    p,
    filename = here::here("analysis", "figures", nameout),
    width = 8,
    height = 6
  )
}


## Plot effects ---------------------------------------------------

for (i in metrics) {
  sub_effect <- simu$effects[simu$effects$ind == i & simu$effects$sad_k == 1, ]

  p <- ggplot(sub_effect) +
    # geom_ribbon(aes(x = m, ymin = lo, ymax = hi, group = as.factor(paste0(rep, effect)),fill =as.factor(effect)),alpha = 0.01)+
    geom_line(aes(
      x = N,
      y = value.med,
      group = as.factor(paste0(rep, part)),
      colour = as.factor(part)
    )) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    facet_grid(
      trait_sd ~ beta_cov,
      labeller = labeller(.rows = label_both, .cols = label_both)
    ) +
    scale_x_log10() +
    # scale_color_gradient2()+
    # scale_color_brewer(palette = palette, drop = FALSE) +
    # scale_fill_brewer(palette = palette, drop = FALSE) +
    labs(
      x = "Individuals sampled (m)",
      y = paste0("effect on ", i)
    )
  nameout <- paste0("simulations_effects_", i, ".png")
  ggsave(
    p,
    filename = here::here("analysis", "figures", nameout),
    width = 8,
    height = 6
  )
}
