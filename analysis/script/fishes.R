## ----setup--------------------------------------------------------------------
N_PERM <- 999 # permutations per assemblage for curves
N_POINTS_CURVE <- 100 # points along the curve
MAX_INDIV <- 10000 # cap individuals used per curve (speed/safety)
SEED <- 42
set.seed(SEED)
# run in 30sec

devtools::load_all()
library(ggplot2)

dat <- read.csv(here::here("analysis", "data", "Fish_Mean.csv"))
#  filter(province == "samoa") %>%
dat <- dat[grep("^samoa_", dat$Site), ]
dat$province <- first_str(dat$Site)
dat$sub_province <- second_str(dat$Site)
dat$site <- last_str(dat$Site)
dat$logT <- log10(dat$Body_Mass)

curves <- list()
for (i in sort(unique(dat$sub_province))) {
  print(i)
  subdat <- dat[dat$sub_province == i, ]
  trasmi <- new_trasm(subdat, "Binomial", "logT", "N")
  out <- data.frame(
    "sub_province" = i,
    summary_trasm(trasmi),
    get_rarefy_curves(
      trasmi,
      n_perm = N_PERM,
      n_points = N_POINTS_CURVE
    ),
    row.names = NULL
  )
  curves[[i]] <- out
}

# merge the curves into a large data.frame
curves_df <- do.call(rbind, curves)

# calculate the effect size per simulation
part_list <- lapply(curves, partition, output = "effect")
part_df <- do.call(rbind, part_list)

# save output
fishes <- list(
  "curves" = curves_df,
  "effects" = part_df
)
# save the simulation
saveRDS(fishes, "analysis/data/fishes_samoa_curves_effects.rds")

# load the previously run simulations
# fishes <- readRDS("analysis/data/fishes_curves_effects.rds")

## Plot curves -----------------------------------------------------------------
for (i in c("S", "rangeT", "CWM", "sumT")) {
  sub_curves <- fishes$curves[fishes$curves$ind == i, ]

  p <- ggplot(
    sub_curves,
    aes(x = N, y = p50, group = scn)
  ) +
    # geom_errorbar(aes(ymin = p10, ymax = p90, y = p50, colour = as.factor(sad_k))) +
    # geom_ribbon(aes(ymin = p10, ymax = p90,y = p50, fill = as.factor(scn)), alpha = 0.1) +
    geom_line(aes(y = p50, col = as.factor(scn)), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_x_log10() +
    labs(x = "Individuals sampled (m)", y = i) +
    facet_wrap(~sub_province)
  nameout <- paste0("samoa_curves_", i, ".png")
  ggsave(
    p,
    filename = here::here("analysis", "figures", nameout),
    width = 8,
    height = 6
  )
}


## Plot effects ---------------------------------------------------

for (i in c("S", "rangeT", "CWM", "sumT")) {
  sub_effect <- fishes$effects[fishes$effects$ind == i, ]

  p <- ggplot(sub_effect) +
    # geom_ribbon(aes(x = m, ymin = lo, ymax = hi, group = as.factor(paste0(rep, effect)),fill =as.factor(effect)),alpha = 0.01)+
    geom_line(aes(
      x = N,
      y = value.med,
      group = as.factor(part),
      colour = as.factor(part)
    )) +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    facet_wrap(~sub_province) +
    scale_x_log10() +
    # scale_color_gradient2()+
    # scale_color_brewer(palette = palette, drop = FALSE) +
    # scale_fill_brewer(palette = palette, drop = FALSE) +
    labs(
      x = "Individuals sampled (m)",
      y = paste0("effect on ", i)
    )
  nameout <- paste0("samoa_effects_", i, ".png")
  ggsave(
    p,
    filename = here::here("analysis", "figures", nameout),
    width = 8,
    height = 6
  )
}
