## to be checked

## data inputs:

**long format** with columns `trait`, `abundance`, (`species`, `sites`, `coordinates` in option)
*advantages*: possibility to integrate intra-specific trait variation
*disadvantage*: not the usual community matrix used in `vegan` or `mobr`
- do we allow NA in abundance? or in traits?
- no NA allowed in species, sites or coordinates when present
- need to have easy function to transform data from vegan or mobr package
- so far only one single trait

**trasm** new R-class 
- Species, Trait, Abundance : in that order
- to be improved: add sites, coordinates and projection system
- add conversion function from mobr and vegan

## Functions

### why scale_traits()
do we need it? does that change the estimates? `SCALE_T`?
flat : for individuals, might be trickier ... should we take unique values of traits? unique combinaison of trait x species?

### do we need `VAR_T` as indicator and ?

### partition over m
covTA_viaTR = EMP - NULL_SAD ? pct_covTA_viaTR_signed = 100 * covTA_viaTR / denom_delta_signed,
    covTA_viaTR = EMP - NULL_SAD, # alternative path

### make_steps() 

- not currently matching the numberof steps because in log/sqrt scale it could be the same integer. Is it an issue? we should clarify it in the documentation
``` r
length(make_steps(m_use, n_points = 100, mode = step_scale))
```
- naming: `make_steps()` is very similar to `seq()` in R, so maybe rename to `seq_scale()`


### compute_curve_from_order()
- renamed to get_curve_from_order
- removed within_sd parameter
- similar to `vegan::rarefy()` function
- to be optimized?



## how to compute rarefaction curves efficiently?
https://vegandevs.github.io/vegan/reference/rarefy.html
https://vegandevs.github.io/vegan/reference/specaccum.html

specaccum
test vs 
sa_all <- vegan::specaccum(wide) #, method="random"


spe_acc_lapply <- function(taxa, site, nseq = 2:min(c(100,nodup(site)/2)), nrep = 100, probs = c(0.25,0.5, 0.75)){
  out <- c()
  usite <- unique(site)
  if (max(nseq)>length(usite)){
    stop("Number of sites lower than max(nseq).")
  }
  out <- lapply(nseq, function(i) {
    ni <- sapply(1:nrep, function(x) nodup(taxa[site%in%sample(usite, i, replace=FALSE)]))
    return(quantile(ni, probs = probs)) 
  })
  out <- data.frame(t(sapply(out,c)))
  names(out) <- paste0("Q",probs)
  row.names(out) <- nseq
  return(out)
}

spe_acc_mclapply <- function(taxa, site, 
                             nseq = 2:min(c(100,nodup(site)/2)), 
                             nrep = 100, probs = c(0.25,0.5, 0.75), 
                             n_cores = 1){
  out <- c()
  usite <- unique(site)
  if (max(nseq)>length(usite)){
    stop("Number of sites lower than max(nseq).")
  }
  out <- parallel::mclapply(nseq, function(i) {
    ni <- sapply(1:nrep, function(x) nodup(taxa[site%in%sample(usite, i, replace=FALSE)]))
    return(quantile(ni, probs = probs)) 
  }, mc.cores = n_cores)
  out <- data.frame(t(sapply(out,c)))
  names(out) <- paste0("Q",probs)
  row.names(out) <- nseq
  return(out)
}

n_cores <- 1 # 8
# compute for all
system.time(
    {sa_all <- spe_acc_mclapply(long$species, long$observationID, 



