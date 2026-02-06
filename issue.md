## to be checked

## data inputs:

**long format** with columns `trait`, `abundance`, `species`, (`sites`, `coordinates` in option)
*advantages*: possibility to integrate intra-specific trait variation
*disadvantage*: not the usual community matrix used in `vegan` or `mobr`
- do we allow NA in abundance? or in traits?
- no NA allowed in species, sites or coordinates when present
- need easy function to transform data from `vegan` and `mobr` package
- so far only one single trait, do we want to include multiple traits?

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
covTA_viaTR = EMP - NULL_SAD ? 
pct_covTA_viaTR_signed = 100 * covTA_viaTR / denom_delta_signed,
    covTA_viaTR = EMP - NULL_SAD, # alternative path

### make_steps() 

- not currently matching the number of steps because in log/sqrt scale it could be the same integer. Is it an issue? we should clarify it in the documentation
``` r
length(make_steps(m_use, n_points = 100, mode = step_scale))
```
- naming: `make_steps()` is very similar to `seq()` in R, so maybe rename to `seq_scale()`


### compute_curve_from_order()
- renamed to get_curve_from_order
- removed within_sd parameter
- similar to `vegan::rarefy()` function
- to be optimized?

### summary_trasw
Should be weigthed indicators (weigthed mean and sd).
same for beta_cov

## how to compute rarefaction curves efficiently?
https://vegandevs.github.io/vegan/reference/rarefy.html
https://vegandevs.github.io/vegan/reference/specaccum.html

specaccum
test vs sa_all <- vegan::specaccum(wide) #, method="random"

