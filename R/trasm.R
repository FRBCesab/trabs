#' Create the 'trasm' object.
#'
#' The 'trasm' object is a trait assemblage dataset in long format
#' that allows the calculation of trait-abundance scaling relationship across scales.
#'
#' @param asm data.frame in long format with columns containing species, traits, and abundance information.
#' @param species_var column name containing species information
#' @param trait_var column name containing trait information
#' @param abundance_var column name containing abundance information
#'
#' @return a "trasm" object which is a `data.frame` with three columns: Sp, Tr, and Ab
#'
#' @export
new_trasm <- function(
  asm,
  species_var = "Sp",
  trait_var = "Tr",
  abundance_var = "Ab"
) {
  # Checking the inputs ------------
  stopifnot("`asm` must be a data.frame." = {
    is.data.frame(asm) & nrow(asm) > 1
  })
  stopifnot(
    "`species_var` must be a single string character matching the name of a column in asm." = {
      is.character(species_var) &
        length(species_var) == 1 &
        species_var %in% names(asm)
    }
  )
  stopifnot(
    "`trait_var` must be a single string character matching the name of a column in asm." = {
      is.character(trait_var) &
        length(trait_var) == 1 &
        trait_var %in% names(asm)
    }
  )
  stopifnot(
    "`abundance_var` must be a single string character matching the name of a column in asm." = {
      is.character(abundance_var) &
        length(abundance_var) == 1 &
        abundance_var %in% names(asm)
    }
  )

  # checking variables
  Ab <- asm[, abundance_var]
  stopifnot(
    "abundances must be positive numeric." = {
      is.numeric(Ab) & all(Ab >= 0)
    }
  )
  stopifnot("abundances should not contains NA, nor be empty" = {
    all(!is.na(Ab)) & sum(Ab) > 0
  })
  Tr <- asm[, trait_var]
  stopifnot("traits must be numeric." = is.numeric(Tr))

  sp <- asm[, species_var]
  stopifnot(
    "species must be character and non-uniform" = {
      is.character(sp) & sum(!duplicated(sp)) > 1
    }
  )

  out <- data.frame(
    "Sp" = sp,
    "Tr" = Tr,
    "Ab" = Ab
  )
  class(out) = c("data.frame", "trasm")
  return(out)
}


#' Calculate a short summary of the trasm data
#'
#' @param x `trasm` object
#' @export
summary_trasm <- function(x) {
  data.frame(
    "S_tot" = nodup(x$Sp[x$Ab > 0]),
    "N_tot" = sum(x$Ab),
    "muT" = mean(x$Tr),
    "sdT" = stats::sd(x$Tr),
    beta_hat = stats::coef(stats::lm(log1p(x$Ab) ~ x$Tr))[2] |>
      unname()
  )
}
