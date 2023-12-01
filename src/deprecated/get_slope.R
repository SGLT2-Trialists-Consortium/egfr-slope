#' NOTE: This function has been deprecated.
#' 
#' Extract eGFR slope results from a mixed model object
#' 
#' @details 
#' `get_slope()` extracts eGFR slope results from a "two slope" mixed model 
#' object using linear combinations of the coefficients. The function is 
#' currently quite cumbersome, as it uses the `multcomp::glht()` function for
#' its calculation of linear combinations. This function requires a contrast 
#' vector to be specified manually. See the repository README for more details.
#' 
#' @param .model_obj A linear model object, i.e. `lme4::lmer()` or
#' `nlme::lme()`.
#' @param .name A character string specifying the type of output you require, 
#' i.e. "Acute slope: Placebo" or "Chronic slope: SGLT2i" etc.
#' @param .contrasts A contrast vector for the calculation of linear
#' combinations. This can be supplied either as a numeric vector or as a
#' string vector specifying the fixed effects required for the calculation.
#' 
#' @return data frame extension of type `tibble`. 

get_slope <- function(.model_obj, .name, .contrasts) {
  
  rnd <- function(x, decimals) {
    # Round numbers correctly to a specified number of decimal places. The 
    # `round()` function in `base` uses the IEC 60559 standard which rounds off 
    # 5 to "to the even digit", which is not what we want
    #
    # Args
    #   x : double (number to round)
    #   decimals : integer (number of decimal places to round to)
    #
    # Returns
    #   tibble
    pos_neg <- sign(x)
    y <- abs(x) * 10^decimals
    y <- y + 0.5 + sqrt(.Machine$double.eps)
    y <- trunc(y)
    y <- y / 10^decimals
    y * pos_neg
  }
  
  if (!is.numeric(.contrasts)) {
    .contrasts <- as.integer(names(fixed.effects(.model_obj)) %in% .contrasts)
  }
  
  mdl <- summary(
    multcomp::glht(.model_obj, linfct = rbind("slope" = .contrasts))
  )
  ci <- confint(mdl)
  res <- 
    tibble::tibble(
      name = .name,
      estimate = rnd(mdl$test$coefficients, decimals = 4),
      lci = rnd(ci$confint[2], decimals = 4),
      uci = rnd(ci$confint[3], decimals = 4),
      se = mdl$test$sigma,
      t_value = mdl$test$tstat,
      p_value = rnd(mdl$test$pvalues[1], decimals = 4)
    ) |> 
    dplyr::mutate(
      dplyr::across(
        c(estimate, lci, uci), \(x) as.character(rnd(x, decimals = 2)),
        .names = "{col}_ch"
      ),
      dplyr::across(
        dplyr::ends_with("_ch"), 
        \(x) dplyr::case_when(
          stringr::str_detect(x, "[.][0-9]$") ~ paste0(x, "0"),
          stringr::str_detect(x, "^[0-9]+$") ~ paste0(x, ".00"),
          TRUE ~ x
        )
      ),
      result = paste0(estimate_ch, " (", lci_ch, " to ", uci_ch, ")"),
      .after = "uci"
    ) |> 
    dplyr::select(-matches("_ch$"))
  message(crayon::black(paste0("Contrast vector supplied: c(", paste(.contrasts, collapse = ", "), ")")))
  return(res)
}