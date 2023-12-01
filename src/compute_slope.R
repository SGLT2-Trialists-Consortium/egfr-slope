#' Compute eGFR slope results from a mixed model object
#' 
#' @details 
#' `compute_slope()` computes eGFR slope results from a "two slope" mixed model 
#' object using linear combinations of the coefficients. The function is at its 
#' core a wrapper for `multcomp::glht()` which calculates linear combinations. 
#' `multcomp::glht()` ordinarily requires a contrast vector to be specified
#' manually, but given this wrapper requires only the names of key variables in
#' in your model. See the repository README for more details.
#' 
#' @param .model_obj A linear model object, i.e. `lme4::lmer()` or
#' `nlme::lme()`.
#' @param .intervention_var A character string specifying the name of your 
#' treatment or intervention. This variable must be binary.
#' @param .time_var A character string specifying the name of your time 
#' variable.
#' @param .spline_var A character string specifying the name of your spline
#' variable.
#' @param .by Optional argument that is a character string specifying the name
#' of your subgroup variable. This must be a factor variable.
#' @param .prop A numeric vector that specifies the proportion of the total
#' slope accounted-for by the chronic slope.
#' @param .output A character string specifying your desired output: acute, 
#' chronic, total, or all three slope types.
#' 
#' @return Data frame extension of type `tibble`. 

compute_slope <- function(.model_obj,
                          .intervention_var,
                          .time_var,
                          .spline_var,
                          .by = NULL,
                          .prop = NULL,
                          .output = "all") {
  
  rnd <- function(x, decimals) {
    # Round numbers correctly to a specified number of decimal places. The 
    # `round()` function in `base` uses the IEC 60559 standard which rounds off 
    # 5 to "to the even digit", which is not what we want.
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
  
  process <- function(.multcomp_obj, 
                      .slope, 
                      .group,
                      .subgroup = NULL, 
                      .subgroup_cat = NULL) {
    # Process output from `multcomp::glht()` to `tibble`.
    #
    # Args
    #   .multcomp_obj : object from `multcomp::glht()`
    #   .slope : character ("Acute" or "Chronic" or "Total")
    #   .group : character ("Control" or "Active" or "Active - Control")
    #   .subgroup : character (Name of subgroup variable)
    #   .subgroup_cat : character (Name of subgroup categories)
    # 
    # Returns
    #   tibble
    ci <- confint(.multcomp_obj)
    res <- 
      tibble::tibble(
        slope = .slope,
        group = .group,
        estimate = rnd(.multcomp_obj$test$coefficients, decimals = 4),
        lci = rnd(ci$confint[2], decimals = 4),
        uci = rnd(ci$confint[3], decimals = 4),
        se = .multcomp_obj$test$sigma,
        t_value = .multcomp_obj$test$tstat,
        p_value = rnd(.multcomp_obj$test$pvalues[1], decimals = 4)
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
      dplyr::select(-tidyselect::matches("_ch$"))
    if (!missing(.subgroup)) {
      res <- res |> 
        mutate(
          subgroup = .subgroup, subgroup_cat = .subgroup_cat,
          .after = "group"
        )
      return(res)
    } else {
      return(res)
    }
  }
  
  # Colour functions
  bblck <- function(.t) crayon::bold(crayon::black(.t))
  bred <- function(.t) crayon::bold(crayon::red(.t))
  blck <- function(.t) crayon::black(.t)
  cyn <- function(.t) crayon::cyan(.t)
  
  # Function to display output header 
  display_header <- function(.output_produced) {
    main <- bblck("Slopes computed: ")
    if (.output_produced == "Acute") {
      message(main, bred(.output_produced))
    } else if (.output_produced == "Chronic") {
      message(main, bred(.output_produced))
    } else if (.output_produced == "Total") {
      message(main, bred(.output_produced))
    } else if (.output_produced == "All") {
      message(main, bred("Acute, Chronic, and Total"))
    }
  }
  
  # Function to display linear combinations and contrast vectors
  display_lc <- function(.title, .fe, .cv, .prop = NULL, .multiply = NULL) {
    message(bblck(.title))
    
    # Add `` around each fixed effect
    .fe <- paste("`", .fe, "`", sep="")
    
    # Add multiplication with proportion
    if (!missing(.multiply)) {
      .fe[.multiply] <- paste0("(`", .prop, "` * ", .fe[.multiply], ")")
    }
    
    # Print linear combinations
    if (length(fe) == 1) {
      message(blck("Linear combination: "), cyn(.fe))
    } else if (length(fe) > 1) {
      message(blck("Linear combination: "), cyn(paste(.fe, collapse = " + ")))
    }
    
    # Edit how contrast vector is displayed if not only 0s and 1s
    .cv <- ifelse(
      .cv > 0 & .cv < 1, as.character(rnd(.cv, 2)), as.character(.cv)
    )
    
    # Print contrast vectors
    message(
      blck(paste0("Contrast vector: c(", paste(.cv, collapse = ", "), ")"))
    )
  }
  
  # Get fixed effects and names of fixed effects
  # Fixed effects
  fe <- nlme::fixed.effects(.model_obj)
  # Fixed effects coefficient names
  nfe <- names(fe)
  
  # Names of interactions required for all calculations
  # Treatment and time
  txi <- names(fe[grepl(paste0("^", .time_var, ":", .intervention_var, "$|^", .intervention_var, ":", .time_var, "$"), nfe)])
  # Treatment and spline
  sxi <- names(fe[grepl(paste0("^", .spline_var, ":", .intervention_var, "$|^", .intervention_var, ":", .spline_var, "$"), nfe)])
  
  if (missing(.by)) {
    # Position of coefficients in the contrast vector that require 
    # multiplication for total slope
    sp <- which(names(fixed.effects(.model_obj)) == .spline_var)
    sxip <- which(names(fixed.effects(.model_obj)) == sxi)
    
    # Acute slope contrast vectors
    acute_0 <- as.integer(nfe %in% .time_var)
    acute_1 <- as.integer(nfe %in% c(.time_var, txi))
    acute_d <- as.integer(nfe %in% c(txi))
    
    # Chronic slope contrast vectors
    chronic_0 <- as.integer(nfe %in% c(.time_var, .spline_var))
    chronic_1 <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi))
    chronic_d <- as.integer(nfe %in% c(txi, sxi))
    
    if (!missing(.prop)) {
      # Total slope contrast vectors
      total_0 <- as.integer(nfe %in% c(.time_var, .spline_var))
      total_0[sp] <- total_0[sp] * .prop
      total_1 <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi))
      total_1[c(sp, sxip)] <- total_1[c(sp, sxip)] * .prop
      total_d <- as.integer(nfe %in% c(txi, sxi))
      total_d[sxip] <- total_d[sxip] * .prop
    }

    # Slope type
    slope <- "Total"
    # Groups
    grp <- c("Control", "Active", "Active - Control")
    
    # Produce results for acute slope only
    if (.output %in% c("Acute", "ACUTE", "acute")) {
      # List of contrast vectors
      vec <- list(acute_0, acute_1, acute_d)
      # Output list
      solo <- list()
      
      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(res, .slope = slope[i], .group = grp[i])
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope)
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Acute slope calculations
      message(bred("Acute slope"))
      display_lc(.title = "Control", .fe = .time_var, .cv = acute_0)
      display_lc(.title = "Active", .fe = c(.time_var, txi), .cv = acute_1)
      display_lc(.title = "Active - Control", .fe = txi, .cv = acute_d)
      return(solo)
    } else if (.output %in% c("Chronic", "CHRONIC", "chronic")) {
      # List of contrast vectors
      vec <- list(chronic_0, chronic_1, chronic_d)
      # Output list
      solo <- list()
      
      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(res, .slope = slope[i], .group = grp[i])
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope)
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Chronic slope calculations
      message(bred("Chronic slope"))
      display_lc(
        .title = "Control", .fe = c(.time_var, .spline_var), .cv = chronic_0
      )
      display_lc(
        .title = "Active", .fe = c(.time_var, .spline_var, txi, sxi), 
        .cv = chronic_1
      )
      display_lc(
        .title = "Active - Control", .fe = c(txi, sxi), .cv = chronic_d
      )
    } else if (.output %in% c("Total", "TOTAL", "total")) {
      # List of contrast vectors
      vec <- list(total_0, total_1, total_d)
      # Output list
      solo <- list()
      
      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(res, .slope = slope[i], .group = grp[i])
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope)
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Total slope calculations
      message(bred("Total slope"))
      display_lc(
        .title = "Control", .fe = c(.time_var, .spline_var), .cv = total_0, 
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = "Active", .fe = c(.time_var, .spline_var, txi, sxi), 
        .cv = total_1, .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = "Active - Control", .fe = c(txi, sxi), .cv = total_d, 
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
    # Produce results for all slope types (acute, chronic, and total)
    } else if (.output %in% c("All", "ALL", "all")) {
      # List of contrast vectors
      vec <- list(
        acute_0, acute_1, acute_d, chronic_0, chronic_1, chronic_d, total_0, 
        total_1, total_d
      )
      # Slope types
      slope <- rep(c("Acute", "Chronic", "Total"), each = 3)
      # Groups
      grp <- rep(c("Control", "Active", "Active - Control"), times = 3)
      all <- list()
      
      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        all[[i]] <- process(res, .slope = slope[i], .group = grp[i])
      }
      all <- dplyr::bind_rows(all)
      display_header(.output_produced = "All")
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Acute slope calculations
      message(bred("Acute slope"))
      display_lc(.title = "Control", .fe = .time_var, .cv = acute_0)
      display_lc(.title = "Active", .fe = c(.time_var, txi), .cv = acute_1)
      display_lc(.title = "Active - Control", .fe = txi, .cv = acute_d)
      # Chronic slope calculations
      message(bred("Chronic slope"))
      display_lc(
        .title = "Control", .fe = c(.time_var, .spline_var), .cv = chronic_0
      )
      display_lc(
        .title = "Active", .fe = c(.time_var, .spline_var, txi, sxi), 
        .cv = chronic_1
      )
      display_lc(
        .title = "Active - Control", .fe = c(txi, sxi), .cv = chronic_d
      )
      # Total slope calculations
      message(bred("Total slope"))
      display_lc(
        .title = "Control", .fe = c(.time_var, .spline_var), .cv = total_0, 
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = "Active", .fe = c(.time_var, .spline_var, txi, sxi), 
        .cv = total_1, .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = "Active - Control", .fe = c(txi, sxi), .cv = total_d, 
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      return(all)
    }
    } else if (!missing(.by) & length(levels(eval(summary(.model_obj)$call$data)[[.by]])) == 2) {
    # Highest level of the subgroup variable (required for regex below)
    ref <- levels(.model_obj@frame[[.by]])[2]
    # Other level
    oth <- levels(.model_obj@frame[[.by]])[1]
    
    # Names of subgroup interactions
    # Treatment and subgroup
    sgxi <- names(fe[grepl(paste0("^", .intervention_var, ":", .by, ref, "$|^", .by, ref, ":", .intervention_var, "$"), nfe)])
    # Time and subgroup
    sgxt <- names(fe[grepl(paste0("^", .time_var, ":", .by, ref, "$|^", .by, ref, ":", .time_var, "$"), nfe)])
    # Spline and subgroup
    sgxs <- names(fe[grepl(paste0("^", .spline_var, ":", .by, ref, "$|^", .by, ref, ":", .time_var, "$"), nfe)])
    # 3-way interactions
    int3 <- names(fe[grepl(":[a-z|A-Z|0-9].+:", nfe)])
    # Time and treatment and subgroup
    sgxtxi <- int3[grepl("time", int3)]
    # Treatment and spline and subgroup
    sgxsxi <- int3[grepl("spline", int3)]
    
    # Position of coefficients in the contrast vector that require 
    # multiplication for total slope
    sp <- which(names(fixed.effects(.model_obj)) == .spline_var)
    sxip <- which(names(fixed.effects(.model_obj)) == sxi)
    sgxsp <- which(names(fixed.effects(.model_obj)) == sgxs)
    sgxsxip <- which(names(fixed.effects(.model_obj)) == sgxsxi)
    
    # Acute, subgroup yes
    acute_0_sgy <- as.integer(nfe %in% .time_var)
    acute_1_sgy <- as.integer(nfe %in% c(.time_var, txi))
    acute_d_sgy <- as.integer(nfe %in% c(txi))
    
    # Acute, subgroup no
    acute_0_sgn <- as.integer(nfe %in% c(.time_var, sgxt))
    acute_1_sgn <- as.integer(nfe %in% c(.time_var, sgxt, sgxi, sgxtxi))
    acute_d_sgn <- as.integer(nfe %in% c(txi, sgxtxi))
    
    # Chronic, subgroup yes
    chronic_0_sgy <- as.integer(nfe %in% c(.time_var, .spline_var))
    chronic_1_sgy <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi))
    chronic_d_sgy <- as.integer(nfe %in% c(txi, sxi))
    
    # Chronic, subgroup no
    chronic_0_sgn <- as.integer(nfe %in% c(.time_var, .spline_var, sgxt, sgxs))
    chronic_1_sgn <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi))
    chronic_d_sgn <- as.integer(nfe %in% c(txi, sxi, sgxtxi, sgxsxi))
    
    if (!missing(.prop)) {
      # Total, subgroup yes
      total_0_sgy <- as.integer(nfe %in% c(.time_var, .spline_var))
      total_0_sgy[sp] <- total_0_sgy[sp] * .prop
      total_1_sgy <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi)) 
      total_1_sgy[c(sp, sxip)] <- total_1_sgy[c(sp, sxip)] * .prop
      total_d_sgy <- as.integer(nfe %in% c(txi, sxi))
      total_d_sgy[sxip] <- total_d_sgy[sxip] * .prop
    
      # Total, subgroup no
      total_0_sgn <- as.integer(nfe %in% c(.time_var, .spline_var, sgxt, sgxs))
      total_0_sgn[c(sp, sgxsp)] <- total_0_sgn[c(sp, sgxsp)] * .prop
      total_1_sgn <- as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi))
      total_1_sgn[c(sp, sxip, sgxsp, sgxsxip)] <- total_1_sgn[c(sp, sxip, sgxsp, sgxsxip)] * .prop
      total_d_sgn <- as.integer(nfe %in% c(txi, sxi, sgxtxi, sgxsxi))
      total_d_sgn[c(sxip, sgxsxip)] <- total_d_sgn[c(sxip, sgxsxip)] * .prop
    }
    
    # Groups
    grp <- rep(c("Control", "Active", "Active - Control"), times = 2)
    subgrp <- rep(.by, times = 6)
    subgrp_cat <- rep(c(oth, ref), each = 3)
    
    # Produce results for acute slope only
    if (.output %in% c("Acute", "ACUTE", "acute")) {
      # List of contrast vectors
      vec <- list(
        acute_0_sgy, acute_1_sgy, acute_d_sgy, acute_0_sgn, acute_1_sgn,
        acute_d_sgn
      )
      # Slope type
      slope <- rep("Acute", times = 6)
      # Output list
      solo <- list()

      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(
          res, .slope = slope[i], .group = grp[i], .subgroup = subgrp[i],
          .subgroup_cat = subgrp_cat[i]
        )
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope[1])
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Acute slope calculations
      message(bred("Acute slope"))
      display_lc(.title = paste0("Control, ", .by, " = ", oth), .fe = .time_var, .cv = acute_0_sgy)
      display_lc(.title = paste0("Active, ", .by, " = ", oth), .fe = c(.time_var, txi), .cv = acute_1_sgy)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", oth), .fe = txi, .cv = acute_d_sgy)
      display_lc(.title = paste0("Control, ", .by, " = ", ref), .fe = c(.time_var, sgxt), .cv = acute_0_sgn)
      display_lc(.title = paste0("Active, ", .by, " = ", ref), .fe = c(.time_var, sgxt, sgxi, sgxtxi), .cv = acute_1_sgn)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", ref), .fe = c(txi, sgxtxi), .cv = acute_d_sgn)
      return(solo)
    } else if (.output %in% c("Chronic", "CHRONIC", "chronic")) {
      # List of contrast vectors
      vec <- list(
        chronic_0_sgy, chronic_1_sgy, chronic_d_sgy, chronic_0_sgn, chronic_1_sgn,
        chronic_d_sgn
      )
      # Slope type
      slope <- rep("Chronic", times = 6)
      # Output list
      solo <- list()

      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(
          res, .slope = slope[i], .group = grp[i], .subgroup = subgrp[i],
          .subgroup_cat = subgrp_cat[i]
        )
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope[1])
      message(bred("Chronic slope"))
      display_lc(.title = paste0("Control, ", .by, " = ", oth), .fe = c(.time_var, .spline_var), .cv = chronic_0_sgy)
      display_lc(.title = paste0("Active, ", .by, " = ", oth), .fe = c(.time_var, .spline_var, txi, sxi), .cv = chronic_1_sgy)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", oth), .fe = c(txi, sxi), .cv = chronic_d_sgy)
      display_lc(.title = paste0("Control, ", .by, " = ", ref), .fe = c(.time_var, .spline_var, sgxt, sgxs), .cv = chronic_0_sgn)
      display_lc(.title = paste0("Active, ", .by, " = ", ref), .fe = c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi), .cv = chronic_1_sgn)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", ref), .fe = c(txi, sxi, sgxtxi, sgxsxi), .cv = chronic_d_sgn)
      return(solo)
    } else if (.output %in% c("Total", "TOTAL", "total")) {
      # List of contrast vectors
      vec <- list(
        total_0_sgy, total_1_sgy, total_d_sgy, total_0_sgn, total_1_sgn,
        total_d_sgn
      )
      # Slope type
      slope <- rep("Total", times = 6)
      # Output list
      solo <- list()

      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        solo[[i]] <- process(
          res, .slope = slope[i], .group = grp[i], .subgroup = subgrp[i],
          .subgroup_cat = subgrp_cat[i]
        )
      }
      solo <- dplyr::bind_rows(solo)
      display_header(.output_produced = slope[1])
      message(bred("Total slope"))
      display_lc(
        .title = paste0("Control, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var), .cv = total_0_sgy,
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = paste0("Active, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, txi, sxi), .cv = total_1_sgy,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = paste0("Active - Placebo, ", .by, " = ", oth),
        .fe = c(txi, sxi), .cv = total_d_sgy,
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = paste0("Control, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, sgxt, sgxs), .cv = total_0_sgn,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = paste0("Active, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi),
        .cv = total_1_sgn, .prop = deparse(substitute(.prop)),
        .multiply = c(2, 4, 6, 8)
      )
      display_lc(
        .title = paste0("Active - Placebo, ", .by, " = ", oth),
        .fe = c(txi, sxi, sgxtxi, sgxsxi), .cv = total_d_sgn,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      return(solo)
    } else if (.output %in% c("All", "ALL", "all")) {
      # List of contrast vectors
      vec <- list(
        acute_0_sgy, acute_1_sgy, acute_d_sgy, acute_0_sgn, acute_1_sgn,
        acute_d_sgn, chronic_0_sgy, chronic_1_sgy, chronic_d_sgy, chronic_0_sgn, 
        chronic_1_sgn, chronic_d_sgn, total_0_sgy, total_1_sgy, total_d_sgy, 
        total_0_sgn, total_1_sgn, total_d_sgn
      )
      # Slope type
      slope <- rep(c("Acute", "Chronic", "Total"), each = 6)
      # Groups
      grp <- rep(c("Control", "Active", "Active - Control"), times = 6)
      subgrp <- rep(.by, times = 18)
      subgrp_cat <- rep(c(oth, oth, oth, ref, ref, ref), times = 3)
      
      # Output list
      all <- list()

      # Loop over contrast vectors to produce results
      for (i in 1:length(vec)) {
        res <- summary(
          multcomp::glht(.model_obj, linfct = rbind("slope" = vec[[i]]))
        )
        # Make tibble
        all[[i]] <- process(
          res, .slope = slope[i], .group = grp[i], .subgroup = subgrp[i],
          .subgroup_cat = subgrp_cat[i]
        )
      }
      all <- dplyr::bind_rows(all)
      display_header(.output_produced = "All")
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
      # Acute slope calculations
      message(bred("Acute slope"))
      display_lc(.title = paste0("Control, ", .by, " = ", oth), .fe = .time_var, .cv = acute_0_sgy)
      display_lc(.title = paste0("Active, ", .by, " = ", oth), .fe = c(.time_var, txi), .cv = acute_1_sgy)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", oth), .fe = txi, .cv = acute_d_sgy)
      display_lc(.title = paste0("Control, ", .by, " = ", ref), .fe = c(.time_var, sgxt), .cv = acute_0_sgn)
      display_lc(.title = paste0("Active, ", .by, " = ", ref), .fe = c(.time_var, sgxt, sgxi, sgxtxi), .cv = acute_1_sgn)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", ref), .fe = c(txi, sgxtxi), .cv = acute_d_sgn)
      # Chronic slope calculations
      message(bred("Chronic slope"))
      display_lc(.title = paste0("Control, ", .by, " = ", oth), .fe = c(.time_var, .spline_var), .cv = chronic_0_sgy)
      display_lc(.title = paste0("Active, ", .by, " = ", oth), .fe = c(.time_var, .spline_var, txi, sxi), .cv = chronic_1_sgy)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", oth), .fe = c(txi, sxi), .cv = chronic_d_sgy)
      display_lc(.title = paste0("Control, ", .by, " = ", ref), .fe = c(.time_var, .spline_var, sgxt, sgxs), .cv = chronic_0_sgn)
      display_lc(.title = paste0("Active, ", .by, " = ", ref), .fe = c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi), .cv = chronic_1_sgn)
      display_lc(.title = paste0("Active - Control, ", .by, " = ", ref), .fe = c(txi, sxi, sgxtxi, sgxsxi), .cv = chronic_d_sgn)
      # Total slope calculations
      message(bred("Total slope"))
      display_lc(
        .title = paste0("Control, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var), .cv = total_0_sgy,
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = paste0("Active, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, txi, sxi), .cv = total_1_sgy,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = paste0("Active - Placebo, ", .by, " = ", oth),
        .fe = c(txi, sxi), .cv = total_d_sgy,
        .prop = deparse(substitute(.prop)), .multiply = 2
      )
      display_lc(
        .title = paste0("Control, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, sgxt, sgxs), .cv = total_0_sgn,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      display_lc(
        .title = paste0("Active, ", .by, " = ", oth),
        .fe = c(.time_var, .spline_var, txi, sxi, sgxt, sgxs, sgxtxi, sgxsxi),
        .cv = total_1_sgn, .prop = deparse(substitute(.prop)),
        .multiply = c(2, 4, 6, 8)
      )
      display_lc(
        .title = paste0("Active - Placebo, ", .by, " = ", oth),
        .fe = c(txi, sxi, sgxtxi, sgxsxi), .cv = total_d_sgn,
        .prop = deparse(substitute(.prop)), .multiply = c(2, 4)
      )
      return(all)
    }
  }
}