#' @title 
#' Compute eGFR Slope Results from a "Two Slope" Mixed Model
#' 
#' @description
#' This function computes acute, chronic, or total eGFR slope estimates (or all
#' three) based on a "two slope" mixed model (e.g., from `lme4::lmer()` or 
#' `nlme::lme()`). It wraps \code{multcomp::glht()} to generate linear 
#' combinations of model coefficients, removing the need to specify contrast 
#' vectors manually. You only need to provide the names of key model variables.
#' 
#' @details
#' - **Acute slope** typically represents the initial, short-term change in eGFR.
#' - **Chronic slope** reflects the longer-term eGFR trajectory post-acute phase.
#' - **Total slope** is a weighted combination of acute and chronic slopes, with
#'   weighting defined by \code{.prop}.
#'   
#' For subgroup analyses, specify a factor variable in \code{.by} to obtain 
#' separate slope estimates by level (e.g., men vs. women, or levels of baseline
#' eGFR). The function will also perform a chi-squared heterogeneity test to 
#' assess differences across subgroups.
#' 
#' @param .model_obj An \emph{fitted} two-slope model object, typically an
#'   object from \code{lme4::lmer()} or \code{nlme::lme()}.
#' @param .time_var A character string specifying the \strong{time} variable in 
#'   \code{.model_obj}, representing the linear term for time.
#' @param .intervention_var A character string specifying the binary treatment 
#'   or intervention variable in \code{.model_obj}.
#' @param .spline_var A character string specifying the \strong{spline} variable
#'   used to model the breakpoint between acute and chronic slope.
#' @param .by (Optional) A character string specifying the name of a factor 
#'   variable for subgroup analyses. When provided, separate slopes are computed 
#'   for each subgroup level, along with a heterogeneity test.
#' @param .prop (Optional) A numeric value indicating the proportion of the 
#'   \strong{total} slope allocated to the chronic phase. This is used to form 
#'   linear combinations for the total slope (when \code{.output = "total"} or 
#'   \code{"all"}).
#' @param .output A character string specifying the slope(s) to return. Valid 
#'   options: \code{"acute"}, \code{"chronic"}, \code{"total"}, or \code{"all"}.
#' 
#' @return A tibble with columns for slope type, group (control, active, or 
#'   difference), estimates, confidence intervals, and p-values. If \code{.by} 
#'   is supplied with multiple subgroup levels, the output also includes 
#'   subgroup identifiers and heterogeneity test results.
#'   
#' @section Warning:
#' This function assumes a specific model structure (e.g., an acute phase 
#' spline and interaction terms for treatment). Ensure the naming of variables
#' in \code{.model_obj} matches the parameters provided here.
#' 
#' 
#' @examples
#' \dontrun{
#' # Fit a two-slope model with lme4 or nlme, then:
#' results <- compute_slope(
#'   .model_obj        = fit,
#'   .time_var         = "time",
#'   .intervention_var = "trt01pn",
#'   .spline_var       = "spline",
#'   .prop             = 0.21
#'   .output           = "all"
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{multcomp::glht}} for details on computing linear contrasts.
#'   \item \code{\link{lme4::lmer}} or \code{\link{nlme::lme}} for model fitting.
#' }
#'
#' @export

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
  
  # Needed later on to match contrast vectors to linear combinations
  match_vars <- function(vars, indicator) {
    # Check if the length of the vars and indicator vectors are the same
    if (length(vars) != length(indicator)) {
      stop("Length of vars and indicator vectors must be the same.")
    }
    
    # Select variables that correspond to 1s in the indicator vector
    selected_vars <- vars[indicator == 1]
    
    return(selected_vars)
  }
  
  # To get levels of factor variable for subgroup analyses
  get_model_frame <- function(.model_obj) {
    if (inherits(.model_obj, "lmerMod")) {
      # For lme4 models
      return(.model_obj@frame)
    } else if (inherits(.model_obj, "lme")) {
      # For lme models
      return(nlme::getData(.model_obj))
    } else {
      stop("Unsupported model type. Please provide an lme4 or lme model object.")
    }
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
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
    } else if (.output_produced == "Chronic") {
      message(main, bred(.output_produced))
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
    } else if (.output_produced == "Total") {
      message(main, bred(.output_produced))
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
    } else if (.output_produced == "All") {
      message(main, bred("Acute, Chronic, and Total"))
      message(blck("Output collated in dataframe extension of type `tibble`."))
      message(blck("--------------------------------------------------------"))
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
  
  # Function to get heterogeneity p-value for difference across subgroups
  get_heterogeneity <- function(.result) {
    
    result <- .result[grep("[-]", .result$group), ]
    
    # Create matrix from the estimates
    bhat <- result$estimate
    
    # Create matrix from the standard errors
    sevec <- result$se
    
    # Matrix multiplication
    smat <- diag(sevec^2)
    
    # Define matrix according to the number of subgroups
    
    # Define the length of the vector
    n <- nrow(result) - 1
    
    # Initialise empty matrix with n rows and n columns
    l <- matrix(0, nrow = n, ncol = n + 1)
    
    # Loop through the positions of the vector
    for (i in 1:n) {
      # Create a vector of 0s of length n
      vec <- c(-1, rep(0, n))
      
      # Set the i-th position to 1
      vec[i + 1] <- 1
      
      # Add the vector as the i-th row of the matrix
      l[i, ] <- vec
    }
    
    # Calculate w
    w <- t(l %*% bhat) %*% solve(l %*% smat %*% t(l)) %*% (l %*% bhat)
    
    # Calculate degrees of freedom (df)
    df <- nrow(l)
    
    # Calculate p-value using the chi-square distribution
    pval <- 1 - pchisq(w, df)
    
    # Assign Chi-square statistic
    chi <- w
    
    # Print the results
    message(blck("--------------------------------------------------------"))
    message(bblck("Heterogeneity test"))
    message(blck("Chi-square statistic: "), bred(chi))
    message(blck("p-value: "), bred(pval))
  }
  
  # Get fixed effects and names of fixed effects
  # Fixed effects
  fe <- nlme::fixed.effects(.model_obj)
  # Fixed effects coefficient names
  nfe <- names(fe)
  
  # Names of interactions required for all calculations
  # Treatment and time
  txi <- names(
    fe[
      grepl(
        paste0(
          "^", .time_var, ":", .intervention_var, "$|^", .intervention_var, 
          ":", .time_var, "$"
        ), 
        nfe
      )
    ]
  )
  # Treatment and spline
  sxi <- names(
    fe[
      grepl(
        paste0(
          "^", .spline_var, ":", .intervention_var, "$|^", .intervention_var, 
          ":", .spline_var, "$"
        ),
        nfe
      )
    ]
  )
  
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

    # Groups
    grp <- c("Control", "Active", "Active - Control")
    
    # Produce results for acute slope only
    if (.output %in% c("Acute", "ACUTE", "acute")) {
      # Slope type
      slope <- rep("Acute", 3)
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
      display_header(.output_produced = slope[1])
      # Acute slope calculations
      message(bred("Acute slope"))
      display_lc(.title = "Control", .fe = .time_var, .cv = acute_0)
      display_lc(.title = "Active", .fe = c(.time_var, txi), .cv = acute_1)
      display_lc(.title = "Active - Control", .fe = txi, .cv = acute_d)
      return(solo)
    } else if (.output %in% c("Chronic", "CHRONIC", "chronic")) {
      # Slope type
      slope <- rep("Chronic", 3)
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
      display_header(.output_produced = slope[1])
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
      return(solo)
    } else if (.output %in% c("Total", "TOTAL", "total")) {
      # Slope type
      slope <- rep("Total", 3)
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
      display_header(.output_produced = slope[1])
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
      return(solo)
    # Produce results for all slope types (acute, chronic, and total)
    } else if (.output %in% c("All", "ALL", "all")) {
      # Slope types
      slope <- rep(c("Acute", "Chronic", "Total"), each = 3)
      # List of contrast vectors
      vec <- list(
        acute_0, acute_1, acute_d, chronic_0, chronic_1, chronic_d, total_0, 
        total_1, total_d
      )
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
  } else if (!missing(.by) & length(levels(eval(summary(.model_obj)$call$data)[[.by]])) >= 2) {
    # Levels of the subgroup variable
    levels_by <- levels(get_model_frame(.model_obj)[[.by]])

    # Function to get interaction names
    get_int <- function(level, var_type) {
      rgx <- paste0(
        "^", var_type, ":", .by, level, "$|^", .by, level, ":", var_type, "$"
      )
      names(fe[grepl(rgx, nfe)])
    }
    
    # List of variable types for obtaining interactions
    var_types <- list(
      intervention = .intervention_var, time = .time_var, spline = .spline_var
    )
    
    # Lists for each interaction type
    list_int <- list()
    
    for (i in seq_along(levels_by)) {
      
      sgxi_val <- get_int(levels_by[i], var_types$intervention)
      sgxt_val <- get_int(levels_by[i], var_types$time)
      sgxs_val <- get_int(levels_by[i], var_types$spline)
      
      # Add to the list of lists only if the values are not NULL or empty
      if ((!is.null(sgxi_val) && length(sgxi_val) > 0) ||
          (!is.null(sgxt_val) && length(sgxt_val) > 0) ||
          (!is.null(sgxs_val) && length(sgxs_val) > 0)) {
        
        list_int[[i]] <- list(
          sgxi = sgxi_val, sgxt = sgxt_val, sgxs = sgxs_val
        )
      }
    }
    
    # Remove NULL values if they exist in the list of interactions
    list_int <- Filter(Negate(is.null), list_int)
    
    # 3-way interactions
    int3 <- names(fe[grepl(":[a-z|A-Z|0-9].+:", nfe)])
    # Time and treatment and subgroup
    sgxtxi <- as.list(int3[grepl(.time_var, int3)])
    # Treatment and spline and subgroup
    sgxsxi <- as.list(int3[grepl(.spline_var, int3)])
    
    # Add 3-way interactions to the main list
    for (i in 1:length(list_int)) {
      list_int[[i]]$sgxtxi <- sgxtxi[[i]]
      list_int[[i]]$sgxsxi <- sgxsxi[[i]]
    }
    
    # Create list for values that contain the .spline_var string (need). This is
    # to then create a list of all the positions of coefficients in the contrast 
    # vector that require multiplication for total slope
    list_spl <- list()
    
    # Loop through each element in list_int
    for (i in seq_along(list_int)) {
      
      # Check each sub-list in list_int (e.g., sgxi, sgxt, sgxs)
      for (sublist_name in names(list_int[[i]])) {
        sublist_values <- list_int[[i]][[sublist_name]]
        
        if (any(grepl(.spline_var, sublist_values))) {
          # Add the sublist to the new list if it contains the name of the 
          # spline variable
          list_spl[[length(list_spl) + 1]] <- sublist_values
        }
      }
    }
    
    # Definitive spline list
    list_spl <- append(list(.spline_var, sxi), list_spl)
    
    # Get the positions of these fixed effects/coefficients
    list_spl_pos <- list()
    for (i in 1:length(list_spl)) {
      val <- which(names(fixed.effects(.model_obj)) == list_spl[[i]])
      list_spl_pos[[i]] <- val
      names(list_spl_pos)[i] <- list_spl[[i]]
    }
    
    # Create contrast vectors for acute slope
    acute_contrasts_list <- list()
    for (i in 1:length(list_int)) {
      acute_0 <- as.integer(nfe %in% c(.time_var, list_int[[i]]$sgxt))
      acute_1 <- as.integer(
        nfe %in% c(.time_var, txi, list_int[[i]]$sgxt, list_int[[i]]$sgxtxi)
      )
      acute_d <- as.integer(nfe %in% c(txi, list_int[[i]]$sgxtxi))
      
      acute <- list(acute_0, acute_1, acute_d)
      names(acute) <- c("acute_0", "acute_1", "acute_d")
      
      acute_contrasts_list[[i]] <- acute
      names(acute_contrasts_list)[i] <- paste0("sg", as.character(i))
    }
    
    # Contrasts for the reference group
    ref_contrasts_list <- list(
      list(
        acute_0 = as.integer(nfe %in% c(.time_var)),
        acute_1 = as.integer(nfe %in% c(.time_var, txi)),
        acute_d = as.integer(nfe %in% c(txi))
      ) 
    )
    names(ref_contrasts_list) <- paste0("sg", length(acute_contrasts_list) + 1)
    
    # Bind all together
    acute_contrasts_list <- append(acute_contrasts_list, ref_contrasts_list)
    
    # Create contrast vectors for chronic slope
    chronic_contrasts_list <- list()
    for (i in 1:length(list_int)) {
      chronic_0 <- as.integer(
        nfe %in% c(
          .time_var, .spline_var, list_int[[i]]$sgxt, list_int[[i]]$sgxs
        )
      )
      chronic_1 <- as.integer(
        nfe %in% c(
          .time_var, .spline_var, txi, sxi, list_int[[i]]$sgxt, 
          list_int[[i]]$sgxs, list_int[[i]]$sgxtxi, list_int[[i]]$sgxsxi
        )
      )
      chronic_d <- as.integer(
        nfe %in% c(txi, sxi, list_int[[i]]$sgxtxi, list_int[[i]]$sgxsxi)
      )
      
      chronic <- list(chronic_0, chronic_1, chronic_d)
      names(chronic) <- c("chronic_0", "chronic_1", "chronic_d")
      
      chronic_contrasts_list[[i]] <- chronic
      names(chronic_contrasts_list)[i] <- paste0("sg", as.character(i))
    }
    
    # Contrasts for the reference group
    ref_contrasts_list <- list(
      list(
        chronic_0 = as.integer(nfe %in% c(.time_var, .spline_var)),
        chronic_1 = as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi)),
        chronic_d = as.integer(nfe %in% c(txi, sxi))
      ) 
    )
    names(ref_contrasts_list) <- paste0(
      "sg", length(chronic_contrasts_list) + 1
    )
    
    # Bind all together
    chronic_contrasts_list <- append(chronic_contrasts_list, ref_contrasts_list)
    
    if (!missing(.prop)) {
      # Create contrast vectors for total slope
      total_contrasts_list <- list()
      for (i in 1:length(list_int)) {
        total_0 <- as.integer(
          nfe %in% c(
            .time_var, .spline_var, list_int[[i]]$sgxt, list_int[[i]]$sgxs
          )
        )
        total_0[
          c(list_spl_pos[[.spline_var]], 
            list_spl_pos[[list_int[[i]]$sgxs]])
        ] <- total_0[
          c(list_spl_pos[[.spline_var]], 
            list_spl_pos[[list_int[[i]]$sgxs]])
        ] * .prop
        total_1 <- as.integer(
          nfe %in% c(
            .time_var, .spline_var, txi, sxi, list_int[[i]]$sgxt, 
            list_int[[i]]$sgxs, list_int[[i]]$sgxtxi, list_int[[i]]$sgxsxi
          )
        )
        total_1[
          c(list_spl_pos[[.spline_var]], 
            list_spl_pos[[sxi]], 
            list_spl_pos[[list_int[[i]]$sgxs]], 
            list_spl_pos[[list_int[[i]]$sgxsxi]])
        ] <- total_1[
          c(list_spl_pos[[.spline_var]], 
            list_spl_pos[[sxi]], 
            list_spl_pos[[list_int[[i]]$sgxs]], 
            list_spl_pos[[list_int[[i]]$sgxsxi]])
        ] * .prop
        total_d <- as.integer(
          nfe %in% c(txi, sxi, list_int[[i]]$sgxtxi, list_int[[i]]$sgxsxi)
        )
        total_d[
          c(list_spl_pos[[sxi]], 
            list_spl_pos[[list_int[[i]]$sgxsxi]])
        ] <- total_d[
          c(list_spl_pos[[sxi]], 
            list_spl_pos[[list_int[[i]]$sgxsxi]])
        ] * .prop
        
        total <- list(total_0, total_1, total_d)
        names(total) <- c("total_0", "total_1", "total_d")
        
        total_contrasts_list[[i]] <- total
        names(total_contrasts_list)[i] <- paste0("sg", as.character(i))
      }
      
      # Contrasts for the reference group
      ref_contrasts_list <- list(
        list(
          total_0 = as.integer(nfe %in% c(.time_var, .spline_var)),
          total_1 = as.integer(nfe %in% c(.time_var, .spline_var, txi, sxi)),
          total_d = as.integer(nfe %in% c(txi, sxi))
        ) 
      )
      # Position of coefficients in the contrast vector that require 
      # multiplication
      sp <- which(names(fixed.effects(.model_obj)) == .spline_var)
      sxip <- which(names(fixed.effects(.model_obj)) == sxi)
      
      # Multiply by `.prop`
      ref_contrasts_list[[1]]$total_0[sp] <- 
        ref_contrasts_list[[1]]$total_0[sp] * .prop
      ref_contrasts_list[[1]]$total_1[c(sp, sxip)] <-
        ref_contrasts_list[[1]]$total_1[c(sp, sxip)] * .prop
      ref_contrasts_list[[1]]$total_d[sxip] <- 
        ref_contrasts_list[[1]]$total_d[sxip] * .prop
      # Rename
      names(ref_contrasts_list) <- 
        paste0("sg", length(total_contrasts_list) + 1)
      
      # Bind all together
      total_contrasts_list <- append(total_contrasts_list, ref_contrasts_list)
    }
    
    # Groups
    reps <- length(levels_by) * 3
    
    grp <- rep(c("Control", "Active", "Active - Control"), times = reps)
    subgrp <- rep(.by, times = reps)
    subgrp_cat <- 
      rep(c(levels_by[-1], levels_by[1]), each = reps / length(levels_by))
    
    # Produce results for acute slope only
    if (.output %in% c("Acute", "ACUTE", "acute")) {
      # List of contrast vectors
      vec <- unlist(acute_contrasts_list, recursive = FALSE)
      # Slope type
      slope <- rep("Acute", times = reps)
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
      message(bred("Acute slope"))
      # Acute slope calculations
      for (i in 1:length(vec)) {
        linear_combination <- match_vars(nfe, vec[[i]])
        display_lc(
          .title = paste0(grp[i], ", ", .by, " = ", subgrp_cat[i]), 
          .fe = linear_combination, .cv = vec[[i]]
        )
      }
      get_heterogeneity(solo)
      return(solo)
    } else if (.output %in% c("Chronic", "CHRONIC", "chronic")) {
      # List of contrast vectors
      vec <- unlist(chronic_contrasts_list, recursive = FALSE)
      # Slope type
      slope <- rep("Chronic", times = reps)
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
      # Chronic slope calculations
      for (i in 1:length(vec)) {
        linear_combination <- match_vars(nfe, vec[[i]])
        display_lc(
          .title = paste0(grp[i], ", ", .by, " = ", subgrp_cat[i]), 
          .fe = linear_combination, .cv = vec[[i]]
        )
      }
      get_heterogeneity(solo)
      return(solo)
    } else if (.output %in% c("Total", "TOTAL", "total")) {
      # List of contrast vectors
      vec <- unlist(total_contrasts_list, recursive = FALSE)
      # Slope type
      slope <- rep("Total", times = reps)
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
      # Total slope calculations
      # Make binary vector to remove all the `.prop` changes so as to index the
      # coefficients and get the linear combinations
      vec_bi <- lapply(vec, function(x) ifelse(x != 0, 1, 0))
      for (i in 1:length(vec)) {
        linear_combination <- match_vars(nfe, vec_bi[[i]])
        # Where to insert multiplication terms into the linear combination
        # (corresponds to all selected coefficients which include `.spline_var`)
        mult <- grep(.spline_var, match_vars(nfe, vec_bi[[i]]))
        display_lc(
          .title = paste0(grp[i], ", ", .by, " = ", subgrp_cat[i]), 
          .fe = linear_combination, .cv = vec[[i]], 
          .prop = deparse(substitute(.prop)), .multiply = mult
        )
      }
      get_heterogeneity(solo)
      return(solo)
    } else if (.output %in% c("All", "ALL", "all")) {
    
      message(
        "The functionality to compute all slopes at once when `.by` is specified has not yet been implemented into this version of the function."
      )
      
    }
  }
}