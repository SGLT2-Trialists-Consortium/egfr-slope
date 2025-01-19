#' @title
#' Generate Synthetic Datasets for eGFR Slope Analysis
#'
#' @description
#' Creates two synthetic datasets, one representing baseline information and the
#' other repeated eGFR measurements over time. These data are generated for 
#' demonstration purposes (e.g., modeling eGFR slope in hypothetical SGLT2 
#' inhibitor trials).
#'
#' @details
#' The baseline dataset mimics a simplified version of an `ADSL`-like table, 
#' while the follow-up dataset contains repeated eGFR measurements, similar to 
#' an `ADLBM` structure. Each dataset includes synthetic variables for treatment 
#' assignment, baseline eGFR, and random follow-up times. Some artificial
#' missingness is introduced to reflect real-world trial data challenges.
#'
#' @param .table A character string indicating which dataset to return. 
#'   Possible values:
#'   \itemize{
#'     \item \code{"baseline"} – returns the baseline dataset (mimicking 
#'       `ADSL`).
#'     \item \code{"follow-up"} – returns the follow-up dataset (mimicking 
#'       repeated eGFR data in `ADLBM`).
#'   }
#'
#' @return A \code{tibble} containing synthetic data.  
#'   \describe{
#'     \item{\code{"baseline"}}{Baseline demographics and eGFR-related variables 
#'       for 5000 simulated subjects.}
#'     \item{\code{"follow-up"}}{Repeated eGFR measurements and associated time 
#'       points (e.g., \code{avisitn}, \code{ady}).}
#'   }
#'
#' @section Notes:
#' \itemize{
#'   \item This function installs \pkg{simstudy} and \pkg{tidyverse} if not 
#'         already present, then uses \code{simstudy} to generate random values 
#'         that approximate marginal distributions observed in typical SGLT2 
#'         inhibitor trials.
#'   \item Subject IDs are prefixed with \code{"id"} for clarity. 
#'   \item Minor amounts of missing data and duplicated rows are introduced to
#'         mimic real clinical trial data structures.
#' }
#'
#' @examples
#' \dontrun{
#' # Generate baseline dataset
#' bl_data <- generate_synthetic_data(.table = "baseline")
#'
#' # Generate follow-up dataset
#' fu_data <- generate_synthetic_data(.table = "follow-up")
#'
#' # Check structure
#' head(bl_data)
#' head(fu_data)
#' }
#'
#' @export

generate_synthetic_data <- function(.table) {
  
  # Set seed
  set.seed(69)
  
  # Install dependencies if not already installed
  libs <- c("simstudy", "tidyverse")
  install.packages(setdiff(libs, rownames(installed.packages())))
  
  # Function to generate unique identifiers
  assign_id <- function(.x) {
    .x <- dplyr::case_when(
      stringr::str_detect(.x, "^[0-9]$")    ~ paste0("id000", .x),
      stringr::str_detect(.x, "^[0-9]{2}$") ~ paste0("id00", .x),
      stringr::str_detect(.x, "^[0-9]{3}$") ~ paste0("id0", .x),
      stringr::str_detect(.x, "^[0-9]{4}$") ~ paste0("id", .x),
      TRUE ~ .x
    )
    return(.x)
  }
  
  # Function simulate eGFR
  simulate_egfr <- function(.def = NULL, .varname, .mean, .var) {
    if (!missing(.def)) {
      .def <- simstudy::defData(
        .def, varname = .varname, dist = "normal", formula = .mean, 
        variance = .var
      )
    } else {
      .def <- simstudy::defData(
        varname = .varname, dist = "normal", formula = .mean, variance = .var
      )
    }
    return(.def)
  }
  
  # Create unique subject identifiers
  id <- seq(from = 1, to = 5000, by = 1) |> as.character() |> assign_id()
  
  # Define simulated baseline variables
  def_bl <-
    # Study strata
    simstudy::defData(
      varname = "strata", formula = "0.3;0.3;0.4", dist = "categorical"
    ) |> 
    # Treatment
    simstudy::defData(
      varname = "trt01pn", formula = "1;1", dist = "trtAssign", 
      variance = "strata"
    ) |> 
    # GLP-1RA subgroups
    simstudy::defData(
      varname = "blglp1", formula = "0.04;0.96", dist = "categorical"
    ) |> 
    # On-treatment flag
    simstudy::defData(
      varname = "trtfl", formula = "0.001;0.999", dist = "categorical"
    ) |> 
    # Baseline eGFR
    simulate_egfr(.varname = "blgfr1", .mean = 37.7, .var = 73.2) |> 
    simulate_egfr(.varname = "blgfr2", .mean = 52.8, .var = 99.8) |> 
    simulate_egfr(.varname = "blgfr3", .mean = 72.0, .var = 179.2)
  
  # Generate baseline data columns
  col_bl <- 
    # `data.table` of 5000 rows
    simstudy::genData(5000, def_bl) |> 
    # Convert to `tibble`
    tibble::tibble() |> 
    dplyr::mutate(
      usubjid = assign_id(as.character(id)), 
      .keep = "unused", .before = 1
    )
  
  # Combine and recode to create baseline data
  bl <- 
    tibble::tibble(usubjid = id) |> 
    dplyr::left_join(col_bl, by = dplyr::join_by(usubjid)) |> 
    dplyr::mutate(
      ittfl = "Y",
      randfl = "Y",
      blglp1 = dplyr::if_else(blglp1 == 1, "Y", "N"),
      trtfl = dplyr::if_else(trtfl == 2, "Y", "N"),
      blgfr = dplyr::case_when(
        strata == 1 ~ blgfr1,
        strata == 2 ~ blgfr2,
        strata == 3 ~ blgfr3
      ), 
      blgfr = as.integer(blgfr),
      strata = dplyr::case_when(
        strata == 1 ~ "Screening eGFR 30 to <45 mL/min/1.73m2",
        strata == 2 ~ "Screening eGFR 45 to <60 mL/min/1.73m2",
        strata == 3 ~ "Screening eGFR 60 to <90 mL/min/1.73m2"
      )
    ) |> 
    dplyr::select(usubjid, trt01pn, randfl, ittfl, trtfl, blgfr, blglp1, strata)
  
  # Create simulated follow-up eGFR variables
  def_fu <-
    # Week 3
    simulate_egfr(.varname = "wk3_0", .mean = 53.4, .var = 43.4) |> 
    simulate_egfr(.varname = "wk3_1", .mean = 51.7, .var = 42.0) |> 
    # Week 13
    simulate_egfr(.varname = "wk13_0", .mean = 54.6, .var = 47.4) |> 
    simulate_egfr(.varname = "wk13_1", .mean = 52.5, .var = 43.2) |> 
    # Week 26
    simulate_egfr(.varname = "wk26_0", .mean = 52.9, .var = 46.7) |> 
    simulate_egfr(.varname = "wk26_1", .mean = 52.0, .var = 44.6) |> 
    # Week 52
    simulate_egfr(.varname = "wk52_0", .mean = 51.0, .var = 52.2) |> 
    simulate_egfr(.varname = "wk52_1", .mean = 51.0, .var = 46.4) |> 
    # Week 78
    simulate_egfr(.varname = "wk78_0", .mean = 48.9, .var = 53.2) |> 
    simulate_egfr(.varname = "wk78_1", .mean = 50.2, .var = 50.1) |> 
    # Week 104
    simulate_egfr(.varname = "wk104_0", .mean = 47.6, .var = 55.3) |> 
    simulate_egfr(.varname = "wk104_1", .mean = 50.0, .var = 51.6) |> 
    # Week 130
    simulate_egfr(.varname = "wk130_0", .mean = 46.4, .var = 52.1) |> 
    simulate_egfr(.varname = "wk130_1", .mean = 49.7, .var = 55.1) |> 
    # Week 156
    simulate_egfr(.varname = "wk156_0", .mean = 49.3, .var = 55.8) |> 
    simulate_egfr(.varname = "wk156_1", .mean = 46.4, .var = 50.6) 
  
  # Generate follow-up data columns
  col_fu <-  
    simstudy::genData(5000, def_fu) |> 
    tibble::tibble() |> 
    dplyr::bind_cols(bl |> dplyr::select(usubjid, trt01pn)) |> 
    dplyr::select(usubjid, trt01pn, dplyr::matches("wk")) |> 
    dplyr::mutate(
      `3` = dplyr::if_else(trt01pn == 1, wk3_1, wk3_0),
      `13` = dplyr::if_else(trt01pn == 1, wk13_1, wk13_0),
      `26` = dplyr::if_else(trt01pn == 1, wk26_1, wk26_0),
      `52` = dplyr::if_else(trt01pn == 1, wk52_1, wk52_0),
      `78` = dplyr::if_else(trt01pn == 1, wk78_1, wk78_0),
      `104` = dplyr::if_else(trt01pn == 1, wk104_1, wk104_0),
      `130` = dplyr::if_else(trt01pn == 1, wk130_1, wk130_0),
      `156` = dplyr::if_else(trt01pn == 1, wk156_1, wk156_0),
      .keep = "unused"
    ) |> 
    tidyr::pivot_longer(
      cols = dplyr::matches("^[0-9]"),
      names_to = "avisitn",
      values_to = "aval"
    ) |>
    dplyr::mutate(
      avisitn = as.double(avisitn),
      avisit = paste("WEEK", avisitn),
      aval = as.integer(aval),
      randfl = "Y",
      ittfl = "Y",
      trtfl = "Y",
      anl01fl = "Y",
      paramcd = "GFRBSCRT",
      param = "GFR from Creatinine Adjusted for BSA (mL/min/1.73m2)",
      ady = as.integer(rnorm(40000, mean = avisitn * 7))
    )
  
  # Records to remove (to simulate missing data)
  rm <- col_fu |> 
    dplyr::slice_sample(prop = 0.1, weight_by = avisitn)
  
  # Duplicate records
  anl01fln <- col_fu |> 
    dplyr::select(-anl01fl) |> 
    dplyr::mutate(aval = aval + 5, ady = ady - 2, anl01fl = NA) |> 
    dplyr::slice_sample(prop = 0.02) 
  
  # Combine with baseline data to produce repeat eGFR data
  fu <- col_fu |> 
    dplyr::anti_join(
      rm,
      by = dplyr::join_by(
        usubjid, avisitn, aval, avisit, randfl, ittfl, trtfl, anl01fl, paramcd, 
        param, ady
      )
    ) |> 
    dplyr::bind_rows(anl01fln) |> 
    dplyr::left_join(
      bl |> dplyr::select(usubjid, base = blgfr), by = dplyr::join_by(usubjid)
    ) |> 
    dplyr::bind_rows(
      bl |>
        dplyr::select(
          usubjid, aval = blgfr, base = blgfr, randfl, ittfl, trtfl
        ) |> 
        dplyr::mutate(anl01fl = "Y", ady = 1, avisitn = 0, avisit = "BASELINE")
    ) |> 
    dplyr::arrange(usubjid, avisitn) |> 
    dplyr::select(
      usubjid, randfl, ittfl, trtfl, anl01fl, aval, base, paramcd, param, ady,
      avisitn, avisit
    ) |> 
    tidyr::fill(c(param, paramcd), .direction = "updown") 
  
    # Recode `trtfl` for individuals not on treatment at baseline
    trtfln <- bl |> dplyr::filter(trtfl == "N") |> dplyr::pull(usubjid)
    
    fu <- fu |> 
      dplyr::mutate(
        trtfl = dplyr::if_else(
          usubjid %in% trtfln, "N", trtfl
        )
      )
  
  if (.table == "baseline") {
    return(bl)
  } else if (.table == "follow-up") {
    return(fu)
  }
} 