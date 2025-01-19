#*******************************************************************************
#
# Project:      Modelling eGFR slope
# Last updated: 19-Jan-2025
# Authors:      Robert Fletcher and Niels Jongs
# Purpose:      Model eGFR slope using synthetic data
# Contact:      rfletcher@georgeinstitute.org.au | n.jongs@umcg.nl
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This reproducible example illustrates how to calculate estimated glomerular 
# filtration rate (eGFR) slope. It demonstrates the required data format and the
# steps needed to replicate the analyses in your own clinical trial data

# For additional details or support, please consult the repository’s README.md 
# or use the contact information provided above

# Please note, this example assumes your dataset includes repeated eGFR
# measurements over time, in addition to any relevant covariates for the 
# analyses


# Install dependencies (if not currently installed) -----------------------

# Library names
libs <- c("glue", "multcomp", "lme4", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(glue)
library(nlme)
library(multcomp)
library(tidyverse)


# Define variables --------------------------------------------------------

# Path to cloned repository
# This assumes you're using RStudio. It extracts the current script's directory
# and removes the "/code" suffix. If you're not using RStudio, specify the path 
# manually
path <- sub("/code$", "", dirname(rstudioapi::getSourceEditorContext()$path))

# Set spline knot point
# In the CREDENCE trial, 21 days (3 weeks) corresponds to the first study visit,
# which marks the initial acute drop in eGFR with SGLT2 inhibition. Adjust this
# value according to the available date in your trial. For example, in DAPA-CKD, 
# the first post-randomisation visit is at 14 days (2 weeks), so `k` would be 
# set to 14
k <- 21


# Source functions --------------------------------------------------------

source(glue::glue("{path}/src/generate_synthetic_data.R"))
source(glue::glue("{path}/src/compute_slope.R"))


# Load data ---------------------------------------------------------------

bl <- generate_synthetic_data(.table = "baseline")
fu <- generate_synthetic_data(.table = "follow-up")

# If you encounter issues running the function above, the two datasets generated 
# by this function are also available in CSV format in the `data` sub-directory 
# of this repository. Use the following functions to load those files:

# bl <- readr::read_csv("{path}/data/synthetic_trial_baseline.csv")
# fu <- readr::read_csv("{path}/data/synthetic_trial_egfr_follow_up.csv")


# Prepare data ------------------------------------------------------------

# Treatment Arms
# Here, we filter for on-treatment individuals only. NOTE: This step is 
# optional. Depending on your study design, you might prefer an 
# intention-to-treat population instead
arm <- bl |> 
  dplyr::filter(trtfl == "Y") |> 
  dplyr::select(usubjid, trt01pn, blglp1, strata) |> 
  dplyr::mutate(blglp1 = dplyr::if_else(blglp1 == "Y", 1, 0))

# Repeat eGFR
# Here, we filter for on-treatment measurements and those flagged for analysis. 
# Adjust these criteria according to your chosen analysis strategy
gfr <- fu |> 
  dplyr::filter(dplyr::if_all(c(trtfl, anl01fl), \(x) x == "Y")) |> 
  dplyr::distinct(usubjid, aval, ady, base, avisitn)

# Join
gfr_c <- arm |> 
  dplyr::left_join(gfr, by = dplyr::join_by(usubjid), multiple = "all") |> 
  dplyr::mutate(
    # Convert days from baseline to years
    time = ady / 365.25,
    # Generate spline term
    spline = dplyr::if_else(time >= k / 365.25, time - k / 365.25, 0)
  )

# Formatting the data for the model is likely the most crucial step. The 
# synthetic data used in this reprex should give you an idea of the required 
# format, but consult the repository README if you're unsure about anything or 
# would like a more detailed explanation


# Fit model (whole population) --------------------------------------------

# Fit mixed effects model with unstructured residual variance-covariance matrix
fit <- lme4::lmer(
  aval ~ base + strata + time * trt01pn + spline * trt01pn - 1 + 
  (time | usubjid), data = gfr_c
)

# The above model is for the overall trial population. Subgroup-specific 
# analyses are detailed further down

# Please see the repository README if you'd like more information on 
# specification of the model parameters


# Extract model results (whole population) --------------------------------

# Define the proportion of the total slope attributed to the chronic slope. In 
# this example, 1095.75 represents the total trial follow-up period (roughly 3 
# years in days). By subtracting 21 days (acute slope), we calculate the
# fraction that remains for the chronic slope
prop <- (1095.75 - k) / 1095.75

# Compute eGFR slope for the whole population
all <- compute_slope(
  .model_obj = fit, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .output = "all"
)


# Fit model and compute slope by binary subgroup variable -----------------

# Recode baseline GLP-1RA use as factor and set `blglp1` == "Yes" as reference
gfr_c <- gfr_c  |> 
  dplyr::mutate(
    blglp1 = factor(blglp1, levels = c(1, 0), labels = c("Yes", "No"))
  )

# Fit mixed effects model with unstructured residual variance-covariance matrix
# this time with adjustment and interactions with `blglp1`
fit_binary <- lme4::lmer(
  aval ~ base + strata + time * trt01pn + spline * trt01pn + time * blglp1 +
  spline * blglp1 + trt01pn * blglp1 + time * trt01pn * blglp1 + 
  spline * trt01pn * blglp1 - 1 + (time | usubjid), data = gfr_c
)

# Compute eGFR slope by binary subgroups (baseline GLP-1RA use)
binary_subgroups <- compute_slope(
  .model_obj = fit_binary, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .by = "blglp1", # `.by` argument specified with subgroup variable
  .output = "total"
)


# Fit model and compute slope by multi-level subgroup variable ------------

# Generate an ordinal categorical variable for baseline eGFR subgroups. The 
# `compute_slope()` function can handle any number of levels, but four are 
# demonstrated here as an example
gfr_c <- gfr_c |> 
  mutate(
    gfr_grp = factor(ntile(base, 4), levels = c(1:4), labels = c(1:4)),
    # Set quartile 4 as the reference level
    gfr_grp = forcats::fct_relevel(gfr_grp, "4")
  )

# Fit mixed effects model with unstructured residual variance-covariance matrix
# with adjustment and interactions with `gfr_grp`
fit_ordinal <- lme4::lmer(
  aval ~ base + strata + time * trt01pn + spline * trt01pn + time * gfr_grp +
  spline * gfr_grp + trt01pn * gfr_grp + time * trt01pn * gfr_grp + 
  spline * trt01pn * gfr_grp - 1 + (time | usubjid), data = gfr_c
)

# Compute eGFR slope by multiple subgroups (baseline eGFR quartiles)
ordinal_subgroups <- compute_slope(
  .model_obj = fit_ordinal, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "gfr_grp", # `.by` argument specified with subgroup variable
  .output = "total"
)