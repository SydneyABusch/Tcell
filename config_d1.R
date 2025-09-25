
# 9/24/25

# assign variables

library(tidyverse)
library(purrr)
library(survival)

here::i_am("config_d1.R")

# files for analysis
data_cleaning_file <- "clean_data_d1.R"
fit_models_file <- "fit_models_d1.R"
effects_file <- "effects_fcn_d1.R"

# COVAIL data and file path
COVAIL_data_name_input <- "/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_data_processed_20250524.csv"

# marker groups (must match groups listed in get_markers fcn below)
# Note: calling these "stim" groups to reference which viral proteins are being stimulated by assays
# avoided calling them "marker groups" so they're not confused with marker families and sub-markers
stims_input <- c("D614G", "BA.4.5")

# trt groups (all_mRNA, mRNA_moderna, mRNA_pfizer, recomb_protein, all_arms)
trt_groups_input <- "all_arms"

# random seed for effects fcn
seed_input = 828

# number of repetitions for Monte Carlo integration
mc_B_input = 10000

# whether or not to run bootstraps
run_BS = TRUE

if (run_BS == TRUE) {
  # number of repetitions for bootstrapping
  bs_B_input = 1000
}

# variable names for weight, covariates, exposure, and outcome info
weights_input = "wt.D15.tcell"
# covariates_input = c("Age", "Sex", "White", "EthnicityHispanic")
covariates_input = c("standardized_risk_score", "FOIstandardized")
exposure_input = "naive"
outcome_time_input = "EventTimePrimaryD15"
outcome_index_input = "EventIndPrimaryD15"

# times to measure outcomes
outcome_time_of_interest_input <- c(seq(from = 0, to = 188, by = 4), 91, 181)

# exposure values
a_input = 1       # naive / no previous infection
aprime_input = 0  # non-naive / previous infection

# total number of "jobs"
total_jobs <- (bs_B_input * length(stims_input) * length(trt_groups_input)) + 
  (length(stims_input) * length(trt_groups_input))
    

# create function that gets markers based off of marker groups
get_markers <- function(stim) {
  
  # markers for D614G analysis
  if (stim == "D614G") {
    markers <- c(
      # m1 - baseline CD4+ T cells
      m1.1 = "Bcd4_IFNg.IL2_BA.4.5.S",
      # m2 - baseline CD8+ T cells
      m2.1 = "Bcd8_IFNg.IL2_BA.4.5.S",
      # m3 - baseline antibodies
      m3.1 = "Bpseudoneutid50_D614G",
      # m4 - day 15 CD4+ T cells
      m4.1 = "Day15cd4_IFNg.IL2_BA.4.5.S",
      # m5 - day 15 CD8+ T cells
      m5.1 = "Day15cd8_IFNg.IL2_BA.4.5.S",
      # m6 - day 15 antibodies
      m6.1 = "Day15pseudoneutid50_BA.4.BA.5"
    )
  }
  
  # markers for BA.4.5 analysis
  if (stim == "BA.4.5") {
    markers <- c(
      # m1 - baseline CD4+ T cells
      m1.1 = "Bcd4_IFNg.IL2_BA.4.5.S",
      # m2 - baseline CD8+ T cells
      m2.1 = "Bcd8_IFNg.IL2_BA.4.5.S",
      # m3 - baseline antibodies
      m3.1 = "Bpseudoneutid50_BA.4.BA.5",
      # m4 - day 15 CD4+ T cells
      m4.1 = "Day15cd4_IFNg.IL2_BA.4.5.S",
      # m5 - day 15 CD8+ T cells
      m5.1 = "Day15cd8_IFNg.IL2_BA.4.5.S",
      # m6 - day 15 antibodies
      m6.1 = "Day15pseudoneutid50_BA.4.BA.5"
    )
  }
  
  return(markers)
}


