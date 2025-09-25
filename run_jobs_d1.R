
# 8/30/25

# read in files
source("config_d1.R")
source(effects_file)

for (job_index in 1:total_jobs) {
  
  # Load data for one jobs
  data_and_fits_list <- readRDS(here::here(paste0(
    "job_output/data_and_fits/data_fits_", job_index, ".rds"
  )))
  
  cat("Running job number", job_index, "\n")
  cat("job_type: ", data_and_fits_list$job_type, ";",
      "stim: ", data_and_fits_list$stim, ";",
      "trt_group: ", data_and_fits_list$trt_group)

  # get natural effects
  effect_parts_results <- get_effects_parts(
    seed = seed_input,
    model_fit_list = data_and_fits_list$fits,
    markers = get_markers(
      stim = data_and_fits_list$stim
    ),
    data = data_and_fits_list$df,
    weights = weights_input,
    covariates = covariates_input,
    mc_B = mc_B_input,
    outcome_time_of_interest = outcome_time_of_interest_input
  )
  
  # save results
  saveRDS(
    effect_parts_results,
    here::here(paste0(
      "job_output/results/job",
      job_index, "_",
      data_and_fits_list$job_type, "_",
      data_and_fits_list$stim, "_", 
      data_and_fits_list$trt_group, ".rds"))
  )
  
} # end job loop


