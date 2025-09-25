
# 8/30/25

# create all data needed for analysis

# create list that contains info needed to get effects (fitted Cox models and GLMs)

# if Cox models do not converge on BS dfs, re-sample until models converge

# naming conventions
  # original = OG
  # bootstrap = BS
  # families: m1, m2, m3, m4, m5, m6
  # sub-markers in family m1: m1.1, m1.2, m1.3, ...

# read in necessary files
source("config_d1.R")
source(data_cleaning_file)
source(fit_models_file)

# get unique combinations of stim and trt groups
stim_trt_combos <- expand.grid(
  stim = stims_input,
  trt_group = trt_groups_input
)

# create empty list to store original dfs
dfs_OG_list <- vector(mode = "list", length = nrow(stim_trt_combos))

# iterate over all stim and trt group combos 
for (i in 1:nrow(stim_trt_combos)) {
  
  # get markers
  markers_input <- get_markers(
    stim = stim_trt_combos$stim[i]
  )
  
  # get dataset that's subset to appropriate trt group and pull variables of interest
  df <- get_data(
    COVAIL_data_name = COVAIL_data_name_input,
    markers = markers_input,
    trt_group = stim_trt_combos$trt_group[i],
    covariates = covariates_input, 
    exposure = exposure_input,
    outcome_time = outcome_time_input,
    outcome_index = outcome_index_input,
    weights = weights_input
  )
  
  # get model fits
  model_fits_OG <- fit_models(
    markers = markers_input, 
    exposure = exposure_input,
    covariates = covariates_input, 
    outcome_time = outcome_time_input, 
    outcome_index = outcome_index_input,
    weights = weights_input,
    a = a_input,
    aprime = aprime_input,
    data = df
  )
  
  # if Cox models converged, continue on
  if (!identical(model_fits_OG, "cox_failed")) {
    
    # save in list
    OG_list_current <- list(
      df = df,
      stim = as.character(stim_trt_combos$stim[i]),
      trt_group = as.character(stim_trt_combos$trt_group[i]),
      job_type = "OG",
      fits = model_fits_OG
    )
    
    # write list to RDS file
    saveRDS(
      OG_list_current,
      here::here(paste0(
        "job_output/data_and_fits/data_fits_", i, ".rds"
      )))
    
    # save in list
    dfs_OG_list[[i]] <- OG_list_current
    
  } else {
    print("Cox models failed on original datasets")
    # Would be nice if this prevented bootstraps from running, but I think code will break here regardless
    # So far, all models have converged on OG datasets
  }
  
}

# if running bootstraps:
if (run_BS == TRUE) {
  
  # iterate through dfs w/ unique stim and trt group combos
  for (i in 1:length(dfs_OG_list)) {
    
    # get original data
    OG_df <- dfs_OG_list[[i]]$df
    
    for (j in 1:bs_B_input) {
      
      # initial seed info
      bs_seed <- j
      retries <- 0
      converged <- FALSE
      
      # repeat until Cox model converges
      while (!converged) {
        
        set.seed(bs_seed)
        
        # draw random bootstrap sample
        boot_rows <- sample(1:nrow(OG_df), replace = TRUE)
        boot_df <- OG_df[boot_rows, ]
        
        # fit models
        model_fits_BS <- fit_models(
          markers = markers_input, 
          exposure = exposure_input,
          covariates = covariates_input, 
          outcome_time = outcome_time_input, 
          outcome_index = outcome_index_input,
          weights = weights_input,
          a = a_input,
          aprime = aprime_input,
          data = boot_df
        )
        
        # If a or a' Cox models fail, shift seed and try again
        if (identical(model_fits_BS, "cox_failed")) {
          
          message("stim/trt: ", i, ". Cox failed at seed ", bs_seed, 
                  ", retrying with new seed...")
          
          # keep track of how many retries are needed
          retries <- retries + 1
          
          # shift seed by number of bs repetitions so no identical seeds are used w/in stim/trt group
          bs_seed <- bs_seed + bs_B_input 
          
        } else {
          
          converged <- TRUE
          
        }
      }
      
      # save in list after convergence
      BS_list_current <- list(
        df = boot_df,
        stim = as.character(dfs_OG_list[[i]]$stim),
        trt_group = as.character(dfs_OG_list[[i]]$trt_group),
        job_type = "BS",
        fits = model_fits_BS,
        seed = bs_seed,
        retries = retries
      )
      
      # write list to RDS file
      saveRDS(
        BS_list_current,
        here::here(paste0(
          "job_output/data_and_fits/data_fits_", nrow(stim_trt_combos) + (i - 1) * bs_B_input + j, ".rds"
        )))
      
    } # end iterating through 1 to # of bootstraps
  } # end iterating over stim and trt combos
  
} # end run bootstrap

