
# 9/5/25

# compile natural effects and bootstrapped RDS files

library(tidyverse)
library(purrr)
library(survminer)

# get path to results
path <- "job_output/results/"

source("config_dX.R")

# List all files from OG datasets
OG_files <- list.files(path, pattern = "OG.*\\.rds$", full.names = TRUE)

# create empty list to store raw results
OG_list <- vector(mode = "list", length = length(OG_files))

for (i in 1:length(OG_files)) {
  
  # get current file name
  file_name <- basename(OG_files[i])
  
  # pull stim and trt group from name
  matches <- str_match(file_name, "OG_(.*?)_(.*?)\\.rds$")
  stim <- matches[2]
  trt_group <- matches[3]
  
  # read in file
  file <- readRDS(OG_files[i])
  
  # store results in list
  OG_list[[i]] <- file
  OG_list[[i]]$stim <- stim
  OG_list[[i]]$trt_group <- trt_group
  
}

# format raw results into wide data frame
# create empty list to store rows
OG_rows <- list()

# iterate over nat_list
for (jobs_index in seq_along(OG_list)) {
  
  # get values
  entry <- OG_list[[jobs_index]]
  stim <- entry$stim
  trt_group <- entry$trt_group
  
  # iterate over outcome times
  for (timepoint in outcome_time_of_interest_input) {
    
    # get name
    results_name <- paste0("results_", timepoint)
    
    # add values
    df <- entry[[results_name]]
    
    # Add identifying columns
    df$stim <- stim
    df$trt_group <- trt_group
    df$outcome_time <- timepoint
    
    # Reorder columns
    df <- df[, c("stim", "trt_group", "outcome_time", 
                 setdiff(names(df), c("stim", "trt_group", "outcome_time")))]
    
    # Append to list
    OG_rows[[length(OG_rows) + 1]] <- df
    
  }
}

# combine into data.frame
OG_data <- bind_rows(OG_rows)

#############################################################################

# List all bs files
BS_files <- list.files(path, pattern = "BS.*\\.rds$", full.names = TRUE)

# create empty list to store raw BS results
BS_list <- list()

for (i in 1:length(BS_files)) {
  
  # get current file name
  file_name <- basename(BS_files[i])
  
  # pull stim and trt group from name
  matches <- str_match(file_name, "job(.*?)_BS_(.*?)_(.*?)\\.rds$")
  job_id <- matches[2]
  stim <- matches[3]
  trt_group <- matches[4]
  
  # read in file
  file <- readRDS(BS_files[i])
  
  # store results in list
  BS_list[[i]] <- file
  BS_list[[i]]$job_id <- job_id
  BS_list[[i]]$stim <- stim
  BS_list[[i]]$trt_group <- trt_group
}

# format raw results into wide data frame
# create empty list to store rows
BS_rows <- list()

# iterate over IIEs_list
for (bs_index in seq_along(BS_list)) {
  
  # get values
  entry <- BS_list[[bs_index]]
  job_id <- entry$job_id
  stim <- entry$stim
  trt_group <- entry$trt_group
  
  # iterate over outcome times
  for (timepoint in outcome_time_of_interest_input) {
    
    # get name
    results_name <- paste0("results_", timepoint)
    
    # add values
    df <- entry[[results_name]]
    
    # Add identifying columns
    df$job_id <- job_id
    df$stim <- stim
    df$trt_group <- trt_group
    df$outcome_time <- timepoint
    
    # Reorder columns
    df <- df[, c("job_id", "stim", "trt_group", "outcome_time",
                 setdiff(names(df), c("job_id", "stim", "trt_group", "outcome_time")))]
    
    # Append to list
    BS_rows[[length(BS_rows) + 1]] <- df
    
  }
}

# combine into data.frame
BS_data <- bind_rows(BS_rows)

################################################################################

# check which iterations did not converge
ids <- as.numeric(unique(BS_data$job_id))
want <- 9:5005
setdiff(want, ids)
length(setdiff(want, ids))
length(setdiff(want, ids)) / length(ids)

nrow(BS_data %>% filter(stim == "D614G", trt_group == "all_arms", outcome_time == "0")) # 1000
nrow(BS_data %>% filter(stim == "BA.4.5", trt_group == "all_arms", outcome_time == "0")) # 1000

################################################################################

# format data from original datasets
effects_data_OG <- OG_data %>%
  mutate(
    total = mean_Y1 - mean_Y2,
    D15_CD4 = mean_Y1 - mean_Y3,
    D15_CD8 = mean_Y3 - mean_Y4,
    D15_Ab = mean_Y4 - mean_Y5,
    B_CD4 = mean_Y5 - mean_Y6,
    B_CD8 = mean_Y6 - mean_Y7,
    B_Ab = mean_Y7 - mean_Y8,
    direct = mean_Y8 - mean_Y2
  ) %>%
  select(stim, trt_group, outcome_time, total, 
         mean_Y1:mean_Y8,
         D15_CD4, D15_CD8, D15_Ab, B_CD4, B_CD8, B_Ab, direct)

# make longer
effects_long_OG <- effects_data_OG %>%
  pivot_longer(
    cols = total:direct,
    names_to = "effect",
    values_to = "value"
  ) %>%
  group_by(stim, trt_group, outcome_time) %>%
  mutate(
    # calculate proportion mediated
    total_val = ifelse(effect %in% c("D15_Ab", "D15_CD8", "D15_CD4", "D15",
                                     "B_Ab", "B_CD8", "B_CD4", "direct"),
                       value[effect == "total"], NA),
    PM = value / total_val,
    log_PM = log10(PM)
  )

# repeat for bootstrapped data
effects_data_BS <- BS_data %>%
  mutate(
    total = mean_Y1 - mean_Y2,
    D15_CD4 = mean_Y1 - mean_Y3,
    D15_CD8 = mean_Y3 - mean_Y4,
    D15_Ab = mean_Y4 - mean_Y5,
    B_CD4 = mean_Y5 - mean_Y6,
    B_CD8 = mean_Y6 - mean_Y7,
    B_Ab = mean_Y7 - mean_Y8,
    direct = mean_Y8 - mean_Y2
  ) %>%
  select(job_id, stim, trt_group, outcome_time, total, 
         mean_Y1:mean_Y8,
         D15_CD4, D15_CD8, D15_Ab, B_CD4, B_CD8, B_Ab, direct)

# make longer
effects_long_BS <- effects_data_BS %>%
  pivot_longer(
    cols = total:direct,
    names_to = "effect",
    values_to = "bs_value"
  ) %>%
  group_by(job_id, stim, trt_group, outcome_time) %>%
  mutate(
    # calculate proportion mediated
    bs_total_val = ifelse(effect %in% c("D15_Ab", "D15_CD8", "D15_CD4", "D15",
                                     "B_Ab", "B_CD8", "B_CD4", "direct"),
                       bs_value[effect == "total"], NA),
    bs_PM = bs_value / bs_total_val
  )

effects_long_both <- effects_long_BS %>%
  left_join(effects_long_OG) %>%
  # drop total values used in PM
  select(!c(bs_total_val, total_val))

################################################################################

# get PM confidence intervals
# 1. raw PM - Wald CI using Delta method
# 2. log PM - Wald CI using Delta method

# duplicate OG data and rename variables to specify original estimates
OG_for_PM <- OG_data
colnames(OG_for_PM) <- c("stim", "trt_group", "outcome_time", 
                         paste0("OG_Y", 1:8))

# duplicate BS data and rename variables to specify bootstrapped estimates
BS_for_PM <- BS_data
colnames(BS_for_PM) <- c("job_id", "stim", "trt_group", "outcome_time", 
                         paste0("BS_Y", 1:8))

# join together wide versions (for each BS job_id, there's one OG estimate per Y)
OG_BS_wide <- BS_for_PM %>%
  left_join(OG_for_PM)

# unique combos of stim, trt group, and time
stim_trt_time_combos <- expand.grid(
  time = outcome_time_of_interest_input,
  stim = stims_input,
  trt_group = trt_groups_input
)

### Define the effect pairs ##
effect_pairs <- list(
  c("Y1","Y3"),
  c("Y3","Y4"),
  c("Y4","Y5"),
  c("Y5","Y6"),
  c("Y6","Y7"),
  c("Y7","Y8"),
  c("Y8","Y2")
)

# create empty list to store results
PM_SE_results_list <- vector(mode = "list", length = nrow(stim_trt_time_combos))

# index <- 26 # Day 100

# iterate through combinations of stim and trt group
for (index in 1:nrow(stim_trt_time_combos)) {
  
  # empty list to store 7 PM SEs
  one_row <- data.frame(
    stim = stim_trt_time_combos$stim[index],
    trt_group = stim_trt_time_combos$trt_group[index],
    outcome_time = stim_trt_time_combos$time[index]
  )  %>%
    mutate(Y1_Y3_SE_raw = NA, Y3_Y4_SE_raw = NA, Y4_Y5_SE_raw = NA,
           Y5_Y6_SE_raw = NA, Y6_Y7_SE_raw = NA, Y7_Y8_SE_raw = NA, Y8_Y2_SE_raw = NA,
           Y1_Y3_SE_log = NA, Y3_Y4_SE_log = NA, Y4_Y5_SE_log = NA,
           Y5_Y6_SE_log = NA, Y6_Y7_SE_log = NA, Y7_Y8_SE_log = NA, Y8_Y2_SE_log = NA)
  
  # get current data (n = 1000)
  dat_current <- OG_BS_wide %>%
    filter(
      stim == stim_trt_time_combos$stim[index],
      trt_group == stim_trt_time_combos$trt_group[index],
      outcome_time == stim_trt_time_combos$time[index]
    )
  
  # raw eqn for PM: PM = (Yi - Yj) / (Y1 - Y2)
  
  # log eqn for PM: log(PM) = log(Yi - Yj) - log(Y1 - Y2)
  
  # total effect
  total <- dat_current$OG_Y1[1] - dat_current$OG_Y2[1]
  
  # values calculated using just Y1 and Y2:
  # raw - d/dYi
  gi_raw = 1 / total
  
  # raw - d/dYj
  gj_raw = -1 / total
  
  # log - d/dY1
  g1_log = -1 / total
  
  # log - d/dY2
  g2_log = 1 / total
  
  # iterate through effect pairs
  for (k in 1:length(effect_pairs)) {
    
    # names for Yi and Yj
    Yi_name <- effect_pairs[[k]][1]
    Yj_name <- effect_pairs[[k]][2]
    
    # values for Yi and Yj
    Yi <- dat_current[[paste0("OG_", Yi_name)]][1]
    Yj <- dat_current[[paste0("OG_", Yj_name)]][1]
    
    # raw - d/dY1
    g1_raw = (Yj - Yi) / (total^2)
    
    # raw - d/dY2
    g2_raw = (Yi - Yj) / (total^2)
    
    # log - d/dYi
    gi_log = 1 / (Yi - Yj)
    
    # log - d/dYj
    gj_log = -1 / (Yi - Yj)
    
    # gradient vectors
    grad_raw <- c(g1_raw, g2_raw, gi_raw, gj_raw)
    
    grad_log <- c(g1_log, g2_log, gi_log, gj_log)
    
    # covariance matrix (same for raw and log PM's)
    Sigma <- cov(data.frame(
      Y1 = dat_current$BS_Y1,
      Y2 = dat_current$BS_Y2,
      Yi = dat_current[[paste0("BS_", Yi_name)]],
      Yj = dat_current[[paste0("BS_", Yj_name)]]
    ), use = "na.or.complete")   
    
    # SE's
    SE_raw <- as.numeric(sqrt(t(grad_raw) %*% Sigma %*% grad_raw))
    
    SE_log <- as.numeric(sqrt(t(grad_log) %*% Sigma %*% grad_log))
    
    # store results for effects
    one_row[[paste0(Yi_name, "_", Yj_name, "_SE_raw")]] <- SE_raw
    
    one_row[[paste0(Yi_name, "_", Yj_name, "_SE_log")]] <- SE_log
    
  }
  
  # store all effects for that row
  PM_SE_results_list[[index]] <- one_row
  
}

# bind rows
PM_SE_results_df <- bind_rows(PM_SE_results_list)

# make longer, rename effects
PM_SE_results_df_long <- PM_SE_results_df %>%
  pivot_longer(
    cols = matches("SE_raw$|SE_log$"),
    names_to = c("effect1", ".value"),
    names_pattern = "(Y\\d+_Y\\d+)_SE_(raw|log)"
  ) %>%
  rename(PM_SE = raw, log_PM_SE = log) %>%
  mutate(effect = case_when(
    effect1 == "Y1_Y3" ~ "D15_CD4",
    effect1 == "Y3_Y4" ~ "D15_CD8",
    effect1 == "Y4_Y5" ~ "D15_Ab",
    effect1 == "Y5_Y6" ~ "B_CD4",
    effect1 == "Y6_Y7" ~ "B_CD8",
    effect1 == "Y7_Y8" ~ "B_Ab",
    effect1 == "Y8_Y2" ~ "direct",
  )) %>%
  select(!effect1)

# join with BS and OG data
effects_both_PM <- effects_long_both %>%
  left_join(PM_SE_results_df_long)  %>%
  mutate(
    type = case_when(
      effect %in% c("mean_Y1", "mean_Y2", "mean_Y3", "mean_Y4",
                    "mean_Y5", 'mean_Y6', "mean_Y7", "mean_Y8") ~ "Y",
      effect == "total" ~ "total",
      effect %in% c("D15_Ab", "D15_CD8", "D15_CD4") ~ "D15",
      effect %in% c("B_Ab", "B_CD8", "B_CD4") ~ "B",
      effect == "direct" ~ "direct"
    ),
    effect_new = case_when(
      effect %in% c("D15_Ab", "B_Ab") ~ "Ab",
      effect %in% c("D15_CD8", "B_CD8") ~ "CD8",
      effect %in% c("D15_CD4", "B_CD4") ~ "CD4",
      effect == "direct" ~ "direct",
      effect %in% c("mean_Y1", "mean_Y2", "total") ~ "total",
      effect %in% c("mean_Y3", "mean_Y4", "mean_Y5", 'mean_Y6', 
                    "mean_Y7", "mean_Y8") ~ "parts"
    )
  )

################################################################################

# summarize data
dat_final <- effects_both_PM %>%
  group_by(stim, trt_group, outcome_time, effect, effect_new, type, value,
           PM, PM_SE, log_PM, log_PM_SE) %>%
  summarize(
    # percentile
    bs_percCI_LB = quantile(bs_value, probs = 0.025, na.rm = TRUE),
    bs_percCI_UB = quantile(bs_value, probs = 0.975, na.rm = TRUE),
    # Wald
    bs_waldCI_LB = value[1] - 1.96 * sd(bs_value, na.rm = TRUE),
    bs_waldCI_UB = value[1] + 1.96 * sd(bs_value, na.rm = TRUE),
    # PM - percentile CI's
    PM_percCI_LB = quantile(bs_PM, probs = 0.025, na.rm = TRUE),
    PM_percCI_UB = quantile(bs_PM, probs = 0.975, na.rm = TRUE),
    .groups = "keep"
  ) %>%
  # PM Wald CI's using Delta method
  mutate(
    # raw PM CI's
    PM_waldCI_LB = PM - (1.96 * PM_SE),
    PM_waldCI_UB = PM + (1.96 * PM_SE),
    # log PM CI's
    PM_waldCI_LB_log = log_PM - (1.96 * log_PM_SE),
    PM_waldCI_UB_log = log_PM + (1.96 * log_PM_SE)
  )

# save as factor
dat_final$effect <- factor(
  dat_final$effect,
  levels = c("total", "mean_Y1", "mean_Y2", "mean_Y3", "mean_Y4", "mean_Y5", 'mean_Y6', 
             "mean_Y7", "mean_Y8", "D15_CD4", "D15_CD8", "D15_Ab",
             "B_CD4", "B_CD8", "B_Ab", "direct")
)

dat_final$effect_new <- factor(
  dat_final$effect_new,
  levels = c("total", "CD4", "CD8", "Ab", "direct", "parts")
)

dat_final$type <- factor(
  dat_final$type,
  levels = c("total", "Y", "B", "D15", "direct")
)

saveRDS(dat_final,
        paste0(path, "/results_summary.rds"))

saveRDS(effects_both_PM, 
        paste0(path, "/results_BS.rds"))








