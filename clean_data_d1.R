
# 8/30/25

# read in COVAIL T cell data
# clean and keep variables needed for natural effects analysis

get_data <- function(
    COVAIL_data_name,
    markers,
    trt_group,
    covariates, 
    exposure,
    outcome_time,
    outcome_index,
    weights
) {
  
  # get raw data
  dat_raw <- read.csv(COVAIL_data_name)

  dat_all <- dat_raw  %>%
    # eligibility criteria
    filter(
      # eligible arms
      !(arm %in% c(3, 16, 17)),
      # keep only ptids that meet inclusion criteria
      ph2.D15.tcell == TRUE
    ) %>%
    # create groupings of trt arms
    mutate(
      all_mRNA = ifelse(arm %in% c(1:2, 4:12), 1, 0),
      mRNA_moderna = ifelse(arm %in% c(1:2, 5:6), 1, 0),
      mRNA_pfizer = ifelse(arm %in% c(7:9, 12), 1, 0),
      recomb_protein = ifelse(arm %in% c(13:15), 1, 0)
    ) %>%
    # keep only relevant variables
    dplyr::select(Ptid, exposure, outcome_time, outcome_index, weights,
           all_mRNA, mRNA_moderna, mRNA_pfizer, recomb_protein,
           all_of(c(covariates, unname(markers))))
  
  # subset by treatment groups
  if (trt_group == "all_mRNA") {dat_final <- dat_all[dat_all$all_mRNA == 1, ]}
  if (trt_group == "mRNA_moderna") {dat_final <- dat_all[dat_all$mRNA_moderna == 1, ]}
  if (trt_group == "mRNA_pfizer") {dat_final <- dat_all[dat_all$mRNA_pfizer == 1, ]}
  if (trt_group == "recomb_protein") {dat_final <- dat_all[dat_all$recomb_protein == 1, ]}
  if (trt_group == "all_arms") {dat_final <- dat_all}
  
  # change names from original variable names to m names
  colnames(dat_final) <- c(
    "Ptid", exposure, outcome_time, outcome_index, weights,
    "all_mRNA", "mRNA_moderna", "mRNA_pfizer", "recomb_protein",
    all_of(c(covariates, names(markers))) 
  )
  
  return(dat_final)
  
}
