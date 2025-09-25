
# 8/21/25

# fcn that estimates natural effects

get_effects_parts <- function(
  seed,
  model_fit_list,
  markers,
  data,
  weights,
  covariates,
  mc_B,
  outcome_time_of_interest
) {
  
  # create list of empty dfs for results
  results_list <- vector(mode = "list", length = length(outcome_time_of_interest))
  
  # each outcome time of interest gets it's own df
  for (k in 1:length(outcome_time_of_interest)) {
    
    results_df <- data.frame(
      i = 1:nrow(data),
      Ptid = data$Ptid,
      weights = data[[weights]],
      # names correspond to formulas in SAP
      Y1 = NA, Y2 = NA, Y3 = NA, Y4 = NA,
      Y5 = NA, Y6 = NA, Y7 = NA, Y8 = NA
    )
    
    results_list[[k]] <- results_df
    
    names(results_list)[k] <- paste0("results_", outcome_time_of_interest[k])
    
  }
  
  # get GLMs (drop Cox models for now)
  glm_models <- keep(model_fit_list, ~ .x$model_type == "glm")
  
  # separate markers into families
  m1 <- markers[grep("^m[1]\\.", names(markers))]
  m2 <- markers[grep("^m[2]\\.", names(markers))]
  m3 <- markers[grep("^m[3]\\.", names(markers))]
  m4 <- markers[grep("^m[4]\\.", names(markers))]
  m5 <- markers[grep("^m[5]\\.", names(markers))]
  m6 <- markers[grep("^m[6]\\.", names(markers))]
  
  # function that makes empty data.frames to store simulation results
  make_empty_df <- function(med_fam){
    df = as.data.frame(matrix(NA, nrow = mc_B, ncol = length(names(med_fam)), 
                              dimnames = list(NULL, names(med_fam))))
    return(df)
  }
  
  # list of simulation results
  m_sim_list <- list(
    m1 = list(
      m1_a = make_empty_df(m1),
      m1_aprime = make_empty_df(m1)
    ),
    m2 = list(
      m2_a_m1_a = make_empty_df(m2),
      m2_aprime_m1_aprime = make_empty_df(m2),
      m2_a_m1_aprime = make_empty_df(m2)
    ),
    m3 = list(
      m3_a_m1_a = make_empty_df(m3),
      m3_aprime_m1_aprime = make_empty_df(m3),
      m3_a_m1_aprime = make_empty_df(m3)
    ),
    m4 = list(
      m4_a = make_empty_df(m4),
      m4_aprime = make_empty_df(m4)
    ),
    m5 = list(
      m5_a_m4_a = make_empty_df(m5),
      m5_aprime_m4_aprime = make_empty_df(m5),
      m5_a_m4_aprime = make_empty_df(m5)
    ),
    m6 = list(
      m6_a_m4_a = make_empty_df(m6),
      m6_aprime_m4_aprime = make_empty_df(m6),
      m6_a_m4_aprime = make_empty_df(m6)
    )
  )
  
  # iterate over each row of the original COVAIL dataset
  for (i in 1:nrow(data)) {
    
    # create df that's covariates for row i repeated B times
    covariates_i <- data[i, ] %>%
      select(all_of(covariates))
    
    covariates_i_rep <- covariates_i[rep(1, mc_B), , drop = FALSE]
    
    ### m1 ###
    # iterate over sub-markers
    for (m1_index in 1:length(m1)) {
      
      # iterate over exposures
      for (A in c("a", "aprime")) {
       
        # get corresponding model and model name
        current_model_name <- paste0(names(m1[m1_index]), "_", A)
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates
        if (m1_index == 1) {
          
          current_data <- covariates_i
          
        } else {
          # for subsequent sub-markers, current data is covariates + previously simulated sub-markers
          pred_meds <- names(m1[1:(m1_index - 1)])
          
          # pull columns from simulated results
          if (A == "a") {
            current_data_simulated <- 
              m_sim_list$m1$m1_a[, colnames(m_sim_list$m1$m1_a) %in% pred_meds, drop = FALSE] 
          }
          if (A == "aprime") {
            current_data_simulated <- 
              m_sim_list$m1$m1_aprime[, colnames(m_sim_list$m1$m1_aprime) %in% pred_meds, drop = FALSE] 
          }
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A == "a") {m_sim_list$m1$m1_a[m1_index] <- m_star}
        if (A == "aprime") {m_sim_list$m1$m1_aprime[m1_index] <- m_star}
      
      } # end m1, a & aprime
    } # end m1
    
    ###########################################################################
    
    ### m2 ###
    for (m2_index in 1:length(m2)) {
      
      for (A_m1_A in c("a_m1_a", "aprime_m1_aprime", "a_m1_aprime")) {
        
        # for M2(a, M1(a)) and M2(a, M1(a')), use "a" model
        if (A_m1_A %in% c("a_m1_a", "a_m1_aprime")) {
          current_model_name <- paste0(names(m2[m2_index]), "_", "a")
        }
        
        # for M2(a', M1(a')), use "aprime" model
        if (A_m1_A == "aprime_m1_aprime") {
          current_model_name <- paste0(names(m2[m2_index]), "_", "aprime")
        }
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates and m1
        if (m2_index == 1) {
          
          if (A_m1_A == "a_m1_a") {
            current_data_simulated <- m_sim_list$m1$m1_a
          }
          if (A_m1_A %in% c("aprime_m1_aprime", "a_m1_aprime")) {
            current_data_simulated <- m_sim_list$m1$m1_aprime
          }
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
          
        } else {
          # for subsequent sub-markers, current data is covariates, m1, m2, and previously simulated sub-markers
          pred_meds <- names(c(m1, m2[1:(m2_index - 1)]))
          
          # pull columns from simulated results
          if (A_m1_A %in% c("a_m1_a", "a_m1_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            
            if (A_m1_A == "a_m1_a") {
              m2_simulated <- m_sim_list$m2$m2_a_m1_a[, colnames(m_sim_list$m2$m2_a_m1_a) %in% pred_meds, drop = FALSE]
            }
            
            if (A_m1_A == "a_m1_aprime") {
              m2_simulated <- m_sim_list$m2$m2_a_m1_aprime[, colnames(m_sim_list$m2$m2_a_m1_aprime) %in% pred_meds, drop = FALSE]
            }
            
          }
          if (A_m1_A == "aprime_m1_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime[, colnames(m_sim_list$m2$m2_aprime_m1_aprime) %in% pred_meds, drop = FALSE]
          }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A_m1_A == "a_m1_a") {m_sim_list$m2$m2_a_m1_a[m2_index] <- m_star}
        if (A_m1_A == "aprime_m1_aprime") {m_sim_list$m2$m2_aprime_m1_aprime[m2_index] <- m_star}
        if (A_m1_A == "a_m1_aprime") {m_sim_list$m2$m2_a_m1_aprime[m2_index] <- m_star}
        
      } # end m2, a_m1_a, aprime_m1_aprime, a_m1_aprime
    } # end m2
    
    ###########################################################################
    
    ### m3 ###
    for (m3_index in 1:length(m3)) {
      
      for (A_m1_A in c("a_m1_a", "aprime_m1_aprime", "a_m1_aprime")) {
        
        # for M3(a, M1(a)) and M3(a, M1(a')), use "a" model
        if (A_m1_A %in% c("a_m1_a", "a_m1_aprime")) {
          current_model_name <- paste0(names(m3[m3_index]), "_", "a")
        }
        
        # for M3(a', M1(a')), use "aprime" model
        if (A_m1_A == "aprime_m1_aprime") {
          current_model_name <- paste0(names(m3[m3_index]), "_", "aprime")
        }
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates and m1
        if (m3_index == 1) {
          
          if (A_m1_A == "a_m1_a") {
            current_data_simulated <- m_sim_list$m1$m1_a
          }
          if (A_m1_A %in% c("aprime_m1_aprime", "a_m1_aprime")) {
            current_data_simulated <- m_sim_list$m1$m1_aprime
          }

          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
          
        } else {
          # for subsequent sub-markers, current data is covariates, m1, m3, and previously simulated sub-markers
          pred_meds <- names(c(m1, m3[1:(m3_index - 1)]))
          
          # pull columns from simulated results
          if (A_m1_A %in% c("a_m1_a", "a_m1_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            
            if (A_m1_A == "a_m1_a") {
              m3_simulated <- m_sim_list$m3$m3_a_m1_a[, colnames(m_sim_list$m3$m3_a_m1_a) %in% pred_meds, drop = FALSE]
            }
            
            if (A_m1_A == "a_m1_aprime") {
              m3_simulated <- m_sim_list$m3$m3_a_m1_aprime[, colnames(m_sim_list$m3$m3_a_m1_aprime) %in% pred_meds, drop = FALSE]
            }
          }
          if (A_m1_A == "aprime_m1_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime[, colnames(m_sim_list$m3$m3_aprime_m1_aprime) %in% pred_meds, drop = FALSE]
          }
          
          current_data_simulated <- cbind(m1_simulated, m3_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A_m1_A == "a_m1_a") {m_sim_list$m3$m3_a_m1_a[m3_index] <- m_star}
        if (A_m1_A == "aprime_m1_aprime") {m_sim_list$m3$m3_aprime_m1_aprime[m3_index] <- m_star}
        if (A_m1_A == "a_m1_aprime") {m_sim_list$m3$m3_a_m1_aprime[m3_index] <- m_star}
        
      } # end m3, a_m1_a, aprime_m1_aprime, a_m1_aprime
    } # end m3
    
    ###########################################################################
    
    ### m4 ###
    for (m4_index in 1:length(m4)) {
      
      for (A in c("a", "aprime")) {
        
        current_model_name <- paste0(names(m4[m4_index]), "_", A)
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates, m1, m2, and m3
        if (m4_index == 1) {
          
          if (A == "a") {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            }
          if (A == "aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated, m3_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
          
        } else {
          # for subsequent sub-markers, current data is covariates, m1, m2, m3, and previously simulated sub-markers
          pred_meds <- names(c(m1, m2, m3, m4[1:(m4_index - 1)]))
          
          # pull columns from simulated results
          if (A == "a") {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            m4_simulated <- m_sim_list$m4$m4_a[, colnames(m_sim_list$m4$m4_a) %in% pred_meds, drop = FALSE]
          }
          if (A == "aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            m4_simulated <- m_sim_list$m4$m4_aprime[, colnames(m_sim_list$m4$m4_aprime) %in% pred_meds, drop = FALSE]
          }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated, 
                                          m3_simulated, m4_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A == "a") {m_sim_list$m4$m4_a[m4_index] <- m_star}
        if (A == "aprime") {m_sim_list$m4$m4_aprime[m4_index] <- m_star}
        
      } # end m4, a & aprime
    } # end m4
    
    ###########################################################################
    
    ### m5 ###
    for (m5_index in 1:length(m5)) {
      
      for (A_m4_A in c("a_m4_a", "aprime_m4_aprime", "a_m4_aprime")) {
        
        # for M5(a, M4(a)) and M5(a_M4(a')), use "a" model
        if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
          current_model_name <- paste0(names(m5[m5_index]), "_", "a")
        }
        
        # for M5(a', M4(a')), use "aprime" model
        if (A_m4_A == "aprime_m4_aprime") {
          current_model_name <- paste0(names(m5[m5_index]), "_", "aprime")
        }
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates, m1, m2, m3, and m4
        if (m5_index == 1) {
          
          if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            
            if (A_m4_A == "a_m4_a") {
              m4_simulated <- m_sim_list$m4$m4_a
            }
            
            if (A_m4_A == "a_m4_aprime") {
              m4_simulated <- m_sim_list$m4$m4_aprime
            }

          }
          
          if (A_m4_A == "aprime_m4_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            m4_simulated <- m_sim_list$m4$m4_aprime
          }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated, 
                                          m3_simulated, m4_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
          
        } else {
          # for subsequent sub-markers, current data is covariates, m1, m2, m3, m4, and previously simulated sub-markers
          pred_meds <- names(c(m1, m2, m3, m4, m5[1:(m5_index - 1)]))
          
          # pull columns from simulated results
          if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            
            if (A_m4_A == "a_m4_a") {
              m4_simulated <- m_sim_list$m4$m4_a
              m5_simulated <- m_sim_list$m5$m5_a_m4_a[, colnames(m_sim_list$m5$m5_a_m4_a) %in% pred_meds, drop = FALSE]
            }
            
            if (A_m4_A == "a_m4_aprime") {
              m4_simulated <- m_sim_list$m4$m4_aprime
              m5_simulated <- m_sim_list$m5$m5_a_m4_aprime[, colnames(m_sim_list$m5$m5_a_m4_aprime) %in% pred_meds, drop = FALSE]
            }
          }
          if (A_m4_A == "aprime_m4_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            m4_simulated <- m_sim_list$m4$m4_aprime
            m5_simulated <- m_sim_list$m5$m5_aprime_m4_aprime[, colnames(m_sim_list$m5$m5_aprime_m4_aprime) %in% pred_meds, drop = FALSE]
          }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated, 
                                          m3_simulated, m4_simulated,
                                          m5_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A_m4_A == "a_m4_a") {m_sim_list$m5$m5_a_m4_a[m5_index] <- m_star}
        if (A_m4_A == "aprime_m4_aprime") {m_sim_list$m5$m5_aprime_m4_aprime[m5_index] <- m_star}
        if (A_m4_A == "a_m4_aprime") {m_sim_list$m5$m5_a_m4_aprime[m5_index] <- m_star}
        
      } # end m5, a_m4_a, aprime_m4_aprime, a_m4_aprime
    } # end m5
    
    ###########################################################################
    
    ### m6 ###
    for (m6_index in 1:length(m6)) {
      
      for (A_m4_A in c("a_m4_a", "aprime_m4_aprime", "a_m4_aprime")) {
        
        # for M6(a, M4(a)) and M6(a_M4(a')), use "a" model
        if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
          current_model_name <- paste0(names(m6[m6_index]), "_", "a")
        }
        
        # for M6(a', M4(a')), use "aprime" model
        if (A_m4_A == "aprime_m4_aprime") {
          current_model_name <- paste0(names(m6[m6_index]), "_", "aprime")
        }
        
        current_model <- glm_models[[current_model_name]]
        
        # for first sub-marker, current data is covariates, m1, m2, m3, and m4
        if (m6_index == 1) {
          
          if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            
            if (A_m4_A == "a_m4_a") {
              m4_simulated <- m_sim_list$m4$m4_a
            }
            
            if (A_m4_A == "a_m4_aprime") {
              m4_simulated <- m_sim_list$m4$m4_aprime
            }
          }
          if (A_m4_A == "aprime_m4_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            m4_simulated <- m_sim_list$m4$m4_aprime
          }

          current_data_simulated <- cbind(m1_simulated, m2_simulated, 
                                          m3_simulated, m4_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
          
        } else {
          # for subsequent sub-markers, current data is covariates, m1, m2, m3, m4, and previously simulated sub-markers
          pred_meds <- names(c(m1, m2, m3, m4, m6[1:(m6_index - 1)]))
          
          # pull columns from simulated results
          if (A_m4_A %in% c("a_m4_a", "a_m4_aprime")) {
            m1_simulated <- m_sim_list$m1$m1_a
            m2_simulated <- m_sim_list$m2$m2_a_m1_a
            m3_simulated <- m_sim_list$m3$m3_a_m1_a
            
            if (A_m4_A == "a_m4_a") {
              m4_simulated <- m_sim_list$m4$m4_a
              m6_simulated <- m_sim_list$m6$m6_a_m4_a[, colnames(m_sim_list$m6$m6_a_m4_a) %in% pred_meds, drop = FALSE]
            }
            
            if (A_m4_A == "a_m4_aprime") {
              m4_simulated <- m_sim_list$m4$m4_aprime
              m6_simulated <- m_sim_list$m6$m6_a_m4_aprime[, colnames(m_sim_list$m6$m6_a_m4_aprime) %in% pred_meds, drop = FALSE]
            }
          }
          
          if (A_m4_A == "aprime_m4_aprime") {
            m1_simulated <- m_sim_list$m1$m1_aprime
            m2_simulated <- m_sim_list$m2$m2_aprime_m1_aprime
            m3_simulated <- m_sim_list$m3$m3_aprime_m1_aprime
            m4_simulated <- m_sim_list$m4$m4_aprime
            m6_simulated <- m_sim_list$m6$m6_aprime_m4_aprime[, colnames(m_sim_list$m6$m6_aprime_m4_aprime) %in% pred_meds, drop = FALSE]
          }
          
          current_data_simulated <- cbind(m1_simulated, m2_simulated, 
                                          m3_simulated, m4_simulated,
                                          m6_simulated)
          
          # add covariates
          current_data <- covariates_i_rep %>%
            cbind(current_data_simulated)
        }
        
        # predict
        mu <- predict(
          current_model$fit,
          current_data,
          type = "response"
        )
        
        # simulate
        set.seed(seed)
        
        m_star <- rnorm(
          n = mc_B,
          mean = mu,
          sd = sd(residuals(current_model$fit))
        )
        
        # save simulated results
        if (A_m4_A == "a_m4_a") {m_sim_list$m6$m6_a_m4_a[m6_index] <- m_star}
        if (A_m4_A == "aprime_m4_aprime") {m_sim_list$m6$m6_aprime_m4_aprime[m6_index] <- m_star}
        if (A_m4_A == "a_m4_aprime") {m_sim_list$m6$m6_a_m4_aprime[m6_index] <- m_star}
        
      } # end m6, a_m4_a, aprime_m4_aprime, a_m4_aprime
    } # end m6
    
    ############################################################################
    
    # iterate over outcome times of interest
    for (k in 1:length(outcome_time_of_interest)) {
      
      # for all versions of Y, use this base data
      base_data <- data.frame(EventIndPrimaryD15 = NA,
                              EventTimePrimaryD15 = outcome_time_of_interest[k],
                              covariates_i_rep)
      
      ##########
      
      # Y1
      data_Y1 <- cbind(
        base_data,
        m_sim_list$m1$m1_a, 
        m_sim_list$m2$m2_a_m1_a,
        m_sim_list$m3$m3_a_m1_a,
        m_sim_list$m4$m4_a,
        m_sim_list$m5$m5_a_m4_a,
        m_sim_list$m6$m6_a_m4_a
      )
      
      # predict using Cox model
      Y1 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y1,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
        ]]$Y1[i] <- mean(Y1)

      ##########
      
      # Y2
      data_Y2 <- cbind(
        base_data,
        m_sim_list$m1$m1_aprime, 
        m_sim_list$m2$m2_aprime_m1_aprime,
        m_sim_list$m3$m3_aprime_m1_aprime,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_aprime_m4_aprime
      )
    
      # predict using Cox model
      Y2 <- 1 - predict(
        model_fit_list$cox_given_aprime$fit,
        newdata = data_Y2,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y2[i] <- mean(Y2)
      
      ##########
      
      # Y3
      data_Y3 <- cbind(
        base_data,
        m_sim_list$m1$m1_a, 
        m_sim_list$m2$m2_a_m1_a,
        m_sim_list$m3$m3_a_m1_a,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_a_m4_aprime,
        m_sim_list$m6$m6_a_m4_aprime
      )
      
      # predict using Cox model
      Y3 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y3,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y3[i] <- mean(Y3)
      
      ##########
      
      # Y4
      data_Y4 <- cbind(
        base_data,
        m_sim_list$m1$m1_a, 
        m_sim_list$m2$m2_a_m1_a,
        m_sim_list$m3$m3_a_m1_a,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_a_m4_aprime
      )
      
      # predict using Cox model
      Y4 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y4,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y4[i] <- mean(Y4)
      
      ##########
      
      # Y5
      data_Y5 <- cbind(
        base_data,
        m_sim_list$m1$m1_a, 
        m_sim_list$m2$m2_a_m1_a,
        m_sim_list$m3$m3_a_m1_a,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_aprime_m4_aprime
      )
      
      # predict using Cox model
      Y5 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y5,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y5[i] <- mean(Y5)
      
      ##########
      
      # Y6
      data_Y6 <- cbind(
        base_data,
        m_sim_list$m1$m1_aprime, 
        m_sim_list$m2$m2_a_m1_aprime,
        m_sim_list$m3$m3_a_m1_aprime,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_aprime_m4_aprime
      )
      
      # predict using Cox model
      Y6 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y6,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y6[i] <- mean(Y6)
      
      ##########
      
      # Y7
      data_Y7 <- cbind(
        base_data,
        m_sim_list$m1$m1_aprime, 
        m_sim_list$m2$m2_aprime_m1_aprime,
        m_sim_list$m3$m3_a_m1_aprime,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_aprime_m4_aprime
      )
      
      # predict using Cox model
      Y7 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y7,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y7[i] <- mean(Y7)
      
      ##########
      
      # Y8
      data_Y8 <- cbind(
        base_data,
        m_sim_list$m1$m1_aprime, 
        m_sim_list$m2$m2_aprime_m1_aprime,
        m_sim_list$m3$m3_aprime_m1_aprime,
        m_sim_list$m4$m4_aprime,
        m_sim_list$m5$m5_aprime_m4_aprime,
        m_sim_list$m6$m6_aprime_m4_aprime
      )
      
      # predict using Cox model
      Y8 <- 1 - predict(
        model_fit_list$cox_given_a$fit,
        newdata = data_Y8,
        type = "survival"
      )
      
      # save results
      results_list[[
        paste0("results_", outcome_time_of_interest[k])
      ]]$Y8[i] <- mean(Y8)
      
    } # end for loop over outcome times
    
  } # end loop over COVAIL data
  
  # take weighted mean over all observations
  results_list_final <- vector(mode = "list", length = length(outcome_time_of_interest))
  
  for (k in 1:length(outcome_time_of_interest)) {
    
    results_list_final[[k]] <- results_list[[k]] %>%
      summarize(
        mean_Y1 = sum((weights) * Y1) / sum(weights),
        mean_Y2 = sum((weights) * Y2) / sum(weights),
        mean_Y3 = sum((weights) * Y3) / sum(weights),
        mean_Y4 = sum((weights) * Y4) / sum(weights),
        mean_Y5 = sum((weights) * Y5) / sum(weights),
        mean_Y6 = sum((weights) * Y6) / sum(weights),
        mean_Y7 = sum((weights) * Y7) / sum(weights),
        mean_Y8 = sum((weights) * Y8) / sum(weights)
      )
    
    names(results_list_final)[k] <- paste0("results_", outcome_time_of_interest[k])
    
  }
  
  # return results
  return(results_list_final)
  
} # end fcn





