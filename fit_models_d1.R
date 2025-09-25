
# 8/30/25

# fit Cox models and GLMS

fit_models <- function(
    markers, 
    exposure,
    covariates, 
    outcome_time, 
    outcome_index,
    weights,
    a,
    aprime,
    data
) {
  
  # create empty list to store model info
  # a and a' versions of Y and m1.X - m6.X sub-markers
  # easier to assign items to list if it's not pre-specified length
  model_fit_list <- list()
  
  ### Cox models ###
  # wrap Cox models in tryCatch so if either a or a' doesn't converge, we exit the fit_models function
  # A = a
  cox_fit_a <- tryCatch(
    {
      coxph(
        formula = as.formula(paste0(
          "Surv(", outcome_time, ",", outcome_index, ") ~", 
          paste(c(covariates, names(markers)), collapse = " + ")
        )),
        data = data[data[[exposure]] == a,],
        weights = data[data[[exposure]] == a, ][[weights]],
        model = TRUE
      )
    },
    error = function(e) {
      message("Cox model (A=a) failed: ", e$message)
      FALSE
    }
  )
  
  cox_fit_aprime <- tryCatch(
    {
      coxph(
        formula = as.formula(paste0(
          "Surv(", outcome_time, ",", outcome_index, ") ~", 
          paste(c(covariates, names(markers)), collapse = " + ")
        )),
        data = data[data[[exposure]] == aprime,],
        weights = data[data[[exposure]] == aprime, ][[weights]],
        model = TRUE
      )
    },
    error = function(e) {
      message("Cox model (A=aprime) failed: ", e$message)
      FALSE
    }
  )
  
  # if both Cox models converge, continue on to save fits and run glm models
  if (!isFALSE(cox_fit_a) && !isFALSE(cox_fit_aprime)) {
    
    # save cox model info and fits in list
    model_fit_list[["cox_given_a"]] <- list(
      model_type = "cox",
      given_exposure = "a",
      fit = cox_fit_a
    )
    
    model_fit_list[["cox_given_aprime"]] <- list(
      model_type = "cox",
      given_exposure = "aprime",
      fit = cox_fit_aprime
    )
    
    ### GLM's ###
    # separate markers into families
    m1 <- markers[grep("^m[1]\\.", names(markers))]
    m2 <- markers[grep("^m[2]\\.", names(markers))]
    m3 <- markers[grep("^m[3]\\.", names(markers))]
    m4 <- markers[grep("^m[4]\\.", names(markers))]
    m5 <- markers[grep("^m[5]\\.", names(markers))]
    m6 <- markers[grep("^m[6]\\.", names(markers))]
    
    # iterate through a and a'
    for (A in c("a", "aprime")) {
      
      ### M1 ###
      for (m1_index in 1:length(m1)) {

        # outcome
        outcome = names(m1[m1_index])
        
        # predictors
        predictors = if (m1_index == 1) {
          # for m1.1, only use covariates
          covariates
        } else {
          # for following sub-markers, use covariates and previous sub-markers
          c(covariates, names(m1[1:(m1_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m1[m1_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M1 sub-markers
      
      ### M2 ###
      for (m2_index in 1:length(m2)) {
        
        # outcome
        outcome = names(m2[m2_index])
        
        # predictors
        predictors = if (m2_index == 1) {
          # for m2.1, only use covariates and m1
          c(covariates, names(m1))
        } else {
          # for following sub-markers, use covariates, m1, and previous sub-markers
          c(covariates, names(m1), names(m2[1:(m2_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m2[m2_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M2 sub-markers
      
      ### M3 ###
      for (m3_index in 1:length(m3)) {
        
        # outcome
        outcome = names(m3[m3_index])
        
        # predictors
        predictors = if (m3_index == 1) {
          # for m3.1, only use covariates and m1
          c(covariates, names(m1))
        } else {
          # for following sub-markers, use covariates, m1, and previous sub-markers
          c(covariates, names(m1), names(m3[1:(m3_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m3[m3_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M3 sub-markers
      
      ### M4 ###
      for (m4_index in 1:length(m4)) {
        
        # outcome
        outcome = names(m4[m4_index])
        
        # predictors
        predictors = if (m4_index == 1) {
          # for m4.1, only use covariates, m1, m2, and m3
          c(covariates, names(m1), names(m2), names(m3))
        } else {
          # for following sub-markers, use covariates, m1, m2, and m3, and previous sub-markers
          c(covariates, names(m1), names(m2), names(m3), names(m4[1:(m4_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m4[m4_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M4 sub-markers
      
      ### M5 ###
      for (m5_index in 1:length(m5)) {
        
        # outcome
        outcome = names(m5[m5_index])
        
        # predictors
        predictors = if (m5_index == 1) {
          # for m5.1, only use covariates, m1, m2, m3, and m4
          c(covariates, names(m1), names(m2), names(m3), names(m4))
        } else {
          # for following sub-markers, use covariates, m1, m2, m3, m4, and previous sub-markers
          c(covariates, names(m1), names(m2), names(m3), names(m4), names(m5[1:(m5_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m5[m5_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M5 sub-markers
      
      ### M6 ###
      for (m6_index in 1:length(m6)) {
        
        # outcome
        outcome = names(m6[m6_index])
        
        # predictors
        predictors = if (m6_index == 1) {
          # for m6.1, only use covariates, m1, m2, m3, and m4
          c(covariates, names(m1), names(m2), names(m3), names(m4))
        } else {
          # for following sub-markers, use covariates, m1, m2, m3, m4, and previous sub-markers
          c(covariates, names(m1), names(m2), names(m3), names(m4), names(m6[1:(m6_index - 1)]))
        }
        
        # save model type, outcome, and predictors in list
        model_fit_list[[paste0(names(m6[m6_index]), "_", A)]] <- list(
          model_type = "glm",
          given_exposure = A,
          outcome = outcome,
          predictors = predictors,
          fit = glm(
            formula = as.formula(paste0(
              outcome, " ~ ", paste(predictors, collapse = " + ")
            )),
            family = gaussian(),
            data = data[data[[exposure]] == get(A),],
            weights = data[data[[exposure]] == get(A), ][[weights]]
          ))
        
      } # end iterate over M6 sub-markers
      
    } # end iterate over a and a'
    
    return(model_fit_list)
    
  } else {
    
    return("cox_failed")
    
  }

} # end fcn
