estimator_ml_for_nuisance <- function(
        dta,
        id_var,
        moderator_var,
        control_var, # vector of variables
        trt_var = "A",
        outcome_var = "Y",
        prob_A_var = "prob_A",
        avail_var = "avail",
        ml_method = c("gam", "lm", "rf", "ranger", "sl.smooth", "sl.all", "zero"),
        control_var_spline = NULL, # the covariates that will be constructed splines for
        cross_fit = FALSE,
        cf_fold = 10,
        weighting_function_var = NULL # omega(t)
) {
    # gam: generalized additive model
    # lm: linear regression
    # rf: random forest
    # ranger: fast implementation of random forest
    # sl: super learner
    #   sl.smooth: sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
    #   sl.all:    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    #                             "SL.ranger", "SL.xgboost", "SL.nnet")
    # zero: set mu_hat_0 and mu_hat_1 to all be 0
    
    ml_method <- match.arg(ml_method)
    dta <- as.data.frame(dta)
    
    if (ml_method == "zero") {
        dta$mu_hat_0 <- dta$mu_hat_1 <- 0
        regfit_a0 <- regfit_a1 <- NULL
    } else {
        if (ml_method == "sl.smooth") {
            sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
        } else if (ml_method == "sl.all") {
            sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
                           "SL.ranger", "SL.xgboost", "SL.nnet")
        }
        
        if (ml_method == "gam") {
            if (is.null(control_var_spline)) {
                control_var_spline <- control_var
                control_var_notspline <- NULL
            } else {
                control_var_notspline <- setdiff(control_var, control_var_spline)
            }
        }
        
        if (cross_fit) {
            # using cross-fitting as implemented in AIPW R package
            # (see their supplementary material for a illustrative diagram)
            
            id_permuted <- sample(unique(dta[[id_var]]))
            id_permuted <- c(id_permuted, 
                             rep(NA, cf_fold - length(id_permuted) %% cf_fold))
            id_folds <- matrix(id_permuted, ncol = cf_fold, byrow = TRUE)
            
            dta_analysis <- data.frame()
            for (k in 1:cf_fold) {
                id_holdout <- na.omit(id_folds[, k])
                id_train <- na.omit(as.vector(id_folds[, -k]))
                
                dta_holdout <- dta[dta[[id_var]] %in% id_holdout, ]
                dta_train <- dta[dta[[id_var]] %in% id_train, ]
                
                dta_train_a0 <- dta_train[dta_train[[trt_var]] == 0, ]
                dta_train_a1 <- dta_train[dta_train[[trt_var]] == 1, ]
                
                if (ml_method == "gam") {
                    if (length(control_var_notspline) == 0) {
                        gam_formula <- as.formula(paste0(outcome_var, " ~ ", paste0("s(", control_var_spline, ", bs = 'cr')", collapse = " + ")))
                    } else {
                        gam_formula <- as.formula(paste0(outcome_var, " ~ ", paste0("s(", control_var_spline, ", bs = 'cr')", collapse = " + "), 
                                                         " + ", paste0(control_var_notspline, collapse = " + ")))
                    }
                    regfit_a0 <- gam(gam_formula, data = dta_train_a0)
                    regfit_a1 <- gam(gam_formula, data = dta_train_a1)
                } else if (ml_method == "lm") {
                    lm_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                    regfit_a0 <- lm(lm_formula, data = dta_train_a0)
                    regfit_a1 <- lm(lm_formula, data = dta_train_a1)
                } else if (ml_method == "rf") {
                    rf_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                    regfit_a0 <- randomForest(rf_formula, data = dta_train_a0)
                    regfit_a1 <- randomForest(rf_formula, data = dta_train_a1)
                } else if (ml_method == "ranger") {
                    rf_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                    regfit_a0 <- ranger(rf_formula, data = dta_train_a0)
                    regfit_a1 <- ranger(rf_formula, data = dta_train_a1)
                } else if (ml_method %in% c("sl.smooth", "sl.all")) {
                    regfit_a0 <- SuperLearner(Y = dta_train_a0 %>% pull(!!outcome_var),
                                              X = dta_train_a0 %>% select(!!control_var),
                                              family = gaussian(),
                                              verbose = FALSE,
                                              SL.library = sl.library)
                    regfit_a1 <- SuperLearner(Y = dta_train_a1 %>% pull(!!outcome_var),
                                              X = dta_train_a1 %>% select(!!control_var),
                                              family = gaussian(),
                                              verbose = FALSE,
                                              SL.library = sl.library)
                }
                
                # predicting eta_hat for the holdout set
                if (ml_method %in% c("sl.smooth", "sl.all")) {
                    newdata_df <- dta_holdout %>% select(!!control_var)
                    dta_holdout <- dta_holdout %>%
                        mutate(mu_hat_0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
                               mu_hat_1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
                } else if (ml_method == "ranger") {
                    dta_holdout <- dta_holdout %>%
                        mutate(mu_hat_0 = predict(regfit_a0, data = dta_holdout, type = "response")$predictions,
                               mu_hat_1 = predict(regfit_a1, data = dta_holdout, type = "response")$predictions)
                } else if (ml_method %in% c("gam", "lm", "rf")) {
                    dta_holdout <- dta_holdout %>%
                        mutate(mu_hat_0 = predict(regfit_a0, newdata = dta_holdout, type = "response"),
                               mu_hat_1 = predict(regfit_a1, newdata = dta_holdout, type = "response"))
                }
                
                dta_analysis <- rbind(dta_analysis, dta_holdout)
            }
            dta <- dta_analysis
            
        } else {
            # not using cross-fitting
            
            dta_a0 <- dta[dta[[trt_var]] == 0, ]
            dta_a1 <- dta[dta[[trt_var]] == 1, ]
            
            if (ml_method == "gam") {
                if (length(control_var_notspline) == 0) {
                    gam_formula <- as.formula(paste0(outcome_var, " ~ ", paste0("s(", control_var_spline, ", bs = 'cr')", collapse = " + ")))
                } else {
                    gam_formula <- as.formula(paste0(outcome_var, " ~ ", paste0("s(", control_var_spline, ", bs = 'cr')", collapse = " + "), 
                                                     " + ", paste0(control_var_notspline, collapse = " + ")))
                }
                regfit_a0 <- gam(gam_formula, data = dta_a0)
                regfit_a1 <- gam(gam_formula, data = dta_a1)
            } else if (ml_method == "lm") {
                lm_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                regfit_a0 <- lm(lm_formula, data = dta_a0)
                regfit_a1 <- lm(lm_formula, data = dta_a1)
            } else if (ml_method == "rf") {
                rf_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                regfit_a0 <- randomForest(rf_formula, data = dta_a0)
                regfit_a1 <- randomForest(rf_formula, data = dta_a1)
            } else if (ml_method == "ranger") {
                rf_formula <- as.formula(paste0(outcome_var, " ~ ", paste0(control_var, collapse = " + ")))
                regfit_a0 <- ranger(rf_formula, data = dta_a0)
                regfit_a1 <- ranger(rf_formula, data = dta_a1)
            } else if (ml_method %in% c("sl.smooth", "sl.all")) {
                regfit_a0 <- SuperLearner(Y = dta_a0 %>% pull(!!outcome_var),
                                          X = dta_a0 %>% select(!!control_var),
                                          family = gaussian(),
                                          verbose = FALSE,
                                          SL.library = sl.library)
                regfit_a1 <- SuperLearner(Y = dta_a1 %>% pull(!!outcome_var),
                                          X = dta_a1 %>% select(!!control_var),
                                          family = gaussian(),
                                          verbose = FALSE,
                                          SL.library = sl.library)
            }
            
            # predicting eta_hat
            if (ml_method %in% c("sl.smooth", "sl.all")) {
                newdata_df <- dta %>% select(!!control_var)
                dta <- dta %>%
                    mutate(mu_hat_0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
                           mu_hat_1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
            } else if (ml_method == "ranger") {
                dta <- dta %>%
                    mutate(mu_hat_0 = predict(regfit_a0, data = dta, type = "response")$predictions,
                           mu_hat_1 = predict(regfit_a1, data = dta, type = "response")$predictions)
            } else if (ml_method %in% c("gam", "lm", "rf")) {
                dta <- dta %>%
                    mutate(mu_hat_0 = predict(regfit_a0, newdata = dta, type = "response"),
                           mu_hat_1 = predict(regfit_a1, newdata = dta, type = "response"))
            }
            
        }
    }
    
    fit <- estimator_core(
        dta = dta,
        id_var = id_var, 
        moderator_var = moderator_var,
        trt_var = trt_var,
        outcome_var = outcome_var,
        mu_hat_1_var = "mu_hat_1",
        mu_hat_0_var = "mu_hat_0",
        prob_A_var = prob_A_var,
        avail_var = avail_var
    )
    
    # return(fit)
    return(c(fit, list(regfit_a0 = regfit_a0, regfit_a1 = regfit_a1)))
}


estimator_core <- function(
        dta,
        id_var,
        moderator_var,
        trt_var = "A",
        outcome_var = "Y",
        mu_hat_1_var = "mu_hat_1",
        mu_hat_0_var = "mu_hat_0",
        prob_A_var = "prob_A",
        avail_var = "avail",
        weighting_function_var = NULL # omega(t)
        ){
    
    total_person_decisionpoint <- nrow(dta)
    sample_size <- length(unique(dta[[id_var]]))
    total_T <- total_person_decisionpoint / sample_size # assuming everyone has the same number of dp
    
    if (is.null(weighting_function_var)) {
        omega <- rep(1/total_T, total_person_decisionpoint)
    } else {
        omega <- dta[[weighting_function_var]]
    }
    
    pt1 <- dta[, prob_A_var] # this is p_t(1|H_t)
    A <- dta[, trt_var]
    Sdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_var] ) )
    Y <- dta[, outcome_var]
    mu_hat_1 <- dta[, mu_hat_1_var]
    mu_hat_0 <- dta[, mu_hat_0_var]
    avail <- dta[, avail_var]
    
    pt0 <- 1 - pt1 # this is p_t(0|H_t)
    pA <- ifelse(A, pt1, pt0) # this is p_t(A_t|H_t)
    
    dim_beta <- ncol(Sdm)
    
    person_first_index <- c(find_change_location(dta[, id_var]), total_person_decisionpoint + 1)
    sample_size <- length(person_first_index) - 1
    
    #### 1. Calculate beta_hat using its analytic form ####
    # from Goodnotes "2025.01.19 - Asymptotic Variance (updated)"
    
    ## compute the front inverse term
    
    sum_SStrans <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (it in 1:total_person_decisionpoint) {
        S_it <- matrix(Sdm[it, ], ncol = 1)
        sum_SStrans <- sum_SStrans + omega[it] * S_it %*% t(S_it)
    }
    avg_SStrans <- sum_SStrans / sample_size
    
    ## compute the second term (involving a residual)
    
    residual <- omega * avail * (-1)^(1-A) / pA * (Y - pt0 * mu_hat_1 - pt1 * mu_hat_0)
    sum_residualS <- matrix(0, nrow = dim_beta, ncol = 1)
    for (it in 1:total_person_decisionpoint) {
        S_it <- matrix(Sdm[it, ], ncol = 1)
        r_it <- as.numeric(residual[it])
        sum_residualS <- sum_residualS + r_it * S_it
    }
    avg_residualS <- sum_residualS / sample_size
    
    ## compute beta_hat
    
    beta_hat <- solve(avg_SStrans) %*% avg_residualS # a single column matrix
    
    #### 2. Calculate standard error of beta_hat ####
    
    ## The derivative (bread) term is avg_SStrans
    
    bread <- avg_SStrans
    
    ## Calculate the cross-product (meat) term
    
    meat <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (i in 1:sample_size) {
        rowid <- person_first_index[i] : (person_first_index[i+1] - 1)
        sum_t_ee_it <- numeric(dim_beta)
        for (it in rowid) {
            S_it <- matrix(Sdm[it, ], ncol = 1)
            r_it <- as.numeric(residual[it]) 
            omega_it <- omega[it]
            
            ee_it <- r_it * S_it - omega_it * S_it %*% t(S_it) %*% beta_hat
            sum_t_ee_it <- sum_t_ee_it + ee_it
        }
        meat <- meat + sum_t_ee_it %*% t(sum_t_ee_it)
    }
    meat <- meat / sample_size
    
    ## Calculate the variance of beta_hat using the sandwich estimator 
    bread_inv <- solve(bread)
    beta_varcov <- bread_inv %*% meat %*% t(bread_inv) / sample_size
    beta_se <- sqrt(diag(beta_varcov))
    
    Snames <- c("Intercept", moderator_var)
    beta_hat <- as.vector(beta_hat)
    names(beta_hat) <- names(beta_se) <- Snames
    colnames(beta_varcov) <- rownames(beta_varcov) <- Snames
    
    #### 3. calculate confidence interval

    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05/2, df = sample_size - dim_beta)
    conf_int_tquantile <- cbind(beta_hat - c * beta_se,
                               beta_hat + c * beta_se)
    colnames(conf_int) <- colnames(conf_int_tquantile) <- c("2.5 %", "97.5 %")
    
    return(list(beta_hat = beta_hat,
                beta_se = beta_se,
                conf_int = conf_int,
                conf_int_tquantile = conf_int_tquantile,
                beta_varcov = beta_varcov))
    
}



#' Find the locations at which the value changes from the previous one
#'
#' @param v a vector that specifies different values
#'
#' @return Return a vector that specifies the indexes of locations
#'         at which the value changes from the previous one
#'
#' @noRd
#' @examples v <- c("a", "a", "b", "c", "c")
#'           find_change_location(v)
find_change_location <- function(v){
    n <- length(v)
    if (n <= 1) {
        stop("The vector need to have length > 1.")
    }
    return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}
