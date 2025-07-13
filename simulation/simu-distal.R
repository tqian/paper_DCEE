library(tidyverse)
library(mgcv)
library(DTRreg)
library(geepack)
library(ranger)
library(MRTAnalysis)

args <- commandArgs(trailingOnly = TRUE)
itask <- as.integer(args[1])

nsim <- 1000 # to change to 1000
nsim_per_task <- 100 # to change to 100
seed <- 1:(nsim / nsim_per_task)

print_every_n_sims <- 10

simulation_design <- expand.grid(
    seed = seed,
    method = c("DCEE-lm", "DCEE-gam", "DCEE-ranger", 
               "DCEE-lm-cf", "DCEE-gam-cf", "DCEE-ranger-cf", 
               "GEE", "SNMM", "WCLS"),
    sample_size = c(30, 50, 100, 200, 300, 400, 500),
    estimand = c("marginal", "moderated")
)

simulation_design <- simulation_design %>%
    filter(!(method == "SNMM" & sample_size < 100))

# # for debugging the cluster run
# nsim <- 10
# nsim_per_task <- 10
# seed <- 1:(nsim / nsim_per_task)
# 
# simulation_design <- expand.grid(
#     estimand = c("marginal", "moderated"),
#     sample_size = c(100),
#     method = c("DCEE-lm", "DCEE-gam", "DCEE-lm-cf", "DCEE-gam-cf", "GEE", "SNMM"),
#     seed = seed
# )
# 
# itask <- 1


source("dgm_distal_outcome_endox_nonlinear_interaction_binaryZ.R")
source("estimator_distal_outcome.R")


estimand <- simulation_design$estimand[itask]
sample_size <- simulation_design$sample_size[itask]
method <- simulation_design$method[itask]
seed <- simulation_design$seed[itask]

print("Conducting simulation for...\n")

print(simulation_design[itask, ])

total_T <- 30
rand_prob_pattern <- "timevar"
availability <- "not-always"

set.seed(seed)

result_collected <- c()

start_time <- Sys.time()
for (isim in 1:nsim_per_task) {
    if (isim %% print_every_n_sims == 0) {
        current_time <- Sys.time()
        hours_diff <- round(difftime(current_time, start_time, units = "hours"), 2)
        cat(paste0("Starting isim: ", isim, "/", nsim_per_task, "; Hours lapsed: ", hours_diff, "\n"))
    }
    dta <- dgm_distal_outcome_endox_nonlinear_interaction_binaryZ(
        sample_size = sample_size, total_T = total_T,
        rand_prob_pattern = rand_prob_pattern,
        availability = availability
    )
    
    if (method == "DCEE-lm" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "lm",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-lm" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "lm",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-lm-cf" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "lm",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-lm-cf" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "lm",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-gam" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "gam",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-gam" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "gam",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-gam-cf" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "gam",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-gam-cf" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "gam",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-ranger" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "ranger",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-ranger" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "ranger",
            cross_fit = FALSE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-ranger-cf" & estimand == "marginal") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = NULL,
            ml_method = "ranger",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "DCEE-ranger-cf" & estimand == "moderated") {
        beta <- estimator_ml_for_nuisance(
            moderator_var = "Z",
            ml_method = "ranger",
            cross_fit = TRUE,
            cf_fold = 5,
            dta = dta, id_var = "userid", control_var = c("X", "Z"),
            trt_var = "A", outcome_var = "Y", prob_A_var = "prob_A",
            avail_var = "avail", control_var_spline = "X"
        )
    } else if (method == "GEE" & estimand == "marginal") {
        fit <- geeglm(Y ~ A + X + Z, data = dta, corstr = "independence", id = dta$userid)
        beta <- list(beta_hat = coef(fit)["A"], beta_se = sqrt(vcov(fit)["A", "A"]))
    } else if (method == "GEE" & estimand == "moderated") {
        fit <- geeglm(Y ~ A * Z + X + Z, data = dta, corstr = "independence", id = dta$userid)
        beta <- list(beta_hat = coef(fit)[c("A", "A:Z")],
                     beta_se = c(sqrt(vcov(fit)["A", "A"]), sqrt(vcov(fit)["A:Z", "A:Z"])))
    } else if (method == "SNMM" & estimand == "marginal") {
        # for specifying the models to be passed into DTRreg()
        generate_DTRreg_models <- function(total_T) {
            # Create blip.mod: List of ~ 1 for each time step
            blip.mod <- lapply(1:total_T, function(t) ~ 1)
            
            # Create treat.mod: Dynamic formula list based on A_t ~ X_t + Z_t
            treat.mod <- lapply(1:total_T, function(t) as.formula(paste0("A_", t, " ~ X_", t, " + Z_", t)))
            
            # Create tf.mod: Dynamic formula list, interaction with A_t starting from the second time step
            tf.mod <- lapply(1:total_T, function(t) {
                if (t == 1) {
                    ~ 1
                } else {
                    as.formula(paste0("~ (X_", t - 1, " + Z_", t - 1, ") * A_", t - 1))
                }
            })
            
            list(blip.mod = blip.mod, treat.mod = treat.mod, tf.mod = tf.mod)
        }
        DTRmodels <- generate_DTRreg_models(total_T)
        dta_wide <- dta %>%
            pivot_wider(
                names_from = dp,
                values_from = c("X", "Z", "avail", "A", "prob_A", "A_lag1")
            )
        beta <- tryCatch({
            # Attempt to run the model fitting code
            suppressMessages({
                mod <- DTRreg(dta_wide$Y, DTRmodels$blip.mod, DTRmodels$treat.mod, 
                              DTRmodels$tf.mod, data = dta_wide,
                              var.est = "sandwich", method = "gest")
            })
            
            # Extract the coefficients and compute beta_hat
            list(beta_hat = mean(unlist(coef(mod))))
            
        }, error = function(e) {
            # In case of any error, return beta_hat as NA
            list(beta_hat = NA)
        })
    } else if (method == "SNMM" & estimand == "moderated") {
        # for specifying the models to be passed into DTRreg()
        generate_DTRreg_models_moderated <- function(total_T) {
            # Create blip.mod: List of ~ 1 for each time step
            blip.mod <- lapply(1:total_T, function(t) as.formula(paste0("~ Z_", t)))
            
            # Create treat.mod: Dynamic formula list based on A_t ~ X_t + Z_t
            treat.mod <- lapply(1:total_T, function(t) as.formula(paste0("A_", t, " ~ X_", t, " + Z_", t)))
            
            # Create tf.mod: Dynamic formula list, interaction with A_t starting from the second time step
            tf.mod <- lapply(1:total_T, function(t) {
                if (t == 1) {
                    ~ 1
                } else {
                    as.formula(paste0("~ (X_", t - 1, " + Z_", t - 1, ") * A_", t - 1))
                }
            })
            
            list(blip.mod = blip.mod, treat.mod = treat.mod, tf.mod = tf.mod)
        }
        DTRmodels <- generate_DTRreg_models_moderated(total_T)
        dta_wide <- dta %>%
            pivot_wider(
                names_from = dp,
                values_from = c("X", "Z", "avail", "A", "prob_A", "A_lag1")
            )
        beta <- tryCatch({
            # Attempt to run the model fitting code
            suppressMessages({
                mod <- DTRreg(dta_wide$Y, DTRmodels$blip.mod, DTRmodels$treat.mod, 
                              DTRmodels$tf.mod, data = dta_wide,
                              var.est = "sandwich", method = "gest")
            })
            
            # Extract the coefficients and compute beta_hat
            list(beta_hat = rowMeans(sapply(coef(mod), rbind)))
            
        }, error = function(e) {
            # In case of any error, return beta_hat as c(NA,NA)
            list(beta_hat = c(NA,NA))
        })
    } else if (method == "WCLS" & estimand == "marginal") {
        fit <- wcls(data = dta,
                    id = "userid",
                    outcome = "Y",
                    treatment = "A",
                    rand_prob = "prob_A",
                    moderator_formula = ~1,
                    control_formula = ~ X + Z,
                    availability = "avail",
                    verbose = FALSE)
        fit_table <- summary(fit)$causal_excursion_effect
        beta_hat <- fit_table[1, "Estimate"]
        names(beta_hat) <- "Intercept"
        beta_se <- fit_table[1, "StdErr"]
        names(beta_se) <- "Intercept"
        conf_int <- matrix(c(fit_table[1, "95% LCL"], fit_table[1, "95% UCL"]), nrow = 1)
        rownames(conf_int) <- "Intercept"
        colnames(conf_int) <- c("2.5 %", "97.5 %")
        conf_int_tquantile <- conf_int
        beta <- list(beta_hat = beta_hat, 
                     beta_se = beta_se,
                     conf_int = conf_int,
                     conf_int_tquantile = conf_int_tquantile)
    } else if (method == "WCLS" & estimand == "moderated") {
        fit <- wcls(data = dta,
                    id = "userid",
                    outcome = "Y",
                    treatment = "A",
                    rand_prob = "prob_A",
                    moderator_formula = ~Z,
                    control_formula = ~ X + Z,
                    availability = "avail",
                    verbose = FALSE)
        fit_table <- summary(fit)$causal_excursion_effect
        beta_hat <- c(fit_table["(Intercept)", "Estimate"], fit_table["Z", "Estimate"])
        names(beta_hat) <- c("Intercept", "Z")
        beta_se <- c(fit_table["(Intercept)", "StdErr"], fit_table["Z", "StdErr"])
        names(beta_se) <- c("Intercept", "Z")
        conf_int <- fit_table[c("(Intercept)", "Z"), c("95% LCL", "95% UCL")]
        rownames(conf_int) <- c("Intercept", "Z")
        colnames(conf_int) <- c("2.5 %", "97.5 %")
        conf_int_tquantile <- conf_int
        beta <- list(beta_hat = beta_hat, 
                     beta_se = beta_se,
                     conf_int = conf_int,
                     conf_int_tquantile = conf_int_tquantile)
    } 
    result_collected <- c(result_collected, list(beta))
}

print("The beta object from the last iteration is:")
print(beta)

dir.create("result", showWarnings = FALSE)
saveRDS(result_collected, file = paste0("result/", itask, ".RDS"))


if (seed == (nsim / nsim_per_task)) {
    # wait till all is finished and collect results
    
    simulation_design_to_collect <- 
        simulation_design[itask, names(simulation_design) != "seed"]
    
    print("Collecting result files for:")
    print(simulation_design_to_collect)
    
    # Assume simulation_design_to_collect is a single-row data frame
    match_row <- simulation_design_to_collect[1, ]
    
    # Find the row indices in simulation_design where the first three columns match
    itask_to_collect <- which(
        simulation_design$estimand   == match_row$estimand &
            simulation_design$sample_size == match_row$sample_size &
            simulation_design$method      == match_row$method
    )
    
    all_result_files <- paste0("result/", itask_to_collect, ".RDS")
    
    while (!all(file.exists(all_result_files))) {
        Sys.sleep(5)
    }
    
    beta_collected <- list()
    for (ifile in all_result_files) {
        beta <- readRDS(ifile)
        beta_collected <- c(beta_collected, beta)
    }
    
    dir.create("collected_result", showWarnings = FALSE)
    saveRDS(beta_collected, file = paste0("collected_result/", 
                                          "estimand=", match_row$estimand,
                                          ",sample_size=", match_row$sample_size,
                                          ",method=", match_row$method,
                                          ".RDS"))
}


if (0) {
    # ----------------------------
    # Step 1. Determine missing file numbers
    # ----------------------------
    
    # Generate the expected file names: "1.RDS", "2.RDS", ..., "800.RDS"
    expected_files <- paste0(1:800, ".RDS")
    
    # List the actual files in the "result" folder
    actual_files <- list.files("result")
    
    # Find which expected files are missing
    missing_files <- setdiff(expected_files, actual_files)
    
    # Extract the numeric part from the missing file names
    missing_numbers <- as.numeric(sub("\\.RDS$", "", missing_files))
    
    # Optional: print missing file numbers
    cat("Missing file numbers (from result folder):\n")
    print(missing_numbers)
    
    
    # ----------------------------
    # Step 2. Compute the corresponding job numbers
    # ----------------------------
    
    start_job <- 35834389  # job number corresponding to 1.RDS
    missing_job_numbers <- start_job + missing_numbers - 1
    
    cat("Corresponding missing job numbers:\n")
    print(missing_job_numbers)
    
    
    # ----------------------------
    # Step 3. Create output_diag folder (if it doesn't already exist)
    # ----------------------------
    
    if (!dir.exists("output_diag")) {
        dir.create("output_diag")
    }
    
    
    # ----------------------------
    # Step 4. Copy the files from "output" to "output_diag"
    # ----------------------------
    
    # Loop over each missing job number and copy the corresponding files.
    for (job in missing_job_numbers) {
        # Construct the file names for the .err and .out files
        err_file <- paste0("distal.", job, ".err")
        out_file <- paste0("distal.", job, ".out")
        
        # Construct the full paths for source and destination files
        src_err <- file.path("output", err_file)
        src_out <- file.path("output", out_file)
        dest_err <- file.path("output_diag", err_file)
        dest_out <- file.path("output_diag", out_file)
        
        # Copy the .err file if it exists; if not, notify
        if (file.exists(src_err)) {
            file.copy(src_err, dest_err)
        } else {
            cat("File not found:", src_err, "\n")
        }
        
        # Copy the .out file if it exists; if not, notify
        if (file.exists(src_out)) {
            file.copy(src_out, dest_out)
        } else {
            cat("File not found:", src_out, "\n")
        }
    }
}


if (0) {
    # Define the folder names
    folders <- c("result", "collected_result", "output", "output_diag")
    
    # Loop over each folder
    for (folder in folders) {
        if (dir.exists(folder)) {
            # List all files in the folder with full paths
            files <- list.files(folder, full.names = TRUE)
            
            if (length(files) > 0) {
                # Remove all files in the folder
                removed <- file.remove(files)
                cat("Removed", sum(removed), "files from", folder, "\n")
            } else {
                cat("No files to remove in", folder, "\n")
            }
        } else {
            cat("Folder", folder, "does not exist.\n")
        }
    }
}