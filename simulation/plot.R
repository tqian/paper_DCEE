rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

# Safely extract a field from an element of a result list.
# If the field is missing (or the element is NA), return NA.
safe_extract <- function(l, field) {
    if (is.null(l) || (is.atomic(l) && length(l) == 1 && is.na(l))) {
        return(NA)
    }
    if (field %in% names(l)) {
        return(l[[field]])
    }
    return(NA)
}

compute_summary <- function(beta_hat, beta_se, conf_int, conf_int_tquantile, true_beta) {
    ## Convert beta_hat (a list of one-number vectors) into a numeric vector:
    bh <- sapply(beta_hat, function(x) as.numeric(x))
    
    ## (1) Bias: average estimate minus true value
    bias <- mean(bh, na.rm = TRUE) - true_beta
    
    ## (2) Standard deviation (empirical variability of beta_hat)
    emp_sd <- sd(bh, na.rm = TRUE)
    
    ## (3) RMSE: square-root of mean squared error
    rmse <- sqrt(mean((bh - true_beta)^2, na.rm = TRUE))
    
    ## (4) Mean estimated standard error (if available)
    if (!is.null(beta_se) && !all(sapply(beta_se, function(x) all(is.na(x))))) {
        bse <- sapply(beta_se, function(x) as.numeric(x))
        mean_se <- mean(bse, na.rm = TRUE)
    } else {
        bse <- NA
        mean_se <- NA
    }
    
    ## (5) Coverage based on sd:
    coverage_sd <- mean((true_beta >= (bh - 1.96 * emp_sd)) & (true_beta <= (bh + 1.96 * emp_sd)), na.rm = TRUE)
    
    ## (6) Coverage based on se:
    if (!is.null(beta_se) && !all(is.na(bse))) {
        coverage_se <- mean((true_beta >= (bh - 1.96 * bse)) & (true_beta <= (bh + 1.96 * bse)), na.rm = TRUE)
    } else {
        coverage_se <- NA
    }
    
    ## (7) Coverage based on conf_int:
    if (!is.null(conf_int) && !all(sapply(conf_int, function(x) all(is.na(x))))) {
        coverage_conf_int <- mean(sapply(conf_int, function(ci) {
            lower <- ci[1]
            upper <- ci[2]
            true_beta >= lower && true_beta <= upper
        }), na.rm = TRUE)
    } else {
        coverage_conf_int <- NA
    }
    
    ## (8) Coverage based on conf_int_tquantile:
    if (!is.null(conf_int_tquantile) && !all(sapply(conf_int_tquantile, function(x) all(is.na(x))))) {
        coverage_conf_int_tquantile <- mean(sapply(conf_int_tquantile, function(ci) {
            lower <- ci[1]
            upper <- ci[2]
            true_beta >= lower && true_beta <= upper
        }), na.rm = TRUE)
    } else {
        coverage_conf_int_tquantile <- NA
    }
    
    ## Return a list with all the summary measures:
    list(bias = bias,
         sd = emp_sd,
         rmse = rmse,
         mean_se = mean_se,
         coverage_sd = coverage_sd,
         coverage_se = coverage_se,
         coverage_conf_int = coverage_conf_int,
         coverage_conf_int_tquantile = coverage_conf_int_tquantile)
}

# Read in the true parameter values.
truth <- read_csv("truth.csv")
true_beta_marginal   <- truth$beta_marginal      # marginal
true_beta0_moderated <- truth$beta0_moderated    # intercept for moderated
true_beta1_moderated <- truth$beta1_moderated    # slope for moderated 

files <- list.files("collected_result", pattern = "\\.RDS$", full.names = TRUE)

# Create empty lists to store the summary results.
marginal_summary_list <- list()
beta0_summary_list <- list()
beta1_summary_list <- list()

for (file in files) {
    # Extract estimand, sample_size, method from the file name
    filename <- basename(file)
    filename_no_ext <- sub("\\.RDS$", "", filename)
    parts <- strsplit(filename_no_ext, ",")[[1]]
    estimand   <- sub("estimand=", "", parts[1])
    sample_size <- sub("sample_size=", "", parts[2])
    method     <- sub("method=", "", parts[3])
    
    cat("Processing file:", filename, "\n")
    cat("  estimand:", estimand, "\n")
    cat("  sample_size:", sample_size, "\n")
    cat("  method:", method, "\n")
    
    # --- Read and process the results in this file ---
    result <- readRDS(file)
    print(names(result[[1]]))
    
    beta_hat           <- lapply(result, function(l) safe_extract(l, "beta_hat"))
    beta_se            <- lapply(result, function(l) safe_extract(l, "beta_se"))
    conf_int           <- lapply(result, function(l) safe_extract(l, "conf_int"))
    conf_int_tquantile <- lapply(result, function(l) safe_extract(l, "conf_int_tquantile"))
    
    if (estimand == "marginal") {
        beta_summary <- compute_summary(beta_hat, beta_se, conf_int, conf_int_tquantile, true_beta_marginal)
        # Add metadata (estimand, sample_size, method)
        beta_summary <- c(estimand = estimand, sample_size = sample_size, method = method, beta_summary)
        marginal_summary_list[[length(marginal_summary_list) + 1]] <- beta_summary
    } else if (estimand == "moderated") {
        # For moderated, extract and summarize beta0 (intercept)
        beta0_hat <- lapply(beta_hat, function(x) x[1])
        beta0_se <- lapply(beta_se, function(x) x[1])
        if (!all(is.na(conf_int))) {
            beta0_conf_int <- lapply(conf_int, function(x) x[1, ])
            beta0_conf_int_tquantile <- lapply(conf_int_tquantile, function(x) x[1, ])
        } else { # all NA
            beta0_conf_int <- conf_int
            beta0_conf_int_tquantile <- conf_int_tquantile
        }
        beta0_summary <- compute_summary(beta0_hat, beta0_se, beta0_conf_int, beta0_conf_int_tquantile, true_beta0_moderated)
        beta0_summary <- c(estimand = estimand, sample_size = sample_size, method = method, param = "beta0", beta0_summary)
        beta0_summary_list[[length(beta0_summary_list) + 1]] <- beta0_summary
        
        # Now extract and summarize beta1 (slope)
        beta1_hat <- lapply(beta_hat, function(x) x[2])
        beta1_se <- lapply(beta_se, function(x) x[2])
        if (!all(is.na(conf_int))) {
            beta1_conf_int <- lapply(conf_int, function(x) x[2, ])
            beta1_conf_int_tquantile <- lapply(conf_int_tquantile, function(x) x[2, ])
        } else { # all NA
            beta1_conf_int <- conf_int
            beta1_conf_int_tquantile <- conf_int_tquantile
        }
        beta1_summary <- compute_summary(beta1_hat, beta1_se, beta1_conf_int, beta1_conf_int_tquantile, true_beta1_moderated)
        beta1_summary <- c(estimand = estimand, sample_size = sample_size, method = method, param = "beta1", beta1_summary)
        beta1_summary_list[[length(beta1_summary_list) + 1]] <- beta1_summary
    }
}

# Convert the summary lists into data frames.
df_marginal <- bind_rows(marginal_summary_list)
df_beta0    <- bind_rows(beta0_summary_list)
df_beta1    <- bind_rows(beta1_summary_list)

# Optionally print the data frames
print(df_marginal)
print(df_beta0)
print(df_beta1)

df_marginal <- df_marginal %>% 
    mutate(coverage = coalesce(coverage_conf_int_tquantile, coverage_se, coverage_sd))
df_beta0 <- df_beta0 %>% 
    mutate(coverage = coalesce(coverage_conf_int_tquantile, coverage_se, coverage_sd))
df_beta1 <- df_beta1 %>% 
    mutate(coverage = coalesce(coverage_conf_int_tquantile, coverage_se, coverage_sd))


# Update method names in all dataframes and ensure factor ordering
df_marginal <- df_marginal %>% 
    filter(!method %in% c("DCEE-lm", "DCEE-lm-cf")) %>% 
    mutate(method = case_when(
        method == "DCEE-gam" ~ "hat_beta",
        method == "DCEE-gam-cf" ~ "tilde_beta",
        method == "SNMM" ~ "SNMM",
        method == "GEE" ~ "GEE",
        TRUE ~ method
    )) %>%
    mutate(method = factor(method, levels = c("hat_beta", "tilde_beta", "SNMM", "GEE")))  # Ensure correct legend order

df_beta0 <- df_beta0 %>% 
    filter(!method %in% c("DCEE-lm", "DCEE-lm-cf")) %>% 
    mutate(method = case_when(
        method == "DCEE-gam" ~ "hat_beta",
        method == "DCEE-gam-cf" ~ "tilde_beta",
        method == "SNMM" ~ "SNMM",
        method == "GEE" ~ "GEE",
        TRUE ~ method
    )) %>%
    mutate(method = factor(method, levels = c("hat_beta", "tilde_beta", "SNMM", "GEE")))

df_beta1 <- df_beta1 %>% 
    filter(!method %in% c("DCEE-lm", "DCEE-lm-cf")) %>% 
    mutate(method = case_when(
        method == "DCEE-gam" ~ "hat_beta",
        method == "DCEE-gam-cf" ~ "tilde_beta",
        method == "SNMM" ~ "SNMM",
        method == "GEE" ~ "GEE",
        TRUE ~ method
    )) %>%
    mutate(method = factor(method, levels = c("hat_beta", "tilde_beta", "SNMM", "GEE")))

my_color_scale_manual <- scale_color_manual(
    values = c(
        "hat_beta" = "#0072B2",    # Blue
        "tilde_beta" = "#E69F00",  # Orange
        "SNMM" = "#009E73",        # Green
        "GEE" = "#CC79A7"          # Purple
    ),
    labels = c(
        "hat_beta" = TeX("$\\hat{\\beta}$"), 
        "tilde_beta" = TeX("$\\tilde{\\beta}$"),  
        "SNMM" = "SNMM", 
        "GEE" = "GEE"
    ),
    name = "Estimator"
)


library(ggplot2)
library(patchwork)
library(latex2exp)

# Bias plot
p_bias_marginal <- ggplot(df_marginal, aes(x = as.numeric(sample_size), y = as.numeric(bias), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_1^*$"), x = "n", y = "Bias") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_bias_beta0 <- ggplot(df_beta0, aes(x = as.numeric(sample_size), y = as.numeric(bias), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_2^*$"), x = "n", y = "Bias") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_bias_beta1 <- ggplot(df_beta1, aes(x = as.numeric(sample_size), y = as.numeric(bias), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_3^*$"), x = "n", y = "Bias") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

bias_plots <- (p_bias_marginal + p_bias_beta0 + p_bias_beta1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

# SD plots
p_sd_marginal <- ggplot(df_marginal, aes(x = as.numeric(sample_size), y = as.numeric(sd), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_1^*$"), x = "n", y = "SD") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_sd_beta0 <- ggplot(df_beta0, aes(x = as.numeric(sample_size), y = as.numeric(sd), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_2^*$"), x = "n", y = "SD") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_sd_beta1 <- ggplot(df_beta1, aes(x = as.numeric(sample_size), y = as.numeric(sd), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_3^*$"), x = "n", y = "SD") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

sd_plots <- (p_sd_marginal + p_sd_beta0 + p_sd_beta1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

# Coverage plots
p_cov_marginal <- ggplot(df_marginal, aes(x = as.numeric(sample_size), y = as.numeric(coverage), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_1^*$"), x = "n", y = "Coverage") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_cov_beta0 <- ggplot(df_beta0, aes(x = as.numeric(sample_size), y = as.numeric(coverage), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_2^*$"), x = "n", y = "Coverage") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

p_cov_beta1 <- ggplot(df_beta1, aes(x = as.numeric(sample_size), y = as.numeric(coverage), color = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(title = TeX("Estimand: $\\beta_3^*$"), x = "n", y = "Coverage") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_color_scale_manual

coverage_plots <- (p_cov_marginal + p_cov_beta0 + p_cov_beta1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

# To display the plots:
bias_plots
sd_plots
coverage_plots

dir.create("figure", showWarnings = FALSE)
path_to_paper <- 'figure'
ggsave(paste0(path_to_paper, "/simulation-bias.pdf"), bias_plots, width = 9, height = 2.5)
ggsave(paste0(path_to_paper, "/simulation-sd.pdf"), sd_plots, width = 9, height = 2.5)
ggsave(paste0(path_to_paper, "/simulation-coverage.pdf"), coverage_plots, width = 9, height = 2.5)