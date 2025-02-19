library(tidyverse)

source("dgm_distal_outcome_endox_nonlinear_interaction_binaryZ.R")

large_sample_size_number <- 1000000

compute_true_marginal_CEE <- function(
        total_T,
        rand_prob_pattern = c("constant", "timevar"),
        availability = c("always", "not-always")
) {
    taut <- numeric(total_T)
    for (t in 1:total_T) {
        cat(t, "")
        dta1 <- dgm_distal_outcome_endox_nonlinear_interaction_binaryZ(
            sample_size = large_sample_size_number, total_T = total_T,
            rand_prob_pattern = rand_prob_pattern,
            availability = availability,
            force_A_dp = t, force_A_value = 1)
        
        dta0 <- dgm_distal_outcome_endox_nonlinear_interaction_binaryZ(
            sample_size = large_sample_size_number, total_T = total_T,
            rand_prob_pattern = rand_prob_pattern,
            availability = availability,
            force_A_dp = t, force_A_value = 0)
        
        taut[t] <- mean(dta1$Y) - mean(dta0$Y)
    }
    
    return(mean(taut))
}


compute_true_CEE_moderatorZ <- function(
        total_T,
        rand_prob_pattern = c("constant", "timevar"),
        availability = c("always", "not-always")
) {
    taut_Z1 <- taut_Z0 <- numeric(total_T)
    for (t in 1:total_T) {
        cat(t, "")
        dta1 <- dgm_distal_outcome_endox_nonlinear_interaction_binaryZ(
            sample_size = large_sample_size_number, total_T = total_T,
            rand_prob_pattern = rand_prob_pattern,
            availability = availability,
            force_A_dp = t, force_A_value = 1)
        
        dta0 <- dgm_distal_outcome_endox_nonlinear_interaction_binaryZ(
            sample_size = large_sample_size_number, total_T = total_T,
            rand_prob_pattern = rand_prob_pattern,
            availability = availability,
            force_A_dp = t, force_A_value = 0)
        
        taut_Z1[t] <- mean(dta1$Y[dta1$dp == t & dta1$Z == 1]) - mean(dta0$Y[dta0$dp == t & dta0$Z == 1])
        taut_Z0[t] <- mean(dta1$Y[dta1$dp == t & dta1$Z == 0]) - mean(dta0$Y[dta0$dp == t & dta0$Z == 0])
    }
    
    CEE <- c(mean(taut_Z0), mean(taut_Z1) - mean(taut_Z0))
    names(CEE) <- c("beta0", "beta1")
    return(CEE)
}

set.seed(123)
truth_marginal <- compute_true_marginal_CEE(
    total_T = 30,
    rand_prob_pattern = "timevar",
    availability = "not-always"
)

print(paste0("truth_marginal: ", truth_marginal))

set.seed(123)
truth_moderated <- compute_true_CEE_moderatorZ(
    total_T = 30,
    rand_prob_pattern = "timevar",
    availability = "not-always"
)

print(paste0("truth_moderated: ", truth_moderated))


# Combine the numbers into a data frame.
beta_df <- data.frame(
    beta_marginal = truth_marginal,
    beta0_moderated = truth_moderated[1],
    beta1_moderated = truth_moderated[2]
)

# Write the data frame to a CSV file.
write_csv(beta_df, file = "truth.csv")
