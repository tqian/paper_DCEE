# Tianchen Qian
# 2024.03.02

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

# publicly available HeartStepsV1 data downloaded from https://github.com/klasnja/HeartStepsV1

jbslot <- read.csv("HeartStepsV1-main/data_files/jbsteps.csv")
gfslot <- read.csv("HeartStepsV1-main/data_files/gfsteps.csv")
users <- read.csv("HeartStepsV1-main/data_files/users.csv")
suggest <- read.csv("HeartStepsV1-main/data_files/suggestions.csv")

CREATE_PLOTS <- FALSE # if TRUE, will make sanity check plots
PRINT_SANITY_CHECK <- TRUE # if TRUE, will print sanity check results

jbslot <- as_tibble(jbslot)

suggest <- as_tibble(suggest)

# remove the baseline covariates (decision.index = NA) -------------------------------------------
jbslot <- jbslot %>% filter(!is.na(decision.index))
# the publicly available dataset has baseline information recorded (i.e. at decision.index = NA)
# we will remove these to keep the form of the datasets to be consistent with our analysis.


# Remove travel dates in suggest
var_from_suggest <- c("user.index", "decision.index.nogap",  "decision.index", "avail", "send",
                      "send.active", "send.sedentary", "sugg.decision.utime", 
                      "dec.location.category", "jbsteps30pre", "jbsteps30pre.zero",
                      "jbsteps30")
suggest_selectedvar <- suggest[, var_from_suggest]

suggest_selectedvar <- suggest_selectedvar %>% 
    mutate(sugg.decision.udate = date(ymd_hms(suggest_selectedvar$sugg.decision.utime)))
tmp_df <- data.frame()

for (i in 1:nrow(users)) {
    
    user_i_dta <- filter(suggest_selectedvar, user.index == users$user.index[i])
    if (users$travel.start[i] != "") {
        nrow_before <- nrow(user_i_dta)
        user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
        nrow_after <- nrow(user_i_dta)
        cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
        cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
                   nrow_before - nrow_after, "\n"))
    }
    tmp_df <- rbind(tmp_df, user_i_dta)
}
suggest_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 278, nrow_after 178, nrow_removed 100
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 255, nrow_after 215, nrow_removed 40
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 212, nrow_after 182, nrow_removed 30
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 215, nrow_after 192, nrow_removed 23
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 265, nrow_after 210, nrow_removed 55
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 214, nrow_after 199, nrow_removed 15
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 309, nrow_after 184, nrow_removed 125

# Construct a new decision.index.nogap variable in suggest dataset, then paste it into jbslot -------------------------------------------
## Issue: Why do we need to do this? b/c somehow there are observations that have decision.index.nogap skipped
## e.g. user 13's orginal decision.index.nogap has labeled wrong at j = 23, 24 (they got skipped)
## Solution: construct a new decision.index.nogap variable

suggest_tdr <- suggest_tdr %>% 
    group_by(user.index) %>% 
    mutate(max.decision.index.nogap = n(),
           decision.index.nogap.new = 0:(unique(max.decision.index.nogap) - 1)) %>% 
    ungroup()
a <- which(suggest_tdr$decision.index.nogap != suggest_tdr$decision.index.nogap.new)
b <- unique(suggest_tdr[a, ]) 

# Paste the newly constructed decision.index.nogap.new into jbslot
jbslot <- suggest_tdr %>% 
    dplyr::select("user.index", "decision.index", "decision.index.nogap.new",
                  "sugg.decision.utime", "sugg.decision.udate") %>%
    right_join(jbslot, by = c("user.index", "decision.index"))


# remove travel dates in jbslot-------------------------------------------

## Issue: At i = 1, j = 203, the suggest dataset has notification sent when they are in travel dates. 
## But the steps are recorded. And the steps recorded time is not in travel dates (the morning they got back). 
## So there are NAs when joining these two datasets. We can remove these observations since their notifications are sent at travel dates.
## Solution: we use the sugg.decision.udate from suggest dataset (to filter the observations in travel dates) to avoid discrepancy between suggest and jbslot.

# travel day is stored in "users" data.frame, with variable names "travel.start" and "travel.end"
tmp_df <- data.frame()
for (i in 1:nrow(users)) {
    
    user_i_dta <- filter(jbslot, user.index == users$user.index[i])
    if (users$travel.start[i] != "") {
        nrow_before <- nrow(user_i_dta)
        user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
        nrow_after <- nrow(user_i_dta)
        cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
        cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
                   nrow_before - nrow_after, "\n"))
    }
    tmp_df <- rbind(tmp_df, user_i_dta)
}
jbslot_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 7965, nrow_after 7542, nrow_removed 423
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 4611, nrow_after 4441, nrow_removed 170
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 7364, nrow_after 7094, nrow_removed 270
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 3389, nrow_after 2866, nrow_removed 523
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 9199, nrow_after 7828, nrow_removed 1371
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 7041, nrow_after 6496, nrow_removed 545
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 5357, nrow_after 5108, nrow_removed 249


##### Check the date range of observations in both data sets

suggest_tdr[suggest_tdr$sugg.decision.utime == "", ]
# There seems to be an erroneous observation being the 0-th decision point for user 4.
# Removing it.

suggest_tdr <- suggest_tdr[!(suggest_tdr$sugg.decision.utime == ""), ]

suggest_tdr <- suggest_tdr %>% 
    mutate(sugg.decision.utime = ymd_hms(sugg.decision.utime))

jbslot_tdr <- jbslot_tdr %>% 
    mutate(steps.utime = ymd_hms(steps.utime),
           steps.utime.local = ymd_hms(steps.utime.local),
           sugg.decision.utime = ymd_hms(sugg.decision.utime))

suggest_timerange <- suggest_tdr %>%
    group_by(user.index) %>% 
    summarize(suggest_earliest = min(sugg.decision.utime),
              suggest_latest = max(sugg.decision.utime))

jbslot_timerange <- jbslot_tdr %>%
    group_by(user.index) %>% 
    summarize(jbslot_earliest = min(steps.utime),
              jbslot_latest = max(steps.utime))

# View(full_join(suggest_timerange, jbslot_timerange) %>% 
#          relocate(jbslot_earliest, .after = suggest_earliest))

# full_join(suggest_timerange, jbslot_timerange) %>% 
#     mutate(start_diff = jbslot_earliest - suggest_earliest,
#            end_diff = jbslot_latest - suggest_latest) %>% 
#     View()

##### Remove out of range data

## Remove jbslot that is earlier than the first decision point
## Remove jbslot that is later than 24 hours following the last decision point

jbslot_tdr <- jbslot_tdr %>% 
    full_join(suggest_timerange)
nrow(jbslot_tdr) # 233542

jbslot_tdr <- jbslot_tdr %>% 
    filter(steps.utime >= suggest_earliest) %>% 
    filter(steps.utime <= (suggest_latest + 3600 * 24))
nrow(jbslot_tdr) # 231691 (removed  0.8% observations from jbslot_tdr)


dta_treatment_and_proximal <- suggest_tdr %>% 
    mutate(rand_prob = 0.6) %>% 
    filter(decision.index.nogap <= 210) %>% 
    mutate(avail = as.numeric(avail == "True"),
           send = as.numeric(send == "True"),
           send.active = as.numeric(send.active == "True"),
           send.sedentary = as.numeric(send.sedentary == "True"),
           home_work = as.numeric(dec.location.category %in% c("home", "work")))

# impute missing jbsteps30 with 0
dta_treatment_and_proximal <- dta_treatment_and_proximal %>% 
    mutate(jbsteps30 = ifelse(is.na(jbsteps30), 0, jbsteps30))

summary(dta_treatment_and_proximal)

dta_treatment_and_proximal <- dta_treatment_and_proximal %>% 
    relocate(rand_prob, .after = send)

## add another intermediate outcome: 
# jbsteps between current decision point and next decision point

dta_jbsteps_between_decision_points <- jbslot_tdr %>% 
    group_by(user.index, decision.index.nogap.new) %>% 
    summarize(jbsteps_between_decision_points = sum(steps)) %>%
    ungroup()

dta_treatment_and_proximal <- dta_treatment_and_proximal %>% 
    full_join(dta_jbsteps_between_decision_points)

saveRDS(dta_treatment_and_proximal, "HeartSteps_dta_treatment_and_proximal.RDS")
write_csv(dta_treatment_and_proximal, "HeartSteps_dta_treatment_and_proximal.csv")
