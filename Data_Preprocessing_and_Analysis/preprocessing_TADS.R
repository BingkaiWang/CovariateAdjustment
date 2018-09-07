# Clean the raw data of trial TADS, which is stored in folder NCT00006286N01MH80008AnneMarie
# The main publication of this trial is 
# Treatment for Adolescents WithDepression Study(TADS) Team(2004).Fluoxetine,cognitive-behavioral therapy, 
# and their combination for adolescents with depres-sion: Treatment for adolescents with depression study
# (tads) randomized con-trolled trial.JAMA 292(7), 807–820.
# 
# set the current directory to be this folder and run the following code
# output a cleaned data set with variable:
#   subject_id: id of each subject
#   age:  age of each subject in month,
#   gender: an indicator with 1 indicating female, 0 indicating male,
#   treatment: "FLX" (fluoxetine), "PBO" (placebo), "CBT" (cognitive-behavioral therapy ) or "COMB" (a combination of FLX and CBT),
#   CDRS_baseline: CDRS-R (Children’s Depression Rating Scale-Revised) score measured at baseline
#   CDRS_12: CDRS-R (Children’s Depression Rating Scale-Revised) score measured at week 12
#   CGI: CGI (Clinical Global Impressions severity score) score measured at baseline
#   CGAS: CGAS (Children’s Global Assessment Scale score) score measured at baseline
#   RADS: Reynolds Adolescent Depression Scale measured at baseline
#   suicide_ideation: Suicidal Ideation Questionnaire-Junior High School Version total score meassured at baseline
#   depression_episode: Current major depressive episode duration measured at baseline
#   comorbidity: indicator of whether co-morbidity exists
#   change_score: CDRS_12 - CDRS_baseline
# 
# Prepocessing method:
# 1. for getting measures at baseline, we use the non-missing value with visit week closest to 0 from visits within 2 weeks of randomization
# 2. for getting measures at week 12, we use the non-missing value with visit week closest to 12 from visits within 2 weeks of randomization
# 3. For missing value in baseline variables, median imputation is used.

library(dplyr)
setwd("../CovariateAdjustment/NCT00006286N01MH80008AnneMarie/")

# Extract age, gender and treatment assigment
tad <- read.table("mastr01.txt", header = T, stringsAsFactors = F)
tad <- select(tad[-1, ], subject_id = src_subject_id, 
              age = interview_age, 
              gender = gender, 
              treatment = txcode) %>% unique

# Extract CDRS-R score (baseline and 12 week score)
tad <- mutate(tad, CDRS_baseline = NA, CDRS_12 = NA)
cdrs <- read.table("cdrsr01.txt", header = T, stringsAsFactors = F)
cdrs <- select(cdrs[-1, ], subject_id = src_subject_id, cdrs_r_b, week = bsit0)
for(i in 2:3) {cdrs[,i] <- as.numeric(cdrs[,i])}
for(i in 1:nrow(tad)){
  temp <- cdrs %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(week < 2) %>% 
    arrange(abs(week))
  tad$CDRS_baseline[i] <- temp$cdrs_r_b[!is.na(temp$cdrs_r_b)][1]
  temp <- cdrs %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(abs(week-12) < 2) %>% 
    arrange(abs(week-12))
  tad$CDRS_12[i] <- temp$cdrs_r_b[!is.na(temp$cdrs_r_b)][1]
}

# Extract CGI score and CGAS score (baseline only) and CGI_12
tad <- mutate(tad, CGI = NA, CGAS = NA, CGI_improvement = NA)
cgi <- read.table("cgi01.txt",  header = T, stringsAsFactors = F)
cgi <- select(cgi[-1, ], subject_id = src_subject_id, CGI = cgi_si, CGAS = cgasrat, 
              CGI_improvement = cgi_sii, week = bsit0)
for(i in 2:5) {cgi[,i] <- as.numeric(cgi[,i])}
for(i in 1:nrow(tad)){
  temp <- cgi %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(week < 2) %>% 
    arrange(abs(week))
  tad$CGI[i] <- temp$CGI[!is.na(temp$CGI)][1]
  tad$CGAS[i] <- temp$CGAS[!is.na(temp$CGAS)][1]
  temp <- cgi %>%
    filter(subject_id == tad$subject_id[i]) %>%
    filter(abs(week - 12) < 2) %>%
    arrange(abs(week - 12))
  tad$CGI_improvement[i] <- temp$CGI_improvement[!is.na(temp$CGI_improvement)][1]
}

# Extract RADS score (baseline only)
tad <- mutate(tad, RADS = NA)
rads <- read.table("rads01.txt",  header = T, stringsAsFactors = F)
rads <- select(rads[-1, ], subject_id = src_subject_id, RADS = rads_scr, week = bsit0)
for(i in 2:3) {rads[,i] <- as.numeric(rads[,i])}
for(i in 1:nrow(tad)){
  temp <- rads %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(week < 2) %>% 
    arrange(abs(week))
  tad$RADS[i] <- temp$RADS[!is.na(temp$RADS)][1]
}

# Extract suicide ideation score (baseline only)
tad <- mutate(tad, suicide_ideation = NA)
suicide <- read.table("siq01.txt",  header = T, stringsAsFactors = F)
suicide <- select(suicide[-1, ], subject_id = src_subject_id, suicide_ideation = siqrisk, week = bsit0)
for(i in 2:3) {suicide[,i] <- as.numeric(suicide[,i])}
for(i in 1:nrow(tad)){
  temp <- suicide %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(week < 2) %>% 
    arrange(abs(week))
  tad$suicide_ideation[i] <- temp$suicide_ideation[!is.na(temp$suicide_ideation)][1]
}

# Extract major depression episode and comorbidity (baseline)
tad <- mutate(tad, depression_episode = NA, comorbidity = NA)
episode <- read.table("ksads_diagnoses01.txt", header = F, stringsAsFactors = F, skip = 2)
episode <- data.frame(subject_id = episode[,4], depression_episode = episode[,82], 
                      comorbidity = episode[,388], week = episode[,59])
for(i in 2:4) {episode[,i] <- as.numeric(episode[,i])}
for(i in 1:nrow(tad)){
  temp <- episode %>% 
    filter(subject_id == tad$subject_id[i]) %>% 
    filter(week < 2) %>% 
    arrange(abs(week))
  tad$depression_episode[i] <- temp$depression_episode[!is.na(temp$depression_episode)][1]
  tad$comorbidity[i] <- temp$comorbidity[!is.na(temp$comorbidity)][1]
}

# Cleaning and median (mode) imputation for baseline variables
tad$gender <- tad$gender == "F"
tad$age <- as.numeric(tad$age)
tad <- mutate(tad, change_score = CDRS_12 - CDRS_baseline, 
              binary_CGI_improvement = (CGI_improvement <= 2))
tad$CDRS_baseline[which(is.na(tad$CDRS_baseline))] <- median(tad$CDRS_baseline, na.rm = T)
tad$CGI[which(is.na(tad$CGI))] <- median(tad$CGI, na.rm = T)
tad$CGAS[which(is.na(tad$CGAS))] <- median(tad$CGAS, na.rm = T)
tad$RADS[which(is.na(tad$RADS))] <- median(tad$RADS, na.rm = T)
tad$suicide_ideation[which(is.na(tad$suicide_ideation))] <- median(tad$suicide_ideation, na.rm = T)
tad$comorbidity[which(is.na(tad$comorbidity))] <- as.numeric(names(table(tad$comorbidity))[1])

#save cleaned data set as TADS.rdata
save(tad, file = "../Data_Preprocessing_and_Analysis/TADS.rdata")
setwd("..")
