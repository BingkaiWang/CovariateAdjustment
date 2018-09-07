# Clean the raw data of trial METS, which is stored in folder NCT00816907
# The main publication of this trial is 
# Jarskog, L. F., R. M. Hamer, D. J. Catellier, D. D. Stewart, L. LaVange, N. Ray, 
# L. H.Golden, J. A. Lieberman, and T. S. Stroup (2013). Metformin for weight loss 
# andmetabolic control in overweight outpatients with schizophrenia and schizoaffective 
# disorder.American Journal of Psychiatry  170(9), 1032–1040.
# 
# set the current directory to be this folder and run the following code
# output a cleaned data set with variable:
#   subject_id: id of each subject
#   age:  age of each subject in month,
#   gender: an indicator with 1 indicating female, 0 indicating male,
#   treatment: "Metformin" or "Placebo"
#   CGI: CGI (Clinical Global Impressions severity score) score measured at baseline
#   weight_baseline: weight (pound) measured at baseline
#   weight_16: weight measured at week 16
#   weightchange: weight_16 - weight_baseline
#   height: height (inch) measured at baseline
#   tobacco: level of using tobacoo
#   alcohol: level of using alcohol
#   drug: level of using drug
#   bmi: BMI index calculated using weight_baseline and height


library(dplyr)
setwd("../CovariateAdjustment/NCT00816907/")

# extract treatment assignment
METS <- read.table("demo01.txt", header = F, stringsAsFactors = F, skip = 2)
METS <- data.frame(subject_id = METS[,4], treatment = as.character(METS[,84])) %>%
  filter(treatment != "")

# extract age, gender and CGI
cgi <- read.table("cgi01.txt", header = T, stringsAsFactors = F)
cgi <- select(cgi[-1, ], subject_id = src_subject_id, age = interview_age,
              gender = gender, CGI = cgi_si, visit = visit) %>%
  filter(visit == "Baseline") %>%
  select(1:4)
for(i in c(1,2,4)) {cgi[,i] <- as.numeric(cgi[,i])}
METS <- inner_join(METS, cgi)

# extract weight_baseline, weight_16, height
METS <- mutate(METS, weight_baseline = NA, weight_16 = NA, height = NA)
weight <- read.table("vitals01.txt", header = T, stringsAsFactors = F)
weight <- select(weight[-1, ], subject_id = src_subject_id, visit = visit,
                 weight = weight_std, height = height_std)
for(i in 1:nrow(METS)){
  temp <- weight %>% filter(subject_id == METS$subject_id[i])
  METS$weight_baseline[i] <- temp$weight[temp$visit == "Baseline"]
  METS$weight_16[i] <- ifelse(sum(temp$visit == "Week 16"), temp$weight[temp$visit == "Week 16"], NA)
  METS$height[i] <- temp$height[temp$visit == "Screening"][1]
}

# extract substance use (tobacco, alcohol, drug)
METS <- mutate(METS, tobacco = NA, alcohol = NA, drug = NA)
substance <- read.table("subuq01.txt", header = F, stringsAsFactors = F, skip = 2)
substance <- data.frame(subject_id = substance[,4], tobacco = substance[,400], 
                        alcohol = substance[,405], drug = substance[,404],
                        visit = substance[, 408]) %>% unique
for(i in 1:nrow(METS)){
  temp <- substance %>% filter(subject_id == METS$subject_id[i]) %>%
    filter(visit == "Baseline")
  METS$tobacco[i] <- temp$tobacco[!is.na(temp$tobacco)][1]
  METS$alcohol[i] <- temp$alcohol[!is.na(temp$alcohol)][1]
  METS$drug[i] <- temp$drug[!is.na(temp$drug)][1]
}

# cleaning
METS$treatment <- as.character(METS$treatment)
METS$weight_baseline <- as.numeric(METS$weight_baseline)
METS$weight_16 <- as.numeric(METS$weight_16)
METS$height <- as.numeric(METS$height)
METS$gender <- METS$gender == "F"
METS <- mutate(METS, weightchange = weight_16 - weight_baseline,
               bmi = (weight_baseline*0.453592)/(height*0.0254)^2)

#save cleaned data set as METS.rdata
save(METS, file = "../Data_Preprocessing_and_Analysis/METS.rdata")
setwd("..")
