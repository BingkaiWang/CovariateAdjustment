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
substance <- read.table("subuq01.txt", header = T, stringsAsFactors = F, skip = 2)
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