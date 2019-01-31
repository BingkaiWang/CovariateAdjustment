library(dplyr)
TADS_result <- read.csv("DataResults/TADS_result.csv")
MCI_result <- read.csv("DataResults/MCI_result.csv")
METS_result <- read.csv("DataResults/METS_result.csv")

summarytable <- rbind(TADS_result, MCI_result, METS_result) %>% select(2, 3, 5, 6, 8, 9)
for(i in c(1,2,4)){ summarytable[,i] <- as.character(summarytable[,i])}
summarytable$treatment_arm[1:3] <- paste0("TADS(", summarytable$treatment_arm[1:3], ")")
summarytable <- data.frame(trial_name = summarytable$treatment_arm,
                           unadj = paste0(summarytable$unadj_est, summarytable$unadj_c.i.),
                           ancova = paste0(summarytable$ancova_est, summarytable$ancova_c.i.),
                           variance_reduction = summarytable$variance_reduction
                           )
summarytable
write.csv(summarytable, file = "DataResults/summarytable.csv")

library(xtable)
xtable(summarytable)
