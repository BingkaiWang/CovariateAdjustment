TADS_result <- read.csv("DataResults/TADS_result.csv")
MCI_result <- read.csv("DataResults/MCI_result.csv")
METS_result <- read.csv("DataResults/METS_result.csv")

summarytable <- rbind(TADS_result, MCI_result, METS_result)
for(i in c(2,6,9)){ summarytable[,i] <- as.character(summarytable[,i])}
summarytable$treatment_arm[1:3] <- paste0("TADS(", summarytable$treatment_arm[1:3], ")")
for(i in c(3,4,7,10:12)){ summarytable[,i] <- round(summarytable[,i], 2)}
summarytable <- data.frame(trial_name = summarytable$treatment_arm,
                           aggregate_imbalance = summarytable$aggregate_imbalance,
                           unadjusted_estimator = paste0(summarytable$unadjusted_estimator, summarytable$unadj_c.i.),
                           adjusted_estimator = paste0(summarytable$adjusted_estimator, summarytable$adj_c.i.),
                           conditional_bias = summarytable$condi_bias,
                           variance_reduction = summarytable$variance_reduction,
                           RE = summarytable$RE
                           )
write.csv(summarytable, file = "DataResults/summarytable.csv")
