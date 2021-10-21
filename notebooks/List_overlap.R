##############################################
# Create merged list of results for Alexey 
###############################################

# Bring in HSMV normal vs HSMV ischemic and SMV normal vs. SMV ischemic 
library(readxl)
HSMVisch_vs_HSMVnormal_results_filtered <- read_excel("~/Desktop/Alexey_results report/excel_files/HSMVisch_vs_HSMVnormal_results_filtered.xlsx", sheet=2)
SMVisch_vs_SMVnormal_results_filtered <- read_excel("~/Desktop/Alexey_results report/excel_files/SMVisch_vs_SMVnormal_results_filtered.xlsx", sheet=2)


# Drop columns in data set that are not needed 

# Drop samples 7, 8, 10, 11, 12, 14 
names(HSMVisch_vs_HSMVnormal_results_filtered)
#HSMVisch_vs_HSMVnormal_results_filtered <- HSMVisch_vs_HSMVnormal_results_filtered[, -c(20:31, 52:57)]
# names(HSMVisch_vs_HSMVnormal_results_filtered)
names(SMVisch_vs_SMVnormal_results_filtered)
#SMVisch_vs_SMVnormal_results_filtered <- SMVisch_vs_SMVnormal_results_filtered[,-c(20:31, 52:57)]
# names(SMVisch_vs_SMVnormal_results_filtered)

# Join same loci 
library(dplyr)
shared_loci <- inner_join(HSMVisch_vs_HSMVnormal_results_filtered, SMVisch_vs_SMVnormal_results_filtered, by="loci", suffix = c(".HSMV", ".SMV"))
head(shared_loci)


# Effect in HSMV but not in SMV
library(dplyr)
HSMV_only <- anti_join(HSMVisch_vs_HSMVnormal_results_filtered, SMVisch_vs_SMVnormal_results_filtered, by="loci")


# Effect in SMV but not in HSMV 
library(dplyr)
SMV_only <- anti_join(SMVisch_vs_SMVnormal_results_filtered, HSMVisch_vs_HSMVnormal_results_filtered, by="loci")


# Write excel files
getwd()
write.xlsx(shared_loci, file = 'both_HSMV_and_SMV.xlsx')
write.xlsx(HSMV_only, file = 'HSMV_effect_only.xlsx')
write.xlsx(SMV_only, file = 'SMV_effect_only.xlsx')
