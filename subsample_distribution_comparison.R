library(dplyr)


acs = read.table("./RealDataAnalysis_CIFTI/acspsw03.txt", header = T, skip = 1)
acs$Subject.ID = stringr::str_remove(acs$Subject.ID.how.it.s.defined.in.lab.project, pattern = "_")
ABCD_with_task_dat <- read.csv("./RealDataAnalysis_CIFTI/nBack_gfactor/ABCD_task_general_gfactor.csv",  as.is=T)
acs_new = left_join(ABCD_with_task_dat, acs, by=c("Subject" = "Subject.ID"))
weight_acs <- acs_new$Imputed.raked.propensity.weight..The.raked.propensity.weight.merges.the.ACS.and.ABCD.data..with.missing.data.imputed...estimates.the.propensity.model..computes.and.scales.trims.the.propensity.weights.and.finally.rakes.the.scaled.weights.to.final.ACS.control.totals.by.age..sex.and.race.ethnicity.
weight_acs <- weight_acs/mean(weight_acs)
weight_final <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/weight_final.rds")
weight_final <- weight_final/mean(weight_final)

all_abcd_general = read.csv("./RealDataAnalysis_CIFTI/nBack/ABCD_task_general.csv", as.is = T)
all_abcd_general = all_abcd_general %>% select(Subject, site_id_l, Age, Gender, HighestParentalEducation, HouseholdIncome, HouseholdMaritalStatus,RaceEthnicity)
# all_abcd_general
tbss_score = read.csv("./RealDataAnalysis_CIFTI/ABCD_lavaan_G_time2.csv")
tbss_score$src_subject_id = stringr::str_remove(tbss_score$subjectkey, pattern = "_")

all_abcd_general$HouseholdIncome[all_abcd_general$HouseholdIncome == "NaN"] = sample(c("[>=50K & <100K]", "[>=100K]", "[<50K]"),
                                                                                     size = sum(all_abcd_general$HouseholdIncome=="NaN"),
                                                                                     replace = T)

all_abcd_general$HighestParentalEducation[all_abcd_general$HighestParentalEducation == "NaN"] = sample(c("HS Diploma/GED",
                                                                                                         "< HS Diploma",
                                                                                                         "Post Graduate Degree",
                                                                                                         "Some College" ,
                                                                                                         "Bachelor"),
                                                                                                       size = sum(all_abcd_general$HighestParentalEducation == "NaN"),
                                                                                                       replace = T)

all_abcd_general$HouseholdMaritalStatus[all_abcd_general$HouseholdMaritalStatus == "NaN"] = sample(c("yes", "no"),
                                                                                                   size = sum(all_abcd_general$HouseholdMaritalStatus=="NaN"),
                                                                                                   replace = T)

all_abcd_general_control <- all_abcd_general %>% select(c(Age, HouseholdIncome, HouseholdMaritalStatus, HighestParentalEducation, Gender, site_id_l, RaceEthnicity))
all_abcd_general_control <- fastDummies::dummy_columns(all_abcd_general_control, c(
  "RaceEthnicity",
  "HighestParentalEducation",
  "HouseholdMaritalStatus",
  "HouseholdIncome",
  "site_id_l",
  "Gender"), 
  remove_selected_columns = T)
all_abcd_general_control$Age10 = as.integer(all_abcd_general_control$Age>=10)
all_abcd_general_control <- all_abcd_general_control %>% select(-Age)
all_abcd_general_control_weight = apply(all_abcd_general_control,2, function(x) x*weight_acs)
# Summarize the demographic characteristics
socio_demographics_summary_abcd_unweighted <- data.frame(
  Variables = colnames(all_abcd_general_control),
  Mean = apply(all_abcd_general_control,2, mean, na.rm = TRUE)
)

socio_demographics_summary_abcd_weighted <- data.frame(
  Variables = colnames(all_abcd_general_control_weight),
  Mean = apply(all_abcd_general_control_weight,2, mean, na.rm = TRUE)
)

Z <- read.csv("./RealDataAnalysis_CIFTI/nBack_gfactor/ABCD_task_general_gfactor.csv",  as.is=T)
Z <- Z %>% select(c(Age, HouseholdIncome, HouseholdMaritalStatus, HighestParentalEducation, Gender, site_id_l, RaceEthnicity))

Z <- fastDummies::dummy_columns(Z, c(
  "RaceEthnicity",
  "HighestParentalEducation",
  "HouseholdMaritalStatus",
  "HouseholdIncome",
  "site_id_l",
  "Gender"), 
  remove_selected_columns = T)

Z$Age10 = as.integer(Z$Age >= 10)
Z <- Z %>% select(-Age)

# Summarize the demographic characteristics
socio_demographics_summary_unweighted <- data.frame(
  Variables = colnames(Z),
  Mean = apply(Z,2, mean, na.rm = TRUE)
)

Z_weighted <- apply(Z,2, function(x) x*weight_final)
# Summarize the demographic characteristics
socio_demographics_summary_weighted <- data.frame(
  Variables = colnames(Z_weighted),
  Mean = apply(Z_weighted,2, mean, na.rm = TRUE)
)

library(dplyr)
socio_demographics_summary_all = socio_demographics_summary_abcd_unweighted %>%
  left_join(socio_demographics_summary_abcd_weighted, by = "Variables", suffix = c(".ABCD", ".ABCD_unweighted")) %>%
  left_join(socio_demographics_summary_unweighted, by = "Variables", suffix = c(".ABCD", ".Subsample_unweighted")) %>%
  left_join(socio_demographics_summary_weighted, by = "Variables", suffix = c(".Subsample_unweighted", ".Subsample_weighted"))
  
colnames(socio_demographics_summary_all) <- c("Variables", "ABCD_unweighted", "ABCD_weighted", "Subsample_unweighted", "Subsample_weighted")
socio_demographics_summary_all[,2:ncol(socio_demographics_summary_all)]  <- 100*socio_demographics_summary_all[,2:ncol(socio_demographics_summary_all)]
print(socio_demographics_summary_all,digits=1)



# # Detail look in tot the income category weights
# income_cat <- c("[<50K]", "[>=50K & <100K]", "[>=100K]")
# weights_income_50K_acs <- weight_acs[all_abcd_general$HouseholdIncome == income_cat[1]]/mean(weight_acs)
# weights_income_50K_100K_acs <- weight_acs[all_abcd_general$HouseholdIncome == income_cat[2]]/mean(weight_acs)
# weights_income_100K_acs <- weight_acs[all_abcd_general$HouseholdIncome == income_cat[3]]/mean(weight_acs)
# 
# 
# weights_income_50K_final <- weight_final[all_abcd_general$HouseholdIncome == income_cat[1]]/mean(weight_final)
# weights_income_50K_100K_final <- weight_final[all_abcd_general$HouseholdIncome == income_cat[2]]/mean(weight_final)
# weights_income_100K_final <- weight_final[all_abcd_general$HouseholdIncome == income_cat[3]]/mean(weight_final)
# 
# Z_income_50K <- Z$`HouseholdIncome_[<50K]`
# Z_income_50K_100K <- Z$`HouseholdIncome_[>=50K & <100K]`
# Z_income_100K <- Z$`HouseholdIncome_[>=100K]`
# 
# Z_weighted_income_50K_acs <- Z_income_50K * weight_acs
# Z_weighted_income_50K_100K_acs <- Z_income_50K_100K * weight_acs
# Z_weighted_income_100K_acs <- Z_income_100K * weight_acs
# 
# Z_weighted_income_50K_final <- Z_income_50K * weight_final
# Z_weighted_income_50K_100K_final <- Z_income_50K_100K * weight_final
# Z_weighted_income_100K_final <- Z_income_100K * weight_final


