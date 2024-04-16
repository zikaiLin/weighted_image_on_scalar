# # Script Description:



library(dplyr)
library(mice)
library(survey)
library(psych)

# --- Load the data --- #
abcd <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/abcd_bl_preprocessed.rds")
# > names(abcd)
# [1] "subjectkey"                   "family_size3"                 "marital_status"              
# [4] "household_income"             "age10"                        "race_ethnicity"              
# [7] "gender"



# --- household size --- #
abcd$household_income <- factor(abcd$household_income,
                                levels = c("Under $50,000", "$50,000 to $100,000", "Over $100,000"),
                                labels = c("<50k", "50k-99k", "100k+"))

abcd <- abcd %>% select(-c(nihtbx_totalcomp_uncorrected))


# Reorder the abcd dataset
abcd <- abcd %>% select(c("subjectkey","age10", "gender", "race_ethnicity", "household_income", "marital_status",
                          "family_size3", "acs_raked_propensity_score"))
abcd$age10 <- factor(abcd$age10, levels = c(0, 1), labels = c("0", "1"))
abcd$family_size3 <- factor(abcd$family_size3, levels = c(0, 1), labels = c("0", "1"))
abcd$complete_case = (complete.cases(abcd))

library(dplyr)
tbss_score = read.csv("./RealDataAnalysis_CIFTI/ABCD_lavaan_G_time2.csv")
tbss_score$src_subject_id = stringr::str_remove(tbss_score$subjectkey, pattern = "_")
abcd_gf <- dplyr::left_join(abcd, tbss_score, by = "subjectkey")
abcd_gf$w_gfactor = (!is.na(abcd_gf$G_lavaan.baseline))

# --- filter imaging subsample --- #
filelist = list.files("./RealDataAnalysis_CIFTI/nBack/baselineYear1Arm1/")
subject_list = stringr::str_remove(filelist, "_zstat11.dtseries.nii")
abcd_gf$nback_included = (abcd_gf$src_subject_id %in% subject_list) # subjects with nback imaging data
abcd_gf$imaging_subsample = (abcd_gf$w_gfactor & abcd_gf$nback_included & abcd_gf$complete_case) 

abcd_gf$complete_case_in_fit = complete.cases(abcd_gf[, c("imaging_subsample","age10", "gender", "race_ethnicity", "household_income", 
                                                          "marital_status", "family_size3", "G_lavaan.baseline", "acs_raked_propensity_score")])

abcd_gf_fit = abcd_gf[abcd_gf$complete_case_in_fit,]
# > nrow(abcd_gf_fit)
# [1] 9907


# --- Propensity score --- #
m1fit2 <- glm(imaging_subsample ~ age10 + race_ethnicity + household_income + marital_status + family_size3 + gender + G_lavaan.baseline, 
              data=abcd_gf_fit, family=binomial(),
                control = glm.control(maxit = 50))


summary(m1fit2)
# inverse propensity score of quality control phrase
abcd_gf_fit$weight_qc  = 1/predict(m1fit2, type="response")
abcd_gf_fit$weight_final = abcd_gf_fit$weight_qc * abcd_gf_fit$acs_raked_propensity_score
hist(abcd_gf_fit$weight_qc)

# Plot the coefficients
library(jtools)
library(dplyr)
library(ggplot2)
jtools::plot_coefs(m1fit2, coefs = c("Age > 10" = "age101", 
                                     "Race/Ethnicity: Black" = "race_ethnicityBlack", 
                                     "Race/Ethnicity: Hispanic" = "race_ethnicityHispanic", 
                                     "Race/Ethnicity: Asian" = "race_ethnicityAsian", 
                                     "Race/Ethnicity: Other" = "race_ethnicityOther", 
                                     "Household Income between 50K-99K" = "household_income50k-99k", 
                                     "Household Income > 100K" = "household_income100k+", 
                                     "Marital Status (Not Married)" = "marital_statusNot Married", 
                                     "Household Size > 3" = "family_size31", 
                                     "Gender (Female)" = "genderF", 
                                     "Baseline g-factor" = "G_lavaan.baseline"),
                   line.size = 2) + theme_minimal(base_size = 20) + ylab("Variable") + xlab("Coefficient Estimate")
abcdsub_gf$propmeth2  = abcd_raw$propmeth2 = NA
abcdsub_gf$propmeth2[abcdsub_gf$w_gfactor] = abcd_raw$propmeth2[abcdsub_gf$w_gfactor] = 1/predict(m1fit2, type="response")

svyabcd_imaging_design <- svydesign(id=~1, strata=NULL, weights=~propmeth2, data=abcdsub_gf[abcdsub_gf$imaging_subsample,])


# Show the weighted distribution of the covariates
svymean(~household_income, svyabcd_imaging_design)
svymean(~household_income, svy_abcd_uw)

svymean(~race_ethnicity, svyabcd_imaging_design)
svymean(~race_ethnicity, svy_abcd_uw)


# Final weight = inversed propensity score * raked weight
abcdsub_gf$weight_final= abcd_raw$weight_final = abcdsub_gf$propmeth2 * abcdsub_gf$acs_raked_propensity_score
abcdsub_gfsub = abcdsub_gf[abcdsub_gf$imaging_subsample == 1,]

# Raked weight
# survey design data for acbd only, set psu =1 and no strata, with weight 
svyabcd_gf1 <- svydesign(id=~1, strata=NULL, weights=~weight_final, data=abcdsub_gfsub)

# rake using population data above 
# pop.age <- data.frame(age10=c(0,1),Freq=c(4074807,4136798)) 
# pop.sex <- data.frame(gender=c("M","F"),Freq=c(4205925,4005680))
# pop.race <- data.frame(race_ethnicity=c("White","Black","Hispanic","Asian","Other"), Freq=c(4305552, 1101297, 1973827, 487673, 343256)) 
# svyabcd_test_4  <- rake(svyabcd_gf1, list(~age10, ~gender, ~race_ethnicity), list(pop.age, pop.sex, pop.race))
# # extract weight 
# svyabcd_test_4$r_weight_final <-(weights(svyabcd_test_4))
# abcdsub_gfsub$r_weight_final <- svyabcd_test_4$r_weight_final


svy_abcd_imaging_uw <- svydesign(id=~1, data=abcdsub_gfsub, na.rm = T)
svy_abcd_imaging <- svydesign(id=~1, data=abcdsub_gfsub, weights=~weight_final, na.rm = T)
print(svymean(~age10, svy_abcd_imaging_uw))
print(svymean(~age10, svy_abcd_imaging))

print(round(svymean(~race_ethnicity, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~race_ethnicity, svy_abcd_imaging)*100,1))

print(round(svymean(~household_income, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~household_income, svy_abcd_imaging)*100,1))

print(round(svymean(~marital_status, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~marital_status, svy_abcd_imaging)*100,1))

print(round(svymean(~gender, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~gender, svy_abcd_imaging)*100,1))

print(round(svymean(~age10, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~age10, svy_abcd_imaging)*100,1))



print(round(svymean(~family_size3, svy_abcd_imaging_uw)*100,1))
print(round(svymean(~family_size3, svy_abcd_imaging)*100,1))
print(round(svymean(~family_size3, svy_abcd_uw)*100,1))
print(round(svymean(~family_size3, svy_abcd)*100,1))


# saveRDS(abcdsub_gf, "./RealDataAnalysis_CIFTI/nBack_gfactor/abcdsub_gf.rds")
# saveRDS(abcdsub_gfsub, "./RealDataAnalysis_CIFTI/nBack_gfactor/abcdsub_gfsub.rds")
# saveRDS(svy_abcd_imaging, "./RealDataAnalysis_CIFTI/nBack_gfactor/svy_abcd_imaging.rds")
# 

