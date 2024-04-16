library(dplyr)

# Preprocessing of the ABCD baseline demographics data
all_abcd_general <- readRDS("./abcd_bl.rds")

all_abcd_general$gender <- as.factor(all_abcd_general$gender)
table(all_abcd_general$race_ethnicity)/nrow(all_abcd_general)


# Race code
# 1: White and ""
# 2: Black or African American
# 3: Asian
# 4: Hispanic
# 5: Other

# Recode the race_ethnicity variable
all_abcd_general$race_ethnicity[all_abcd_general$race_ethnicity == "1" | all_abcd_general$race_ethnicity == ""] <- "White"
all_abcd_general$race_ethnicity[all_abcd_general$race_ethnicity == "2"] <- "Black"
all_abcd_general$race_ethnicity[all_abcd_general$race_ethnicity == "3"] <- "Hispanic"
all_abcd_general$race_ethnicity[all_abcd_general$race_ethnicity == "4"] <- "Asian"
all_abcd_general$race_ethnicity[all_abcd_general$race_ethnicity == "5"] <- "Other"

# Recode the Age variable and binarize with Age >= 10 and Age < 10
all_abcd_general$age <- as.integer(all_abcd_general$interview_age.x)/12
all_abcd_general$age10 <- as.integer(all_abcd_general$age >= 10)

# Recode the gender and drop labels
all_abcd_general$gender[all_abcd_general$gender == ""] <- NA  # Drop the ""
all_abcd_general$gender <- droplevels(all_abcd_general$gender)

# Recode the Household Income variable

# 1 to 6: Under $50,000
# 7 to 8: $50,000 to $100,000
# 9 to 10: Over $100,000
# other: unknown or NA

all_abcd_general$demo_comb_income_v2[all_abcd_general$demo_comb_income_v2 %in% c("1", "2", "3", "4", "5", "6")] <- "Under $50,000"
all_abcd_general$demo_comb_income_v2[all_abcd_general$demo_comb_income_v2 %in% c("7", "8")] <- "$50,000 to $100,000"
all_abcd_general$demo_comb_income_v2[all_abcd_general$demo_comb_income_v2 %in% c("9", "10")] <- "Over $100,000"
all_abcd_general$demo_comb_income_v2[all_abcd_general$demo_comb_income_v2 %in% c("999", "777", "")] <- NA
all_abcd_general$household_income <- as.factor(all_abcd_general$demo_comb_income_v2)

# Recode the Marital variable

# 1 and 6: Married
# Other

all_abcd_general$demo_prnt_marital_v2[all_abcd_general$demo_prnt_marital_v2 != 1 & all_abcd_general$demo_prnt_marital_v2 != 6] <- "Not Married"
all_abcd_general$demo_prnt_marital_v2[all_abcd_general$demo_prnt_marital_v2 == 1 | all_abcd_general$demo_prnt_marital_v2 == 6] <- "Married"
all_abcd_general$marital_status <- as.factor(all_abcd_general$demo_prnt_marital_v2)

# Recode the household size variable
# Binarized to family size >= 3 and family size < 3
all_abcd_general$demo_roster_v2 <- as.integer(all_abcd_general$demo_roster_v2)
all_abcd_general$family_size3 <- as.integer(all_abcd_general$demo_roster_v2 >= 4)



# Select subjectkey, gender, marrital status, race, family size, household income, and age
all_abcd_general <- all_abcd_general %>% select(subjectkey, family_size3, marital_status, household_income, age10, race_ethnicity, gender, nihtbx_totalcomp_uncorrected, acs_raked_propensity_score)
all_abcd_general$acs_raked_propensity_score = as.numeric(all_abcd_general$acs_raked_propensity_score)
saveRDS(all_abcd_general, "./RealDataAnalysis_CIFTI/nBack_gfactor/abcd_bl_preprocessed.rds")

