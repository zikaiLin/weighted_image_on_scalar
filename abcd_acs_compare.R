acs <- readRDS("./RealDataAnalysis_CIFTI/nBack/ACS2011_15_I.Rds")
abcd <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/abcd_bl_preprocessed.rds")
# > names(abcd)
# [1] "subjectkey"                   "family_size3"                 "marital_status"              
# [4] "household_income"             "age10"                        "race_ethnicity"              
# [7] "gender"



# > names(acs)
# [1] "PWGTP"     "age"       "sex"       "region"    "race_eth"  "faminc"    "famtype"  
# [8] "hhsize"    "allmult"   "id_redcap" "site_name" "fesabcd"   "acsflag"  

# --- Race --- #

abcd$race_ethnicity <- factor(abcd$race_ethnicity)
# > unique(abcd$race_ethnicity)
# [1] "White"    "Asian"    "Black"    "Other"    "Hispanic"
acs$race_eth <- factor(acs$race_eth, levels = c(1,2,3,4,5,6,8,9),
 labels = c("White", "Black", "Hispanic", "Asian" ,"AIAN", "NHPI", "Other" ,"MULTIPLE")) 
table (acs$race_eth, exclude=NULL)

# Combine categories and refactor
acs$race_eth[acs$race_eth %in% c("AIAN", "NHPI", "Other")] <- "Asian"
acs$race_eth[acs$race_eth %in% c("MULTIPLE")] <- "Other"

# Remove redundant level 
acs$race_eth <- droplevels(acs$race_eth)
table (acs$race_eth, exclude=NULL)/nrow(acs)

# --- age --- #
abcd$age <- ifelse(abcd$age10 == 1, 10, 9)
abcd$age <- factor(abcd$age, levels = c(9, 10), labels = c("9", "10"))
table(abcd$age, exclude = NULL)/nrow(abcd)
acs$age <- factor(acs$age, levels = c(9,10), labels = c("0", "1"))

# --- Sex --- #
abcd$gender[abcd$gender == ""] <- NA  # Drop the ""
abcd$gender <- droplevels(abcd$gender)
table(abcd$gender, exclude = NULL)/nrow(abcd)

acs$gender <- factor(acs$sex, levels = c(1,2), labels = c("M", "F")) 
table (acs$gender, exclude=NULL) 

# --- Household income --- #
acs$faminc <- factor(acs$faminc, levels = c(1:6,888,999,777), 
 labels = c("<25k", "25k-49k", "50k-74k", "75k-99k", "100k-199k", "200k+", "<50k", "50k-99k", "100k+")) 
table (acs$faminc, exclude=NULL)/nrow(acs)

acs$faminc[acs$faminc %in% c("<25k", "25k-49k", "<50k")] <- "<50k"
acs$faminc[acs$faminc %in% c("50k-74k", "75k-99k", "50k-99k")] <- "50k-99k"
acs$faminc[acs$faminc %in% c("100k-199k", "200k+", "100k+")] <- "100k+"
acs$faminc <- droplevels(acs$faminc)
table (acs$faminc, exclude=NULL)/nrow(acs)

abcd$household_income <- factor(abcd$household_income,
                                levels = c("Under $50,000", "$50,000 to $100,000", "Over $100,000"),
                                labels = c("<50k", "50k-99k", "100k+"))


# --- household size --- #
acs$hhsize  <- factor(acs$hhsize, levels = c(1,4,5,6,7,888), 
                          labels = c( "2-3 Persons", 
                          "4 Persons" ,
                          "5 Persons",
                           "6 Persons",
                           "7 or More Persons",
                           "4 or More Persons" )) 

# Binarized to family size >= 4 and family size < 4
acs$hhsize[acs$hhsize %in% c("2-3 Persons")] <- factor("2-3 Persons")
acs$hhsize[acs$hhsize %in% c("4 Persons", "5 Persons", "6 Persons", "7 or More Persons")] <- factor("4 or More Persons")
acs$hhsize <- droplevels(acs$hhsize)
acs$hhsize <- factor(acs$hhsize, levels = c("2-3 Persons", "4 or More Persons"),
    labels = c(0, 1))



table (acs$hhsize, exclude=NULL) 
table(acs$hhsize, exclude=NULL)/nrow(acs)

# --- Family type --- #
acs$famtype  <- factor(acs$famtype, levels = c(1:2), labels = c("Married", "Not Married")) 
table (acs$famtype , exclude=NULL) 
table(acs$famtype, exclude=NULL)/nrow(acs)


##################### Combine the two datasets #####################
acs <- acs %>% select(-c(allmult, id_redcap, site_name, fesabcd, region, sex))
acs$acs_raked_propensity_score <- NA

# acs <- acs[,1:8]
abcd <- abcd %>% select(-c(nihtbx_totalcomp_uncorrected, subjectkey, age))
abcd$acsflag <- 0
abcd$PWGTP <- 1
# > names(acs)
#  [1] "age"                        "sex"                        "race_eth"                   "faminc"                    
#  [5] "famtype"                    "hhsize"                     "acsflag"                    "race_ethc"                 
#  [9] "gender"                     "acs_raked_propensity_score"

# Rename variables of acs to match abcd
names(acs) <- c("PWGTP","age10", "race_ethnicity",  "household_income", "marital_status", "family_size3",
                "acsflag","gender", "acs_raked_propensity_score")

# Reorder the abcd dataset
abcd <- abcd %>% select(c("PWGTP","age10", "gender", "race_ethnicity", "household_income", "marital_status", "family_size3", "acs_raked_propensity_score", "acsflag"))
abcd$age10 <- factor(abcd$age10, levels = c(0, 1), labels = c("0", "1"))
abcd$family_size3 <- factor(abcd$family_size3, levels = c(0, 1), labels = c("0", "1"))
names(abcd)

##################### Survey weight mean #####################
abcd_i <- mice(abcd, maxit = 0)
abcd_i <- complete(abcd_i, 1)


# Combine the two datasets
all_samples_acs_abcd = rbind(acs, abcd_i)

names(all_samples_acs_abcd) 
head(all_samples_acs_abcd)
summary(all_samples_acs_abcd)

all_samples_acs_abcd$acsflag <- relevel(factor(all_samples_acs_abcd$acsflag), ref="1")
all_samples_acs_abcd$abcdflag <- factor(as.integer(all_samples_acs_abcd$acsflag)-1)

m1fit1 <- glm(abcdflag ~ age10 + race_ethnicity + household_income + marital_status + family_size3 + gender, 
              data=all_samples_acs_abcd, family=quasibinomial(), weights=PWGTP,
              control = glm.control(maxit = 50))
summary(m1fit1)

all_samples_acs_abcd$propmeth1 = predict(m1fit1, type="response") 
summary(all_samples_acs_abcd$PWGTP, exclude = NULL)
summary(all_samples_acs_abcd$propmeth1) 

# assign different weight values depending on acsflag status 
all_samples_acs_abcd$pwgtmeth1=1/all_samples_acs_abcd$propmeth1
summary(all_samples_acs_abcd$pwgtmeth1) 


describeBy(all_samples_acs_abcd$pwgtmeth1, all_samples_acs_abcd$abcdflag)
quantile(all_samples_acs_abcd$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))


# acsflag==0
abcdsub <- subset(all_samples_acs_abcd,acsflag==0)
quantile(abcdsub$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))


# Trim weight at .02 and .98 percentiles for acsflag==0 only, note these are manual entry from quantiles above
all_samples_acs_abcd$pwgtmeth1 [all_samples_acs_abcd$abcdflag==1 & all_samples_acs_abcd$pwgtmeth1 < 170] <- 170
all_samples_acs_abcd$pwgtmeth1 [all_samples_acs_abcd$abcdflag==1 & all_samples_acs_abcd$pwgtmeth1 > 1321] <- 1321

describeBy(all_samples_acs_abcd$pwgtmeth1, all_samples_acs_abcd$acsflag)

# Adjust the weights for acsflag==0
all_samples_acs_abcd$pwgtmeth1 <- ifelse(all_samples_acs_abcd$acsflag %in% c(0), 
(all_samples_acs_abcd$pwgtmeth1*8211605/8072664), all_samples_acs_abcd$pwgtmeth1) 

describeBy(all_samples_acs_abcd$pwgtmeth1, all_samples_acs_abcd$acsflag) 

# ---------------------------------------------------- #
library(survey)

abcdsub <- subset(all_samples_acs_abcd,acsflag==0)
acssub <- subset(all_samples_acs_abcd,acsflag==1)

# set survey design with single psu for each respondent and weight, no strata or FPC 
svy_acs <- svydesign(id=~1, weights = ~PWGTP, data=acssub)
svy_abcd <- svydesign(id=~1, weights = ~acs_raked_propensity_score, data=abcdsub)
svy_abcd_uw <- svydesign(id=~1, data=abcdsub)

print(svymean(~age10, svy_abcd))
print(svymean(~age10, svy_acs))



print(svymean(~household_income, svy_abcd))
print(svymean(~household_income, svy_abcd_uw))
print(svymean(~household_income, svy_acs))

print(svymean(~race_ethnicity, svy_abcd))
print(svymean(~race_ethnicity, svy_abcd_uw))
print(svymean(~race_ethnicity, svy_acs))

# ---------------------------------------------------- #
#                         Rake                        #
# ---------------------------------------------------- #

pop.age <- data.frame(age10=c(0,1),Freq=c(4074807,4136798)) 
pop.sex <- data.frame(gender=c("M","F"),Freq=c(4205925,4005680))
pop.race <- data.frame(race_ethnicity=c("White","Black","Hispanic","Asian","Other"), Freq=c(4305552, 1101297, 1973827, 487673, 343256)) 
pop.age
pop.sex
pop.race


# survey design data for acbd only, set psu =1 and no strata, with weight 
svyabcd1 <- svydesign(id=~1, strata=NULL, weights=~pwgtmeth1, data=abcdsub)

# rake using population data above 
svyabcd_test_3  <- rake(svyabcd1, list(~race_ethnicity), list(pop.race ))
# extract weight 
svyabcd_test_3$rpwgtmeth1 <-(weights(svyabcd_test_3))
# add to abcdsub data frame 
abcdsub$rpwgtmeth1 <- svyabcd_test_3$rpwgtmeth1 
names(abcdsub)

# check raked totals 
svytable (~age10, svyabcd_test_3) 
svytable (~gender, svyabcd_test_3)
svytable (~race_ethnicity, svyabcd_test_3)

svymean (~race_ethnicity, svyabcd_test_3)
svymean (~race_ethnicity, svyabcd_test_3)


# check with raked weight rpwgtmeth1 
wtd.table(abcdsub$age10, weights=abcdsub$rpwgtmeth1) 
wtd.table(abcdsub$gender, weights=abcdsub$rpwgtmeth1)
wtd.table(abcdsub$race_ethnicity, weights=abcdsub$rpwgtmeth1)




# load psych package for use of descriptive tools 
library(psych)
describeBy(all_samples_acs_abcd$pwgtmeth1, all_samples_acs_abcd$acsflag)

svy_abcd <- svydesign(id=~1, data=abcd_i, weights=~acs_raked_propensity_score, na.rm = T)

svy_abcd_uw <- svydesign(id=~1, data=abcd_i, na.rm = T)
svy_acs <- svydesign(id=~1, data=acs, na.rm = T)
# Use survey mean to calculate the mean of the survey weighted data
print(svymean(~age10, svy_abcd))
print(svymean(~gender, svy_abcd))
print(svymean(~household_income, svy_abcd))
print(svymean(~household_income, svy_acs))
print(svymean(~household_income, svy_abcd_uw))


print(svymean(~race_ethnicity, svy_abcd_uw))
print(svymean(~race_ethnicity, svy_abcd))
print(svymean(~race_ethnicity, svy_acs))

print(svymean(~marital_status, svy_abcd_uw))
print(svymean(~marital_status, svy_abcd))
print(svymean(~marital_status, svy_acs))


print(svymean(~household_income, svy_abcd_uw))
print(svymean(~household_income, svy_abcd))
print(svymean(~household_income, svy_acs))

# ---------------------------- Imaging Reweight --------------------- #
abcdsub$subjectkey = abcd$subjectkey

library(dplyr)
tbss_score = read.csv("./RealDataAnalysis_CIFTI/ABCD_lavaan_G_time2.csv")
tbss_score$src_subject_id = stringr::str_remove(tbss_score$subjectkey, pattern = "_")
abcdsub_gf <- dplyr::left_join(abcdsub, tbss_score, by = "subjectkey")
abcdsub_gf$w_gfactor = (!is.na(abcdsub_gf$G_lavaan.baseline))



filelist = list.files("./RealDataAnalysis_CIFTI/nBack/baselineYear1Arm1/")
subject_list = stringr::str_remove(filelist, "_zstat11.dtseries.nii")
abcdsub_gf$nback_included = (abcdsub_gf$src_subject_id %in% subject_list)
abcdsub_gf$imaging_subsample = (abcdsub_gf$w_gfactor & abcdsub_gf$nback_included)

m1fit2 <- glm(imaging_subsample ~ age10 + race_ethnicity + household_income + marital_status + family_size3 + gender, 
              data=abcdsub_gf, family=quasibinomial(),
              control = glm.control(maxit = 50))
summary(m1fit2)
abcdsub_gf$propmeth2 = predict(m1fit2, type="response")
abcdsub_gf$weight_final=1/abcdsub_gf$propmeth2 * abcdsub_gf$rpwgtmeth1

# Raked weight
# survey design data for acbd only, set psu =1 and no strata, with weight 
svyabcd_gf1 <- svydesign(id=~1, strata=NULL, weights=~weight_final, data=abcdsub_gf)

# rake using population data above 
svyabcd_test_4  <- rake(svyabcd_gf1, list(~age10, ~gender, ~race_ethnicity), list(pop.age, pop.sex, pop.race))
# extract weight 
svyabcd_test_4$r_weight_final <-(weights(svyabcd_test_4))
abcdsub_gf$r_weight_final <- svyabcd_test_4$r_weight_final


describeBy(abcdsub_gf$weight_final, abcdsub_gf$imaging_subsample)

svy_abcd_imaging_uw <- svydesign(id=~1, data=abcdsub_gf, subset = imaging_subsample == 1, na.rm = T)
svy_abcd_imaging <- svydesign(id=~1, data=abcdsub_gf, weights=~r_weight_final, subset = imaging_subsample == 1, na.rm = T)
print(svymean(~age10, svy_abcd_imaging_uw))
print(svymean(~age10, svy_abcd_imaging))

print(svymean(~household_income, svy_abcd_imaging_uw))
print(svymean(~household_income, svy_abcd_imaging))


saveRDS(abcdsub_gf, "./RealDataAnalysis_CIFTI/nBack_gfactor/abcdsub_gf.rds")
