# Weighted Image-on-scalar Regression Analyses of Large Scale Neuroimaging Data


# Data Preprocessing

- ACS2011_15_I.Rds (See https://github.com/ABCD-STUDY/abcd_acs_raked_propensity): - This is the ACS 2011-2015 data set with the demographic (demo) and socio-economic (SES) variables that are used as the benchmark for the propensity weight estimation. The "I" suffix indicates that the missing data in this input file are already imputed with a single multivariate imputation. This file will be used as the input for the Step 2 program below. It serves as our permanent benchmark/reference for the population. The data file, "ACS2011_15_I.Rds" is created directly from the American Community Survey 2011-2015 Public Use Microsample File. Some variables were recoded to make them consistent with ABCD demographic and SES constructs but the ACS data is in the public domain. The link to the ACS user guide and best citation is: https://www2.census.gov/programssurveys/acs/tech_docs/pums/ACS2011_2015_PUMS_README.pdf. In addition to some variable recoding, missing data in the original ACS 2011-2015 PUMS file was singly imputed using the SAS V9.4 MI procedure.

- 
abcd_bl_preprocessed.rds (not provided)
