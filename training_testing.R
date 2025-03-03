library(dplyr)
library(ciftiTools)
library(Rcpp)
library(dplyr)

# library(ciftiTools)
# ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosx64/wb_command')
# std_img = read_cifti("./Simulation_Cifti/cifti_template.dtseries.nii")

# Rcpp::sourceCpp("./RealDataAnalyis/sd_est_svcm_main_only.cpp")

setwd("/Volumes/ExternalDisk/Dropbox (Personal)/U-Mich/Research/WIR_ABCD/")

sourceCpp("./RealDataAnalysis_CIFTI/population_SGD_interaction.cpp")
sourceCpp("./scripts/weighted_mvr.cpp")

ABCD_task_general_gfactor <- read.csv("./RealDataAnalysis_CIFTI/nBack_gfactor/ABCD_task_general_gfactor.csv", as.is = T)
X <- ABCD_task_general_gfactor$G_lavaan.baseline
basisFuns_L <- readRDS("./RealDataAnalysis_CIFTI/basis_umat_L.rds")
basisFuns_R <- readRDS("./RealDataAnalysis_CIFTI/basis_umat_R.rds")
yu_L <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nback_imaging_L_basis.rds")
yu_R <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nBack_imaging_R_basis.rds")
y_L <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nback_imaging_L.rds")
y_R <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nback_imaging_R.rds")
weight_dat <- ABCD_task_general_gfactor$weight_final
weight_dat <- weight_dat/mean(weight_dat)



# Z <- read.csv("./RealDataAnalysis_CIFTI/nBack_gfactor/ABCD_task_general_gfactor.csv", as.is = T)
abcdsub_gf <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/abcdsub_gfsub.rds")
abcdsub_gf$age10 <- as.integer(abcdsub_gf$age10)
Z <- abcdsub_gf %>% select(c(age10, site_id_l.baseline, gender, race_ethnicity, marital_status, family_size3, household_income))
Z <- fastDummies::dummy_columns(Z, c(
  "race_ethnicity",
  "family_size3",
  "marital_status",
  "household_income",
  "site_id_l.baseline",
  "gender"), 
  remove_first_dummy = T,
  remove_selected_columns = T)
library(Rcpp)
library(dplyr)

inclusion_probability <- 1/ABCD_task_general_gfactor$r_weight_final
n_rep = 100
sample_prob = 0.9
n = length(X)
rmse_unweighted_L = matrix(0, n_rep, nrow(basisFuns_L))
rmse_weighted_R = matrix(0, n_rep, nrow(basisFuns_R))
rmse_unweighted_R = matrix(0, n_rep, nrow(basisFuns_R))
rmse_weighted_L = matrix(0, n_rep, nrow(basisFuns_L))

# ------------------------------------------------------------------------------------- #


for(i in 1:n_rep){
  if (i %% 10 == 0) cat(i, "\n")
  training_sample_idx = sample(1:n, size = ceiling(n*sample_prob), replace = F, prob = inclusion_probability)
    test_sample_idx = setdiff(1:n, training_sample_idx)
  X_training_sample = X[training_sample_idx]
  Z_training_sample = Z[training_sample_idx,]
  Z_test_sample = Z[test_sample_idx,]
  X_test_sample = X[test_sample_idx]
  
    yu_L_training_sample = yu_L[training_sample_idx,]
    yu_R_training_sample = yu_R[training_sample_idx,]
    yu_L_test_sample = yu_L[test_sample_idx,]
    yu_R_test_sample = yu_R[test_sample_idx,]
    y_L_test_sample = y_L[test_sample_idx,]
    y_R_test_sample = y_R[test_sample_idx,]
    weight_dat_training_sample = weight_dat[training_sample_idx]
    
  EM_unweighted_L <- weighted_mvr(X = cbind(1,as.matrix(X_training_sample)),
                                  Y = yu_L_training_sample,
                                  w = rep(1, length(weight_dat_training_sample)),
                                  const_var = 0)
  
  EM_weighted_L <- weighted_mvr(X = cbind(1,as.matrix(X_training_sample)),
                                Y = yu_L_training_sample,
                                w = weight_dat_training_sample,
                                const_var = 0)
  
  EM_unweighted_R <- weighted_mvr(X = cbind(1,as.matrix(X_training_sample)),
                                  Y = yu_R_training_sample,
                                  w = rep(1, length(weight_dat_training_sample)),
                                  const_var = 0)
  
  EM_weighted_R <- weighted_mvr(X = cbind(1,as.matrix(X_training_sample)),
                                Y = yu_R_training_sample,
                                w = weight_dat_training_sample,
                                const_var = 0)
  
  
    # Evaluate the prediction
    y_L_test_pred_uw = cbind(1,as.matrix(X_test_sample)) %*% (EM_unweighted_L$beta %*% t(basisFuns_L))
    y_R_test_pred_uw = cbind(1,as.matrix(X_test_sample)) %*% (EM_unweighted_R$beta %*% t(basisFuns_R))

    y_L_test_pred_w = cbind(1,as.matrix(X_test_sample)) %*% (EM_weighted_L$beta %*% t(basisFuns_L))
    y_R_test_pred_w = cbind(1,as.matrix(X_test_sample)) %*% (EM_weighted_R$beta %*% t(basisFuns_R))

    rmse_unweighted_L[i,] = sqrt(colMeans((y_L_test_sample - y_L_test_pred_uw)^2))
    rmse_unweighted_R[i,] = sqrt(colMeans((y_R_test_sample - y_R_test_pred_uw)^2))
    rmse_weighted_L[i,] = sqrt(colMeans((y_L_test_sample - y_L_test_pred_w)^2))
    rmse_weighted_R[i,] = sqrt(colMeans((y_R_test_sample - y_R_test_pred_w)^2))

}

# Use ggpubr to plot the results (boxplot)
rmse_unweighted_L_mean = apply(rmse_unweighted_L, 2, mean)
rmse_unweighted_R_mean = apply(rmse_unweighted_R, 2, mean)
rmse_weighted_L_mean = apply(rmse_weighted_L, 2, mean)
rmse_weighted_R_mean = apply(rmse_weighted_R, 2, mean)

# Create a data frame with the mean RMSE values and network names
rmse_data <- data.frame(
  Method_Hemisphere = c(rep("Unweighted Left", length(rmse_unweighted_L_mean)),
              rep("Unweighted Right", length(rmse_unweighted_R_mean)),
              rep("Weighted Left", length(rmse_weighted_L_mean)),
              rep("Weighted Right", length(rmse_weighted_R_mean))),
  RMSE = c(rmse_unweighted_L_mean, rmse_unweighted_R_mean, rmse_weighted_L_mean, rmse_weighted_R_mean)
)

library(ggpubr)
ggboxplot(rmse_data, x = "Method_Hemisphere", y = "RMSE", 
          fill = "Method_Hemisphere", palette = "jco",
          title = "RMSE Comparison",
          caption = "Boxplot showing RMSE comparison between different methods")


cat("Without control variable adjusted: ", "\n")
cat("Left (unweighted) RMSE: ", mean(rmse_unweighted_L_mean),  "\n")
cat("Left (weighted) RMSE: ", mean(rmse_weighted_L_mean),  "\n")
cat("Right (unweighted) RMSE: ", mean(rmse_unweighted_R_mean),  "\n")
cat("Right (weighted) RMSE: ", mean(rmse_weighted_R_mean),  "\n")


est_unweighted_corrected <- readcii("./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_unweighted_corrected.dtseries.nii")
est_weighted_corrected <- readcii("./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_weighted_corrected.dtseries.nii")

idx_unweighted_L_nonzero = which(est_unweighted_corrected$data$cortex_left[,1] != 0)
idx_unweighted_R_nonzero = which(est_unweighted_corrected$data$cortex_right[,1] != 0)
idx_weighted_L_nonzero = which(est_weighted_corrected$data$cortex_left[,1] != 0)
idx_weighted_R_nonzero = which(est_weighted_corrected$data$cortex_right[,1] != 0)


# Save the results
saveRDS(list(rmse_unweighted_L_mean = rmse_unweighted_L_mean, rmse_unweighted_R_mean = rmse_unweighted_R_mean, 
             rmse_weighted_L_mean = rmse_weighted_L_mean, rmse_weighted_R_mean = rmse_weighted_R_mean, 
             idx_unweighted_L_nonzero = idx_unweighted_L_nonzero, idx_unweighted_R_nonzero = idx_unweighted_R_nonzero, 
             idx_weighted_L_nonzero = idx_weighted_L_nonzero, idx_weighted_R_nonzero = idx_weighted_R_nonzero,
             est_unweighted_corrected = est_unweighted_corrected, est_weighted_corrected = est_weighted_corrected), 
              "./RealDataAnalysis_CIFTI/nBack_gfactor/rmse_results_without_Z.rds")




