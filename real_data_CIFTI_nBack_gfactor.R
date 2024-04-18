library(Rcpp)
library(dplyr)
library(ciftiTools)
library(ggplot2)

ciftiTools.setOption('wb_path', './workbench_macos/bin_macosx64/wb_command')
std_img = read_cifti("./RealDataAnalysis_CIFTI/cifti_template.dtseries.nii")

# Rcpp::sourceCpp("./RealDataAnalyis/sd_est_svcm_main_only.cpp")
sourceCpp("./scripts/weighted_mvr.cpp")


ABCD_task_general_gfactor <- read.csv("./RealDataAnalysis_CIFTI/nBack_gfactor/ABCD_task_general_gfactor.csv", as.is = T)
X <- ABCD_task_general_gfactor$G_lavaan.baseline
basisFuns_L <- readRDS("./RealDataAnalysis_CIFTI/basis_umat_L.rds")
basisFuns_R <- readRDS("./RealDataAnalysis_CIFTI/basis_umat_R.rds")
yu_L <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nback_imaging_L_basis.rds")
yu_R <- readRDS("./RealDataAnalysis_CIFTI/nBack_gfactor/nBack_imaging_R_basis.rds")

# Normalize the weights
weight_dat <- ABCD_task_general_gfactor$weight_final
weight_dat <- weight_dat/mean(weight_dat)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Weighted image-on-scalar regression
EM_unweighted_L <- weighted_mvr(X = cbind(1,as.matrix(X)),
                                Y = yu_L,
                                w = rep(1, length(weight_dat)),
                                const_var = 0)

EM_weighted_L <- weighted_mvr(X = cbind(1,as.matrix(X)),
                              Y = yu_L,
                              w = weight_dat,
                              const_var = 0)

EM_unweighted_R <- weighted_mvr(X = cbind(1,as.matrix(X)),
                                Y = yu_R,
                                w = rep(1, length(weight_dat)),
                                const_var = 0)

EM_weighted_R <- weighted_mvr(X = cbind(1,as.matrix(X)),
                              Y = yu_R,
                              w = weight_dat,
                              const_var = 0)

## Standard error estimation
sd_est_unweighted_L <- sqrt(transform_vcov(vcov_basis = EM_unweighted_L$H_sd[1:2,],
                                           basisFuns = basisFuns_L))
sd_est_weighted_L <- sqrt(transform_vcov(vcov_basis = EM_weighted_L$H_sd[1:2,],
                                         basisFuns = basisFuns_L))

sd_est_unweighted_R <- sqrt(transform_vcov(vcov_basis = EM_unweighted_R$H_sd[1:2,],
                                           basisFuns = basisFuns_R))
sd_est_weighted_R <- sqrt(transform_vcov(vcov_basis = EM_weighted_R$H_sd[1:2,],
                                         basisFuns = basisFuns_R))


# Save the SE results
sd_rda_data_L = data.frame(SD_intercept_unweighted = sd_est_unweighted_L[1,],
                           SD_intercept_weighted = sd_est_weighted_L[1,],
                           SD_ext_unweighted = sd_est_unweighted_L[2,],
                           SD_ext_weighted = sd_est_weighted_L[2,])

sd_rda_data_R = data.frame(SD_intercept_unweighted = sd_est_unweighted_R[1,],
                           SD_intercept_weighted = sd_est_weighted_R[1,],
                           SD_ext_unweighted = sd_est_unweighted_R[2,],
                           SD_ext_weighted = sd_est_weighted_R[2,])

write.csv(sd_rda_data_L, "./RealDataAnalysis_CIFTI/nBack_gfactor/standard_deviation_L.csv", row.names = F)
write.csv(sd_rda_data_R, "./RealDataAnalysis_CIFTI/nBack_gfactor/standard_deviation_R.csv", row.names = F)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Transform the basis coefficients to the original space
est_unweighted_L = EM_unweighted_L$beta[1:2,] %*% t(basisFuns_L)
est_weighted_L = EM_weighted_L$beta[1:2,] %*% t(basisFuns_L)
est_unweighted_R = EM_unweighted_R$beta[1:2,] %*% t(basisFuns_R)
est_weighted_R = EM_weighted_R$beta[1:2,] %*% t(basisFuns_R)


# Save the estimated coefficients
est_rda_data_L = data.frame(est_intercept_unweighted = est_unweighted_L[1,],
                            est_intercept_weighted = est_weighted_L[1,],
                            est_ext_unweighted = est_unweighted_L[2,],
                            est_ext_weighted = est_weighted_L[2,])
est_rda_data_R = data.frame(est_intercept_unweighted = est_unweighted_R[1,],
                            est_intercept_weighted = est_weighted_R[1,],
                            est_ext_unweighted = est_unweighted_R[2,],
                            est_ext_weighted = est_weighted_R[2,])

write.csv(est_rda_data_L, "./RealDataAnalysis_CIFTI/nBack_gfactor/coef_est_L.csv", row.names = F)
write.csv(est_rda_data_R, "./RealDataAnalysis_CIFTI/nBack_gfactor/coef_est_R.csv", row.names = F)

std_img$data$cortex_left[,1] = est_unweighted_L[1,]
std_img$data$cortex_right[,1] = est_unweighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_intercept_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = est_weighted_L[1,]
std_img$data$cortex_right[,1] = est_weighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_intercept_weighted.dtseries.nii")

std_img$data$cortex_left[,1] = est_unweighted_L[2,]
std_img$data$cortex_right[,1] = est_unweighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = est_weighted_L[2,]
std_img$data$cortex_right[,1] = est_weighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_weighted.dtseries.nii")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Z-statistics
zstat_unweighted_L = est_unweighted_L/sd_est_unweighted_L
zstat_weighted_L = est_weighted_L/sd_est_weighted_L
zstat_unweighted_R = est_unweighted_R/sd_est_unweighted_R
zstat_weighted_R = est_weighted_R/sd_est_weighted_R

std_img$data$cortex_left[,1] = zstat_unweighted_L[1,]
std_img$data$cortex_right[,1] = zstat_unweighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/zstat_intercept_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = zstat_weighted_L[1,]
std_img$data$cortex_right[,1] = zstat_weighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/zstat_intercept_weighted.dtseries.nii")

std_img$data$cortex_left[,1] = zstat_unweighted_L[2,]
std_img$data$cortex_right[,1] = zstat_unweighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/zstat_slope_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = zstat_weighted_L[2,]
std_img$data$cortex_right[,1] = zstat_weighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/zstat_slope_weighted.dtseries.nii")



pval_unweighted_L = (1-pnorm(abs((est_unweighted_L/sd_est_unweighted_L)),0,1) )*2 # tail z-test
pval_weighted_L = (1-pnorm(abs((est_weighted_L/sd_est_weighted_L)),0,1) )*2 # tail z-test
pval_unweighted_R = (1-pnorm(abs((est_unweighted_R/sd_est_unweighted_R)),0,1) )*2 # tail z-test
pval_weighted_R = (1-pnorm(abs((est_weighted_R/sd_est_weighted_R)),0,1) )*2 # tail z-test

std_img$data$cortex_left[,1] = sd_est_unweighted_L[1,]
std_img$data$cortex_right[,1] = sd_est_unweighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/sd_est_intercept_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = sd_est_weighted_L[1,]
std_img$data$cortex_right[,1] = sd_est_weighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/sd_est_intercept_weighted.dtseries.nii")

std_img$data$cortex_left[,1] = sd_est_unweighted_L[2,]
std_img$data$cortex_right[,1] = sd_est_unweighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/sd_est_slope_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = sd_est_weighted_L[2,]
std_img$data$cortex_right[,1] = sd_est_weighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/sd_est_slope_weighted.dtseries.nii")

## P values adjustment BH
pval_unweighted_L_adjusted = pval_unweighted_L
pval_weighted_L_adjusted = pval_weighted_L
pval_unweighted_R_adjusted = pval_unweighted_R
pval_weighted_R_adjusted = pval_weighted_R

for (j in 1:nrow(pval_unweighted_L)){
  pval_unweighted_L_adjusted[j,] = p.adjust(pval_unweighted_L[j,], method = "holm")
  pval_weighted_L_adjusted[j,] = p.adjust(pval_weighted_L[j,], method = "holm")
  pval_unweighted_R_adjusted[j,] = p.adjust(pval_unweighted_R[j,], method = "holm")
  pval_weighted_R_adjusted[j,] = p.adjust(pval_weighted_R[j,], method = "holm")
}

std_img$data$cortex_left[,1] = pval_unweighted_L_adjusted[1,]
std_img$data$cortex_right[,1] = pval_unweighted_R_adjusted[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_intercept_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = pval_weighted_L_adjusted[1,]
std_img$data$cortex_right[,1] = pval_weighted_R_adjusted[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_intercept_weighted.dtseries.nii")

std_img$data$cortex_left[,1] = pval_unweighted_L_adjusted[2,]
std_img$data$cortex_right[,1] = pval_unweighted_R_adjusted[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_slope_unweighted.dtseries.nii")

std_img$data$cortex_left[,1] = pval_weighted_L_adjusted[2,]
std_img$data$cortex_right[,1] = pval_weighted_R_adjusted[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_slope_weighted.dtseries.nii")

# Save the original p-values
std_img$data$cortex_left[,1] = pval_unweighted_L[1,]
std_img$data$cortex_right[,1] = pval_unweighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_intercept_unweighted_original.dtseries.nii")

std_img$data$cortex_left[,1] = pval_weighted_L[1,]
std_img$data$cortex_right[,1] = pval_weighted_R[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_intercept_weighted_original.dtseries.nii")

std_img$data$cortex_left[,1] = pval_unweighted_L[2,]
std_img$data$cortex_right[,1] = pval_unweighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_slope_unweighted_original.dtseries.nii")

std_img$data$cortex_left[,1] = pval_weighted_L[2,]
std_img$data$cortex_right[,1] = pval_weighted_R[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/pval_slope_weighted_original.dtseries.nii")


# ---------------------------------  pvalue correction --------------------------------- #

est_unweighted_L_corrected = est_unweighted_L
est_weighted_L_corrected = est_weighted_L
est_unweighted_R_corrected = est_unweighted_R
est_weighted_R_corrected = est_weighted_R

# Threshold for p-value
pval_thres = 0.01
est_unweighted_L_corrected[pval_unweighted_L_adjusted > pval_thres] = 0
est_weighted_L_corrected[pval_weighted_L_adjusted > pval_thres] = 0
est_unweighted_R_corrected[pval_unweighted_R_adjusted > pval_thres] = 0
est_weighted_R_corrected[pval_weighted_R_adjusted > pval_thres] = 0

std_img$data$cortex_left[,1] = est_unweighted_L_corrected[1,]
std_img$data$cortex_right[,1] = est_unweighted_R_corrected[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_intercept_unweighted_corrected.dtseries.nii")

std_img$data$cortex_left[,1] = est_weighted_L_corrected[1,]
std_img$data$cortex_right[,1] = est_weighted_R_corrected[1,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_intercept_weighted_corrected.dtseries.nii")

std_img$data$cortex_left[,1] = est_unweighted_L_corrected[2,]
std_img$data$cortex_right[,1] = est_unweighted_R_corrected[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_unweighted_corrected.dtseries.nii")

std_img$data$cortex_left[,1] = est_weighted_L_corrected[2,]
std_img$data$cortex_right[,1] = est_weighted_R_corrected[2,]
write_cifti(std_img, "./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_weighted_corrected.dtseries.nii")



# ---------------------- Gordon 333 Parcellation ---------------------- #
extract_text_between_underscores <- function(row_names) {
  # Define a regular expression pattern to match text between underscores
  pattern <- "_(.*?)_"
  
  # Initialize an empty vector to store the extracted text
  extracted_text <- character(length(row_names))
  
  # Loop through each row name and extract the text
  for (i in 1:length(row_names)) {
    # Use regmatches and regexpr to find and extract the text between underscores
    matches <- regmatches(row_names[i], gregexpr(pattern, row_names[i]))
    
    # Check if any matches were found
    if (length(matches[[1]]) > 0) {
      # Store the first match (text between the first pair of underscores)
      extracted_text[i] <- matches[[1]][1]
    } else {
      # Store NA if no match was found
      extracted_text[i] <- NA
    }
  }
  
  # Return the vector of extracted text with underscores removed
  extracted_text <- gsub("_", "", extracted_text)
  
  return(extracted_text)
}

# Read the Gordon 333 parcellation
gordon_atlas = read_cifti("./Parcellations/GordonParcellation/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii")
gordon_atlas_labels = gordon_atlas$meta$cifti$labels$`333Cort+SubCort`

est_slope_weighted_corrected <- read_cifti("./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_weighted_corrected.dtseries.nii")
est_slope_unweighted_corrected <- read_cifti("./RealDataAnalysis_CIFTI/nBack_gfactor/est_slope_unweighted_corrected.dtseries.nii")
weight_L_sig_regions_ids = unique(gordon_atlas$data$cortex_left[est_slope_weighted_corrected$data$cortex_left[,1] != 0])
weight_R_sig_regions_ids = unique(gordon_atlas$data$cortex_right[est_slope_weighted_corrected$data$cortex_right[,1] != 0])
unweight_L_sig_regions_ids = unique(gordon_atlas$data$cortex_left[est_slope_unweighted_corrected$data$cortex_left[,1] != 0])
unweight_R_sig_regions_ids = unique(gordon_atlas$data$cortex_right[est_slope_unweighted_corrected$data$cortex_right[,1] != 0])

weight_sig_regions_ids = union(weight_L_sig_regions_ids, weight_R_sig_regions_ids) # Obtain the union of significant regions
unweight_sig_regions_ids = union(unweight_L_sig_regions_ids, unweight_R_sig_regions_ids) # Obtain the union of significant regions
sig_regions_weighted = gordon_atlas_labels[(gordon_atlas_labels$Key %in% weight_sig_regions_ids),] # Extract the significant regions
sig_network_weighted = unique(extract_text_between_underscores(row_name = rownames(sig_regions_weighted))) # Extract the network names
sig_regions_unweighted = gordon_atlas_labels[(gordon_atlas_labels$Key %in% unweight_sig_regions_ids),]
sig_network_unweighted = unique(extract_text_between_underscores(row_name = rownames(sig_regions_unweighted)))

# Print out the significant networks
cat("Networks in weighted but not in unweighted: ", setdiff(sig_network_weighted, sig_network_unweighted), "\n")
cat("Networks in unweighted but not in weighted: ", setdiff(sig_network_unweighted, sig_network_weighted), "\n")
cat("Networks both in unweighted: ", intersect(sig_network_unweighted, sig_network_weighted), "\n")


###################################################################################################################

# Function to summarize the differences between two sets of counts
summarize_difference <- function(weighted_counts, unweighted_counts, ids) {
  # Ensure input vectors are of the same length
  if (length(weighted_counts) != length(unweighted_counts) || length(weighted_counts) != length(ids)) {
    stop("Input vectors must be of the same length.")
  }

  # Calculate the differences
  differences <- unweighted_counts - weighted_counts

  # Combine IDs and differences
  results <- data.frame(ID = ids, Difference = differences)

  # Sort results by differences in descending order
  sorted_results <- results[order(-results$Difference), ]

  return(sorted_results)
}


# Count the number of significant vertices in each region
region_vertices_count = function(parc, est_slope_weighted_corrected, est_slope_unweighted_corrected, gordon_atlas_labels){
  library(dplyr)
  if(sum(parc$meta$cortex$medial_wall_mask$left) > nrow(parc$data$cortex_left)){
    parc$data$cortex_left = parc$data$cortex_left[est_slope_unweighted_corrected$meta$cortex$medial_wall_mask$left,,drop=F]
  }
  if(sum(parc$meta$cortex$medial_wall_mask$right) > nrow(parc$data$cortex_right)){
    parc$data$cortex_right = parc$data$cortex_right[est_slope_unweighted_corrected$meta$cortex$medial_wall_mask$right,,drop=F]  
  }
  
  n_parc = length(unique(parc$meta$cifti$labels[[1]][,"Key"]))-1
  # print(unique(parc$meta$cifti$labels$parcels[,"Key"]))
  res = data.frame(region_id = 1:n_parc, unweighted = rep(NA, n_parc), weighted = rep(NA, n_parc))
  
  for(i in 1:n_parc){
    # reset sample cii
    # sample_cii = est_slope_weighted_corrected
    # sample_cii$data$cortex_left[,1] = 0
    # sample_cii$data$cortex_right[,1] = 0
    # sample_cii_uw = sample_cii
    # sample_cii_w = sample_cii
    # 
    region_i_ids_L = (parc$data$cortex_left[,1] == i)
    region_i_ids_R = (parc$data$cortex_right[,1] == i)
    # 
    # 
    # # Write CII files
    # sample_cii_w$data$cortex_left[region_i_ids_L,1] = est_slope_weighted_corrected$data$cortex_left[region_i_ids_L,1]
    # sample_cii_w$data$cortex_right[region_i_ids_R,1] = est_slope_weighted_corrected$data$cortex_right[region_i_ids_R,1]
    # 
    # sample_cii_uw$data$cortex_left[region_i_ids_L,1] = est_slope_unweighted_corrected$data$cortex_left[region_i_ids_L,1]
    # sample_cii_uw$data$cortex_right[region_i_ids_R,1] = est_slope_unweighted_corrected$data$cortex_right[region_i_ids_R,1]
    # 
    # writecii(sample_cii_w, sprintf("./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/region_%d_w", i))
    # writecii(sample_cii_uw, sprintf("./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/region_%d_uw", i))
    # 
    # Update summary table
    
    # sum of significant voxel (weighted and unweighted)
    weighted_L_sig = (est_slope_weighted_corrected$data$cortex_left[region_i_ids_L,1] != 0)
    weighted_R_sig = (est_slope_weighted_corrected$data$cortex_right[region_i_ids_R,1] != 0)
    
    unweighted_L_sig = (est_slope_unweighted_corrected$data$cortex_left[region_i_ids_L,1] != 0)
    unweighted_R_sig = (est_slope_unweighted_corrected$data$cortex_right[region_i_ids_R,1] != 0)
    
    
    res[i, "weighted"] <- sum(est_slope_weighted_corrected$data$cortex_left[region_i_ids_L,1] != 0) + sum(est_slope_weighted_corrected$data$cortex_right[region_i_ids_R,1] != 0)
    res[i, "unweighted"] <- sum(est_slope_unweighted_corrected$data$cortex_left[region_i_ids_L,1] != 0) + sum(est_slope_unweighted_corrected$data$cortex_right[region_i_ids_R,1] != 0)
    res[i, "num_diff_vertices"] <- (sum(region_i_ids_L) - sum(weighted_L_sig == unweighted_L_sig)) + (sum(region_i_ids_R) - sum(weighted_R_sig == unweighted_R_sig))
    res[i, "pct_diff_vertices"] <- ((sum(region_i_ids_L) - sum(weighted_L_sig == unweighted_L_sig)) + (sum(region_i_ids_R) - sum(weighted_R_sig == unweighted_R_sig)))/(sum(region_i_ids_L) + sum(region_i_ids_R))
    res[i, "region_name"] <- rownames(parc$meta$cifti$labels[[1]])[parc$meta$cifti$labels[[1]][,"Key"] == i]
  }
  
  
  saveRDS(res, "./RealDataAnalysis_CIFTI/nBack_gfactor/region_vertices_count.rds")
  return(res)
}
gordon_atlas = read_cifti("./Parcellations/GordonParcellation/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii")
gordon_atlas_labels = gordon_atlas$meta$cifti$labels$`333Cort+SubCort`
# Schaefer_100 = load_parc("Schaefer_100")
# region_vertex_tab = region_vertices_count(Schaefer_100,est_slope_weighted_corrected, est_slope_unweighted_corrected)
region_vertex_tab = region_vertices_count(gordon_atlas,est_slope_weighted_corrected, est_slope_unweighted_corrected)
region_vertex_tab %>% arrange(desc(pct_diff_vertices)) %>% filter(num_diff_vertices >= 30)
region_vertex_tab %>% arrange(desc(pct_diff_vertices)) 
write.csv(region_vertex_tab, "./RealDataAnalysis_CIFTI/nBack_gfactor/region_vertex_tab_summary.csv", row.names = F)
#####################################
library(grid)
simpleColormap <- function(n) {
  # Interpolate between #F4F1DE and #E07A5F
  colors <- colorRampPalette(c("#F4F1DE", "#E07A5F"))(n)
  return(colors)
}
# Use the simpleColormap function to generate the colors
n <- 100  # Choose the number of color gradations you want
colormap <- simpleColormap(n)


# 
# default_ids <- c(48:88)
# FrontoParietal_ids <- c(89:120)
plotNetwork <- function(self_defined_ids, est_slope_unweighted_corrected, 
                        est_slope_weighted_corrected, parc, colormap, network_name = NULL, basedir = "./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/visualize_sig_regions/") {
  
  network_to_plot <- self_defined_ids
  parc$data$cortex_left = parc$data$cortex_left[est_slope_unweighted_corrected$meta$cortex$medial_wall_mask$left,,drop=F]
  parc$data$cortex_right = parc$data$cortex_right[est_slope_unweighted_corrected$meta$cortex$medial_wall_mask$right,,drop=F]
  
  cii_to_plot_uw <- est_slope_unweighted_corrected
  cii_to_plot_w <- est_slope_weighted_corrected
  network_ids_L <- (parc$data$cortex_left[,1] %in% network_to_plot)
  network_ids_R <- (parc$data$cortex_right[,1] %in% network_to_plot)
  
  # cii_to_plot_uw$data$cortex_left[cii_to_plot_uw$data$cortex_left != 0, 1] <- 1
  # cii_to_plot_uw$data$cortex_right[cii_to_plot_uw$data$cortex_right != 0, 1] <- 1
  cii_to_plot_uw$data$cortex_left[!network_ids_L,1] <- NA
  cii_to_plot_uw$data$cortex_right[!network_ids_R,1] <- NA
  
  # cii_to_plot_w$data$cortex_left[cii_to_plot_w$data$cortex_left != 0,1] <- 1
  # cii_to_plot_w$data$cortex_right[cii_to_plot_w$data$cortex_right != 0,1] <- 1
  cii_to_plot_w$data$cortex_left[!network_ids_L,1] <- NA
  cii_to_plot_w$data$cortex_right[!network_ids_R,1] <- NA
  
  if(is.null(network_name)){
    plot(cii_to_plot_uw, colors = colormap, zlim = c(-0.03,0.03),
         # NA_color = "#3D405B",
         fname = sprintf("%sregion_%d_unweighted", basedir, self_defined_ids))
    plot(cii_to_plot_w , colors = colormap, zlim = c(-0.03,0.03),
         # NA_color = "#3D405B",
         fname = sprintf("%sregion_%d_weighted", basedir, self_defined_ids))
  }else{
    plot(cii_to_plot_uw, colors = colormap, zlim = c(-0.03,0.03),
         # NA_color = "#3D405B",
         fname = sprintf("%sregion_%s_unweighted", basedir, self_defined_ids,network_name))
    plot(cii_to_plot_w , colors = colormap, zlim = c(-0.03,0.03),
         # NA_color = "#3D405B", 
         fname = sprintf("%sregion_%s_weighted", basedir, self_defined_ids,network_name))
    
  }
}


#####################################
library(grid)

# Define the color template
color_template <- c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#e2e2e2", "#e1a692", "#de6e56", "#e14b31", "#c23728")

# Number of color gradations you want
n <- 100

# Create a continuous colormap
continuous_colormap <- colorRampPalette(color_template)(n)

cii_to_plot_uw <- est_slope_unweighted_corrected
cii_to_plot_w <- est_slope_weighted_corrected


plot(cii_to_plot_uw, colors = continuous_colormap, zlim = c(-0.6,0.6),
     # NA_color = "#F4F1DE",
     fname = sprintf("%sregion_all_unweighted_raw", "./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/visualize_sig_regions/"))
plot(cii_to_plot_w , colors = continuous_colormap, zlim = c(-0.6,0.6),
     # NA_color = "#F4F1DE",
     fname = sprintf("%sregion_all_weighted_raw", "./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/visualize_sig_regions/"  ))

plot(cii_to_plot_uw, colors = continuous_colormap, zlim = c(-0.6,0.6))
plot(cii_to_plot_w , colors = continuous_colormap, zlim = c(-0.6,0.6))


###### Binary Mask ###########
network_ids_L <- unique(gordon_atlas$data$cortex_left[,1])
network_ids_R <- unique(gordon_atlas$data$cortex_right[,1])

cii_to_plot_uw$data$cortex_left[cii_to_plot_uw$data$cortex_left != 0, 1] <- 1
cii_to_plot_uw$data$cortex_right[cii_to_plot_uw$data$cortex_right != 0, 1] <- 1
cii_to_plot_uw$data$cortex_left[cii_to_plot_uw$data$cortex_left == 0,1] <- NA
cii_to_plot_uw$data$cortex_right[cii_to_plot_uw$data$cortex_right == 0,1] <- NA

cii_to_plot_w$data$cortex_left[cii_to_plot_w$data$cortex_left != 0,1] <- 1
cii_to_plot_w$data$cortex_right[cii_to_plot_w$data$cortex_right != 0,1] <- 1
cii_to_plot_w$data$cortex_left[cii_to_plot_w$data$cortex_left == 0,1] <- NA
cii_to_plot_w$data$cortex_right[cii_to_plot_w$data$cortex_right == 0,1] <- NA

plot(cii_to_plot_uw, colors = continuous_colormap, zlim = c(0,1),
     # NA_color = "#F4F1DE",
     fname = sprintf("%sregion_all_unweighted_raw_mask", "./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/visualize_sig_regions/"))
plot(cii_to_plot_w , colors = continuous_colormap, zlim = c(0,1), 
     # NA_color = "#F4F1DE",
     fname = sprintf("%sregion_all_weighted_raw_mask", "./RealDataAnalysis_CIFTI/nBack_gfactor/regionwise_sig_plot/visualize_sig_regions/"  ))

plot(cii_to_plot_uw, colors = continuous_colormap, zlim = c(0,1))
plot(cii_to_plot_w , colors = continuous_colormap, zlim = c(0,1))

###### Plot test statistics ###########
cii_to_plot_uw$data$cortex_left[,1] = zstat_unweighted_L[2,]
cii_to_plot_uw$data$cortex_right[,1] = zstat_unweighted_R[2,]
cii_to_plot_w$data$cortex_left[,1] = zstat_weighted_L[2,]
cii_to_plot_w$data$cortex_right[,1] = zstat_weighted_R[2,]

plot(cii_to_plot_uw, colors = BayesGPfit::GP.create.cols(), zlim = c(-12, 12),
 # NA_color = "#F4F1DE",
  fname = sprintf("%szstat_unweighted",
   "./RealDataAnalysis_CIFTI/nBack_gfactor/"))

plot(cii_to_plot_w, colors = BayesGPfit::GP.create.cols(), zlim =c(-12, 12),
  # NA_color = "#F4F1DE",
    fname = sprintf("%szstat_weighted",
    "./RealDataAnalysis_CIFTI/nBack_gfactor/"))



############## Plot difference between weighted and unwegithed ###################
diff_uw_w_cii = cii_to_plot_uw
diff_uw_w_cii$data$cortex_left[, 1] <- NA
diff_uw_w_cii$data$cortex_right[, 1] <- NA
diff_uw_w_cii$data$cortex_left[(est_slope_unweighted_corrected$data$cortex_left != 0) & (est_slope_weighted_corrected$data$cortex_left == 0), 1] <- 1
diff_uw_w_cii$data$cortex_left[(est_slope_unweighted_corrected$data$cortex_left == 0) & (est_slope_weighted_corrected$data$cortex_left != 0), 1] <- 1
diff_uw_w_cii$data$cortex_right[(est_slope_unweighted_corrected$data$cortex_right != 0) & (est_slope_weighted_corrected$data$cortex_right == 0), 1] <- 1
diff_uw_w_cii$data$cortex_right[(est_slope_unweighted_corrected$data$cortex_right == 0) & (est_slope_weighted_corrected$data$cortex_right != 0), 1] <- 1

plot(diff_uw_w_cii , colors = continuous_colormap, zlim = c(0,1), legend_embed = F,
     fname = "./RealDataAnalysis_CIFTI/nBack_gfactor/plot_figures/difference_without_Z")

