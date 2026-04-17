#!/usr/bin/env Rscript

# This script models gene detection rates across sample clusters:
#
# (1) Loads gene-level feature data (long format) and sample-to-cluster assignments
#     generated in previous pipeline steps
#
# (2) Assigns each sample to a cluster and computes gene detection rates
#     within each cluster (fraction of samples in which a gene is detected)
#
# (3) Fits a one-dimensional Gaussian Mixture Model (GMM) to the distribution
#     of detection rates to identify threshold(s) separating low- and high-detection genes
#
# (4) Visualizes the fitted mixture model and detected threshold(s)
#
# (5) Saves the fitted model and the resulting plot for downstream analysis
#
# Input:
# - Long-format feature table (raw_ft_long.tsv)
# - Sample-to-cluster mapping (sample_kit_cluster_map.tsv)
#
# Output:
# - GMM model object (GMM_det_rate.RDS)
# - Detection rate distribution plot with thresholds (dpGMM_detection_rate_plot.pdf)
#
# Note:
# Input data are not included in the repository and must be provided by the user.

library(data.table)
library(ggplot2)
library(httpgd)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)


dir <- "/path/to/results/dir" # directory containing sudirectories from Variant_post_processing and Variant_to_gene steps
dir.create(file.path(dir, "Results"), showWarnings = FALSE)

# Load cluster assignments and long-format feature table from previous step
ft_long <- fread(file.path(dir, "Results", "raw_ft_long.tsv"))
clusters <- fread( file.path(dir, "Results", "sample_kit_cluster_map.tsv") )

# add clustering 
ft_long$Group <- as.factor(clusters$Cluster[match(ft_long$Sample, clusters$Sample)])

# compute detecion rate for genes in each cluster
Genes <- unique(ft_long$Gene)
n_total <- unique(ft_long[, .(Sample, Group)])[ , .N, by = Group]
setnames(n_total, "N", "Total")

det_rate <- ft_long %>% group_by(Group) %>% 
                        count(Gene) %>% 
                        mutate(Detection_rate = n / n_total$Total[match(Group, n_total$Group)]) %>% 
                        ungroup()

# Fit Gaussian mixture model to detection rate distribution
library(dpGMM)
custom.settings <- GMM_1D_opts

mix_test <- runGMM(det_rate$Detection_rate, opts = custom.settings)
saveRDS(mix_test, file.path(dir, "Results", "GMM_det_rate.RDS"))


thr <- mix_test[["threshold"]]
p <- mix_test[["fig"]] + 
  geom_text(
    data = data.frame(x = thr, y = 15), 
    aes(x = x, y = y, label = round(x, 3)),
    angle = 90, vjust = -0.5, hjust = 0, color = "red", size = 3
  ) + labs(x = "Detection rate")


pdf(file.path(dir, "Results", "Figures", "dpGMM_detection_rate_plot.pdf"), width = 6, height = 4)
print(p)
dev.off()
