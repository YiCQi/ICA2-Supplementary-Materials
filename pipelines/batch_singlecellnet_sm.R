setwd("/mnt/volume1/qyc/")
source("auto-annot/pipelines/cross_validation.R")
source("auto-annot/pipelines/evaluate.R")
source("auto-annot/pipelines/scripts/single_cell_net_sm.R")

pbmc_10x_v2_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/10Xv2/10Xv2_pbmc1.csv"
pbmc_10x_v2_dir <- dirname(pbmc_10x_v2_path)
pbmc_10x_v3_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.csv"
pbmc_10x_v3_dir <- dirname(pbmc_10x_v3_path)
pbmc_cel_seq_path <- "/mnt/volume1/qyc/data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.csv"
pbmc_cel_seq_dir <- dirname(pbmc_cel_seq_path)
lung_cancer_10x_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/10x_5cl/10x_5cl_data.csv"
lung_cancer_10x_dir <- dirname(lung_cancer_10x_path)
lung_cancer_cel_seq_path <- "/mnt/volume1/qyc/data/Intra-dataset/CellBench/CelSeq2_5cl/CelSeq2_5cl_data.csv"
lung_cancer_cel_seq_dir <- dirname(lung_cancer_cel_seq_path)
pancreas_cel_seq_path <- "/mnt/volume1/qyc/data/Intra-dataset/Pancreatic_data/Muraro/Filtered_Muraro_HumanPancreas_data.csv"
pancreas_cel_seq_dir <- dirname(pancreas_cel_seq_path)
pancreas_smart_seq_path <- "/mnt/volume1/qyc/data/Intra-dataset/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.csv"
pancreas_smart_seq_dir <- dirname(pancreas_smart_seq_path)
brain_allen_10x_path <- "/mnt/volume1/qyc/data/Intra-dataset/AMB/Filtered_mouse_allen_brain_data.csv"
brain_allen_10x_dir <- dirname(brain_allen_10x_path)

results_dir <- "/mnt/volume1/qyc/auto-annot/results"

sample_list <- list(
  pbmc_10x_v2_path,
  pbmc_10x_v3_path,
  pbmc_cel_seq_path,
  lung_cancer_10x_path,
  lung_cancer_cel_seq_path,
  pancreas_cel_seq_path,
  pancreas_smart_seq_path,
  brain_allen_10x_path
)

dir_list <- list(
  pbmc_10x_v2_dir,
  pbmc_10x_v3_dir,
  pbmc_cel_seq_dir,
  lung_cancer_10x_dir,
  lung_cancer_cel_seq_dir,
  pancreas_cel_seq_dir,
  pancreas_smart_seq_dir,
  brain_allen_10x_dir
)

for (i in 1:length(sample_list)) {
  sample_path <- sample_list[[i]]
  dir_path <- dir_list[[i]]

  labels_path <- paste0(dir_path, "/Labels.csv")
  if (!file.exists(labels_path)) {
    labels_path <- paste0(tools::file_path_sans_ext(sample_path), "Labels.csv")
  }

  # Cross_Validation(labels_path, 1, dir_path)

  sample_dir <- paste0(results_dir, "/singleCellNet_sm/", sub("\\.csv$", "", basename(sample_path)), "/")
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  run_singleCellNet(sample_path, labels_path, paste0(dir_path, "/CV_folds.RData"), sample_dir)

  # result <- evaluate(paste0(results_dir, "/singleCellNet/singleCellNet_True_Labels.csv"), paste0(results_dir, "/singleCellNet/singleCellNet_Pred_Labels.csv"))
  # result <- evaluate(paste0(results_dir, "/SingleR/SingleR_True_Labels.csv"), paste0(results_dir, "/SingleR/SingleR_Pred_Labels.csv"))
  # print(result)
}