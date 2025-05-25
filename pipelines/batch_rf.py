import os
from scripts.RF import run_RF

# 设置根路径和脚本路径（假设已导入 run_singleR 函数）
base_dir = "/mnt/volume1/qyc"
results_dir = os.path.join(base_dir, "auto-annot/results")

# 样本路径与目录列表
sample_list = [
    # os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv2/10Xv2_pbmc1.csv"),
    # os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.csv"),
    # os.path.join(base_dir, "data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.csv"),
    # os.path.join(base_dir, "data/Intra-dataset/CellBench/10x_5cl/10x_5cl_data.csv"),
    # os.path.join(base_dir, "data/Intra-dataset/CellBench/CelSeq2_5cl/CelSeq2_5cl_data.csv"),
    # os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Muraro/Filtered_Muraro_HumanPancreas_data.csv"),
    # os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.csv"),
    os.path.join(base_dir, "data/Intra-dataset/AMB/Filtered_mouse_allen_brain_data.csv"),
]

dir_list = [os.path.dirname(path) for path in sample_list]

# 主循环
for sample_path, dir_path in zip(sample_list, dir_list):
    # 构造标签路径
    labels_path = os.path.join(dir_path, "Labels.csv")
    if not os.path.exists(labels_path):
        labels_path = os.path.splitext(sample_path)[0] + "Labels.csv"

    # sample 子目录
    sample_name = os.path.basename(sample_path).replace(".csv", "")
    sample_dir = os.path.join(results_dir, "RF", sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    folds_path = os.path.join(dir_path, "CV_folds.RData")

    # 运行 SVM
    run_RF(sample_path, labels_path, folds_path, sample_dir)
