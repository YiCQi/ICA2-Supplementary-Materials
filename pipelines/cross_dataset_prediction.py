import os
import numpy as np
import pandas as pd
import time as tm
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from scvi.dataset import CsvDataset
from scvi.models import SCANVI
from scvi.inference import SemiSupervisedTrainer
import torch

# Set random seeds
torch.manual_seed(666)
np.random.seed(666)

# Define paths
base_dir = "/mnt/volume1/qyc"
# train_data_path = os.path.join(base_dir, "data/Intra-dataset/CellBench/10x_5cl/10x_5cl_data.csv")
# train_labels_path = os.path.join(base_dir, "data/Intra-dataset/CellBench/10x_5cl/Labels.csv")
# test_data_path = os.path.join(base_dir, "data/Intra-dataset/CellBench/CelSeq2_5cl/CelSeq2_5cl_data.csv")
# test_labels_path = os.path.join(base_dir, "data/Intra-dataset/CellBench/CelSeq2_5cl/Labels.csv")
# train_data_path = os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.csv")
# train_labels_path = os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Segerstolpe/Labels.csv")
# test_data_path = os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Muraro/Filtered_Muraro_HumanPancreas_data.csv")
# test_labels_path = os.path.join(base_dir, "data/Intra-dataset/Pancreatic_data/Muraro/Labels.csv")
train_data_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.csv")
train_labels_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1Labels.csv")
test_data_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.csv")
test_labels_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1Labels.csv")

# Create output directories
results_dir = os.path.join(base_dir, "auto-annot/result_across")
rf_dir = os.path.join(results_dir, "RF")
svm_dir = os.path.join(results_dir, "SVM")
scvi_dir = os.path.join(results_dir, "scVI")

for dir_path in [rf_dir, svm_dir, scvi_dir]:
    os.makedirs(dir_path, exist_ok=True)

def preprocess_data(data, labels=None):
    # Ensure all values are non-negative
    data[data < 0] = 0

    # Remove all-zero genes and cells
    data = data.loc[(data.sum(axis=1) > 0), (data.sum(axis=0) > 0)]
    
    # If labels are provided, ensure they match the data
    if labels is not None:
        common_cells = list(set(data.index).intersection(set(labels.index)))
        data = data.loc[common_cells]
        labels = labels.loc[common_cells]
        return data, labels

    return data

def get_common_genes(train_data, test_data):
    """Get common genes between train and test data"""
    train_genes = set(train_data.columns)
    test_genes = set(test_data.columns)
    common_genes = list(train_genes.intersection(test_genes))
    print(f"Number of common genes: {len(common_genes)}")
    return common_genes

def run_rf_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, output_dir):
    # Read data
    train_data = pd.read_csv(train_data_path, index_col=0)
    train_labels = pd.read_csv(train_labels_path, index_col=0)
    test_data = pd.read_csv(test_data_path, index_col=0)
    test_labels = pd.read_csv(test_labels_path, index_col=0)
    
    # Get common genes
    common_genes = get_common_genes(train_data, test_data)
    
    # Select common genes
    train_data = train_data[common_genes]
    test_data = test_data[common_genes]
    
    # Preprocess data and align labels
    train_data, train_labels = preprocess_data(train_data, train_labels)
    test_data, test_labels = preprocess_data(test_data, test_labels)
    
    # Normalize (CPM then log1p)
    train_data = train_data.div(train_data.sum(axis=1), axis=0) * 1e6
    train_data = np.log1p(train_data)
    test_data = test_data.div(test_data.sum(axis=1), axis=0) * 1e6
    test_data = np.log1p(test_data)
    
    # Initialize classifier
    classifier = RandomForestClassifier(n_estimators=50, random_state=666)
    
    # Training phase
    start_time = tm.time()
    classifier.fit(train_data, train_labels['x'])
    training_time = tm.time() - start_time
    
    # Testing phase
    start_time = tm.time()
    predicted = classifier.predict(test_data)
    testing_time = tm.time() - start_time
    
    # Save results
    pd.DataFrame(test_labels['x']).to_csv(os.path.join(output_dir, 'RF_True_Labels.csv'), index=False)
    pd.DataFrame(predicted).to_csv(os.path.join(output_dir, 'RF_Pred_Labels.csv'), index=False)
    pd.DataFrame([training_time]).to_csv(os.path.join(output_dir, 'RF_Training_Time.csv'), index=False)
    pd.DataFrame([testing_time]).to_csv(os.path.join(output_dir, 'RF_Testing_Time.csv'), index=False)

def run_svm_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, output_dir):
    # Read data
    train_data = pd.read_csv(train_data_path, index_col=0)
    train_labels = pd.read_csv(train_labels_path, index_col=0)
    test_data = pd.read_csv(test_data_path, index_col=0)
    test_labels = pd.read_csv(test_labels_path, index_col=0)
    
    # Get common genes
    common_genes = get_common_genes(train_data, test_data)
    
    # Select common genes
    train_data = train_data[common_genes]
    test_data = test_data[common_genes]
    
    # Preprocess data and align labels
    train_data, train_labels = preprocess_data(train_data, train_labels)
    test_data, test_labels = preprocess_data(test_data, test_labels)
    
    # Normalize (CPM then log1p)
    train_data = train_data.div(train_data.sum(axis=1), axis=0) * 1e6
    train_data = np.log1p(train_data)
    test_data = test_data.div(test_data.sum(axis=1), axis=0) * 1e6
    test_data = np.log1p(test_data)
    
    # Initialize classifier
    classifier = LinearSVC(random_state=666)
    
    # Training phase
    start_time = tm.time()
    classifier.fit(train_data, train_labels['x'])
    training_time = tm.time() - start_time
    
    # Testing phase
    start_time = tm.time()
    predicted = classifier.predict(test_data)
    testing_time = tm.time() - start_time
    
    # Save results
    pd.DataFrame(test_labels['x']).to_csv(os.path.join(output_dir, 'SVM_True_Labels.csv'), index=False)
    pd.DataFrame(predicted).to_csv(os.path.join(output_dir, 'SVM_Pred_Labels.csv'), index=False)
    pd.DataFrame([training_time]).to_csv(os.path.join(output_dir, 'SVM_Training_Time.csv'), index=False)
    pd.DataFrame([testing_time]).to_csv(os.path.join(output_dir, 'SVM_Testing_Time.csv'), index=False)

def run_scvi_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, output_dir):
    # Read data
    train_data = pd.read_csv(train_data_path, index_col=0)
    train_labels = pd.read_csv(train_labels_path, index_col=0)
    test_data = pd.read_csv(test_data_path, index_col=0)
    test_labels = pd.read_csv(test_labels_path, index_col=0)
    
    # Get common genes
    common_genes = get_common_genes(train_data, test_data)
    
    # Select common genes
    train_data = train_data[common_genes]
    test_data = test_data[common_genes]
    
    # Preprocess data and align labels
    train_data, train_labels = preprocess_data(train_data, train_labels)
    test_data, test_labels = preprocess_data(test_data, test_labels)
    
    # Normalize (CPM then log1p)
    train_data = train_data.div(train_data.sum(axis=1), axis=0) * 1e6
    train_data = np.log1p(train_data)
    test_data = test_data.div(test_data.sum(axis=1), axis=0) * 1e6
    test_data = np.log1p(test_data)
    
    # Save data for scVI
    train_data.to_csv(os.path.join(output_dir, 'Data_scvi_train.csv'))
    train_labels.to_csv(os.path.join(output_dir, 'Labels_scvi_train.csv'))
    test_data.to_csv(os.path.join(output_dir, 'Data_scvi_test.csv'))
    test_labels.to_csv(os.path.join(output_dir, 'Labels_scvi_test.csv'))
    
    # Create datasets
    train_dataset = CsvDataset('Data_scvi_train.csv', save_path=output_dir, 
                             sep=",", labels_file="Labels_scvi_train.csv", gene_by_cell=False)
    test_dataset = CsvDataset('Data_scvi_test.csv', save_path=output_dir, 
                            sep=",", labels_file="Labels_scvi_test.csv", gene_by_cell=False)
    
    # Initialize model
    scanvi = SCANVI(train_dataset.nb_genes, train_dataset.n_batches, train_dataset.n_labels)
    trainer = SemiSupervisedTrainer(scanvi, train_dataset, frequency=5)
    
    # Training phase
    start_time = tm.time()
    trainer.train(n_epochs=200)
    training_time = tm.time() - start_time
    
    # Testing phase
    start_time = tm.time()
    y_true, y_pred = trainer.unlabelled_set.compute_predictions()
    testing_time = tm.time() - start_time
    
    # Save results
    pd.DataFrame(y_true).to_csv(os.path.join(output_dir, 'scVI_True_Labels.csv'), index=False)
    pd.DataFrame(y_pred).to_csv(os.path.join(output_dir, 'scVI_Pred_Labels.csv'), index=False)
    pd.DataFrame([training_time]).to_csv(os.path.join(output_dir, 'scVI_Training_Time.csv'), index=False)
    pd.DataFrame([testing_time]).to_csv(os.path.join(output_dir, 'scVI_Testing_Time.csv'), index=False)

# Run predictions
print("Running RF cross-dataset prediction...")
run_rf_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, rf_dir)

print("Running SVM cross-dataset prediction...")
run_svm_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, svm_dir)

# print("Running scVI cross-dataset prediction...")
# run_scvi_cross(train_data_path, train_labels_path, test_data_path, test_labels_path, scvi_dir) 