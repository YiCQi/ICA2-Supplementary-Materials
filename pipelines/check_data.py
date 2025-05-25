import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define paths
base_dir = "/mnt/volume1/qyc"
train_data_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.csv")
train_labels_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1Labels.csv")
test_data_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.csv")
test_labels_path = os.path.join(base_dir, "data/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1Labels.csv")

# Read data
train_data = pd.read_csv(train_data_path, index_col=0)
train_labels = pd.read_csv(train_labels_path, index_col=0)
test_data = pd.read_csv(test_data_path, index_col=0)
test_labels = pd.read_csv(test_labels_path, index_col=0)

# Print basic information
print("Training data shape:", train_data.shape)
print("Test data shape:", test_data.shape)
print("\nTraining labels distribution:")
print(train_labels['x'].value_counts())
print("\nTest labels distribution:")
print(test_labels['x'].value_counts())

# Check for common genes
common_genes = set(train_data.columns).intersection(set(test_data.columns))
print("\nNumber of common genes:", len(common_genes))

# Check for zero values
print("\nPercentage of zeros in training data:", (train_data == 0).sum().sum() / (train_data.shape[0] * train_data.shape[1]))
print("Percentage of zeros in test data:", (test_data == 0).sum().sum() / (test_data.shape[0] * test_data.shape[1]))

# Check for negative values
print("\nNumber of negative values in training data:", (train_data < 0).sum().sum())
print("Number of negative values in test data:", (test_data < 0).sum().sum())

# Plot label distributions
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
train_labels['x'].value_counts().plot(kind='bar')
plt.title('Training Labels Distribution')
plt.xticks(rotation=45)

plt.subplot(1, 2, 2)
test_labels['x'].value_counts().plot(kind='bar')
plt.title('Test Labels Distribution')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('label_distributions.png')
plt.close()

# Check if test labels are a subset of training labels
train_unique_labels = set(train_labels['x'].unique())
test_unique_labels = set(test_labels['x'].unique())
print("\nLabels in training but not in test:", train_unique_labels - test_unique_labels)
print("Labels in test but not in training:", test_unique_labels - train_unique_labels)

# Check if the data is already normalized
print("\nTraining data statistics:")
print(train_data.describe())
print("\nTest data statistics:")
print(test_data.describe()) 