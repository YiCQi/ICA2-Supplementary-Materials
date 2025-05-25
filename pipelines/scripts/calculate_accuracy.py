import pandas as pd

# Read the CSV files
pred_labels = pd.read_csv('/mnt/volume1/qyc/auto-annot/results/RF/10Xv2_pbmc1/RF_Pred_Labels.csv')
true_labels = pd.read_csv('/mnt/volume1/qyc/auto-annot/results/RF/10Xv2_pbmc1/RF_True_Labels.csv')

# Calculate the number of matching labels
# col1, from row2
matches = (pred_labels.iloc[1:, 0] == true_labels.iloc[1:, 0]).sum()
total = len(pred_labels)

# Calculate accuracy
accuracy = (matches / total) * 100

print(f"Total number of cells: {total}")
print(f"Number of correct predictions: {matches}")
print(f"Accuracy: {accuracy:.2f}%") 