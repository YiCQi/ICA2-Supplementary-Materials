import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

# Define tools and result directory
result_dir = "/mnt/volume1/qyc/auto-annot/result_across/"
tools = ["RF", "singleCellNet", "SingleR", "SVM"]
datasets = ["VISp"]  # Add more datasets if needed

# Load data
results = {}
for dataset in datasets:
    dataset_results = {}
    pred_labels_list = []
    true_labels_list = []
    time_list = []
    
    for tool in tools:
        pred_labels_list.append(pd.read_csv(f"{result_dir}{dataset}/{tool}/{tool}_Pred_Labels.csv"))
        true_labels_list.append(pd.read_csv(f"{result_dir}{dataset}/{tool}/{tool}_True_Labels.csv"))
        time_list.append(pd.read_csv(f"{result_dir}{dataset}/{tool}/{tool}_Total_Time.csv").iloc[0,0])

    # Function to calculate metrics
    def calculate_metrics(y_true, y_pred):
        # Calculate per-class F1 scores
        f1_scores = f1_score(y_true, y_pred, average=None)
        
        metrics = {
            'Accuracy': accuracy_score(y_true, y_pred),
            'Median F1': np.median(f1_scores),  # Changed from Macro F1 to Median F1
            'Weighted F1': f1_score(y_true, y_pred, average='weighted'),
            'Macro Precision': precision_score(y_true, y_pred, average='macro'),
            'Macro Recall': recall_score(y_true, y_pred, average='macro'),
            'ALM': np.mean(y_true == y_pred),  # Average Label Match
            'MTG': np.mean([len(set(y_true[y_pred == label])) for label in np.unique(y_pred)])  # Mean Top Gene
        }
        return metrics

    # Calculate metrics for each tool
    for i, tool in enumerate(tools):
        y_true = true_labels_list[i].iloc[:, 0]
        y_pred = pred_labels_list[i].iloc[:, 0]
        dataset_results[tool] = calculate_metrics(y_true, y_pred)
    
    results[dataset] = dataset_results

# Create DataFrames with results for each dataset
for dataset in datasets:
    results_df = pd.DataFrame(results[dataset]).T
    results_df['Time (s)'] = time_list
    
    # Display results
    print(f"\nPerformance Metrics for {dataset} Dataset:")
    print(results_df.round(4))
    
    # Plot confusion matrices
    plt.figure(figsize=(15, 10))
    for i, tool in enumerate(tools):
        plt.subplot(2, 2, i+1)
        cm = confusion_matrix(true_labels_list[i].iloc[:, 0], pred_labels_list[i].iloc[:, 0])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
        plt.title(f'Confusion Matrix - {tool}')
        plt.xlabel('Predicted')
        plt.ylabel('True')
    plt.tight_layout()
    plt.savefig(f'{dataset}_confusion_matrices.png')
    plt.close()
    
    # Plot performance metrics
    plt.figure(figsize=(12, 6))
    metrics_to_plot = ['Accuracy', 'Median F1', 'Weighted F1', 'Macro Precision', 'Macro Recall', 'ALM', 'MTG']
    results_df[metrics_to_plot].plot(kind='bar')
    plt.title(f'Performance Metrics Comparison - {dataset}')
    plt.xlabel('Tools')
    plt.ylabel('Score')
    plt.xticks(rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{dataset}_performance_metrics.png')
    plt.close()
    
    # Save results to CSV
    results_df.to_csv(f'{dataset}_performance_metrics_results.csv') 