from scvi.dataset import CsvDataset
import os
import numpy as np
import pandas as pd
from scvi.models import SCANVI
from scvi.inference import SemiSupervisedTrainer
import time as tm
import rpy2.robjects as robjects
import torch

torch.manual_seed(666)
np.random.seed(666)

def preprocess_data(data):
    # Ensure all values are non-negative
    data[data < 0] = 0

    # Remove all-zero genes and cells
    data = data.loc[(data.sum(axis=1) > 0), (data.sum(axis=0) > 0)]

    # Normalize (CPM then log1p)
    data = data.div(data.sum(axis=1), axis=0) * 1e6
    data = np.log1p(data)

    return data

def run_scVI(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run scVI
    Wrapper script to run scVI on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
    
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1 
    
    # Convert R list to Python list first, then to numpy array
    test_ind = [np.array(x, dtype='int') for x in robjects.r['Test_Idx']]
    train_ind = [np.array(x, dtype='int') for x in robjects.r['Train_Idx']]

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    data = preprocess_data(data)
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep] 
    
    # Convert string labels to integer labels
    unique_labels = labels['x'].unique()
    label_to_int = {label: i for i, label in enumerate(unique_labels)}
    labels['x'] = labels['x'].map(label_to_int)
    
    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')
        
    os.chdir(OutputDir)
    
    if (NumGenes == 0):
        #save labels as csv file with header and index column
        labels.to_csv('Labels_scvi.csv')
        data.to_csv('Data_scvi.csv')    
        
        train = CsvDataset('Data_scvi.csv', save_path = "", sep = ",", labels_file = "Labels_scvi.csv", gene_by_cell = False)
        
        ## this semisupervised trainer automatically uses a part of the input data for training and a part for testing
        scanvi = SCANVI(train.nb_genes, train.n_batches, train.n_labels)
        trainer_scanvi = SemiSupervisedTrainer(scanvi, train, frequency=5)
    
    n_epochs = 200
    
    truelab = []
    pred = []
    tr_time = []
    ts_time = []
    
    for i in range(np.squeeze(nfolds)):
        test_ind_i = test_ind[i] - 1
        train_ind_i = train_ind[i] - 1
        
        if (NumGenes > 0):
            feat_to_use = features.iloc[0:NumGenes,i]
            data2 = data.iloc[:,feat_to_use]
            data2 = preprocess_data(data2)
            
            labels.to_csv('Labels_scvi.csv')
            data2.to_csv('Data_scvi.csv')    
            
            train = CsvDataset('Data_scvi.csv', save_path = "", sep = ",", labels_file = "Labels_scvi.csv", gene_by_cell = False, new_n_genes = False)
            
            ## this semisupervised trainer automatically uses a part of the input data for training and a part for testing
            scanvi = SCANVI(train.nb_genes, train.n_batches, train.n_labels)
            trainer_scanvi = SemiSupervisedTrainer(scanvi, train, frequency=5)

        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=train_ind_i, shuffle = False)
        trainer_scanvi.labelled_set.to_monitor = ['accuracy']  # Only monitor accuracy
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=test_ind_i, shuffle = False)
        trainer_scanvi.unlabelled_set.to_monitor = ['accuracy']  # Only monitor accuracy
        
        # Add epsilon to ensure positive scale values
        epsilon = 1e-6
        scanvi.epsilon = epsilon
    
        start = tm.time()
        trainer_scanvi.train(n_epochs)
        tr_time.append(tm.time()-start)
    
        ## labels of test set are in y_pred
        ## labels are returned in numbers, should be mapped back to the real labels
        ## indices are permutated
        start = tm.time()
        y_true, y_pred = trainer_scanvi.unlabelled_set.compute_predictions()
        ts_time.append(tm.time()-start)
        
        # Convert predictions back to original labels
        int_to_label = {i: label for label, i in label_to_int.items()}
        y_true = [int_to_label[int(x)] for x in y_true]
        y_pred = [int_to_label[int(x)] for x in y_pred]
        
        truelab.extend(y_true)
        pred.extend(y_pred)
    
    #write results
    
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    
    if (NumGenes == 0):  
        truelab.to_csv("scVI_True_Labels.csv", index = False)
        pred.to_csv("scVI_Pred_Labels.csv", index = False)
        tr_time.to_csv("scVI_Training_Time.csv", index = False)
        ts_time.to_csv("scVI_Testing_Time.csv", index = False)
    else:
        truelab.to_csv("scVI_" + str(NumGenes) + "_True_Labels.csv", index = False)
        pred.to_csv("scVI_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        tr_time.to_csv("scVI_" + str(NumGenes) + "_Training_Time.csv", index = False)
        ts_time.to_csv("scVI_" + str(NumGenes) + "_Testing_Time.csv", index = False)