import os
import pandas as pd

# base_dir = '/mnt/volume1/qyc/auto-annot/results/'
base_dir = '/mnt/volume1/qyc/auto-annot/result_across/'

for root, dirs, files in os.walk(base_dir):
    # Find testing and training time files
    testing_files = [f for f in files if f.endswith('Testing_Time.csv')]
    
    for testing_file in testing_files:
        # Get the corresponding training file
        prefix = testing_file.split('Testing_Time.csv')[0]
        training_file = prefix + 'Training_Time.csv'
        
        if training_file in files:
            # Full paths to the files
            testing_path = os.path.join(root, testing_file)
            training_path = os.path.join(root, training_file)
            
            # Read the files
            test_df = pd.read_csv(testing_path, header=None)
            train_df = pd.read_csv(training_path, header=None)
            
            # # Sum the values from line 2-6 (assuming 0-based index, lines 1-5)
            # total_values = test_df.iloc[1:6, 0].astype(float) + train_df.iloc[1:6, 0].astype(float)

            # Sum the values from line 2-6 (assuming 0-based index, lines 1)
            total_values = test_df.iloc[1, 0].astype(float) + train_df.iloc[1, 0].astype(float)
            
            # Create new dataframe with header from test file and summed values
            total_df = pd.DataFrame({
                # 0: [test_df.iloc[0, 0]] + total_values.tolist()
                0: [test_df.iloc[0, 0]] + total_values
            })
            
            # Create output filename
            output_file = os.path.join(root, f"{prefix}Total_Time.csv")
            
            # Save the new file without index and header
            total_df.to_csv(output_file, index=False, header=False)
            
            print(f"Created: {output_file}")

print("Processing complete!")