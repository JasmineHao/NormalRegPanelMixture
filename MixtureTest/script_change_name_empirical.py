
# %%
import os
import pandas as pd

# Specify the folder containing the CSV files
folder_path = 'empirical_test'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv') & file.startswith('result') ]

csv_files_stats = [file for file in os.listdir(folder_path) if file.endswith('.csv') & file.startswith('statistics') ]

# %%
for file in csv_files_stats:
    if not (file.startswith('statistics_chile_mLshare_') | file.startswith('statistics_chile_mYshare_')):
        new_name = 'statistics_chile_mYshare_' + file.split('statistics_chile_', 1)[1]
        # print(file, new_name)
        os.rename(os.path.join(folder_path, file), os.path.join(folder_path, new_name))

for file in csv_files:
    if not (file.startswith('result_chile_mLshare_') | file.startswith('result_chile_mYshare_')):
        new_name = 'result_chile_mYshare_' + file.split('result_chile_', 1)[1]
        # print(file, new_name)
        os.rename(os.path.join(folder_path, file), os.path.join(folder_path, new_name))

# %%
