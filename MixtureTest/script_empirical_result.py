# %%
import os
import pandas as pd

# Specify the folder containing the CSV files
folder_path = 'empirical_test'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv') & file.startswith('result') ]

# %%
# Read each CSV file into a DataFrame and store them in a list
dataframes_jp = pd.DataFrame()
dataframes_cl = pd.DataFrame()
for csv_file in csv_files:
    file_path = os.path.join(folder_path, csv_file)
    _, country, model = csv_file.strip('.csv').split('_')
    df = pd.read_csv(file_path, header=1, index_col=0)
    df = df[df.columns[1:]]    
    df_row = df.loc[['AIC','BIC','LR'],:].stack().to_frame().T
    df_row['country'] = country.capitalize()
    df_row['model'] = model.capitalize()
    if country == 'japan':
        dataframes_jp = pd.concat([dataframes_jp, df_row])
    else:
        dataframes_cl = pd.concat([dataframes_cl, df_row])



combined = pd.concat([dataframes_jp.set_index('model'), dataframes_cl.set_index('model')],axis=1)
combined.to_csv('result_empirical.csv')
# %%
