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
    csv_file_name = csv_file.strip('.csv').split('_')
    country = csv_file_name[1]
    model = ' '.join(csv_file_name[2:])
    
    df = pd.read_csv(file_path, header=1, index_col=0)
    df = df[df.columns[1:]]    
    df_row = df.loc[['AIC','BIC','LR'],:].stack().to_frame().T
    df_row['country'] = country.capitalize()
    df_row['model'] = model.capitalize()
    if country == 'japan':
        dataframes_jp = pd.concat([dataframes_jp, df_row])
    else:
        dataframes_cl = pd.concat([dataframes_cl, df_row])


# %%
# combined = pd.concat([dataframes_jp.set_index('model'), dataframes_cl.set_index('model')],axis=1)
combined = dataframes_cl.set_index('model')
combined = combined.replace(-1,'10+')
model_order= ['Plain', 'Plain mixture', 'K', 'Kl', 'Kl spline',  'K mixture', 'Kl mixture', 'Kl mixture spline', 'Ar1 plain', 'Ar1 k', 'Ar1 kl', 'Ar1 plain mixture']

combined = combined.loc[model_order,:]
combined.index = ['No Covariate Normal', 'No Covariate 2-component mixture', 'K Normal', 'Kl Normal', 'Kl spline Normal',  'K 2-component mixture', 'Kl 2-component mixture', 'Kl spline 2-component mixture', 'Ar1 Normal', 'Ar1 k Normal', 'Ar1 kl Normal', 'Ar1 2-component mixture']

print(combined)
# %%
combined.to_csv('result_empirical.csv')
# %%

model_not_included = ['Plain mixture']