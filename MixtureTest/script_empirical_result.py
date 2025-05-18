# %%
"""
This script processes empirical test results from CSV files and generates a consolidated output.
The script performs the following steps:
1. Imports necessary libraries.
2. Specifies the folder containing the CSV files.
3. Retrieves a list of CSV files that start with 'result' and 'statistics' in the specified folder.
4. Reads each 'result' CSV file into a DataFrame, processes the data, and concatenates it into separate DataFrames for Japan and another country.
5. Reads each 'statistics' CSV file into a DataFrame and concatenates it into separate DataFrames for Japan and another country.
6. Processes the concatenated DataFrames to generate a final output DataFrame.
7. Replaces specific values in the final output DataFrame.
8. Reorders and renames the models in the final output DataFrame.
9. Prints the final output DataFrame.
10. Saves the final output DataFrame to a CSV file named 'result_empirical.csv'.
Variables:
- folder_path: The path to the folder containing the CSV files.
- csv_files: List of CSV files that start with 'result'.
- csv_files_stats: List of CSV files that start with 'statistics'.
- dataframes_jp: DataFrame to store processed data for Japan from 'result' CSV files.
- dataframes_cl: DataFrame to store processed data for another country from 'result' CSV files.
- dataframes_jp_stats: DataFrame to store processed data for Japan from 'statistics' CSV files.
- dataframes_cl_stats: DataFrame to store processed data for another country from 'statistics' CSV files.
- result_output: Final output DataFrame containing processed results.
- model_order: List of model names in the desired order.
- model_name: List of model names to be used as index in the final output DataFrame.
"""
import os
import pandas as pd

# Specify the folder containing the CSV files
folder_path = 'empirical_test'

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv') & file.startswith('result') ]

csv_files_stats = [file for file in os.listdir(folder_path) if file.endswith('.csv') & file.startswith('statistics') ]


# %%
# Read each CSV file into a DataFrame and store them in a list
dataframes_jp = pd.DataFrame()
dataframes_cl = pd.DataFrame()
for csv_file in csv_files:
    file_path = os.path.join(folder_path, csv_file)
    csv_file_name = csv_file.strip('.csv').split('_')
    country = csv_file_name[1]
    y_variable = csv_file_name[2]
    model = ' '.join(csv_file_name[3:-1])
    T_length = csv_file_name[-1]
    
    df = pd.read_csv(file_path, header=1, index_col=0)
    df = df[df.columns[1:]]    
    stats_to_select = ['AIC', 'BIC', 'LR', 'LR 90', 'LR 99']
    available_stats = [stat for stat in stats_to_select if stat in df.index]
    df_row = df.loc[available_stats, :].stack().to_frame().T
    df_row['country'] = country.capitalize()
    df_row['model'] = model.capitalize()
    df_row['T'] = T_length
    df_row['Y'] = y_variable
    # if T_length == '3':
    if country == 'japan':
        dataframes_jp = pd.concat([dataframes_jp, df_row])
    else:
        dataframes_cl = pd.concat([dataframes_cl, df_row])

result_output = dataframes_cl.set_index('model')
result_output = result_output.replace(-1,'10+')

# %%
# Read each Statistics CSV file into a DataFrame and store them in a list
dataframes_jp_stats = pd.DataFrame()
dataframes_cl_stats = pd.DataFrame()
for csv_file in csv_files_stats:
    file_path = os.path.join(folder_path, csv_file)
    csv_file_name = csv_file.strip('.csv').split('_')
    country = csv_file_name[1]
    y_variable = csv_file_name[2]

    model = ' '.join(csv_file_name[3:-1])
    T_length = csv_file_name[-1]
    df = pd.read_csv(file_path, header=0, index_col=0)
    df['Model'] = model.capitalize()
    df['T'] = T_length
    df['Y'] = y_variable
    if country == 'japan':
        dataframes_jp_stats = pd.concat([dataframes_jp_stats, df])
    else:
        dataframes_cl_stats = pd.concat([dataframes_cl_stats, df])

result_combined = dataframes_cl_stats
# %%
# Output test result


mapping = {
    'Nonpar': 'Nonparametric',
    'Plain': 'No Covariate Normal',
    'Plain mixture': 'No Covariate 2-Component Mixture',
    'Plain mixture 3': 'No Covariate 3-Component Mixture',
    'Plain mixture 4': 'No Covariate 4-Component Mixture',
    'Plain ciiu': 'No Covariate Normal control CIIU4',
    'Plain mixture ciiu': 'No Covariate 2-Component Mixture control CIIU4',
    'Plain mixture 3 ciiu': 'No Covariate 3-Component Mixture control CIIU4',
    'Plain mixture 4 ciiu': 'No Covariate 4-Component Mixture control CIIU4',
    'K': 'K Normal',
    'K mixture': 'K 2-Component Mixture',
    'K mixture 3': 'K 3-Component Mixture',
    'K mixture 4': 'K 4-Component Mixture',
    'Kmshare': 'K Material Share Normal',
    'Kmshare mixture': 'K Material Share 2-Component Mixture',
    'Kmshare mixture 3': 'K Material Share 3-Component Mixture',
    'Kmshare mixture 4': 'K Material Share 4-Component Mixture',
    'Kmshare ciiu': 'K Material Share Normal control CIIU4',
    'Kmshare mixture ciiu': 'K Material Share 2-Component Mixture control CIIU4',
    'Kmshare mixture 3 ciiu': 'K Material Share 3-Component Mixture control CIIU4',
    'Kmshare mixture 4 ciiu': 'K Material Share 4-Component Mixture control CIIU4',
    'K ciiu': 'K Normal control CIIU4',
    'K mixture ciiu': 'K 2-Component Mixture control CIIU4',
    'K mixture 3 ciiu': 'K 3-Component Mixture control CIIU4',
    'K mixture 4 ciiu': 'K 4-Component Mixture control CIIU4',
    'Ar1 plain': 'AR1 No Covariate Normal',
    'Ar1 plain mixture': 'AR1 No Covariate 2-Component Mixture',
    'Ar1 plain mixture 3': 'AR1 No Covariate 3-Component Mixture',
    'Ar1 k': 'AR1 K Normal',
    'Ar1 k mixture': 'AR1 K 2-Component Mixture',
    'Ar1 k mixture 3': 'AR1 K 3-Component Mixture',
    'Ar1 k ciiu': 'AR1 K Normal control CIIU4',
    'Ar1 k mixture ciiu': 'AR1 K 2-Component Mixture control CIIU4',
    'Ar1 k mixture 3 ciiu': 'AR1 K 3-Component Mixture control CIIU4',
    'Ar1 kmshare ciiu': 'AR1 kmshare Normal control CIIU4',
    'Ar1 kmshare mixture ciiu': 'AR1 kmshare 2-Component Mixture control CIIU4',
    'Ar1 kmshare mixture 3 ciiu': 'AR1 kmshare 3-Component Mixture control CIIU4'

}

result_output = result_output.loc[mapping.keys(),:]
result_output['Model'] = result_output.index
result_output['Model'] = result_output['Model'].replace(mapping)
# result_output = result_output.sort_values('T')
result_output.to_csv('result_empirical.csv')

# %%
# Filter data for LR and LR Crit statistics with models in `model_order`
LR_stat = dataframes_cl_stats[
    (dataframes_cl_stats['Stat'] == 'LR') &
    (dataframes_cl_stats['Model'].isin(model_order))  # Using .isin() for readability
]

LR_crit_stat = dataframes_cl_stats[
    (dataframes_cl_stats['Stat'] == 'LR Crit') &
    (dataframes_cl_stats['Model'].isin(model_order))  # Using .isin() for readability
]

# Function to reformat the critical values string
def reformat_crit(crit_str):
    """Reformats the 'crit_str'. Converts 'inf, inf, inf' to '-' or rounds numeric values."""
    if crit_str == 'inf, inf, inf':
        return '-'
    else:
        crit_list = crit_str.split(',')
        return ', '.join([f'{round(float(each), 1)}' for each in crit_list])

# Select columns that start with 'M='
M_columns = [col for col in LR_crit_stat.columns if col.startswith('M=')]

# Apply the reformat_crit function to the relevant columns
LR_crit_stat[M_columns] = LR_crit_stat[M_columns].applymap(reformat_crit)

# Define the industry names
ind_names = ['Fabricated metal products', 'Food products', 'Textiles']

# Combine LR_stat and LR_crit_stat, then sort by 'Industry' and 'Model'
LR_stat_report = pd.concat([LR_stat, LR_crit_stat]).sort_values(['Industry', 'Model'])

# Iterate over models and industries to adjust LR statistics
result_output_T = result_output[result_output['T'] == '3']
for each_model in model_order:
    for each_ind in ind_names:
        # Retrieve the chosen number for the current model and industry
        number_chosen = result_output_T.replace('10+', -1).loc[each_model, ('LR', each_ind)]
        
        # If '10+' was replaced with -1, set number_chosen to 10
        if number_chosen == -1:
            number_chosen = 10
        
        # Update LR_stat_report for columns where M > number_chosen
        for m_col in M_columns:
            if int(m_col.strip('M=')) > number_chosen:
                LR_stat_report.loc[
                    (LR_stat_report['Industry'] == each_ind) &
                    (LR_stat_report['Model'] == each_model) &
                    (LR_stat_report['Stat'] == 'LR'),
                    m_col
                ] = '-'

# Export the final DataFrame to a CSV file
LR_stat_report.to_csv("LR_stat_report.csv")

# %%
