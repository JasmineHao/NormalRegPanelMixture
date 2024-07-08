# %%
import pandas as pd
import os,re
import numpy as np

def parse_stars(x):
    # Extract stars from the LRT value
    stars = re.findall(r"\^{(\*+)}", x)
    return len(stars[0]) if stars else 0  # Return the number of stars


def find_best_model_per_metric(sub_df):
     # Apply star parsing on LRT and convert AIC/BIC to float
    sub_df.loc['LRT'] = sub_df.loc['LRT'].apply(parse_stars)
    sub_df.loc['AIC':'BIC'] = sub_df.loc['AIC':'BIC'].applymap(lambda x: float(re.sub(r"[^\d.-]", "", x)))

    # LRT: Find the first column with less than 2 stars
    if sum(sub_df.loc['LRT'] < 2) > 0:
        best_lrt_model = sub_df.columns[(sub_df.loc['LRT'] < 2).argmax()]
    else:
        best_lrt_model = 'M=10'
    # AIC & BIC: Find the minimum point before the values start increasing
    def find_optimal_model(series):
        differences = np.diff(series.values)  # Calculate the first derivative
        min_point = np.where(differences > 0)[0]  # Find where the difference becomes positive
        
        if len(min_point) > 0:
            return series.index[min_point[0]]  # Return the last model before increase
        else:
            return series.index[-1]  # Return the global minimum if no increase is found

    best_aic_model = find_optimal_model(sub_df.loc['AIC'])
    best_bic_model = find_optimal_model(sub_df.loc['BIC'])

    return {
        'Best LRT Model': best_lrt_model,
        'Best AIC Model': best_aic_model,
        'Best BIC Model': best_bic_model
    }


def process_country_data(files_path, country_prefix):
    all_results = []
    
    # List all files for the specified country
    industry_files = [f for f in os.listdir(files_path) if f.startswith(country_prefix) and f.endswith(".csv")]
    
    
    # Process each file (each file represents one industry)
    for file in industry_files:
        
        model_name = (file.strip(country_prefix).strip("_|.csv"))
        df = pd.read_csv(os.path.join(files_path, file))
        df = df.loc[df['M=1']!='0',]
        df = df.set_index("Unnamed: 0")
        
        # Assuming the structure is rows for LRT, AIC, BIC and columns for models M=1 to M=10
        for i in range(3):
            sub_df = df.iloc[((i)*3):((i+1)*3)]
            
            industry_name = sub_df.index[0]
            sub_df.index = ['LRT', 'AIC', 'BIC']
            
            if 'AR1' in file:
                sub_df = sub_df.iloc[:,0:5]
            else:
                sub_df = sub_df.iloc[:,0:10]
                    
            best_models = find_best_model_per_metric(sub_df)
            # Append results including the industry name
            all_results.append({
            'Model': model_name,
            'Industry': industry_name.strip("1"),
            **best_models
             })
    # Convert results to DataFrame
    result_df = pd.DataFrame(all_results)
    return(result_df)
    # Output to Excel


# %%
# Example usage
country_prefix = 'Japan'  # Adjust based on the file naming convention
country_prefix = 'Chile' #['Chile', 'Japan']
files_path = r'C:\Users\JasmineHao\Desktop\EmpiricalTestResults'

files_path = "/Users/haoyu/Desktop/Empirical"
# %%

for country_prefix in ['Chile', 'Japan']:
    result_df = process_country_data(files_path, country_prefix)
    melted_df = result_df.melt(id_vars=['Model', 'Industry'], var_name='Metric', value_name='Best Model')
    pivot_df = pd.pivot_table(melted_df, index=['Industry', 'Metric'], columns='Model', values='Best Model', aggfunc='first').sort_index(level=1).T
    output_file = os.path.join(files_path, f"{country_prefix}_best_models.xlsx")
    pivot_df.to_excel(output_file)

    print(f"Results saved to {output_file}")


# %%