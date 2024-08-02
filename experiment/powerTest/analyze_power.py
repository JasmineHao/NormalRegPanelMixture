# %%
import pandas as pd

# Load the CSV data
file_path = '/Users/haoyu/Documents/GitHub/NormalRegPanelMixture/results/powerTest/powerTestM2SimPLR.csv'

file_path = '/Users/haoyu/Documents/GitHub/NormalRegPanelMixture/results/powerTest/powerTestM2_bootstrap.csv'

data = pd.read_csv(file_path)

# Display the first few rows of the data to understand its structure
print(data.head())

# Check the columns and unique values
print("Columns:", data.columns)
print("Unique V5 values:", data['V5'].unique())

for each_v5 in data['V5'].unique():
    # Pivot the data to create a multi-indexed DataFrame
    pivot_data = data[data['V5']==each_v5].pivot_table(
        index=['V2', 'V3'],   # Row index for (N, T) combinations
        columns=['V4', 'V6'],  # Column index for (C, G) combinations, separate panel for each V5
        values='V1',  # Assuming 'value' is the column for test results
        aggfunc='mean'  # Use mean in case of duplicate entries, adjust as needed
    )

    # Display the pivot table
    print(each_v5)

    print(pivot_data)


# %%
