import pandas as pd

# Create a sample DataFrame
data = {'Name': ['Alice', 'Bob', 'Charlie'],
        'Age': [25, 30, 35],
        'City': ['New York', 'San Francisco', 'Los Angeles']}
df = pd.DataFrame(data)

# Write the DataFrame to a text file
# Open the file manually
file = open('output.txt', 'w')

# Write the DataFrame to the file
file.write(df.to_string(index=False))  # `index=False` to exclude row indices

# Close the file
file.close()