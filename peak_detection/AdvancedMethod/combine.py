import pandas as pd
import os
import re

# Directory containing CSV files
input_dir = 'P:/Projects/strain/dockerstrain/output/longit/results'
output_file = '...'

# Initialize an empty list to store dataframes
global_rows = []

# Iterate through CSV files in the directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('.csv'):
        # Extract the 7-digit number and instance (V2 or V3) from the file name
        match = re.search(r'^(\d{7})_(V[2|3])\.csv$', file_name)
        if not match:
            continue  # Skip files without the required naming format
        digital_number = match.group(1)
        instance = match.group(2)

        # Read the CSV file
        file_path = os.path.join(input_dir, file_name)
        df = pd.read_csv(file_path)

        # Filter rows where signal_name is 'Global'
        global_df = df[df['signal_name'] == 'Global'].copy()

        # Add the 'digital_number' and 'instance' columns
        global_df.insert(0, 'digital_number', digital_number)  # Add as the first column
        global_df.insert(1, 'instance', instance)  # Add as the second column

        # Append the result to the list
        global_rows.append(global_df)

# Combine all filtered rows into a single DataFrame
combined_global_rows = pd.concat(global_rows, ignore_index=True)

# Save the result to a CSV file
combined_global_rows.to_csv(output_file, index=False)

print(f"Extracted 'Global' rows saved to {output_file}")
