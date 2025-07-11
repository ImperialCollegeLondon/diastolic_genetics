import pandas as pd
from pathlib import Path

# Define input and output file paths (modify these paths as needed)
input_file1 = Path("path/to/first_input.csv")  # Path to the first input file (e.g., radial strain data)
input_file2 = Path("path/to/second_input.csv")  # Path to the second input file (e.g., longitudinal strain data)
output_file = Path("path/to/output.csv")  # Path where the merged output file will be saved


# Function to read CSV safely
def read_csv_safe(filepath):
    try:
        return pd.read_csv(filepath)
    except FileNotFoundError:
        print(f"Error: File not found - {filepath}")
        exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File is empty - {filepath}")
        exit(1)
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        exit(1)

# Read input CSV files
df1 = read_csv_safe(input_file1)
df2 = read_csv_safe(input_file2)

# Select and rename required columns
columns_needed = ["digital_number", "instance", "der_peak_value"]
df1 = df1[columns_needed].rename(columns={"der_peak_value": "der_peak_value_1"})
df2 = df2[columns_needed].rename(columns={"der_peak_value": "der_peak_value_2"})

# Merge the dataframes on "digital_number" and "instance"
merged_df = pd.merge(df1, df2, on=["digital_number", "instance"], how="inner")

# Save the merged dataframe
merged_df.to_csv(output_file, index=False)

print(f"Output file saved to: {output_file}")

