import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
from sklearn.metrics import mean_squared_error
# Read the CSV files
file1_path = '//isd_netapp/cardiac$/Majid/strain/peak-detect-strain/Data/Strain/PDSR/81K_PDSR_radial_longit_Final.csv'
file2_path = '//isd_netapp/cardiac$/UKBB_40616/Phenotypes/48K_PDSR.csv'
# Output file path for differences
output_diff_path = "differencesradi.csv"

# Read the CSV files
df1 = pd.read_csv(file1_path)
df2 = pd.read_csv(file2_path)
# Filter the first file to include only rows where "instance" is "V2"
df1_filtered = df1[df1["instance"] == "V2"]
# Merge the two files on "eid_40616" to find matching IDs
merged_df = pd.merge(
    df1_filtered[["eid_40616", "PDSR_radial"]],  # Take relevant columns from file1    longit
    df2[["eid_40616", "PDSR_radial"]],          # Take relevant columns from file2
    on="eid_40616",
    suffixes=("_81K_PDSR", "_48k")
)
# Remove rows with invalid or missing numeric values ("-" or "NA")
merged_df.replace(["-", "NA"], np.nan, inplace=True)  # Replace invalid values with NaN
merged_df.dropna(inplace=True)  # Drop rows with NaN values
# Convert columns to numeric to ensure proper calculations
merged_df["PDSR_radial_81K_PDSR"] = pd.to_numeric(merged_df["PDSR_radial_81K_PDSR"])
merged_df["PDSR_radial_48k"] = pd.to_numeric(merged_df["PDSR_radial_48k"])
# Calculate the absolute difference
merged_df["absolute_difference"] = abs(merged_df["PDSR_radial_81K_PDSR"] - merged_df["PDSR_radial_48k"])
# Calculate the difference
merged_df["difference"] = merged_df["PDSR_radial_81K_PDSR"] - merged_df["PDSR_radial_48k"]




# Save IDs and differences where the values differ in a CSV file
diff_df = merged_df[merged_df["difference"] != 0][["eid_40616", "difference"]]
diff_df.to_csv(output_diff_path, index=False)
# Calculate metrics
median_abs_diff = merged_df["absolute_difference"].median()
max_abs_diff = merged_df["absolute_difference"].max()
mse = ((merged_df["PDSR_radial_81K_PDSR"] - merged_df["PDSR_radial_48k"]) ** 2).mean()
# Print metrics
print(f"Median Absolute difference: {median_abs_diff}")
print(f"Maximum Absolute difference: {max_abs_diff}")
print(f"Mean Squared Error: {mse}")
# Plot the differences
plt.figure(figsize=(10, 6))
plt.scatter(
    merged_df["PDSR_radial_81K_PDSR"],
    merged_df["PDSR_radial_48k"],
    color="blue",
    alpha=0.6,
    label="Values"
)
plt.plot(
    [merged_df["PDSR_radial_81K_PDSR"].min(), merged_df["PDSR_radial_81K_PDSR"].max()],
    [merged_df["PDSR_radial_81K_PDSR"].min(), merged_df["PDSR_radial_81K_PDSR"].max()],
    color="red",
    linestyle="--",
    label="y = x"
)
plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48k_PDSR)")
plt.title("Comparison of PDSR_radial Values")
plt.legend()
plt.grid(True)




# Highlight points where values differ
diff_points = merged_df[merged_df["difference"] != 0]
plt.scatter(diff_points["PDSR_radial_81K_PDSR"], diff_points["PDSR_radial_48k"], c="red", label="Differences")

# Add a legend
plt.legend()

# Show the plot
plt.show()

print(f"IDs with differences and their values are saved to: {output_diff_path}")

# Filter out IDs with differences for the second plot
filtered_df = merged_df[~merged_df["eid_40616"].isin(diff_df["eid_40616"])]

# Plot 2: Scatter plot excluding IDs from differences.csv
plt.figure(figsize=(10, 6))
plt.scatter(filtered_df["PDSR_radial_81K_PDSR"], filtered_df["PDSR_radial_48k"], c="green", label="Filtered Matching IDs")
plt.axline((0, 0), slope=1, color="red", linestyle="--", label="y=x (Perfect Match)")
plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48k_PDST)")
plt.title("Comparison of PDSR_radial Values (Filtered)")
plt.legend()
plt.grid(True)
plt.show()
=====================


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error
from pathlib import Path

# Define input and output file paths (modify these paths as needed)
file1_path = Path("path/to/first_input.csv")  # Path to the first dataset (e.g., 81K PDSR radial & longit)
file2_path = Path("path/to/second_input.csv")  # Path to the second dataset (e.g., 48K PDSR dataset)
output_diff_path = Path("path/to/output_differences.csv")  # Output file for differences

# Function to read CSV safely
def read_csv_safe(filepath):
    """Reads a CSV file safely with error handling."""
    try:
        return pd.read_csv(filepath)
    except FileNotFoundError:
        print(f" Error: File not found - {filepath}")
        exit(1)
    except pd.errors.EmptyDataError:
        print(f" Error: File is empty - {filepath}")
        exit(1)
    except Exception as e:
        print(f" Error reading file {filepath}: {e}")
        exit(1)

# Read input CSV files
df1 = read_csv_safe(file1_path)
df2 = read_csv_safe(file2_path)

# Filter df1 to include only rows where "instance" is "V2"
df1_filtered = df1[df1["instance"] == "V2"]

# Merge the two files on "eid_40616" to match subject IDs
merged_df = pd.merge(
    df1_filtered[["eid_40616", "PDSR_radial"]],  # Select relevant columns from first dataset
    df2[["eid_40616", "PDSR_radial"]],           # Select relevant columns from second dataset
    on="eid_40616",
    suffixes=("_81K_PDSR", "_48k")
)

# Handle missing or invalid values
merged_df.replace(["-", "NA"], np.nan, inplace=True)
merged_df.dropna(inplace=True)

# Convert columns to numeric for calculations
merged_df["PDSR_radial_81K_PDSR"] = pd.to_numeric(merged_df["PDSR_radial_81K_PDSR"])
merged_df["PDSR_radial_48k"] = pd.to_numeric(merged_df["PDSR_radial_48k"])

# Compute absolute and relative differences
merged_df["absolute_difference"] = abs(merged_df["PDSR_radial_81K_PDSR"] - merged_df["PDSR_radial_48k"])
merged_df["difference"] = merged_df["PDSR_radial_81K_PDSR"] - merged_df["PDSR_radial_48k"]

# Save differences to a CSV file
diff_df = merged_df[merged_df["difference"] != 0][["eid_40616", "difference"]]
diff_df.to_csv(output_diff_path, index=False)

# Compute error metrics
median_abs_diff = merged_df["absolute_difference"].median()
max_abs_diff = merged_df["absolute_difference"].max()
mse = mean_squared_error(merged_df["PDSR_radial_81K_PDSR"], merged_df["PDSR_radial_48k"])

# Print summary statistics
print(f" Median Absolute Difference: {median_abs_diff}")
print(f" Maximum Absolute Difference: {max_abs_diff}")
print(f" Mean Squared Error: {mse}")

# Plot comparison of PDSR_radial values
plt.figure(figsize=(10, 6))
plt.scatter(
    merged_df["PDSR_radial_81K_PDSR"], 
    merged_df["PDSR_radial_48k"], 
    color="blue", alpha=0.6, label="Values"
)

# Add reference y = x line
plt.plot(
    [merged_df["PDSR_radial_81K_PDSR"].min(), merged_df["PDSR_radial_81K_PDSR"].max()],
    [merged_df["PDSR_radial_81K_PDSR"].min(), merged_df["PDSR_radial_81K_PDSR"].max()],
    color="red", linestyle="--", label="y = x"
)

plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48K_PDSR)")
plt.title("Comparison of PDSR_radial Values")
plt.legend()
plt.grid(True)
plt.show()

# Highlight points where values differ
plt.figure(figsize=(10, 6))
diff_points = merged_df[merged_df["difference"] != 0]
plt.scatter(diff_points["PDSR_radial_81K_PDSR"], diff_points["PDSR_radial_48k"], c="red", label="Differences")

plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48K_PDSR)")
plt.title("Highlighted Differences in PDSR_radial")
plt.legend()
plt.grid(True)
plt.show()

# Exclude IDs with differences for a filtered plot
filtered_df = merged_df[~merged_df["eid_40616"].isin(diff_df["eid_40616"])]

plt.figure(figsize=(10, 6))
plt.scatter(filtered_df["PDSR_radial_81K_PDSR"], filtered_df["PDSR_radial_48k"], c="green", label="Filtered Matching IDs")
plt.axline((0, 0), slope=1, color="red", linestyle="--", label="y = x (Perfect Match)")
plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48K_PDSR)")
plt.title("Comparison of PDSR_radial Values (Filtered)")
plt.legend()
plt.grid(True)
plt.show()

print(f"Differences saved in: {output_diff_path}")
