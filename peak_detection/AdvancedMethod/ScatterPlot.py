import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error
from pathlib import Path

# Define input and output file paths (modify these paths as needed)
file1_path = Path("path/to/A_dataset.csv")  # Path to the first dataset (e.g., 81K PDSR radial & longit)
file2_path = Path("path/to/B_dataset.csv")  # Path to the second dataset (e.g., 48K PDSR dataset)
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

# Merge the two files on "eid_40" to match subject IDs
merged_df = pd.merge(
    df1_filtered[["eid_40", "PDSR_radial"]],  # Select relevant columns from first dataset
    df2[["eid_40", "PDSR_radial"]],           # Select relevant columns from second dataset
    on="eid_40",
    suffixes=("_81K_PDSR", "_48k")
)

# Handle missing or invalid values
merged_df.replace(["-", "NA"], np.nan, inplace=True)
merged_df.dropna(inplace=True)

# Convert columns to numeric for calculations
merged_df["A_dataset"] = pd.to_numeric(merged_df["A_dataset"])
merged_df["B_dataset"] = pd.to_numeric(merged_df["B_dataset"])

# Compute absolute and relative differences
merged_df["absolute_difference"] = abs(merged_df["A_dataset"] - merged_df["B_dataset"])
merged_df["difference"] = merged_df["A_dataset"] - merged_df["B_dataset"]

# Save differences to a CSV file
diff_df = merged_df[merged_df["difference"] != 0][["eid_40", "difference"]]
diff_df.to_csv(output_diff_path, index=False)

# Compute error metrics
median_abs_diff = merged_df["absolute_difference"].median()
max_abs_diff = merged_df["absolute_difference"].max()
mse = mean_squared_error(merged_df["A_dataset"], merged_df["B_dataset"])

# Print summary statistics
print(f" Median Absolute Difference: {median_abs_diff}")
print(f" Maximum Absolute Difference: {max_abs_diff}")
print(f" Mean Squared Error: {mse}")

# Plot comparison of PDSR_radial values
plt.figure(figsize=(10, 6))
plt.scatter(
    merged_df["A_dataset"], 
    merged_df["B_dataset"], 
    color="blue", alpha=0.6, label="Values"
)

# Add reference y = x line
plt.plot(
    [merged_df["A_dataset"].min(), merged_df["A_dataset"].max()],
    [merged_df["A_dataset"].min(), merged_df["A_dataset"].max()],
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
plt.scatter(diff_points["A_dataset"], diff_points["B_dataset"], c="red", label="Differences")

plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48K_PDSR)")
plt.title("Highlighted Differences in PDSR_radial")
plt.legend()
plt.grid(True)
plt.show()

# Exclude IDs with differences for a filtered plot
filtered_df = merged_df[~merged_df["eid_40"].isin(diff_df["eid_40"])]

plt.figure(figsize=(10, 6))
plt.scatter(filtered_df["A_dataset"], filtered_df["B_dataset"], c="green", label="Filtered Matching IDs")
plt.axline((0, 0), slope=1, color="red", linestyle="--", label="y = x (Perfect Match)")
plt.xlabel("PDSR_radial (81K_PDSR)")
plt.ylabel("PDSR_radial (48K_PDSR)")
plt.title("Comparison of PDSR_radial Values (Filtered)")
plt.legend()
plt.grid(True)
plt.show()

print(f"Differences saved in: {output_diff_path}")
