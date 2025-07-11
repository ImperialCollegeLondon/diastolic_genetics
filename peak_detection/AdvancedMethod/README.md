# PDSR Data Analysis and Pipeline

## Overview

This repository contains a pipeline for calculating **Peak Derivative Strain Rate (PDSR)** from cardiac strain data. The pipeline processes strain data, detects peaks using **Generalized Additive Models (GAMs)**, and extracts PDSR values for further analysis.

The pipeline includes multiple scripts for **data processing, integration, and visualization**, ensuring a structured approach to cardiac strain analysis.

---

## **Pipeline Overview**

### **1. Peak Detection and Data Processing**
- **`main_detect_peaks_ukb.py`**: The core script for processing strain data and detecting peaks using GAMs.
  - Reads input strain data from CSV files.
  - Uses Expectile GAMs for peak detection.
  - Extracts peak derivative strain rate (PDSR) values.

### **2. Data Integration**
- **`combine.py`**: Aggregates extracted PDSR values from multiple CSV files into a single dataset.
  - Reads PDSR values from individual subject files.
  - Filters relevant data (e.g., `Global` signal).
  - Saves the consolidated dataset for further analysis.

- **`Combine_longi_radial.py`**: Merges **radial** and **longitudinal** strain datasets into one file.
  - Reads **radial** and **longitudinal** strain results.
  - Aligns data based on `digital_number` and `instance`.
  - Produces a unified dataset for comprehensive PDSR analysis.

### **3. Data Comparison and Visualization**
- **`ScatterPlot.py`**: Compares **PDSR values** from two datasets and generates visualizations.
  - Reads two PDSR datasets (e.g.,  A_dataset vs. B_dataset).
  - Merges data based on `eid` (patient ID).
  - Filters relevant data and calculates differences.
  - Generates **scatter plots** and **comparison metrics**.

---

## **Installation**
Ensure you have Python and the required libraries installed:

```bash
pip install pandas matplotlib scipy numpy sklearn
```

---

### **Python Environment**
The pipeline was developed and tested using:
- **Python 3.9**
- The following packages (as defined in `requirements.txt`):
  - `pandas==2.2.3`
  - `numpy==1.26.4`
  - `matplotlib==3.9.4`
  - `scipy==1.11.4`
  - `scikit-learn==1.6.1`
  - `pygam==0.9.1`
  - `tqdm==4.67.1`
## **Usage Instructions**

### **1. Running the Peak Detection Script**
To detect peaks in strain data:

```bash
python main_detect_peaks_ukb.py --signal radial --input_dir ./data --output_dir ./results
```

**Arguments Explained:**
- `--signal`: Type of strain data (`radial`, `circum`, `longit`).
- `--input_dir`: Path to input CSV files.
- `--output_dir`: Folder where results will be saved.

---

### **2. Combining Extracted PDSR Data**
To merge PDSR values into a single file:

```bash
python combine.py
```

- Make sure to update `input_dir` and `output_file` in the script.

---

### **3. Merging Radial & Longitudinal Data**
To combine radial and longitudinal PDSR datasets:

```bash
python Combine_longi_radial.py
```

- Ensure the correct input file paths are set in the script.

---

### **4. Comparing Datasets and Generating Scatter Plots**
To compare PDSR values from two datasets:

```bash
python ScatterPlot.py
```

- Updates the script with paths to the **A_dataset** and **B_dataset**.
- The script generates **scatter plots** for visualization.

---

## **File Descriptions**

### **1. main_detect_peaks_ukb.py**
- Detects **strain peaks** and calculates **PDSR values**.
- Uses **Expectile GAMs** for smooth signal processing.
- Saves results in CSV format.

### **2. combine.py**
- Collects and consolidates PDSR values from individual **CSV files**.
- Filters **Global** strain data.
- Outputs a **single merged CSV file**.

### **3. Combine_longi_radial.py**
- Merges **radial** and **longitudinal** PDSR datasets.
- Ensures matching **`digital_number`** and **`instance`** values.
- Saves the final dataset in CSV format.

### **4. ScatterPlot.py**
- Compares PDSR values between two datasets (e.g., **A_dataset vs. B_dataset** subjects).
- Filters based on **matching IDs** (`eid_40616`).
- Computes **differences, median error, and mean squared error (MSE)**.
- Generates **scatter plots**.

---

## **Data Output Structure**
```bash
output/
├── results/               # Processed PDSR values
│   ├── subject1_PDSR.csv
│   ├── subject2_PDSR.csv
│   ├── combined_PDSR.csv  # Final merged file
│   ├── merged_radial_longi.csv  # Combined radial & longitudinal data
│   ├── comparison_plot.png  # Scatter plot output
└── logs/                  # Logs and metadata
```

---

## **Contributors**
- **Majid Vafaeezadeh** – Developer of the pipeline (20-01-2025)
- Contributions are welcome! Feel free to submit **issues or pull requests**.

---





