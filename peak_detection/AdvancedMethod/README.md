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

---

## Running the Script

To run the script, execute it from the command line or terminal with the required arguments:

```python
 python main_detect_peaks_ukb.py --signal --input_dir --output_dir 
```


Arguments Explained:
- `--signal`: The type of strain data to process ("radial", "circum", "longit").
- `--input_dir`: The directory containing subject strain data.
- `--output_dir`: The directory where results will be saved.

# Records
Created by Majid Vafaeezadeh, 20-01-2025
