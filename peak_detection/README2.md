# Peak Detection in Cardiac Strain Data

## Overview
This repository provides a pipeline for detecting peaks in cardiac strain data from the UK Biobank (UKBB). The repository now includes two different methods for peak detection:

1. **Original Peak Detection Method** (Baseline)
2. **New Peak Detection Method using Expectile GAMs** (Advanced)

Both methods process strain data and extract key features used in downstream analyses.

---

## Methods

### 1. Original Peak Detection Method
This method is based on the `peakdetect.py` script and identifies peaks using basic peak detection techniques. It relies on predefined thresholding and simple peak searching algorithms.

## Usage
### Running the Original Peak Detection Method
```
Enter in main.py and please add the position of the corresponding paths:

```
# 1) Input data folder
path_strain    = ""

# 2) Phenotype data csv - you must to have it
path_ages      = "./phenotype_ukbb_13k.csv"

# 3) Output folder where you will save the plots 
path_save      = "./plots_path_folder"

# 4) Output folder where you final csv will be save - usually where you have your plot folders 
path_f_out     = "./strain_path_folder"

# 5) True/False if you want to generate the plots for checking. 
generate_plots = True 

# 6) Distance to look ahead from a peak candidate to determine if it is the actual peak - default 5 or 10 ( you must to try the optimal)
look_ahead     = 5
```

Then execute it:

```
python main.py
```
### 2. New Peak Detection Method (PDSR-GAM)
This new approach enhances the peak detection process by leveraging **Expectile Generalized Additive Models (GAMs)** to smooth the signal and accurately locate peaks in the **first derivative** of the strain signal.

#### **Key Features of the New Method**
- Uses **pyGAM** to model and smooth the strain signal.
- Identifies **both strain peaks and first derivative peaks**.
- Robust against signal noise and artifacts.
- Optimized for large-scale processing with parallelization.

#### **Pipeline Overview**
1. **Read input data**: Loads strain signals from `.csv` files.
2. **Preprocessing**: Cleans missing values and normalizes the signal.
3. **Peak Detection with Expectile GAMs**:
   - Fits a GAM model to smooth the signal.
   - Detects the highest peak in the **strain signal**.
   - Detects the **first derivative peak**, which is critical for calculating PDSR values.
4. **Output and Visualization**:
   - Saves extracted peak values in `.csv` format.
   - Generates plots for quality control.

---

## Installation
To use this repository, install the required dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

### Running the New Peak Detection Method
```bash
python main_detect_peaks_ukb.py --signal radial --input_dir ./data --output_dir ./results
```

#### **Command-Line Arguments**

| Argument       | Description |
|---------------|-------------|
| `--signal`    | Type of strain signal (`radial`, `circum`, `longit`) |
| `--input_dir` | Directory containing strain data files |
| `--output_dir`| Directory to save results |
| `--do_plots`  | Generate plots for peak detection (optional) |
---
## File Structure
```
peak_detection/
├── peakdetect.py  # Original peak detection script
├── main.py        # Original pipeline
├── main_detect_peaks_ukb.py  # New PDSR-GAM method
├── requirements.txt  # Dependencies
├── README.md      # Documentation
```
## Results
The output includes:
- `results/strain_peaks.csv`: Contains detected peak values for each subject.
- `results/plots/`: Folder with visualizations of detected peaks.

Example output format:
```
subject_id, strain_peak_value, derivative_peak_value
1000018, 1.75, -6.66
1000079, 1.91, -7.21
```

---

## Contributions
- Nicolo Savioli (Imperial College London)
- Majid Vafaeezadeh - Added PDSR-GAM method
This repository was updated to include the **PDSR-GAM** method for improved peak detection. Contributions are welcome! Please open an issue or submit a pull request for any improvements.
