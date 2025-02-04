# PDSR Data,and Analysis Pipeline

## Code

The PDSR values were calculated using a Python pipeline, which processes the strain data to detect peaks and derivatives. The main_detect_peaks_ukb.py code is the main script.

## Pipeline Overview

The processing follows these key steps:

1. Reading input data:
   - Load strain data from CSV files.
   - Convert time series strain values to a numerical format, handling missing values (NA, -).
2. Peak Detection Algorithm:
   - Detect strain peaks using the pyGAM package.
   - Fit Expectile Generalized Additive Models (GAM) to smooth signals.
   - Locate the derivative peaks to determine PDSR values.
3. Saving Result:
   - Save processed data into output CSV files.

## Running the Script

To run the script, execute it from the command line or terminal with the required arguments:

```python
 python main_detect_peaks_ukb.py --signal --input_dir --output_dir 
```


Arguments Explained:
- `--signal`: The type of strain data to process ("radial", "circum", "longit").
- `--input_dir`: The directory containing subject strain data.
- `--output_dir`: The directory where results will be saved.

## Supplementary Scripts: 
**PDSR Data Integration Process**

The combine.py script is used to aggregate and process the extracted PDSR values from multiple CSV files generated during the peak detection process. The script selects relevant data from the results and compiles it into a single consolidated CSV file for further analysis.

**PDSR Data Merging Process**

Combine_longi_radial.py: This script combines the extracted PDSR values from two separate datasets—radial and longitudinal strain measurements—into a single final file for comprehensive analysis.


**PDSR Comparison and Analysis Process**

ScatterPlot.py: This script compares PDSR_radial values between two datasets—one with 81K instances and another with 48K instances—to evaluate differences and generate visual insights.

# Records
Created by Majid Vafaeezadeh, 20-01-2025
