# Peak detection

## Usage

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


