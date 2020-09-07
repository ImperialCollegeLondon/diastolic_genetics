#################################################################
# Nicolo Savioli, PhD -  Imperial College London, 15/01/2020    #
#################################################################

import os 
from   shutil import copyfile
import shutil
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from   peakdetect import peakdetect
import collections
import matplotlib.pyplot as plt
import csv
from   tqdm import tqdm
  
def plots(list_data,list_data_dt,pathSave,type_data,strain_data):
    fig            = plt.figure()
    ax             = fig.add_subplot(111)
    ax.set_title   (str(type_data) + " segmental strain and strain rate")
    ax.plot        (list_data, '-',  label="strain",color='r')
    ax.plot        (list_data_dt, '-',  label="strain rate",color='g')
    plt.plot       (strain_data[0][0] ,strain_data[0][1] , "b+")
    plt.plot       (strain_data[1][0] ,strain_data[1][1] , "b+")
    ax.set_ylabel  ('Strain')
    ax.set_xlabel  ('Number of cardiac phases')
    ax.legend      (loc='lower right')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label        = collections.OrderedDict(zip(labels, handles))
    plt.legend     (by_label.values(), by_label.keys())
    fig.savefig    (os.path.join  (pathSave, type_data+".jpg"))
    plt.close('all')

def CorrelationPlot(csv_path):
    import seaborn as sns
    df = pd.read_csv(csv_path, sep=',', header=0) 
    sns.set(style="white", color_codes=True)
    sns.jointplot(x=df["Age"], y=df["Longitudinal strain"], kind='kde', color="skyblue")
    plt.show()
    # Usual boxplot
    ax = sns.boxplot(x=df["Age"], y=df["Longitudinal strain"], data=df)
    # Add jitter with the swarmplot function.
    ax = sns.swarmplot(x=df["Age"],y=df["Longitudinal strain"], data=df, color="grey")
    plt.show()
    sns.set(style="white", color_codes=True)
    sns.jointplot(x=df["Age"], y=df["Radial strain"], kind='kde', color="skyblue")
    plt.show()
    # Usual boxplot
    ax = sns.boxplot(x=df["Age"], y=df["Radial strain"], data=df)
    # Add jitter with the swarmplot function.
    ax = sns.swarmplot(x=df["Age"],y=df["Radial strain"], data=df, color="grey")
    plt.show()

def makeFolder(dirName):
    if not os.path.exists(dirName):
        os.mkdir(dirName)

def get_strain(seqs,look_ahead):
    seq_list             = list(seqs)
    dt_seq_list          = np.diff(seq_list)
    basline_seq          = np.mean(seq_list)
    ####################################
    max_strain,\
    min_strain           = peakdetect(seq_list,    lookahead=look_ahead)
    ####################################
    max_dt_strain,\
    min_dt_strain        = peakdetect(dt_seq_list, lookahead=look_ahead)
    ####################################
    xm_max_strain        = [p[0] for p in max_strain]
    ym_max_strain        = [p[1] for p in max_strain]
    xm_max_dt_strain     = [p[0] for p in max_dt_strain]
    ym_max_dt_strain     = [p[1] for p in max_dt_strain]
    ####################################
    xm_min_strain        = [p[0] for p in min_strain]
    ym_min_strain        = [p[1] for p in min_strain]
    xm_min_dt_strain     = [p[0] for p in min_dt_strain]
    ym_min_dt_strain     = [p[1] for p in min_dt_strain]
    ####################################
    strain_coordinates = []
    final_strain       = 0
    if basline_seq>0:
        xm_max_strain = xm_max_strain[0]
        ym_max_strain = ym_max_strain[0]
        x_final = y_final = 0 
        if len(xm_min_dt_strain)>1:
            index = min(range(len(xm_min_dt_strain)), key=lambda i: abs(xm_min_dt_strain[i]-xm_max_strain))
            x_final = xm_min_dt_strain[index]
            y_final = ym_min_dt_strain[index]
        else:
            x_final  = xm_min_dt_strain
            y_final  = ym_min_dt_strain

        strain_coordinates = [[xm_max_strain,ym_max_strain],[x_final,y_final]]
        final_strain = y_final
    else:
        xm_min_strain = xm_min_strain[0]
        ym_min_strain = ym_min_strain[0]
        x_final = y_final = 0 
        if len(xm_max_dt_strain)>1:
            index = min(range(len(xm_max_dt_strain)), key=lambda i: abs(xm_max_dt_strain[i]-xm_min_strain))
            x_final = xm_max_dt_strain[index]
            y_final = ym_max_dt_strain[index]
        else:
            x_final   = xm_max_dt_strain
            y_final   = ym_max_dt_strain
        strain_coordinates = [[xm_min_strain,ym_min_strain],[x_final,y_final]]
        final_strain  = y_final
    return final_strain,seq_list,dt_seq_list,strain_coordinates

def get_global_strain(path_file,type_data,EID_number,look_ahead):
    data             = pd.read_csv(path_file,index_col=[0])
    rows_cont        = 1 
    out_final_strain = 0
    out_seq_list     = []
    out_dt_seq_list  = [] 
    for row in data.values:
        # global segment 
        if rows_cont == len(data.values):
            final_strain,\
            seq_list,\
            dt_seq_list,\
            strain_coordinates     = get_strain(row,look_ahead)
            out_final_strain       = final_strain
            out_seq_list           = seq_list
            out_dt_seq_list        = dt_seq_list 
            out_strain_coordinates = strain_coordinates
        rows_cont += 1
    return out_final_strain,out_seq_list,\
           out_dt_seq_list,out_strain_coordinates

def get_data(strain_path,csv_path,save_path,\
             path_f_out,generate_plots,look_ahead):
    df    = pd.read_csv(csv_path, sep=',', header=0) 
    raw_data = {}
    eid_list = []
    age_list = []
    longitudinal_strain_list = []
    radial_strain_list = []
    for folder in tqdm(os.listdir(strain_path)):
        try:
            match = df["Encoded_ID"].astype(str).str.contains(folder)
            age_np   = df["Age"][match].to_numpy()
            if True not in match or len(age_np) ==[]:
                continue 
            #print("\n ... EID number: " + str(folder))
            makeFolder(os.path.join(save_path,str(folder)))
            age      = age_np[0]
            longit   = os.path.join(strain_path,folder,"cine_2d_strain_la_4ch_longit.csv")
            radial   = os.path.join(strain_path,folder,"cine_2d_strain_sa_radial.csv")
            global_longit_strain,\
            global_longit_list,\
            global_longit_dev_list,\
            global_longit_strain_coordinates = get_global_strain(longit,"longit",str(folder),look_ahead)
            if generate_plots:
                plots(global_longit_list,global_longit_dev_list,os.path.join(save_path,str(folder)),"longitudinal",global_longit_strain_coordinates)
            global_radial_strain,\
            global_radial_list,\
            global_radial_dev_list,\
            global_radial_strain_coordinates = get_global_strain(radial,"radial",str(folder),look_ahead)
            if generate_plots:
                plots(global_radial_list,global_radial_dev_list,os.path.join(save_path,str(folder)),"radial",global_radial_strain_coordinates)
            df_longit = pd.DataFrame(global_longit_dev_list)
            df_longit.to_csv(os.path.join(save_path,str(folder),"derived_global_strain_longit"+".csv"), index=False, header=False)
            df_radial = pd.DataFrame(global_radial_dev_list)
            df_radial.to_csv(os.path.join(save_path,str(folder),"derived_global_strain_radial"+".csv"), index=False, header=False)
            final_radial_strain = final_longit_strain =  0 
            if isinstance(global_radial_strain, list):
                final_radial_strain = global_radial_strain[0]
            else:
                final_radial_strain = global_radial_strain
            if isinstance(global_longit_strain, list):
                final_longit_strain = global_longit_strain[0]
            else:
                final_longit_strain = global_longit_strain
            df_all    = pd.DataFrame([str(folder),final_radial_strain,final_longit_strain])
            df_all.to_csv(os.path.join(save_path,str(folder),"strain"+".csv"), index=False, header=False)
            eid_list.append(folder)
            age_list.append(age)
            longitudinal_strain_list.append(final_longit_strain)
            radial_strain_list.append(final_radial_strain)
            raw_data["eid"]                 = eid_list
            raw_data["Age"]                 = age_list
            raw_data["long_PDSR"]           = longitudinal_strain_list
            raw_data["radial_PDSR"]         = radial_strain_list
            df_out                          = pd.DataFrame(raw_data, columns = ['eid', 'Age', 'long_PDSR', 'radial_PDSR'])
            df_out.to_csv(os.path.join(path_f_out,'table.csv'))
        except:
            shutil.rmtree(os.path.join(save_path,str(folder)))
            continue

def PostProcessingData(csv_path_strain):
    df         = pd.read_csv(csv_path_strain, sep=',', header=0) 
    max_strain = np.max(df["Radial strain"].to_numpy())
    match      = df["Radial strain"].astype(str).str.contains(str(max_strain))
    p_name     = df["eid"][match].to_numpy()[0]

def plot_from_folder(data_path,data_path_strain):
    #path           = os.path.join(data_path,"derived_global_strain_radial.csv")
    df_no          = pd.read_csv(data_path_strain, sep=',', header=0) 
    #df             = pd.read_csv(path, sep=',', header=0) 
    #list_data_dt   = df.to_numpy()
    data           = df_no.columns.to_numpy()
    data_list = []
    for i in data:
        data_list.append(float(i))    
    final_strain,seq_list,dt_seq_list,strain_coordinates = get_strain(data_list)    
    fig            = plt.figure()
    ax             = fig.add_subplot(111)
    ax.set_title  ("Wrong")
    ax.plot        (data_list, '-',  label="strain",color='r')
    ax.plot        (dt_seq_list, '-',  label="strain derived",color='g')
    plt.plot      (strain_coordinates[0][0] ,strain_coordinates[0][1] , "b+")
    plt.plot      (strain_coordinates[1][0] ,strain_coordinates[1][1] , "b+")
    ax.set_ylabel  ('Strain')
    ax.set_xlabel  ('Time')
    ax.legend      (loc='lower right')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label        = collections.OrderedDict(zip(labels, handles))
    plt.legend     (by_label.values(), by_label.keys())
    fig.savefig    ("strain.jpg")

if __name__ == "__main__":
    print("\n ~~ StrainPeak UKBB 0.1 ~~ \n")
    # 1) Input data folder
    path_strain    = "/mnt/storage/home/nsavioli/cardiac/UKBB_40616/13k_UKBB_NIFTI"
    # 2) Phenotype data csv - you must to have it
    path_ages      = "/mnt/storage/home/nsavioli/cardiac/UKBB_New_Data/phenotype_ukbb_13k.csv"
    # 3) Output folder where you will save the plots 
    path_save      = "/mnt/storage/home/nsavioli/cardiac/UKBB_New_Data/strain/plots"
    # 4) Output folder where you final csv will be save - usually where you have your plot folders 
    path_f_out     = "/mnt/storage/home/nsavioli/cardiac/UKBB_New_Data/strain"
    # 5) True/False if you want to generate the plots for checking. 
    generate_plots = True 
    # 6) Distance to look ahead from a peak candidate to determine if it is the actual peak - default 5 or 10 ( you must to try the optimal)
    look_ahead     = 5
    ##############################################################################
    get_data(path_strain,path_ages,path_save,path_f_out,generate_plots,look_ahead)