
# 4DSegment-2.0

> Semi Deep learning pipline for 4D heart analysis 

The code in this repository implements the new version of 4Dsegment, a pipeline for carrying out deep learning 4D segmentation, non-rigid co-registration, 4D mesh generation and motion tracking using raw grey-scale cardiac MRI data in NIfTI format.

The new implementation can process the motion from 4D segmentation made with deep learning technique and increase the total speedup of 30% respect the previous pipeline.
The 4D segmentation is made with a direct mapping between low grayscale segmentation to 3D segmentation for all cardiac phases. The coregistration core is still MIRTK (https://github.com/BioMedIA/MIRTK)  

![](img/1.gif)

*The figure shows the epicardium as well with endocardium and right ventricle surface coregister with the motion script of 4DSegment2.0*

![](img/2.gif)

*The figure shows the motion surface of the entire heart coregister with the motion script of 4DSegment2.0*

![](img/3.gif)

*The figure shows the 4D motion coregistration processed with the new super-resolution deep learning core of 4DSegment2.0*

## Installation

### Install Docker

For Windows 10 Pro first install [Docker](https://www.docker.com/docker-windows). Windows 10 Home users will require [Docker toolbox](https://docs.docker.com/toolbox/toolbox_install_windows/).

Ensure you have the C drive selected as a [shared drive](https://docs.docker.com/docker-for-windows/) in Docker settings (or in VirtualBox on W10 Home).

To visualise the segmentations download [ITKsnap](http://www.itksnap.org/pmwiki/pmwiki.php).

### NVIDIA driver - (only for GPU support)

If you are in a university/company, do this step with the IT support of your organization.

Download the latest NVIDIA driver version 440.33 from [NVIDIA DRIVER](https://www.nvidia.com/Download/index.aspx?lang=it).

### Download 4DSegment-2.0 Docker image

In W10 open _PowerShell_ from the Windows search box (`Win` + `X` then `I`), in macOS navigate Finder > Applications > Utilities > Terminal, or in Linux any terminal can be used. Then download the pre-compiled image:
   
```
docker pull cardiacimperial/4d_segment-2.0-cuda:latest

docker images
```

should show `cardiacimperial/4d_segment-2.0-cuda:latest` on the list of images on your local system


## Download the 4DSegment-2.0  code

### Linux - server/cloud

If you have a Linux system (recommended) you just need to:

-  Enter in your shell with ssh/PuTTY or Ctrl+Alt+T (local) 
-  Git the code with the following command: git clone https://github.com/ImperialCollegeLondon/4DSegment2.0.git


### Windows 

If you have a windows system 

- Go to https://github.com/ImperialCollegeLondon/4DSegment2.0
- Green button: Clone or Download and download the software 
- unzip 4DSegment2.0.zip folder 


## Set DMACS_docker.py file

Now enter in your code 4DSegment2.0 folder and open DMACS_docker.py - you just need to set:


```
# 1. Number of CPUs cores in your server.

coreNo        =  4

# 2. Select the device where running the neural networks (CPU/GPU) 

device        =  "GPU"

# 3. if option 2 is GPU select the device number of your GPUs (ex id_gpu = "0,1,2") 

id_gpu        =  "0"

# 4. Path to the test set directory (your_data_folder) under which images are organised in subdirectories for each subject.

test_dir      = "/cardiac/your_data_folder"
``` 

## Run 4DSegment-2.0 Docker image

### Linux - server/cloud

If you have a Linux system (recommended) you just need to:

-  Enter in your home direcotry
-  Then command pwd for discovering your local path  \<folder-path\>
-  If your data are not in a cardiac folder, please create a folder called cardiac and insert all your_data_folder - with the patient data inside
-  Then sudo nvidia-docker run -it --rm -v  /folder-path/4DSegment2.0:/workspace/code -v folder-path/cardiac:/cardiac cardiacimperial/4d_segment-2.0-cuda:latest /bin/bash


### Windows

If you have a Windows system you just need to:

- Enter in the folder that you had download before and unzip 
- Note the path to the folder on your windows desktop eg /c/Users/home/Desktop/4Dsegment-2.0, your \<folder-path\> is c/Users/home/Desktop/  
- Open a terminal window (Command Prompt or PowerShell, but not PowerShell ISE)
- If your data are not in a cardiac folder, please create a folder called cardiac and insert all your_data_folder - with the patient data inside
- Then sudo nvidia-docker run -it --rm -v  /folder-path/4DSegment2.0:/workspace/code -v folder-path/cardiac:/cardiac cardiacimperial/4d_segment-2.0-cuda:latest /bin/bash


## Inside 4DSegment-2.0 Docker


You are now inside the workspace where you will find the code folder

```
cd code
```

Then: 

```
python DMACS_docker.py
```


## Outputs from the pipeline

Once the pipeline is finished, under the root directory of each subject, you have eight nifti files, i.e., lvsa_.nii.gz, LVSA_seg_ED.nii.gz, and LVSA_seg_ES.nii.gz. However, the lvsa_.nii.gz is the preprocessing 4D raw data LVSA.nii. While LVSA_seg_ED.nii.gz and LVSA_seg_ES.nii.gz are low-resolution segmentations and seg_lvsa_SR_ED.nii.gz and seg_lvsa_SR_ES.nii.gz the high-resolution one. Note that these segmentations are smooth, high-resolution bi-ventricular three-dimensional models.

You also have meshes (txt files) for left and right ventricles at ED and ES under the root directory. For example, lv_myoed_curvature.txt records the curvature of each vertex on the myocardium of left ventricle at ED. lv_myoed_wallthickness.txt records the wall thickness of each vertex on the epicardium of the left ventricle at ED. lv_myoed_signeddistances.txt records the sign distance of each vertex on the epicardium of the left ventricle at ED, by referring to a template. lv_myoed_curvature.txt, lv_myoes_wallthickness.txt, and lv_myoes_signeddistances.txt have the same meanings for left ventricle at ES. There are also counterparts for the right ventricle at ED and ES.

In addition, the pipeline also produces the folders of dofs, segs, sizes, tmps, vtks, 4D_rview, and motion under the root directory. Apart from the motion folder where u can finally find the txt/vtk motion mesh, the files in other folders are intermediate results, which may not be useful for sequential analysis. In the motion folder, you have 50 UKBB computational both ENDO, EPI, and RV meshes (i.e vtk and txt files) for a complete cardiac cycle. In each of the 50 meshes, only spatial locations of vertices are recorded. Vertex spatial position (x, y, and z) on the same row in the txt files corresponds to the same anatomical location across the cardiac cycle.
Finally, if you want to see the rview cardiac motion you can use the file inside of the 4d_rview folder where 4Dimg.nii is the 4D grayscale file while 4Dseg.nii is the 4D segmentation.


## Release History

* 0.1.1 - with nvidia cuda support - 15/06/2020 - Nicolo Savioli, PhD

## Meta

Distributed under the GNU GENERAL PUBLIC LICENSE license. See ``LICENSE`` for more information.

Nicolo Savioli, PhD - Jinming Duan, PhD - Antonio De-Marvao, MD,PhD - Declan D'Oregon, MD,PhD
