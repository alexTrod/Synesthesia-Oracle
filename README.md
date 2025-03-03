# Synesthesia-Oracle

## Predicting Synesthetic Experience from Resting State fMRI: A graph-based Deep Learning approach

Synesthesia is a phenomenon that causes sensory crossovers, such as tasting colors or seeing sounds. While Magnetic Resonance Imaging (MRI) and functional MRI (fMRI) are great modalities to map functional and structural properties of the brain; currently, the detection of synesthetic experience exclusively based on behavioral tasks. In this paper, we construct various models that can predict the presence of synesthesia from resting state fMRI recordings and aim to facilitate the identification of neurological biomarkers underlying this condition. To this end, we employ Graph Neural Networks (GNNs) as it is the most natural choice for the fMRI data. We also integrate Kolmogorov Arnold Networks (KANs) to explore its potential for interpreting the results.

This is the code repository for the corresponding conference paper.

## Dataset

In order to replicate this research project you need to download the MSMall_data_all.mat file from the link: https://osf.io/v7px5/. Next, you need to run the code (https://osf.io/7gevx for the original, a slightly modified version is present in the current repository) in matlab to get the files fullcor.mat and parcor.mat.

For matlab, you also need the folder fslnets (http://www.fmrib.ox.ac.uk/~steve/ftp/fslnets.tar.gz). Unzip it and add it to the current matlab path.

## Methods

The following jupyter notebooks are provided for the respective tasks:

synesthesia 1.ipynb -> data preprocessing

synesthesia 2.ipynb -> vanilla gcn, braingnn, no kan

synesthesia 3.ipynb -> braingnn with kan

## Required Files and Repositories

The anaconda environment is exported in the file: environment.yml
To create an identical environment run: 
conda env create -f environment.yml
conda activate syn_oracle

Please clone the relevant original github repositories:

BrainGNN:   https://github.com/xxlya/BrainGNN_Pytorch

KAN:        https://github.com/KindXiaoming/pykan

GraphKAN:   https://github.com/WillHua127/GraphKAN-Graph-Kolmogorov-Arnold-Networks

Folder tree required:

```

syn_oracle
│   braingnn_plus_kan.py
│   BrainGNN_with_KAN.py
│   environment.yml             (conda environment with all required packages)
│   README.md
│   synaesthesia 1.ipynb        (jupyter notebook for data preprocessing)
│   synaesthesia 2.ipynb        (jupyter notebook for vanilla gcn, braingnn, no kan)
│   synaesthesia 3.ipynb        (jupyter notebook for braingnn with kan and kan stand alone)
│
├───BrainGNN_Pytorch        (https://github.com/xxlya/BrainGNN_Pytorch) a slightly modified version is provided
│
├───data
│       fullcor.mat             
│       MSMall_data_all.mat     (https://osf.io/v7px5/)
│       parcor.mat
│
├───gnn_model           (where the model weights are saved)
│
├───GraphKAN            (https://github.com/WillHua127/GraphKAN-Graph-Kolmogorov-Arnold-Networks)
│
├───log                 (where the tensorboardX logs are saved)
│
├───matlab
│   │   CISC_IDs.xlsx                       (https://osf.io/bsahz)
│   │   fullcor.mat                         (generated by matlab)
│   │   MSMall_data_all.mat                 (https://osf.io/v7px5)
│   │   parcel_names_N360.xlsx              (https://osf.io/znh32)
│   │   parcor.mat                          (generated by matlab)
│   │   wcs_partialcorr_finalsample.m       (https://osf.io/7gevx)
│   │   partialcorr_mod.m                   (slightly modified version of wcs_partialcorr_finalsample.m)
│   │
│   └───fslnets                             (http://www.fmrib.ox.ac.uk/~steve/ftp/fslnets.tar.gz)
│       
│
└───pykan       (https://github.com/KindXiaoming/pykan)

```