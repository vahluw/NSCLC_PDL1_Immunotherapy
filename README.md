# Explainable machine learning to predict treatment response in advanced non-small-cell lung cancer

## Why use this?

Immune checkpoint inhibitors (ICIs) have demonstrated promise in the treatment of various cancers. Single-drug ICI therapy (IO monotherapy) that targets programmed death ligand-1 (PD-L1) is the standard of care in patients with advanced non-small-cell lung cancer (NSCLC) with PD-L1 expression ≥ 50%. We sought to find out if a machine learning (ML) algorithm can perform better as a predictive biomarker than PD-L1 alone. Using a real-world, nationwide electronic health record (EHR)-derived de-identified database of 38048 patients with advanced NSCLC, we trained binary prediction algorithms to predict likelihood of 12-month progression-free survival (12-month PFS) and 12-month overall survival (12-month OS) from initiation of first-line therapy. We evaluated the algorithms by calculating the area under the receiver operator curve (AUC) on the test set. We plotted Kaplan-Meier curves and fit Cox survival models comparing survival between patients who were classified as “low-risk” for 12-month disease progression or 12-month mortality versus those classified as “high-risk.” 

In order for a clinician to use this algorithm, we envision that they would input the necessary variables into a graphical user interface that is either hosted on a third-party smartphone application or online. The clinician would indicate which first-line therapy they are considering (i.e., IO monotherapy, doublet chemotherapy, combination therapy, etc.) and the algorithm would predict a likelihood of 12-month progression free survival and 12-month overall survival by outputting a continuous value between 0 and 1. Using the respective probability outputted by the software, the clinician can determine if they wish to continue with the therapy in question or if they should consider a different treatment. 


## Documentation

Please cite our work as follows:

Ahluwalia VS, Parikh, RB. Explainable machine learning to predict treatment response in advanced non-small-cell lung cancer. Accepted at JCO Clin Cancer Inform. 2024.

### Contents of Repository
    1. ml_model_OS.pickle  --> Machine learning model designed to predict likelihood of 12-month mortality for any first-line therapy type.
    2. ml_model_OS_io_monotherapy_only.pickle --> Machine learning model designed to predict likelihood of 12-month mortality for first-line IO monotherapy only.
    3. ml_model_PFS.pickle --> Machine learning model designed to predict likelihood of 12-month disease progression for any first-line therapy type.
    4. ml_model_PFS_io_monotherapy.pickle --> Machine learning model designed to predict likelihood of 12-month disease progression for first-line IO monotherapy only.
    5. NSCLC.py --> Code used to extract data from raw datasets, perform machine learning model training, and validate model.
    6. read_in.py --> Code used to create Shapley plots and output dataset to readable csv format.

The necessary libraries that must be installed to run this code include copy, pickle, sys, time, datetime, decimal, matplotlib, numpy, pandas, shap, xgboost, and sklearn. After installing necessary dependencies, users can use the prediction powers of each of the four models described in our paper.


The data that support the findings of this study have been originated by Flatiron Health, Inc. Requests for data sharing by license or by permission for the specific purpose of replicating results in this manuscript can be submitted to dataaccess@flatiron.com. 


### Installation and Running ML Models
```bash
(base) $> pip install sklearn                           # Install necessary dependencies if needed
(base) $> pip install xgboost
(base) $> pip install shap
(base) $> pip install datetime

