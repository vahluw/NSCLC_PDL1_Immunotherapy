import copy
import pickle
import sys
import time
from datetime import date
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import shap
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.utils import class_weight
from sklearn.utils import shuffle
import keras
from keras import layers, callbacks, regularizers

ts = np.load('/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/ablate/test_set_365_ablation_1000101.npy')
pred_mort = np.load('/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/ablate/y_pred_365_ablation_1000101_xgb_mort_0.593_365_ablation_1000101.npy')
pred_prog = np.load('/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/ablate/y_pred_365_ablation_1000101_xgb_prog_0.681_365_ablation_1000101.npy')
pred_mort = np.expand_dims(pred_mort, axis=1)
pred_prog = np.expand_dims(pred_prog, axis=1)

headers_all_mort=  ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                             "ALK+", "EGFR+", "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported",
                             "First-Line Combination Therapy", "First-Line Chemotherapy", "Non-First-Line Chemotherapy",
                       "Anti-ALK Drug", "Anti-EGFR Drug","Anti-BRAF Drug", "Anti-ROS1 Drug", "Anti-RAS Drug",
                       "Other First-Line Therapy", "Clinical Study Drug Used", "Bevacizumab Used",
                       "Three Or More Chemotherapy Drugs",  "Carboplatin Monotherapy", "Cisplatin Monotherapy", "Pembrolizumab Used",
                       "TRK Inhibitor", "MET Inhibitor", "Days from Advanced Diagnosis to Treatment", "Therapy Year",
                       "Renal Failure", "Chronic Kidney Disease", "General Renal Disease",
                        "Prior Kidney Transplant", "Cirrhosis", "Hepatitis", "Prior Liver Transplant", "Connective Tissue Disease",
                        "Scleroderma", "Lupus", "Rheumatoid Arthritis", "Interstitial Lung Disease", "Diabetes",
                        "Bone Metastases", "Brain Metastases", "Other CNS Metastases", "Digestive System Metastases",
                       "Adrenal Metastases", "Unspecified Metastases", "SDH", "Glucocorticoid Use Prior to Treatment",
                         "Anti-Infective Use Prior to Treatment", "Albumin", "Creatinine", "Bilirubin", "AST", "ALT",
                     "Male", "White", "Asian", "Other Race", "Hispanic Race", "Black", 'WI Residence', 'MN Residence', 'IN Residence',
                            'VA Residence', 'PR Residence', 'DC Residence', 'UT Residence', 'ID Residence', 'MO Residence',
                       'CT Residence', 'NH Residence', 'CA Residence', 'AR Residence', 'NV Residence', 'DE Residence',
                       'MD Residence', 'TN Residence', 'AL Residence', 'NJ Residence', 'PA Residence', 'NY Residence',
                       'NE Residence', 'WA Residence', 'WV Residence', 'AZ Residence', 'LA Residence', 'OR Residence',
                       'OK Residence', 'TX Residence', 'CO Residence', 'IA Residence', 'MS Residence', 'RI Residence',
                       'OH Residence', 'SC Residence', 'GA Residence', 'MI Residence', 'NC Residence', 'ME Residence',
                            'FL Residence', 'IL Residence', 'NM Residence', 'HI Residence', 'KS Residence', 'KY Residence',
                       'MA Residence', "Academic Medical Center", "ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3",
                       "ECOG 4", 'Stage 0','Stage I', 'Stage IA', 'Stage IA1', 'Stage IA2', 'Stage IA3', 'Stage IB',
                       'Stage II', 'Stage IIA', 'Stage IIB', 'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC',
                       'Stage IV', 'Stage IVA', 'Stage IVB', 'Occult',  "Squamous Cell Carcinoma",
                            "Nonsquamous Cell Carcinoma", "Never Smoker", "Previous Smoker", "First-Line Nivolumab Monotherapy",
                       "First-Line Pembrolizumab Monotherapy", "First-Line Cemiplimab Monotherapy",
                       "First-Line Atezolizumab Monotherapy",  "First-Line Durvalumab Monotherapy", "First-Line Ipilimumab/Nivolumab"]
headers_all_prog =  ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                             "ALK+", "EGFR+", "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported",
                             "First-Line Combination Therapy", "First-Line Chemotherapy", "Non-First-Line Chemotherapy",
                       "Anti-ALK Drug", "Anti-EGFR Drug","Anti-BRAF Drug", "Anti-ROS1 Drug", "Anti-RAS Drug",
                       "Other First-Line Therapy", "Clinical Study Drug Used", "Bevacizumab Used",
                       "Three Or More Chemotherapy Drugs",  "Carboplatin Monotherapy", "Cisplatin Monotherapy", "Pembrolizumab Used",
                       "TRK Inhibitor", "MET Inhibitor", "Days from Advanced Diagnosis to Treatment", "Therapy Year",
                       "Renal Failure", "Chronic Kidney Disease", "General Renal Disease",
                        "Prior Kidney Transplant", "Cirrhosis", "Hepatitis", "Prior Liver Transplant", "Connective Tissue Disease",
                        "Scleroderma", "Lupus", "Rheumatoid Arthritis", "Interstitial Lung Disease", "Diabetes",
                        "Bone Metastases", "Brain Metastases", "Other CNS Metastases", "Digestive System Metastases",
                       "Adrenal Metastases", "Unspecified Metastases", "SDH", "Glucocorticoid Use Prior to Treatment",
                         "Anti-Infective Use Prior to Treatment", "Albumin", "Creatinine", "Bilirubin", "AST", "ALT",
                     "Male", "White", "Asian", "Other Race", "Hispanic Race", "Black", 'WI Residence', 'MN Residence', 'IN Residence',
                            'VA Residence', 'PR Residence', 'DC Residence', 'UT Residence', 'ID Residence', 'MO Residence',
                       'CT Residence', 'NH Residence', 'CA Residence', 'AR Residence', 'NV Residence', 'DE Residence',
                       'MD Residence', 'TN Residence', 'AL Residence', 'NJ Residence', 'PA Residence', 'NY Residence',
                       'NE Residence', 'WA Residence', 'WV Residence', 'AZ Residence', 'LA Residence', 'OR Residence',
                       'OK Residence', 'TX Residence', 'CO Residence', 'IA Residence', 'MS Residence', 'RI Residence',
                       'OH Residence', 'SC Residence', 'GA Residence', 'MI Residence', 'NC Residence', 'ME Residence',
                            'FL Residence', 'IL Residence', 'NM Residence', 'HI Residence', 'KS Residence', 'KY Residence',
                       'MA Residence', "Academic Medical Center", "ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3",
                       "ECOG 4", 'Stage 0','Stage I', 'Stage IA', 'Stage IA1', 'Stage IA2', 'Stage IA3', 'Stage IB',
                       'Stage II', 'Stage IIA', 'Stage IIB', 'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC',
                       'Stage IV', 'Stage IVA', 'Stage IVB', 'Occult',  "Squamous Cell Carcinoma",
                            "Nonsquamous Cell Carcinoma", "Never Smoker", "Previous Smoker", "First-Line Nivolumab Monotherapy",
                       "First-Line Pembrolizumab Monotherapy", "First-Line Cemiplimab Monotherapy",
                       "First-Line Atezolizumab Monotherapy",  "First-Line Durvalumab Monotherapy", "First-Line Ipilimumab/Nivolumab"]

for remove_item in ['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'Diagnosis Year', 'Albumin', 'ALK+', 'EGFR+']:
    headers_all_prog.remove(remove_item)

for remove_item in ['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'PDL1 Reported', 'Albumin', 'ALK+', 'EGFR+']:
    headers_all_mort.remove(remove_item)

concat_prog = np.concatenate((ts[:, 76:80],  pred_prog, pred_mort), axis=1)
df_prog = pd.DataFrame(data=concat_prog)
df_prog.columns = ['y_prog', 'prog_days', 'mort_days','y_mort', 'pred_prog', 'pred_mort']
df_prog.to_csv('ablate_roc.csv')

