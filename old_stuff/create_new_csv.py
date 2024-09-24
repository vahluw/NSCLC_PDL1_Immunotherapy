import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import shap
import matplotlib.pyplot as plt
import pickle


limit = "365"
extender = "100000"
folder = 'research_letter/'
all_data = np.array(np.load(folder + 'whole_dataset_0_100000.npy', allow_pickle=True))[:, :-5]
y  = np.array(np.load(folder + 'whole_dataset_0_100000.npy', allow_pickle=True))[:, -5:]

headers = ["Physician ID", "Practice ID", "Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                             "ALK+", "EGFR+", "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported",  "HER2/ERBB2+", "MET+",
                        "RET+", "NTRK1+", "NTRK2+", "NTRK3+",
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
                  "Female", "Male", "White", "Asian", "Other Race", "Hispanic Race", "Black", 'WI Residence', 'MN Residence', 'IN Residence',
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
                       "First-Line Atezolizumab Monotherapy",  "First-Line Durvalumab Monotherapy", "First-Line Ipilimumab/Nivolumab",
           "Progression Outcome", "Progression Days", "Mortality Days", "Mortality Outcome", "Time to Censor"]

all_data = pd.DataFrame(data=all_data)
categorical_indices = [5, 6, 8, 17, 18, 19, 20, 21, 35]
X_static_df_final = pd.get_dummies(all_data, columns=categorical_indices, drop_first=True, dtype=int)
X_static_df_final = pd.concat((X_static_df_final, pd.DataFrame(data=y)),axis=1)
X_static_df_final.columns = headers
X_static_df_final.to_csv(folder + 'all_data.csv')
