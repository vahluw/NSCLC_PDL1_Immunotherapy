# This is a sample Python script.

# Press ^R to execute it or replace it with your code.
# Press Double ? to search everywhere for classes, files, tool windows, actions, and settings.
#import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve

limit = "365"
extender = "1001"
all_dataset = np.array(np.load('whole_dataset_' + limit +'_' + extender + '.npy'))
test_dataset = np.array(np.load('test_set_' + limit + '_' + extender + '.npy'))
prog_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_1001_hgb_prog_1.00_365_1001.npy')), axis=1)
#mort_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_100_hgbmort_0.74_365_100.npy')), axis=1)

all_test_data = np.concatenate((test_dataset, prog_gb_preds), axis=1)
headers_test_set = [ "physicianID", "practiceID",  "diag_year", "age_at_diagnosis", "birth_year", "gender", "race",
                    "ethnicity", "state", "other_no_insurance","workers_comp","self_pay","pt_assistance",
                    "other_gov_insurance","medicare", "medicaid", "commercial_health_plan", "practice_type",  "ecog", "stage",
                    "histology", "smoking_status",  "ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given",
                   "io_mono_used", "combo_therapy", "first_line_chemo", "secondary_chemo_drug", "alk_drug", "egfr_drug",
                    "braf_drug", "ros1_drug", "ras_drug", "other_first_line_therapy", "clinical_study_drug", "bev_used",
                   "three_plus_chemo_drugs",  "carboplatin_only", "cisplatin_only", "pembrolizumab_used", "trk_inhibitor",
                        "met_drug", "days_from_dx_to_tx", "therapy_year", "kidney_failure", "chronic_kidney_disease", "renal_disease",
                    "kidney_transplant", "cirrhosis", "hepatitis", "liver_transplant", "connective_tissue",
                    "scleroderma", "lupus", "rheumatoid_arthritis", "interstitial_lung_disease", "diabetes",
                    "bone_mets", "brain_mets", "cns_mets", "digestive_mets", "adrenal_mets", "unspecified_mets","steroid",
                     "abx", "albumin", "creatinine", "bilirubin", "ast", "alt", "progression_outcome",
                     "progression_days", "mortality_days", "mortality_outcome", "censor_days", "prog_gb_preds"]



data = pd.DataFrame(data=all_test_data)

data.columns = headers_test_set
data.to_csv('test_set_' + limit + '_' + extender + '.csv')


fpr, tpr, thresholds = roc_curve(data['progression_outcome'], prog_gb_preds)
# get the best threshold
J = tpr - fpr
ix = np.argmax(J)
best_thresh = thresholds[ix]
print('Best Threshold For Temporal =%f' % (best_thresh))

'''
fpr, tpr, thresholds = roc_curve(data['progression_outcome'], mort_gb_preds)
# get the best threshold
J = tpr - fpr
ix = np.argmax(J)
best_thresh = thresholds[ix]
print('Best Threshold For Temporal =%f' % (best_thresh))
'''

in_test_set = []
test_dataset_list = test_dataset.tolist()
for i in range(all_dataset.shape[0]):
    row = all_dataset[i].tolist()
    if row in test_dataset_list:
        in_test_set.append(1)
    else:
        in_test_set.append(0)

in_test_set = np.array(in_test_set)
in_test_set = np.expand_dims(in_test_set, axis=1)
headers = ["in_test_set", "physicianID", "practiceID",  "diag_year", "age_at_diagnosis", "birth_year", "gender", "race",
                    "ethnicity", "state", "other_no_insurance","workers_comp","self_pay","pt_assistance",
                    "other_gov_insurance","medicare", "medicaid", "commercial_health_plan", "practice_type",  "ecog", "stage",
                    "histology", "smoking_status",  "ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given",
                   "io_mono_used", "combo_therapy", "first_line_chemo", "secondary_chemo_drug", "alk_drug", "egfr_drug",
                    "braf_drug", "ros1_drug", "ras_drug", "other_first_line_therapy", "clinical_study_drug", "bev_used",
                   "three_plus_chemo_drugs",  "carboplatin_only", "cisplatin_only", "pembrolizumab_used", "trk_inhibitor",
                        "met_drug", "days_from_dx_to_tx", "therapy_year", "kidney_failure", "chronic_kidney_disease", "renal_disease",
                    "kidney_transplant", "cirrhosis", "hepatitis", "liver_transplant", "connective_tissue",
                    "scleroderma", "lupus", "rheumatoid_arthritis", "interstitial_lung_disease", "diabetes",
                    "bone_mets", "brain_mets", "cns_mets", "digestive_mets", "adrenal_mets", "unspecified_mets",
           "steroid", "abx", "albumin", "creatinine", "bilirubin", "ast", "alt",
                    "progression_outcome",  "progression_days", "mortality_days", "mortality_outcome", "censor_days"]


all_data = np.concatenate((in_test_set, all_dataset), axis=1)
data = pd.DataFrame(data=all_data)
data.columns = headers
data.to_csv('all_data_' + limit + '_' + extender + '.csv')

