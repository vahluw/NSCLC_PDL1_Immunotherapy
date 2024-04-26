# This is a sample Python script.

# Press ^R to execute it or replace it with your code.
# Press Double ? to search everywhere for classes, files, tool windows, actions, and settings.
#import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve

new_therapy = 1

all_dataset = np.array(np.load('whole_dataset_progression_0.005_365_01.npy'))
test_dataset = np.array(np.load('test_set_progression_0.005_365_01.npy'))
#gb_preds = np.expand_dims(np.array(np.load('y_pred_gb_0.81_0110111.npy')), axis=1)
#xgb_preds = np.expand_dims(np.array(np.load('y_pred_xgb_0.81_0110111.npy')), axis=1)
#rf_preds = np.expand_dims(np.array(np.load('y_pred_rf_0.80_0110111.npy')), axis=1)
#knn_preds = np.expand_dims(np.array(np.load('y_pred_knn_0.73_0110111.npy')), axis=1)
#logistic_preds = np.expand_dims(np.array(np.load('y_pred_logistic_0.63_0110111.npy')), axis=1)
all_test_data = test_dataset #np.concatenate((test_dataset, gb_preds, xgb_preds, rf_preds, knn_preds, logistic_preds), axis=1)
if new_therapy == 1:
    headers_test_set = ["days_from_dx_to_tx", "contra1", "contra2", "contra3", "contra4", "contra5", "contra6", "bilirubin", "creatinine", "AST", "ALT", "diag_year",
                    "age_at_diagnosis", "birth_year", "gender", "race", "ethnicity", "state", "other_no_insurance",
                    "workers_comp","self_pay","pt_assistance","other_gov_insurance","medicare", "medicaid",
                    "commercial_health_plan", "practice_ID", "practice_type", "physician_ID", "histology",
                "stage", "smoking_status","ecog",  "ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given",
                  "io_mono", "io_mono_used", "combo_therapy", "first_line_chemo", "secondary_chemo_drug", "other_therapy",
                    "alk_drug", "egfr_drug", "braf_drug", "ros1_drug", "ras_drug", "other_first_line_therapy", "no_first_line",
                    "progression_12",  "progression_days", "censor_days"]#, "gb_preds", "xgb_preds", "rf_preds", "knn_preds", "logistic_preds"]
else:
    headers_test_set = ["contra1", "contra2", "contra3", "contra4", "contra5", "contra6", "bilirubin", "diag_year", "AST", "ALT", "diag_year",
                    "age_at_diagnosis", "birth_year", "gender", "race", "ethnicity", "state", "other_no_insurance",
                    "workers_comp","self_pay","pt_assistance","other_gov_insurance","medicare", "medicaid",
                    "commercial_health_plan", "practice_ID", "practice_type", "physician_ID", "histology",
                "stage", "smoking_status","ecog",  "ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given",
                   "io_mono", "combo_therapy", "chemo",  "other_therapy", "io_mono_used",
                    "progression_12",  "progression_days", "censor_days", "gb_preds", "xgb_preds", "rf_preds", "knn_preds", "logistic_preds"]


data = pd.DataFrame(data=all_test_data)

data.columns = headers_test_set
data.to_csv('test_set.csv')

'''
fpr, tpr, thresholds = roc_curve(data['progression_12'], gb_preds)
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
headers = ["in_test_set", "days_from_dx_to_tx", "contra1", "contra2", "contra3", "contra4", "contra5", "contra6", "bilirubin", "creatinine", "AST", "ALT", "diag_year",
                    "age_at_diagnosis", "birth_year", "gender", "race", "ethnicity", "state", "other_no_insurance",
                    "workers_comp","self_pay","pt_assistance","other_gov_insurance","medicare", "medicaid",
                    "commercial_health_plan", "practice_ID", "practice_type", "physician_ID", "histology",
                "stage", "smoking_status","ecog",  "ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given",
                  "io_mono", "io_mono_used", "combo_therapy", "first_line_chemo", "secondary_chemo_drug", "other_therapy",
                    "alk_drug", "egfr_drug", "braf_drug", "ros1_drug", "ras_drug", "other_first_line_therapy", "no_first_line",
                    "progression_12",  "progression_days", "censor_days"]

all_data = np.concatenate((in_test_set, all_dataset), axis=1)
data = pd.DataFrame(data=all_data)
data.columns = headers
data.to_csv('all_data.csv')

