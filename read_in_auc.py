import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import shap
import matplotlib.pyplot as plt
import pickle

limit = "365"
extender = "1000"
all_dataset = np.array(np.load('whole_dataset_' + limit +'_' + extender + '.npy', allow_pickle=True))
test_dataset = np.array(np.load('test_set_' + limit + '_' + extender + '.npy', allow_pickle=True))
file = open('clf_xgb_prog_365_1000.pickle', 'rb')
with open("clf_hgb_prog_365_1000.pickle", 'rb') as f:
    best_estimator = pickle.load(f)
X_train_df = pd.read_csv('clf_xgb_prog_365_1000_train.csv')
X_test_df = pd.read_csv('clf_xgb_prog_365_1000_test.csv')
print(best_estimator)
explainer = shap.TreeExplainer(best_estimator, X_train_df)

final_shap_values = explainer.shap_values(X_test_df.values)
shap.summary_plot(final_shap_values, X_test_df.values, max_display=20, feature_names=X_test_df.headers, show=False)
plt.savefig('shap_summary_plot_ci_.png')
plt.close()

shap.summary_plot(final_shap_values, X_test_df.values, max_display=20, feature_names=X_test_df.headers, plot_type='bar', show=False)
plt.savefig('shap_summary_plot_bar.png')
plt.close()
exit(0)

prog_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_1000_xgb_prog_0.70_365_1000.npy')), axis=1)
#mort_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_1000_xgb_mort_0.74_365_1000.npy')), axis=1)

all_test_data = np.concatenate((test_dataset, prog_gb_preds), axis=1)

data = pd.DataFrame(data=all_test_data)

data.columns = headers_test_set
data.to_csv('test_set_' + limit + '_' + extender + '.csv')


fpr, tpr, thresholds = roc_curve(data['progression_outcome'], prog_gb_preds)
# get the best threshold
J = tpr - fpr
ix = np.argmax(J)
best_thresh = thresholds[ix]
print('Best Threshold For Temporal =%f' % (best_thresh))


fpr, tpr, thresholds = roc_curve(data['progression_outcome'], mort_gb_preds)
# get the best threshold
J = tpr - fpr
ix = np.argmax(J)
best_thresh = thresholds[ix]
print('Best Threshold For Temporal =%f' % (best_thresh))


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

