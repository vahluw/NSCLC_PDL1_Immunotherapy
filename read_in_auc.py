import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import shap
import matplotlib.pyplot as plt
import pickle


limit = "365"
extender = "10000"
folder = 'clf_xgb/'
io_extenders = ['1', '0']
for io_extender in io_extenders:
    test_data = np.array(np.load(folder + 'test_set_365_10000' + io_extender + '.npy', allow_pickle=True))[:, -5:]

    if io_extender == '0':
        headers = ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
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
                       "Adrenal Metastases", "Unspecified Metastases", "Glucocorticoid Use Prior to Treatment",
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
                       "First-Line Atezolizumab Monotherapy",  "First-Line Durvalumab Monotherapy", "First-Line Ipilimumab/Nivolumab"]
    else:
        headers = ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance",
                       "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                             "ALK+", "EGFR+", "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported",
                       "Days from Advanced Diagnosis to Treatment", "Therapy Year",
                       "Renal Failure", "Chronic Kidney Disease", "General Renal Disease",
                        "Prior Kidney Transplant", "Cirrhosis", "Hepatitis", "Prior Liver Transplant", "Connective Tissue Disease",
                        "Scleroderma", "Lupus", "Rheumatoid Arthritis", "Interstitial Lung Disease", "Diabetes",
                        "Bone Metastases", "Brain Metastases", "Other CNS Metastases", "Digestive System Metastases",
                       "Adrenal Metastases", "Unspecified Metastases", "Glucocorticoid Use Prior to Treatment",
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
                            "Nonsquamous Cell Carcinoma", "Never Smoker", "Previous Smoker"]

    all_dataset = pd.read_csv(folder + 'clf_xgb_prog_365_10000' + io_extender + '_train.csv')
    all_dataset = all_dataset.drop('Unnamed: 0', axis=1)
    test_dataset = pd.read_csv(folder + 'clf_xgb_prog_365_10000' + io_extender + '_test.csv')
    test_dataset = test_dataset.drop('Unnamed: 0', axis=1)

    for objective in ["prog", "mort"]:
        file = open(folder + 'clf_xgb_' + objective + '_365_10000' + io_extender + '.pickle', 'rb')

        with open(folder + 'clf_xgb_' + objective + '_365_10000' + io_extender + '.pickle', 'rb') as f:
            grid_search = pickle.load(f)

        best_estimator = grid_search.best_estimator_
        X_train_df = pd.read_csv(folder + 'clf_xgb_prog_365_10000' + io_extender + '_train.csv')
        X_test_df = pd.read_csv(folder + 'clf_xgb_prog_365_10000' + io_extender + '_test.csv')
        X_test_df = X_test_df.drop('Unnamed: 0', axis=1)
        X_train_df = X_train_df.drop('Unnamed: 0', axis=1)
        print(X_train_df.shape)
        print(best_estimator)
        explainer = shap.TreeExplainer(best_estimator, X_train_df)

        final_shap_values = explainer.shap_values(X_test_df.values)
        shap.summary_plot(final_shap_values, X_test_df.values, max_display=10, feature_names=headers, show=False)
        plt.savefig(folder + 'shap_summary_plot_ci_' + objective + '.png')
        plt.close()

        shap.summary_plot(final_shap_values, X_test_df.values, max_display=10, feature_names=headers, plot_type='bar', show=False)
        plt.savefig(folder + 'shap_summary_plot_bar_' + objective + '.png')
        plt.close()

    if io_extender == '0':
        prog_gb_preds = np.expand_dims(np.array(np.load(folder + 'y_pred_365_100000_xgb_prog_0.70_365_100000.npy')), axis=1)
        mort_gb_preds = np.expand_dims(np.array(np.load(folder + 'y_pred_365_100000_xgb_mort_0.74_365_100000.npy')), axis=1)
    else:
        prog_gb_preds = np.expand_dims(np.array(np.load(folder + 'y_pred_365_100001_xgb_prog_0.66_365_100001.npy')), axis=1)
        mort_gb_preds = np.expand_dims(np.array(np.load(folder + 'y_pred_365_100001_xgb_mort_0.66_365_100001.npy')), axis=1)

    all_test_data = np.concatenate((test_dataset.values, test_data, prog_gb_preds, mort_gb_preds), axis=1)
    print(all_test_data.shape)
    print(len(headers))
    test_set_columns = headers + ["progression_outcome",  "progression_days", "mortality_days", "mortality_outcome",
                                           "censor_days", "prog_preds", "mort_preds"]
    print(len(test_set_columns))

    all_test_data = pd.DataFrame(data=all_test_data)
    all_test_data.columns = test_set_columns
    all_test_data.to_csv('test_set_' + limit + '_' + extender + '.csv')

    progression_true_outcomes = all_test_data['progression_outcome'].values.astype(int)
    mortality_true_outcomes = all_test_data['progression_outcome'].values.astype(int)

    fpr, tpr, thresholds = roc_curve(progression_true_outcomes, prog_gb_preds)
    # get the best threshold
    J = tpr - fpr
    ix = np.argmax(J)
    best_thresh = thresholds[ix]
    print('Best Threshold For Temporal =%f' % (best_thresh))


    fpr, tpr, thresholds = roc_curve(mortality_true_outcomes, mort_gb_preds)
    # get the best threshold
    J = tpr - fpr
    ix = np.argmax(J)
    best_thresh = thresholds[ix]
    print('Best Threshold For Temporal =%f' % (best_thresh))

    del all_dataset
    in_test_set = []
    test_dataset = all_test_data.values
    whole_dataset = pd.read_csv(folder + 'all_x_static_365_10000' + io_extender + '.csv').values
    test_dataset_list = test_dataset.tolist()
    for i in range(whole_dataset.shape[0]):
        row = whole_dataset[i].tolist()
        if row in test_dataset_list:
            in_test_set.append(1)
        else:
            in_test_set.append(0)

    in_test_set = np.array(in_test_set)
    in_test_set = np.expand_dims(in_test_set, axis=1)
    headers_all = ["physicianid", "practiceid"] + test_set_columns[:-2]

    whole_dataset_pickle = np.array(np.load(folder + 'whole_dataset_365_10000' + io_extender + '.npy', allow_pickle=True))
    whole_dataset_labels  = whole_dataset_pickle[:, -5:]
    whole_dataset_ids  = whole_dataset_pickle[:, 0:2]

    all_data = np.concatenate((whole_dataset_ids, whole_dataset[:, 1:], whole_dataset_labels), axis=1)
    data = pd.DataFrame(data=all_data)
    data.columns = headers_all
    data.to_csv('all_data_' + limit + '_' + extender + '.csv')

