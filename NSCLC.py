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

no_progression_date = date(2024, 12, 31)
dir_path1 = '/origdata/Parikh_Flatirons/updateJuly2024/'
dir_path2 = '/Users/vahluw/Downloads/NSCLC_Updated/'

headers_all =  ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                             "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported", "HER2/ERBB2+", "MET+",
                        "RET+", "NTRK1+", "NTRK2+", "NTRK3+", "NTRK_Other",
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

headers_to_drop = ["First-Line Combination Therapy", "First-Line Chemotherapy", "Non-First-Line Chemotherapy",
                       "Anti-ALK Drug", "Anti-EGFR Drug","Anti-BRAF Drug", "Anti-ROS1 Drug", "Anti-RAS Drug",
                       "Other First-Line Therapy", "Clinical Study Drug Used", "Bevacizumab Used",
                       "Three Or More Chemotherapy Drugs",  "Carboplatin Monotherapy", "Cisplatin Monotherapy", "Pembrolizumab Used",
                       "TRK Inhibitor", "MET Inhibitor",  "First-Line Nivolumab Monotherapy",
                       "First-Line Pembrolizumab Monotherapy", "First-Line Cemiplimab Monotherapy",
                       "First-Line Atezolizumab Monotherapy",  "First-Line Durvalumab Monotherapy", "First-Line Ipilimumab/Nivolumab"]

def perform_grid_search(param_grid, clf, X_train, y_train, X_test, y_test, filename_extender, type='rf', weights=None, io_only=0):

    #shap.initjs()
    X_train_df = pd.DataFrame(data=X_train)
    X_test_df = pd.DataFrame(data=X_test)
    X_train_df = X_train_df.astype(float)
    X_test_df = X_test_df.astype(float)
    if io_only == 0:
        headers_long = ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance", "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                              "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported", "HER2/ERBB2+", "MET+",
                        "RET+", "NTRK1+", "NTRK2+", "NTRK3+", "NTRK_Other",
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
    else:
        headers_long = ["Diagnosis Year", "Age At Diagnosis", "Birth Year",  "Hispanic Ethnicity", "No Insurance",
                       "Worker's Compensation ", "Self-Pay", "Patient Assistance Program",
                       "Other Governmental Insurance", "Medicare", "Medicaid", "Commercial Health Plan",
                              "KRAS+",  "ROS1+", "BRAF+", "PDL1", "PDL1 Reported", "HER2/ERBB2+", "MET+",
                        "RET+", "NTRK1+", "NTRK2+", "NTRK3+", "NTRK_Other",
                       "Days from Advanced Diagnosis to Treatment", "Therapy Year",
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
                            "Nonsquamous Cell Carcinoma", "Never Smoker", "Previous Smoker"]
    print(len(headers_long))

    # Initialize GridSearchCV
    grid_search = GridSearchCV(estimator=clf, param_grid=param_grid, cv=5, scoring='roc_auc', n_jobs=-1)

    # Perform grid search
    if type=='xgb':
        grid_search.fit(X_train_df, y_train, sample_weight=weights)
    else:
        grid_search.fit(X_train_df, y_train)

    # Best parameters
    best_params = grid_search.best_params_
    print(f"Best parameters: {best_params}")

    # Best estimator
    best_estimator = grid_search.best_estimator_
    print(f"Best estimator: {best_estimator}")

    # Predictions
    y_pred_prob = best_estimator.predict_proba(X_test)[:,1]
    auc_final = roc_auc_score(y_test_final, y_pred_prob)
    auc_dec = Decimal(str(auc_final))
    rounded_num = round(auc_dec, 3)
    print("AUC: ", str(rounded_num))
    np.save('y_pred_' + filename_extender + '_' + type + '_' + str(rounded_num) + '_' + file_name_extender + '.npy', y_pred_prob)
    y_pred = best_estimator.predict(X_test_df)

    # Evaluate model
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Accuracy: {accuracy:.4f}")
    with open("clf_" + type + '_' + filename_extender + ".pickle", 'wb') as f:
        pickle.dump(grid_search, f)

    if use_dynamic == 0 and exclude_diagnoses==1:
        X_train_df.headers = headers_long
        X_test_df.headers = headers_long

    X_train_df.to_csv("clf_" + type + '_' + filename_extender + "_train.csv")
    X_test_df.to_csv("clf_" + type + '_' + filename_extender + "_test.csv")
    return

def convert_date_to_iso(orig_date):
    ind1 = orig_date.find('/')
    mo = int(orig_date[:ind1])
    next = orig_date[ind1+1:]
    ind2 = next.find('/')
    day = int(next[:ind2])
    year = int(next[ind2+1:])
    if year <= 24:
        year = 2000 + year
    else:
        year = 1900 + year

    return date(year, mo, day)

def update_last_note(patientID, dictionary, new_date):
    if isinstance(new_date, str):
        try:
            modified_date = date.fromisoformat(new_date)
        except:
            modified_date = convert_date_to_iso(new_date)
    else:
        modified_date = new_date

    if modified_date >= no_progression_date:
        return dictionary

    if patientID not in dictionary:
        dictionary[patientID] = modified_date
        return dictionary
    else:
        orig = dictionary[patientID]
        if modified_date > orig:
            dictionary[patientID] = modified_date
        return dictionary

def perform_imputation_df(df, columns, type='mode'):
    print(df)
    # Replace 0 values with NaN in the specified columns
    df[columns] = df[columns].replace(0, pd.NA)

    if type == 'mode':
        # Perform mode imputation
        return df.fillna(df.mode().iloc[0])
    else:
        return df.fillna(df.mean())

def find_key_by_value(dictionary, value):
    for key, val in dictionary.items():
        if isinstance(key, float) and np.isnan(key):
            return value
        if val == value:
            return value
    return None  # Return None if the value is not found in the dictionary

def create_set(data, start_at_zero=False):
    items = set()

    for i in range(len(data)):
        items.add(data[i])

    entry_to_number = dict()

    if start_at_zero:
        counter = 0
    else:
        counter = 1

    for item in items:
        entry_to_number[item] = counter
        counter += 1

    return entry_to_number

def get_meds_value(med_dose_temp, med_units_temp, count_temp, previous_val):
    if isinstance(med_dose_temp, str) and ('%' in med_dose_temp or 'L' in med_dose_temp or 'H' in med_dose_temp):
        return float(med_dose_temp[:-1])
    elif count_temp >= 1 and (med_dose_temp == 0.0 or med_dose_temp == 'nan' or np.isnan(med_dose_temp) or med_dose_temp < 0):
        return previous_val
    elif count_temp == 0 and (med_dose_temp == 0.0 or med_dose_temp == 'nan' or np.isnan(med_dose_temp) or med_dose_temp < 0):
        return 0.0
    else:
        if med_units_temp == 'g':
            return med_dose_temp * 1000.0
        elif med_units_temp == 'ug' or med_units_temp == 'ug/hr' or med_units_temp == 'mcg' or med_units_temp == 'mcg/kg/min':
            return med_dose_temp/1000.0
        else:
            return med_dose_temp * 1.0


def get_labs_value(lab_value_temp, count_temp, previous_val):
    try:
        if isinstance(lab_value_temp, str) and ('%' in lab_value_temp or 'L' in lab_value_temp or 'H' in lab_value_temp or '+' in lab_value_temp):
            return float(lab_value_temp[:-1])
        elif isinstance(lab_value_temp, str) and lab_value_temp == "Negative":
            return 0.0
        elif isinstance(lab_value_temp, str) and lab_value_temp == "Positive":
            return 1.0
        elif count_temp >= 1 and (lab_value_temp == 0.0 or lab_value_temp == 'nan' or np.isnan(lab_value_temp)):
            return previous_val
        elif count_temp == 0 and (lab_value_temp == 0.0 or lab_value_temp == 'nan' or np.isnan(lab_value_temp)):
            return 0.0
        elif lab_value_temp < 0 and count_temp >= 1:
            return previous_val
        elif lab_value_temp < 0 and count_temp == 0:
            return 0.0
        else:
            return lab_value_temp
    except:
        if count_temp >= 1:
            return previous_val
        else:
            return 0.0


def get_vitals_value(value_temp, units_temp, count_temp, previous_val):

    try:
        if isinstance(value_temp, str) and "nan" in value_temp:
            if count_temp == 0:
                return 0.0
            else:
                return previous_val
        elif isinstance(value_temp, str) and ('%' in value_temp or '.' in value_temp):
            try:
                value_temp = float(value_temp)
                if value_temp < 0 and count_temp == 0:
                    return 0.0
                elif value_temp < 0 and count_temp >= 1:
                    return previous_val
                elif units_temp == 'oz':                     # Convert ounces to kg
                    return value_temp * 0.0283
                elif units_temp == 'lb' or units_temp == 'lbs':   # Convert lbs to kg
                    return value_temp * 0.454
                elif units_temp == 'in':                     # Convert inches to cm
                    return value_temp * 2.54
                else:
                    return value_temp * 1.0
            except:
                value_temp = float(value_temp[:-1])
                if value_temp < 0 and count_temp == 0:
                    return 0.0
                elif value_temp < 0 and count_temp >= 1:
                    return previous_val
                elif units_temp == 'oz':                     # Convert ounces to kg
                    return value_temp * 0.0283
                elif units_temp == 'lb' or units_temp == 'lbs':   # Convert lbs to kg
                    return value_temp * 0.454
                elif units_temp == 'in':                     # Convert inches to cm
                    return value_temp * 2.54
                else:
                    return value_temp * 1.0
        elif isinstance(value_temp, str) and 'L' in value_temp:
            return float(value_temp[:-1])
        elif isinstance(value_temp, float) and np.isnan(value_temp) and count_temp >= 1:
            return previous_val
        elif isinstance(value_temp, float) and np.isnan(value_temp) and count_temp == 0:
            return 0.0
        elif isinstance(value_temp, str) and value_temp == 'nan' and previous_val != 0:
            return previous_val
        elif isinstance(value_temp, str) and value_temp == 'nan':
            return 0.0
        else:
            value_temp = float(value_temp)
            if value_temp < 0 and count_temp == 0:
                return 0.0
            elif value_temp < 0 and count_temp >= 1:
                return previous_val
            elif units_temp == 'oz':                     # Convert ounces to kg
                return value_temp * 0.0283
            elif units_temp == 'lb' or units_temp == 'lbs':   # Convert lbs to kg
                return value_temp * 0.454
            elif units_temp == 'in':                     # Convert inches to cm
                return value_temp * 2.54
            else:
                return value_temp * 1.0

    except:
        return 0.0


if __name__ == '__main__':
    min_time, lr, dir, exclude_diagnoses, tx_interval, use_dl, tx_start, use_dynamic, use_dx = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9])
    use_imputation = int(sys.argv[10])
    include_dynamic = int(sys.argv[11])
    io_only = int(sys.argv[12])
    exclude_mutations = int(sys.argv[13])
    starting_year = int(sys.argv[14])
    ablation = int(sys.argv[15])
    patientID_to_censor_date = dict()

    if dir == 0:
        dir_path = dir_path1
    else:
        dir_path = dir_path2
    print(len(headers_to_drop))
    diagnosis_dates = pd.read_csv(dir_path+'Enhanced_AdvancedNSCLC.csv')
    patientID_to_advanced_diagnosis_date = dict()

    for i in range(len(diagnosis_dates['PatientID'])):
        patient_ID, diag_date  = diagnosis_dates['PatientID'][i], diagnosis_dates['AdvancedDiagnosisDate'][i]

        try:
            final_date = date.fromisoformat(diag_date)
            patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, final_date)
            patientID_to_advanced_diagnosis_date[patient_ID] = final_date
        except:
            continue

    sdh = pd.read_csv(dir_path+'SocialDeterminantsOfHealth.csv')
    patientID_to_sdh = dict()

    for i in range(len(sdh['PatientID'])):
        patient_ID, sdh_individual  = sdh['PatientID'][i], sdh['SESIndex2015_2019'][i]

        try:
            if "Lowest" in sdh_individual:
                patientID_to_sdh[patient_ID] = 1
            elif "Highest" in sdh_individual:
                patientID_to_sdh[patient_ID] = 5
            else:
                patientID_to_sdh[patient_ID] = int(sdh_individual)
        except:
            patientID_to_sdh[patient_ID] = 0

    practiceIDs = pd.read_csv(dir_path+'Practice.csv')
    practiceID_to_number = create_set(practiceIDs['PracticeID'])
    physician_ID_to_number = create_set(practiceIDs['PrimaryPhysicianID'])
    _nan = float("nan")

    race_to_number = {'': 0, _nan: 0, 'nan': 0, float("nan"): 0, 'White': 1,  'Asian': 2, 'Other Race': 3, 'Hispanic or Latino': 4,
                      'Black or African American': 5}
    histologies = {'Squamous cell carcinoma': 1, 'Non-squamous cell carcinoma': 2, 'NSCLC histology NOS': 0}
    gender_to_number = {'': 0, float("nan"): 0, 'F': 1, 'M': 2}
    stages = {'Group stage is not reported': 0, 'Stage 0': 1, 'Stage I': 2, 'Stage IA': 3, 'Stage IA1': 4,
              'Stage IA2': 5, 'Stage IA3': 6, 'Stage IB': 7, 'Stage II': 8, 'Stage IIA': 9, 'Stage IIB': 10,
              'Stage III': 11, 'Stage IIIA': 12, 'Stage IIIB': 13, 'Stage IIIC': 14, 'Stage IV': 15, 'Stage IVA': 16,
              'Stage IVB': 17,    'Occult': 18}
    smoking_statuses = {'Unknown/Not documented': 0, 'No history of smoking': 1, 'History of smoking': 2}
    practice_type_to_number = {'ACADEMIC': 2, 'COMMUNITY': 1}

    practiceIDs_data = dict()

    for i in range(len(practiceIDs['PatientID'])):

        patientID_current, practiceID_current, practice_type_current, physician_current = practiceIDs['PatientID'][i], \
                    practiceIDs['PracticeID'][i], practiceIDs['PracticeType'][i], practiceIDs['PrimaryPhysicianID'][i]
        practiceIDs_data.update({patientID_current:np.zeros(3,)})
        practiceIDs_data[patientID_current][0] = practiceID_to_number[practiceID_current]
        practiceIDs_data[patientID_current][1] = practice_type_to_number[practice_type_current]
        practiceIDs_data[patientID_current][2] = physician_ID_to_number[physician_current]

    insurance = pd.read_csv(dir_path + 'Insurance.csv')
    patientID_to_insurance = dict()

    for i in range(len(insurance['PatientID'])):
        patientID, insurance_type, start_date, end_date, practiceID_current = insurance['PatientID'][i], insurance['PayerCategory'][i], insurance['StartDate'][i], insurance['EndDate'][i], insurance['PracticeID'][i]
        diagnosis_date = patientID_to_advanced_diagnosis_date[patientID]
        practiceID_true = int(practiceIDs_data[patientID][0])
        practiceID_current_set = practiceID_to_number[practiceID_current]
        if practiceID_current_set != practiceID_true:
            continue

        start_date_final = date.today()
        try:
            start_date_final = date.fromisoformat(start_date)
            patientID_to_censor_date = update_last_note(patientID, patientID_to_censor_date, start_date_final)
        except:
            try:
                end_date_final = date.fromisoformat(end_date)
                patientID_to_censor_date = update_last_note(patientID, patientID_to_censor_date, end_date_final)
            except:
                continue

        try:
            end_date_final = date.fromisoformat(end_date)
            patientID_to_censor_date = update_last_note(patientID, patientID_to_censor_date, end_date_final)
        except:
            end_date_final = date.today()

        if 'Workers Compensation' == insurance_type:
            insuranceID = 1
        elif 'Self Pay' == insurance_type:
            insuranceID = 2
        elif 'Patient Assistance Program' == insurance_type:
            insuranceID = 3
        elif 'Other Government Program' == insurance_type:
            insuranceID = 4
        elif 'Medicare' == insurance_type:
            insuranceID = 5
        elif 'Medicaid' == insurance_type:
            insuranceID = 6
        elif 'Commercial Health Plan' == insurance_type:
            insuranceID = 7
        else:
            insuranceID = 0

        all_insurances = [0, 0, 0, 0, 0, 0, 0, 0]

        if diagnosis_date >= start_date_final and diagnosis_date <= end_date_final:
            all_insurances[insuranceID] = 1

        if patientID not in patientID_to_insurance:
            patientID_to_insurance[patientID] = all_insurances
        else:
            patientID_to_insurance[patientID][insuranceID] = 1

    patientIDs_demos = dict()
    patientIDs_cancer_types = dict()

    cancer_types = pd.read_csv(dir_path + 'Enhanced_AdvancedNSCLC.csv')

    for i in range(len(cancer_types['PatientID'])):
        [patientID, histology, stage, smoking_status] = [cancer_types['PatientID'][i], cancer_types['Histology'][i],
                                                 cancer_types['GroupStage'][i], cancer_types['SmokingStatus'][i]]
        patientIDs_cancer_types[patientID] = np.array([stages[stage], histologies[histology], smoking_statuses[smoking_status]])

    demographics = pd.read_csv(dir_path+'Demographics.csv')
    patients_IDs_to_number = create_set(demographics['PatientID'])
    state_to_number = {float("nan"): 0, 'WI': 1, 'MN': 2, 'IN': 3, 'VA': 4, 'PR': 5, 'DC': 6, 'UT': 7, 'ID': 8, 'MO': 9,
                       'CT': 10, 'NH': 11, 'CA': 12, 'AR': 13, 'NV': 14, 'DE': 15, 'MD': 16, 'TN': 17, 'AL': 18,
                       'NJ': 19, 'PA': 20, 'NY': 21, 'NE': 22, 'WA': 23, 'WV': 24, 'AZ': 25, 'LA': 26, 'OR': 27,
                       'OK': 28, 'TX': 29, 'CO': 30, 'IA': 31, 'MS': 32, 'RI': 33, 'OH': 34, 'SC': 35, 'GA': 36,
                       'MI': 37, 'NC': 38, 'ME': 39, 'FL': 40, 'IL': 41, 'NM': 42, 'HI': 43, 'KS': 44, 'KY': 45, 'MA': 46}

    physicianID_nan = -1

    for key, val in physician_ID_to_number.items():
        if isinstance(key, float) and np.isnan(key):
            physicianID_nan = val
            break

    for i in range(len(demographics['PatientID'])):
        [patientID, birth_year, gender, race, state, ethnicity] = [demographics['PatientID'][i], demographics['BirthYear'][i],
            demographics['Gender'][i], demographics['Race'][i], demographics['State'][i], demographics['Ethnicity'][i]]
        try:
            gender_ = gender_to_number[gender]
        except:
            gender_ = 0

        try:
            race_ = race_to_number[race]
        except:
            race_ = 0

        try:
            state_ = state_to_number[state]
        except:
            state_ = 0
        try:
            if ethnicity == 'Hispanic or Latino':
                ethnicity_ = 1
            else:
                ethnicity_ = 0
        except:
            ethnicity_ = 0

        patientIDs_demos[patientID] = np.array([birth_year, gender_, race_, ethnicity_, state_])

    ECOGs = pd.read_csv(dir_path + 'BaselineECOG.csv')
    patientID_to_ecog = dict()

    # Remember that you added 1 to ECOG values so that missing values could be 0
    for i in range(len(ECOGs['PatientID'])):
        patient_ID, ecog = ECOGs['PatientID'][i], ECOGs["ECOGValue"][i]
        try:
            patientID_to_ecog[patient_ID] = int(ecog) + 1
        except:
            patientID_to_ecog[patient_ID] = 0

    # Line of therapy (first-line mono IO vs. chemo vs combo)
    line_of_therapy = pd.read_csv(dir_path +'/LineOfTherapy.csv')
    patientID_to_therapyline = dict()

    approved_first_line_immunomonotherapies = ["nivolumab", "pembrolizumab", "cemiplimab", "atezolizumab", "durvalumab", "ipilimumab,nivolumab", "nivolumab,ipilimumab"]

    for i in range(len(line_of_therapy['PatientID'])):

        patient_ID, line_number, line_name, end_date, start_date, is_maintenance = line_of_therapy['PatientID'][i], line_of_therapy["LineNumber"][i], \
                                  line_of_therapy["LineName"][i], line_of_therapy["EndDate"][i], line_of_therapy["StartDate"][i], line_of_therapy["IsMaintenanceTherapy"][i]
        is_maintenance = bool(is_maintenance)
        patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, end_date)

        if is_maintenance:
            del patient_ID
            continue

        try:
            start_date = convert_date_to_iso(start_date)
            end_date = convert_date_to_iso(end_date)
        except:
            start_date = date.fromisoformat(start_date)
            end_date = date.fromisoformat(end_date)

        if patient_ID not in patientID_to_therapyline:
            patientID_to_therapyline[patient_ID] = [(line_number, line_name.lower(), end_date, start_date)]
        else:
            patientID_to_therapyline[patient_ID].append((line_number, line_name.lower(), end_date, start_date))

    biomarkers = pd.read_csv(dir_path +'Enhanced_AdvNSCLCBiomarkers.csv')
    patientID_to_biomarkers = dict()
    staining = {"nan": -1.0, float('nan'): -1.0, "": -1.0, "0%": 0.0, "< 1%": 0.0, "1%": 0.01, "2% - 4%": 0.02,
                "5% - 9%": 0.05,"10% - 19%": 0.10, "20% - 29%": 0.2, "30% - 39%": 0.3,"40% - 49%": 0.4,
                "50% - 59%": 0.5, "60% - 69%": 0.6, "70% - 79%": 0.7, "80% - 89%": 0.8, "90% - 99%": 0.9, "100%": 1.0}

    for i in range(len(biomarkers['PatientID'])):
        patient_ID, biomarker_name, cell_type, biomarker_status, expression_level, sample_type, staining_intensity=\
              biomarkers['PatientID'][i], biomarkers["BiomarkerName"][i], biomarkers["CellType"][i], \
                  biomarkers["BiomarkerStatus"][i], biomarkers["ExpressionLevel"][i], biomarkers["SampleType"][i], \
                  biomarkers["PercentStaining"][i]

        if patient_ID not in patientID_to_biomarkers:
            patientID_to_biomarkers[patient_ID] = {"ALK": 0, "EGFR": 0, "KRAS": 0, "ROS1": 0, "BRAF": 0, "PDL1": -1.0,
                                                   "PDL1_given": 0, "HER2/ERBB2": 0, "MET": 0, "RET": 0, "NTRK1": 0,
                                                   "NTRK2": 0, "NTRK3": 0, "NTRK_other": 0}

        if biomarker_name == "PDL1":
            test_type = biomarkers["TestType"][i]
            if patientID_to_biomarkers[patient_ID]["PDL1"] != -1.0:
                continue
            if cell_type != "Tumor cells":
                continue
            elif test_type != "IHC":
                continue
            elif biomarker_status == "Unsuccessful/indeterminate test":
                continue
            else:
                try:
                    perc = staining[staining_intensity]
                except:
                    perc = -1.0
                if perc == -1.0:
                    if expression_level == "PD-L1 low expression":
                        perc = 0.01
                    elif expression_level == "PD-L1 high expression":
                        perc = 0.5
                    elif "PD-L1 negative/not detected" == biomarker_status:
                        perc = 0.0
                    else:
                        continue
                if perc >= 0.0:
                    patientID_to_biomarkers[patient_ID]["PDL1"] = perc
                else:
                    continue
        else:
            if ("rearrangement present" in biomarker_status or "positive" in biomarker_status or
                    "Rearrangement present" in biomarker_status):

                if "NTRK1" in biomarker_name:
                    biomarker_name = "NTRK1"
                elif "NTRK2" in biomarker_name:
                    biomarker_name = "NTRK2"
                elif "NTRK3" in biomarker_name:
                    biomarker_name = "NTRK3"
                elif "NTRK" in biomarker_name:
                    biomarker_name = "NTRK_other"

                patientID_to_biomarkers[patient_ID][biomarker_name] = 1
            else:
                continue

    for patientID in patientID_to_biomarkers:
        pdl1 = patientID_to_biomarkers[patientID]["PDL1"]
        if pdl1 >= 0.0:
            patientID_to_biomarkers[patientID]["PDL1_given"] = 1
        else:
            patientID_to_biomarkers[patientID]["PDL1"] = 0.0

    diagnosis = pd.read_csv(dir_path + 'Diagnosis.csv')
    diagnosis_values = list(map(lambda x: x.lower(), diagnosis["DiagnosisDescription"].values))
    all_diagnoses_to_numbers = create_set(diagnosis_values, start_at_zero=True)
    patientIDs_diagnoses = dict()
    patientID_to_contraindications = dict()
    patientID_to_important_diagnoses = dict()
    dynamic_variables = dict()

    for i in range(len(diagnosis['PatientID'])):

        patientID_current = diagnosis['PatientID'][i]
        diagnosis_name = diagnosis['DiagnosisDescription'][i].lower()
        diagnosis_date = diagnosis['DiagnosisDate'][i]

        if isinstance(diagnosis_date, float) or '1800' in diagnosis_date or '1899' in diagnosis_date:
            continue

        try:
            dx_date_temp = time.strptime(diagnosis_date, "%Y-%m-%d")
        except:
            dx_date_temp = time.strptime(diagnosis_date, "%m/%d/%y")

        dx_date = date(dx_date_temp.tm_year, dx_date_temp.tm_mon, dx_date_temp.tm_mday)
        adv_dx_date = patientID_to_advanced_diagnosis_date[patientID_current]
        patientID_to_censor_date = update_last_note(patientID_current, patientID_to_censor_date, dx_date)

        if adv_dx_date < dx_date:
            continue

        index_diagnosis = all_diagnoses_to_numbers[diagnosis_name]
        possible_contraindications_kidney = ["kidney failure", "chronic kidney disease",
                                             "renal disease", "kidney transplant"]
        possible_contraindications_liver = ["cirrhosis", "hepatitis", "liver transplant"]
        possible_contraindications_ctd = ["connective tissue", "scleroderma", "lupus", "rheumatoid arthritis",
                                          "interstitial lung disease"]

        comorbidities = [list(), list(), list(), list(), list(), list()]

        if patientID_current in patientIDs_diagnoses:
            patientIDs_diagnoses[patientID_current][index_diagnosis] = 1
        else:
            patientIDs_diagnoses[patientID_current] = np.zeros(len(all_diagnoses_to_numbers))
            patientIDs_diagnoses[patientID_current][index_diagnosis] = 1

        if patientID_current not in patientID_to_important_diagnoses:
            patientID_to_important_diagnoses[patientID_current] = [0] * 19

        for i in range(len(possible_contraindications_kidney)):
            kidney_contra = possible_contraindications_kidney[i]
            if kidney_contra in diagnosis_name or diagnosis_name == kidney_contra:
                patientID_to_important_diagnoses[patientID_current][i] = 1

        for i in range(len(possible_contraindications_liver)):
            liver_contra = possible_contraindications_liver[i]
            if liver_contra in diagnosis_name or diagnosis_name == liver_contra:
                patientID_to_important_diagnoses[patientID_current][4+i] = 1

        for i in range(len(possible_contraindications_ctd)):
            ctd_contra = possible_contraindications_ctd[i]
            if ctd_contra in diagnosis_name or diagnosis_name == ctd_contra:
                patientID_to_important_diagnoses[patientID_current][7+i] = 1

        if "diabetes" in diagnosis_name:
            patientID_to_important_diagnoses[patientID_current][-7] = 1

        if "secondary malignant neoplasm" in diagnosis_name:
            name_len = len("secondary malignant neoplasm of ")
            organ = diagnosis_name[name_len:]
            if "bone" in organ or "bone marrow" in organ:
                patientID_to_important_diagnoses[patientID_current][-6] = 1
            elif "brain" in organ:
                patientID_to_important_diagnoses[patientID_current][-5] = 1
            elif "nervous system" in organ or "spinal cord" in organ:
                patientID_to_important_diagnoses[patientID_current][-4] = 1
            elif "retroperitoneum" in organ or "liver" in organ or "digestive" in organ:
                patientID_to_important_diagnoses[patientID_current][-3] = 1
            elif "adrenal" in organ:
                patientID_to_important_diagnoses[patientID_current][-2] = 1
            else:
                patientID_to_important_diagnoses[patientID_current][-1] = 1

    y = []
    static_variables = dict()
    total_meds_vitals_labs = 0
    len_static = 0
    total_len_static = 0
    chemotherapy_drugs = ["cisplatin", "carboplatin"]
    first_line_mono_io_to_index = {"nivolumab": 1, "pembrolizumab": 2, "atezolizumab": 4,
                                   "durvalumab": 5, "ipilimumab,nivolumab": 6, "nivolumab,ipilimumab": 6, "cemiplimab": 3,
                                   "cemiplimab-rwlc": 3}
    first_line_chemotherapy_drugs = ["cisplatin", "carboplatin"]
    egfr_drugs = ["osimertinib", "erlotinib", "gefitinib", "afatinib"]
    alk_drugs = ["alectinib", "ceritinib", "brigatinib"]
    ros1_drugs = ["crizotinib", "entrectinib", "repotrectinib"]
    braf_drugs = ["dabrafenib,trametinib", "binimetinib,encorafenib", "dabrafenib", "trametinib", "binimetinib", "encorafenib"]
    ras_drugs = ["sotorasib", "adagrasib"]
    secondary_chemo_drugs = ["vinorelbine", "docetaxel", "gemcitabine", "pemetrexed", "paclitaxel",
                             "paclitaxel protein-bound"]
    trk_inhibitors = ["imatinib", "dasatinib", "nilotinib", "bosutinib"]
    met_inhibitors = ["capmatinib", "tepotinib"]
    patientID_to_first_line_start_date = dict()

    for patientID in patientID_to_advanced_diagnosis_date:
        x_demos = patientIDs_demos[patientID]

        age_diagnosis_stats = [0, 0]

        if patientID in patientID_to_advanced_diagnosis_date:
            diagnosis_year = patientID_to_advanced_diagnosis_date[patientID].year
            age_at_diagnosis = diagnosis_year - x_demos[0]
            age_diagnosis_stats = [diagnosis_year, age_at_diagnosis]

        insurance_patient = [0, 0, 0, 0, 0, 0, 0, 0]
        if patientID in patientID_to_insurance:
            insurance_patient = patientID_to_insurance[patientID]

        x_practice = [0, 0, 0]
        if patientID in practiceIDs_data:
            x_practice = practiceIDs_data[patientID]

        cancer_vec = [0, 0, 0]
        if patientID in patientIDs_cancer_types:
            cancer_vec = patientIDs_cancer_types[patientID]

        ecog_ = [0]
        if patientID in patientID_to_ecog:
            ecog_ = [patientID_to_ecog[patientID]]

        all_biomarkers = []
        if patientID in patientID_to_biomarkers:

            for biomarker in ["ALK", "EGFR", "KRAS", "ROS1", "BRAF", "PDL1", "PDL1_given", "HER2/ERBB2", "MET", "RET", "NTRK1",
                                                   "NTRK2", "NTRK3", "NTRK_other"]:
                all_biomarkers.append(patientID_to_biomarkers[patientID][biomarker])
        else:
            all_biomarkers = [0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0]

        if exclude_mutations == 1:
            if all_biomarkers[0] == 1 or all_biomarkers[1] == 1:
                continue

        all_biomarkers = all_biomarkers[2:]

        temp_combo_list = []
        combo_therapy = 0
        first_line_chemo = 0
        io_mono = 0
        io_mono_used = 0
        egfr_drug = 0
        alk_drug = 0
        ros1_drug = 0
        braf_drug = 0
        ras_drug = 0
        secondary_chemo_drug = 0
        other_first_line_therapy = 0
        days_from_dx_to_tx = 0
        clinical_study_drug = 0
        bevacizumab_used = 0
        three_plus_chemo_drugs = 0
        carboplatin_only = 0
        cisplatin_only = 0
        pembrolizumab_used = 0
        trk_inhibitor = 0
        met_drug = 0

        if patientID in patientID_to_therapyline:
            diagnosis_date = patientID_to_advanced_diagnosis_date[patientID]
            all_therapies_used = patientID_to_therapyline[patientID]

            # Loop through all therapies this patient was prescribed
            for i in range(len(all_therapies_used)):
                (therapy_line, therapy_name, end_date, start_date) = all_therapies_used[i]

                # Only examine first-line therapies that had at least one dose given after advanced diagnosis date
                if therapy_line == 1:
                    some_other_immuno = 0
                    patientID_to_first_line_start_date[patientID] = start_date

                    if start_date < diagnosis_date:
                        continue

                    else:
                        days_from_dx_to_tx = (start_date-diagnosis_date).days
                    if "clinical study drug" in therapy_name:
                        clinical_study_drug = 1

                    if "pembrolizumab" in therapy_name:
                        pembrolizumab_used = 1

                    if "bevacizumab" in therapy_name or "bevacizumab-awwb" in therapy_name:
                        bevacizumab_used = 1

                    if therapy_name in approved_first_line_immunomonotherapies:
                        io_mono = 1
                        io_mono_used = first_line_mono_io_to_index[therapy_name]

                    elif therapy_name == "carboplatin":
                        carboplatin_only = 1
                    elif therapy_name == "cisplatin":
                        cisplatin_only = 1
                    elif therapy_name in secondary_chemo_drugs:
                        secondary_chemo_drug = 1
                    elif therapy_name in egfr_drugs:
                        egfr_drug = 1
                    elif therapy_name in alk_drugs:
                        alk_drug = 1
                    elif therapy_name in ros1_drugs:
                        ros1_drug = 1
                    elif therapy_name in braf_drugs:
                        braf_drug = 1
                    elif therapy_name in ras_drugs:
                        ras_drug = 1
                    elif therapy_name in trk_inhibitors:
                        trk_inhibitor = 1
                    elif therapy_name in met_inhibitors:
                        met_drug = 1
                    else:
                        temp_combo = [0, 0, 0, 0]
                        all_meds_at_instance = therapy_name.split(',')
                        for j in range(len(all_meds_at_instance)):
                            current_med = all_meds_at_instance[j]

                            if current_med in egfr_drugs:
                                egfr_drug = 1
                            elif current_med in alk_drugs:
                                alk_drug = 1
                            elif current_med in ros1_drugs:
                                ros1_drug = 1
                            elif current_med in braf_drugs:
                                braf_drug = 1
                            elif current_med in ras_drugs:
                                ras_drug = 1
                            elif current_med in trk_inhibitors:
                                trk_inhibitor = 1
                            elif current_med in met_inhibitors:
                                met_drug = 1
                            elif current_med in approved_first_line_immunomonotherapies:
                                temp_combo[0] = 1
                            elif current_med in first_line_chemotherapy_drugs:
                                temp_combo[1] = 1
                            elif current_med in secondary_chemo_drugs:
                                if temp_combo[1] == 1 and temp_combo[2] == 1:
                                    three_plus_chemo_drugs = 1
                                temp_combo[2] = 1
                            else:
                                temp_combo[3] = 1

                        if clinical_study_drug == 1:
                            clinical_study_drug = 1
                            other_first_line_therapy = 1
                        # Chemotherapy only: No IO monotherapy but cisplatin/carboplatin (first-line)
                        elif temp_combo[0] == 0 and temp_combo[1] == 1  and temp_combo[2] == 1:
                            first_line_chemo = 1

                        # Combination therapy: IO monotherapy plus cisplatin/carboplatin (first-line) +/- other drug(s)
                        elif temp_combo[0] == 1 and temp_combo[1] == 1:
                            combo_therapy = 1

                        # Other therapy: No IO monotherapy and no first-line cisplatin/carboplatin
                        else:
                            other_first_line_therapy = 1

        if patientID not in patientID_to_first_line_start_date and use_dx == 0:
            continue

        if days_from_dx_to_tx > tx_interval and use_dx == 0:
            continue

        if patientID in patientID_to_important_diagnoses:
            all_important_dx = patientID_to_important_diagnoses[patientID]
        else:
            all_important_dx = [0] * 19

        if patientID in patientID_to_sdh:
            sdh_final = patientID_to_sdh[patientID]
        else:
            sdh_final = 0

        if use_dx == 0:
            therapy_year = patientID_to_first_line_start_date[patientID].year
            therapy_info = [io_mono_used, combo_therapy, first_line_chemo, secondary_chemo_drug, alk_drug, egfr_drug,
                        braf_drug, ros1_drug, ras_drug, other_first_line_therapy, clinical_study_drug, bevacizumab_used,
                        three_plus_chemo_drugs, carboplatin_only, cisplatin_only, pembrolizumab_used, trk_inhibitor,
                        met_drug, days_from_dx_to_tx, therapy_year]

            x_demos_no_diagnoses = np.concatenate(([x_practice[2], x_practice[0]], age_diagnosis_stats, x_demos,
                                               insurance_patient, [x_practice[1]], ecog_, cancer_vec, all_biomarkers,
                                               therapy_info, all_important_dx, [sdh_final]))
        else:
            x_demos_no_diagnoses = np.concatenate(([x_practice[2], x_practice[0]], age_diagnosis_stats, x_demos,
                                               insurance_patient, [x_practice[1]], ecog_, cancer_vec, all_biomarkers,
                                               all_important_dx, [sdh_final]))

        static_variables[patientID] = x_demos_no_diagnoses

    patientID_to_steroids = dict()
    patientID_to_abx = dict()
    patientID_to_albumin = dict()
    patientID_to_bmi = dict()
    patientID_to_creatinine = dict()
    patientID_to_bilirubin = dict()
    patientID_to_ALT = dict()
    patientID_to_AST = dict()

    if use_dynamic:
        patientID_to_labs = dict()
        labs = pd.read_csv(dir_path + 'Lab.csv')

        for i in range(len(labs['PatientID'])):
            patient_ID, date_, test_name, test_result = labs['PatientID'][i], labs['ResultDate'][i], labs['LabComponent'][i], labs['TestResult'][i]
            easy_lab_name = labs['TestBaseName'][i]

            if not isinstance(date_, str) or test_result == 'nan':
                continue

            try:

                try:
                    labs_date_temp = time.strptime(date_, "%Y-%m-%d")
                except:
                    labs_date_temp = time.strptime(date_, "%m/%d/%y")

                labs_date = date(labs_date_temp.tm_year, labs_date_temp.tm_mon, labs_date_temp.tm_mday)
                adv_dx_date = patientID_to_advanced_diagnosis_date[patient_ID]
                patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, labs_date)

                if use_dx == 0:
                    if patient_ID not in patientID_to_first_line_start_date:
                        continue
                    else:
                        tx_init_date = patientID_to_first_line_start_date[patient_ID]
                        starting_date = tx_init_date
                else:
                    starting_date = adv_dx_date

                if starting_date < labs_date:
                    continue

                test_name = test_name.lower()

                if patient_ID not in patientID_to_labs:
                    patientID_to_labs[patient_ID] = dict()
                    patientID_to_labs[patient_ID][labs_date] = {test_name:test_result}
                elif labs_date not in patientID_to_labs[patient_ID]:
                    patientID_to_labs[patient_ID][labs_date] = {test_name:test_result}
                else:
                    patientID_to_labs[patient_ID][labs_date][test_name] = test_result

                if "Test not performed" in test_result or "Negative" in test_result or "Pending" in test_result or "Note_Comment" in test_result or test_result=="":
                    continue
                elif '<' in test_result or '>' in test_result:
                    test_result = test_result[1:]
                elif '%' in test_result or 'L' in test_result or 'H' in test_result or '+' in test_result:
                    test_result = test_result[:-1]
                test_result = float(test_result)

                if test_result <= 0.0:
                    continue

                if test_name == "albumin, serum" and labs_date <= starting_date:
                    if test_result > 10:
                        continue

                    if patient_ID not in patientID_to_albumin:
                        patientID_to_albumin[patient_ID] = (labs_date, test_result)
                    else:
                        (orig_date, orig_test_result) = patientID_to_albumin[patient_ID]
                        if orig_date < labs_date:
                            patientID_to_albumin[patient_ID] = (labs_date, test_result)

                if test_name == "creatinine, serum" and labs_date <= starting_date:

                    if patient_ID not in patientID_to_creatinine:
                        patientID_to_creatinine[patient_ID] = (labs_date, test_result)
                    else:
                        (orig_date, orig_test_result) = patientID_to_creatinine[patient_ID]
                        if orig_date < labs_date:
                            patientID_to_creatinine[patient_ID] = (labs_date, test_result)

                if 'bilirubin' in test_name and 'total' in test_name and labs_date <= starting_date:
                    if test_result > 50:
                        continue

                    if patient_ID not in patientID_to_bilirubin:
                        patientID_to_bilirubin[patient_ID] = (labs_date, test_result)
                    else:
                        (orig_date, orig_test_result) = patientID_to_bilirubin[patient_ID]
                        if orig_date < labs_date:
                            patientID_to_bilirubin[patient_ID] = (labs_date, test_result)

                if 'alanine aminotransferase' in test_name  and labs_date <= starting_date:
                    if patient_ID not in patientID_to_ALT:
                        patientID_to_ALT[patient_ID] = (labs_date, test_result)
                    else:
                        (orig_date, orig_test_result) = patientID_to_ALT[patient_ID]
                        if orig_date < labs_date:
                            patientID_to_ALT[patient_ID] = (labs_date, test_result)

                if 'aspartate aminotransferase' in test_name and labs_date <= starting_date:
                    if patient_ID not in patientID_to_AST:
                        patientID_to_AST[patient_ID] = (labs_date, test_result)
                    else:
                        (orig_date, orig_test_result) = patientID_to_AST[patient_ID]
                        if orig_date < labs_date:
                            patientID_to_AST[patient_ID] = (labs_date, test_result)
            except:
                continue

        patientID_to_med_administration = dict()
        med_administration = pd.read_csv(dir_path + 'MedicationAdministration.csv')

        for i in range(len(med_administration['PatientID'])):
            patient_ID, date_, drug_name, drug_dose, dose_units = med_administration['PatientID'][i], med_administration['AdministeredDate'][i], med_administration['CommonDrugName'][i], med_administration['AdministeredAmount'][i], med_administration['AdministeredUnits'][i]
            drug_type = med_administration['DetailedDrugCategory'][i]

            if not isinstance(date_, str):
                continue

            try:
                meds_date_temp = time.strptime(date_, "%Y-%m-%d")
            except:
                meds_date_temp = time.strptime(date_, "%m/%d/%y")

            drug_name = drug_name.lower()
            meds_date = date(meds_date_temp.tm_year, meds_date_temp.tm_mon, meds_date_temp.tm_mday)
            adv_dx_date = patientID_to_advanced_diagnosis_date[patient_ID]
            patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, meds_date)

            if use_dx == 0:
                if patient_ID not in patientID_to_first_line_start_date:
                    continue
                else:
                    tx_init_date = patientID_to_first_line_start_date[patient_ID]
                    starting_date = tx_init_date
            else:
                starting_date = adv_dx_date

            if starting_date < meds_date:
                continue

            if patient_ID not in patientID_to_med_administration:
                patientID_to_med_administration[patient_ID] = dict()
                patientID_to_med_administration[patient_ID][meds_date] = {drug_name: (drug_dose, dose_units)}
            elif meds_date not in patientID_to_med_administration[patient_ID]:
                patientID_to_med_administration[patient_ID][meds_date] = {drug_name: (drug_dose, dose_units)}
            else:
                patientID_to_med_administration[patient_ID][meds_date][drug_name] = (drug_dose, dose_units)

            if drug_type == "glucocorticoid" and 0 <= (starting_date - meds_date).days <= 28:
                    patientID_to_steroids[patient_ID] = 1
            if drug_type == "anti-infective" and 0 <= (starting_date - meds_date).days <= 28:
                patientID_to_abx[patient_ID] = 1

        vitals = pd.read_csv(dir_path + 'Vitals.csv')
        patientID_to_vitals = dict()

        for i in range(len(vitals['PatientID'])):

            patient_ID, date_, test_name, test_result, test_units = vitals['PatientID'][i], vitals['TestDate'][i], \
                                                       vitals['LabComponent'][i], vitals['TestResult'][i], vitals['TestUnits'][i]

            if not isinstance(date_, str):
                continue

            test_name = test_name.lower()

            try:
                vitals_date_temp = time.strptime(date_, "%Y-%m-%d")
            except:
                vitals_date_temp = time.strptime(date_, "%m/%d/%y")

            vitals_date = date(vitals_date_temp.tm_year, vitals_date_temp.tm_mon, vitals_date_temp.tm_mday)
            adv_dx_date = patientID_to_advanced_diagnosis_date[patient_ID]
            patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, vitals_date)

            if use_dx == 0:
                if patient_ID not in patientID_to_first_line_start_date:
                    continue
                else:
                    tx_init_date = patientID_to_first_line_start_date[patient_ID]
                    starting_date = tx_init_date
            else:
                starting_date = adv_dx_date

            if starting_date < vitals_date:
                continue

            if patient_ID not in patientID_to_vitals:
                patientID_to_vitals[patient_ID] = dict()
                patientID_to_vitals[patient_ID][vitals_date] = {test_name: (test_result, test_units)}
            elif vitals_date not in patientID_to_vitals[patient_ID]:
                patientID_to_vitals[patient_ID][vitals_date] = {test_name: (test_result, test_units)}
            else:
                patientID_to_vitals[patient_ID][vitals_date][test_name] = (test_result, test_units)

        med_admin_set, vitals_set, labs_set = [], [], []
        all_dates_for_all_meds_vitals_labs = dict()

        for patientID in patientID_to_med_administration:
            all_meds_all_dates = patientID_to_med_administration[patientID]
            for med_date, all_meds_single_date in sorted(all_meds_all_dates.items()):

                if patientID not in all_dates_for_all_meds_vitals_labs:
                    all_dates_for_all_meds_vitals_labs[patientID] = set()
                all_dates_for_all_meds_vitals_labs[patientID].add(med_date)

                for drug in all_meds_single_date:
                    med_admin_set.append(drug)

        for patientID in patientID_to_vitals:
            all_vitals_all_dates = patientID_to_vitals[patientID]
            for vitals_date, all_vitals_single_date in sorted(all_vitals_all_dates.items()):

                if patientID not in all_dates_for_all_meds_vitals_labs:
                    all_dates_for_all_meds_vitals_labs[patientID] = set()
                all_dates_for_all_meds_vitals_labs[patientID].add(vitals_date)

                for vital in all_vitals_single_date:
                    vitals_set.append(vital)

        for patientID in patientID_to_labs:
            all_labs_all_dates = patientID_to_labs[patientID]
            for labs_date, all_labs_single_date in sorted(all_labs_all_dates.items()):

                if patientID not in all_dates_for_all_meds_vitals_labs:
                    all_dates_for_all_meds_vitals_labs[patientID] = set()
                all_dates_for_all_meds_vitals_labs[patientID].add(labs_date)

                for lab in all_labs_single_date:
                    labs_set.append(lab)

        med_admin_set_final, vitals_set_final, labs_set_final = create_set(med_admin_set, start_at_zero=True), create_set(vitals_set, start_at_zero=True), create_set(labs_set, start_at_zero=True)

        total_num_meds = len(med_admin_set_final)
        total_num_vitals = len(vitals_set_final)
        total_num_labs = len(labs_set_final)
        total_meds_vitals_labs = total_num_labs + total_num_meds + total_num_vitals
        total_vitals_labs = total_num_vitals + total_num_labs

        print("Number of unique meds: ", total_num_meds)
        print("Number of unique vitals: ", total_num_vitals)
        print("Number of unique labs: ", total_num_labs)
        print("Dynamic length: ", total_meds_vitals_labs)
        patientID_to_height_weight = dict()

        # Go through each dict of patientID --> all dates for which a vital, med, or lab was recorded to add info to dynamic variables for that respective patient
        for patientID, dates_set in all_dates_for_all_meds_vitals_labs.items():
            if patientID in patientID_to_med_administration:
                all_meds_for_patient = patientID_to_med_administration[patientID]
            else:
                all_meds_for_patient = []

            if patientID in patientID_to_vitals:
                all_vitals_for_patient = patientID_to_vitals[patientID]
            else:
                all_vitals_for_patient = []

            if patientID in patientID_to_labs:
                all_labs_for_patient = patientID_to_labs[patientID]
            else:
                all_labs_for_patient = []

            dynamic_holder_all_dates_single_patient = []
            current_count = 0   # counter for current timestep in patient's array

            for current_date in sorted(dates_set):
                dynamic_holder_single_date = np.zeros(total_meds_vitals_labs, dtype='float64')    # value holder

                if current_date in all_vitals_for_patient:                      # If vitals were taken on this date, add them
                    date_specific_vitals = all_vitals_for_patient[current_date]
                    for vital, value_units in sorted(date_specific_vitals.items()):     # Iterate through all vitals on date
                        value, units = value_units[0], value_units[1]
                        index = vitals_set_final[vital]                         # New index in dynamic array
                        if current_count > 0:
                            final_value = get_vitals_value(value, units, current_count, dynamic_holder_all_dates_single_patient[current_count-1][index])
                        else:
                            final_value = get_vitals_value(value, units, current_count, 0)

                        if "body height" in vital and final_value > 0:
                            if patientID not in patientID_to_height_weight:
                                patientID_to_height_weight[patientID] = ([], [])
                            patientID_to_height_weight[patientID][0].append(final_value)#[current_date, final_value])

                        if "body weight" in vital and final_value > 0:
                            if patientID not in patientID_to_height_weight:
                                patientID_to_height_weight[patientID] = ([], [])
                            patientID_to_height_weight[patientID][1].append(final_value)#[current_date, final_value])

                        dynamic_holder_single_date[index] = final_value

                if current_date in all_labs_for_patient:                        # If labs were taken on this date, add them
                    date_specific_labs = all_labs_for_patient[current_date]
                    for lab, lab_value in sorted(date_specific_labs.items()):   # Iterate through all labs on date
                        index2 = labs_set_final[lab] + total_num_vitals         # New index in dynamic array
                        if current_count > 0:
                            final_value = get_labs_value(lab_value, current_count, dynamic_holder_all_dates_single_patient[current_count-1][index2])
                        else:
                            final_value = get_labs_value(lab_value, current_count, 0)
                        dynamic_holder_single_date[index2] = final_value

                if current_date in all_meds_for_patient:                        # If meds were taken on this date, add them
                    date_specific_meds = all_meds_for_patient[current_date]
                    for med_name, (med_dose, med_units) in sorted(date_specific_meds.items()):  # Iterate through all meds on date
                        index3 = med_admin_set_final[med_name] + total_num_vitals + total_num_labs   # New index in dynamic array
                        if current_count > 0:
                            final_value = get_meds_value(med_dose, med_units, current_count, dynamic_holder_all_dates_single_patient[current_count-1][index3])
                        else:
                            final_value = get_meds_value(med_dose, med_units, current_count, 0)
                        dynamic_holder_single_date[index3] = final_value

                dynamic_holder_all_dates_single_patient.append(dynamic_holder_single_date)
                current_count += 1

            dynamic_holder_all_dates_single_patient_array = np.array(dynamic_holder_all_dates_single_patient)
            total_rows = dynamic_holder_all_dates_single_patient_array.shape[0]
            for row in range(1, total_rows, 1):
                for col in range(total_meds_vitals_labs):
                    if dynamic_holder_all_dates_single_patient_array[row][col] == 0.0 and dynamic_holder_all_dates_single_patient_array[row-1][col] != 0.0:
                        dynamic_holder_all_dates_single_patient_array[row][col] = dynamic_holder_all_dates_single_patient_array[row-1][col]

            for row in range(total_rows-2, -1, -1):
                for col in range(total_meds_vitals_labs):
                    if dynamic_holder_all_dates_single_patient_array[row][col] == 0.0 and dynamic_holder_all_dates_single_patient_array[row+1][col] != 0.0:
                        dynamic_holder_all_dates_single_patient_array[row][col] = dynamic_holder_all_dates_single_patient_array[row+1][col]

            dynamic_variables[patientID] = dynamic_holder_all_dates_single_patient_array

    patientID_to_mixed_pseudo_date = dict()
    progressions = pd.read_csv(dir_path+'Enhanced_AdvNSCLC_Progression.csv')
    pd1s_progression = []
    pd1s_progression_dict = dict()
    patientID_to_progression = dict()
    location = 0

    for i in range(len(progressions['PatientID'])):

        patient_ID, rads, path, clinical, mixed, pseudo, progression_date = progressions['PatientID'][i], progressions['IsRadiographicEvidence'][i],\
                                progressions['IsPathologicEvidence'][i], progressions['IsClinicalAssessmentOnly'][i], \
                                progressions['IsMixedResponse'][i], progressions['IsPseudoprogressionMentioned'][i], progressions['ProgressionDate'][i]
        diagnosis_date = patientID_to_advanced_diagnosis_date[patient_ID]

        try:
            progression_date_final = date.fromisoformat(progression_date)
        except:
            if patient_ID not in patientID_to_progression:
                patientID_to_progression[patient_ID] = no_progression_date
            continue

        if rads == 'Yes':
            progression_positive_rads = 1
        else:
            progression_positive_rads = 0

        if path == 'Yes':
            progression_positive_path = 1
        else:
            progression_positive_path = 0

        if clinical == 'Yes':
            progression_positive_clinical = 1
        else:
            progression_positive_clinical = 0

        if mixed == 'Yes':
            progression_positive_mixed = 1
        else:
            progression_positive_mixed = 0

        if pseudo == 'Yes':
            progression_positive_pseudo = 1
        else:
            progression_positive_pseudo = 0

        if progression_positive_mixed == 1 or progression_positive_pseudo == 1:
            if patient_ID not in patientID_to_mixed_pseudo_date:
                patientID_to_mixed_pseudo_date[patient_ID] = progression_date_final
            else:
                curr_date = patientID_to_mixed_pseudo_date[patient_ID]
                if curr_date > progression_date_final:
                    patientID_to_mixed_pseudo_date[patient_ID] = progression_date_final

        patientID_to_censor_date = update_last_note(patient_ID, patientID_to_censor_date, progression_date_final)
        progression_bool = (progression_positive_rads or progression_positive_path or progression_positive_clinical)
        adv_dx_date = patientID_to_advanced_diagnosis_date[patient_ID]

        if use_dx == 0:
            if patient_ID not in patientID_to_first_line_start_date:
                continue
            else:
                tx_init_date = patientID_to_first_line_start_date[patient_ID]
                starting_date = tx_init_date
        else:
            starting_date = adv_dx_date

        if starting_date > progression_date_final:
            continue

        if patient_ID not in patientID_to_progression:
            if progression_bool:
                patientID_to_progression[patient_ID] = progression_date_final
            else:
                patientID_to_progression[patient_ID] = no_progression_date

        elif patient_ID in patientID_to_progression and progression_bool:
            orig = patientID_to_progression[patient_ID]
            if no_progression_date > orig > progression_date_final:
                patientID_to_progression[patient_ID] = progression_date_final
        else:
            continue

    X_dynamic, X_static, y = [], [], []
    mortality_list_y = []
    max_features = -1
    max_time_steps = -1

    mortality_dates = pd.read_csv(dir_path+'Enhanced_Mortality_V2.csv')
    mortality_dict = dict()
    for i in range(len(mortality_dates['PatientID'])):
        patientID, date_of_death = mortality_dates['PatientID'][i], mortality_dates['DateOfDeath'][i]
        try:
            mortality_date = date(int(date_of_death[0:4]), int(date_of_death[5:7]), 1)
            mortality_dict[patientID] = mortality_date
            patientID_to_censor_date = update_last_note(patientID, patientID_to_censor_date, mortality_date)
        except:
            mortality_dict[patientID] = patientID_to_censor_date[patientID]

    patientIDs_used = []

    num_pts_mixed_pseudo_to_full_prog = 0
    mixed_pseudo_full_prog_array = []

    for patientID, vals in static_variables.items():
        adv_dx_date = patientID_to_advanced_diagnosis_date[patientID]

        if patientID not in patientID_to_censor_date:
            continue

        last_final_recorded_date_records = patientID_to_censor_date[patientID]

        if last_final_recorded_date_records == no_progression_date or last_final_recorded_date_records > no_progression_date:
            print("Error")

        if use_dx == 0:
            tx_start_date = patientID_to_first_line_start_date[patientID]
            starting_date = tx_start_date
        else:
            starting_date = adv_dx_date

        time_to_censor = (last_final_recorded_date_records - starting_date).days

        if patientID in patientID_to_progression:
            progression_date = patientID_to_progression[patientID]

            if time_to_censor < 0:
                print(patientID)
                print(time_to_censor)
                print(last_final_recorded_date_records)
                print(starting_date)

            if progression_date == no_progression_date:
                progression = 0
            else:
                progression = (progression_date - starting_date).days

            if time_to_censor < min_time and progression_date == no_progression_date:
                continue

        else:
            continue

        if patientID in mortality_dict:
            mortality_days = (mortality_dict[patientID]-starting_date).days
        else:
            mortality_days = 0

        steroid_use, abx_use, albumin, creatinine, bilirubin, ast, alt = 0, 0, 0, 0, 0, 0, 0
        if patientID in patientID_to_steroids:
            steroid_use = 1
        if patientID in patientID_to_abx:
            abx_use = 1
        if patientID in patientID_to_albumin:
            albumin = patientID_to_albumin[patientID][1]
        if patientID in patientID_to_creatinine:
            creatinine = patientID_to_creatinine[patientID][1]
        if patientID in patientID_to_AST:
            ast = patientID_to_AST[patientID][1]
        if patientID in patientID_to_ALT:
            alt = patientID_to_ALT[patientID][1]
        if patientID in patientID_to_bilirubin:
            bilirubin = patientID_to_bilirubin[patientID][1]

        imp_predictors = [steroid_use, abx_use, albumin, creatinine, bilirubin, ast, alt]
        x_demos_no_diagnoses = np.concatenate((vals, imp_predictors))
        len_static = len(x_demos_no_diagnoses)
        location = x_demos_no_diagnoses.shape[0]

        x_diagnosis = np.zeros(len(all_diagnoses_to_numbers))
        if patientID in patientIDs_diagnoses:
            x_diagnosis = patientIDs_diagnoses[patientID]

        if exclude_diagnoses:
            x_static_final = x_demos_no_diagnoses
        else:
            x_static_final = np.concatenate((x_demos_no_diagnoses, x_diagnosis))
        total_len_static = x_static_final.shape[0]

        if include_dynamic:
            if patientID in dynamic_variables:
                dynamic_data = (dynamic_variables[patientID][-1]).flatten()
            else:
                dynamic_data = np.zeros(total_meds_vitals_labs)

            final_vals = np.concatenate((x_static_final, dynamic_data))
        else:
            final_vals = x_static_final

        if time_to_censor < 0:
            continue

        if patientID in patientID_to_progression and patientID in patientID_to_mixed_pseudo_date:
            mixed_pseudo_date = patientID_to_mixed_pseudo_date[patientID]
            if progression_date < no_progression_date and mixed_pseudo_date < progression_date:
                num_pts_mixed_pseudo_to_full_prog += 1
                interval = (progression_date-mixed_pseudo_date).days
                mixed_pseudo_full_prog_array.append([patientID, interval])

        prog_bool = int(min_time >= progression > 0)
        mort_bool = int(min_time >= mortality_days > 0)

        if time_to_censor < min_time and (prog_bool == 0):
            continue

        patientIDs_used.append(patientID)
        X_static.append(final_vals)
        y.append([int(min_time >= progression > 0), progression, mortality_days, int(min_time >= mortality_days > 0), time_to_censor])

    num_patients = len(X_static)
    print(num_patients)

    final_arr_df = pd.DataFrame(data=mixed_pseudo_full_prog_array)
    final_arr_df.to_csv('mixed_pseudo.csv')

    adv_dx_first_line_arr = []
    for patientID in patientID_to_first_line_start_date:
        start_date = patientID_to_first_line_start_date[patientID]
        if patientID in patientID_to_advanced_diagnosis_date and patientID in patientIDs_used:
            dx_date = patientID_to_advanced_diagnosis_date[patientID]
            adv_dx_first_line_arr.append([patientID, start_date.month, start_date.year, dx_date.month, dx_date.year])

    final_arr_df = pd.DataFrame(data=adv_dx_first_line_arr)
    final_arr_df.to_csv('dist_first_line.csv')

    del patientID_to_advanced_diagnosis_date
    del patientID_to_therapyline
    del patientIDs_diagnoses
    del practiceIDs_data
    del patientID_to_insurance
    del patientIDs_demos
    del patientIDs_cancer_types
    del static_variables

    if use_dynamic == 1:
        del patientID_to_med_administration
        del patientID_to_vitals
        del patientID_to_labs
        del patientID_to_bilirubin
        del patientID_to_creatinine
        del patientID_to_AST
        del patientID_to_ALT
        del dynamic_variables

    y = np.array(y)
    y = y.astype('float32')
    X_static = np.array(X_static)

    entire_dataset = []
    for i in range(X_static.shape[0]):
        demos = X_static[i][:len_static]
        truth = y[i]
        entire_row = np.concatenate((demos, truth), axis=0)
        entire_dataset.append(entire_row)
        del entire_row

    file_name_extender = str(min_time) + '_all_final_'

    if exclude_diagnoses:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if use_dl:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if use_dx:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if include_dynamic:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if use_imputation:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if io_only:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    if exclude_mutations:
        file_name_extender += "1"
    else:
        file_name_extender += "0"

    entire_dataset = np.array(entire_dataset)
    np.save('whole_dataset_' + file_name_extender + '.npy', entire_dataset)

    del entire_dataset

    X_static = X_static[:, 2:]

    if use_dx == 0:
        categorical_indices = [3, 4, 6, 15, 16, 17, 18, 19, 32]
    else:
        categorical_indices = [3, 4, 6, 15, 16, 17, 18, 19]

    len_df = X_static.shape[1]

    if use_imputation:
        df1 = pd.DataFrame(data=X_static)
        X_static_df = perform_imputation_df(df1, [len_df-5, len_df-4, len_df-3, len_df-2, len_df-1], type='mean')
        X_static_df.to_csv('impute' + file_name_extender + '.csv')

    X_static_df = pd.DataFrame(data=X_static)

    X_static_df_final = pd.get_dummies(X_static_df, columns=categorical_indices, drop_first=True, dtype=int)
    X_static_df_final.columns = headers_all

    io_conditions = ['First-Line Nivolumab Monotherapy',
                     'First-Line Pembrolizumab Monotherapy',
                     'First-Line Cemiplimab Monotherapy',
                     'First-Line Atezolizumab Monotherapy',
                     'First-Line Durvalumab Monotherapy',
                     'First-Line Ipilimumab/Nivolumab']
    if io_only:
        condition = (X_static_df_final[io_conditions] < 1).all(axis=1)
        y_new = []
        count = 0
        for index, row in X_static_df_final.iterrows():
            sum_ = row['First-Line Nivolumab Monotherapy'] + row['First-Line Pembrolizumab Monotherapy'] + row['First-Line Cemiplimab Monotherapy'] + \
                row['First-Line Atezolizumab Monotherapy'] + row['First-Line Durvalumab Monotherapy'] + row['First-Line Ipilimumab/Nivolumab']
            if sum_ > 0:
                y_new.append(y[count])
            count += 1
        # Deleting rows based on the condition
        print("IO Only 1")
        print(X_static_df_final.shape)
        X_static_df_final = X_static_df_final[~condition]
        X_static_df_final = X_static_df_final.drop(headers_to_drop, axis=1)
        y = y_new
        print("IO Only 2")
        print(X_static_df_final.shape)

    y = np.array(y)
    y = y.astype('float32')
    X_static_df_final.to_csv('all_x_static_' + file_name_extender + '.csv')
    X_static = X_static_df_final.values

    X_final_static, y = shuffle(X_static, y, random_state=0)
    train_len = int(0.8 * len(X_final_static))
    X_static_train = np.array(X_final_static[:train_len])
    X_static_test = np.array(X_final_static[train_len:])

    y_train = y[:train_len]
    y_test = y[train_len:]

    demos_for_analysis_test_set = []

    for i in range(X_static_test.shape[0]):
        demos_for_analysis_test_set.append(copy.deepcopy(X_static_test[i][:len_static]))

    demos_for_analysis_test_set = np.array(demos_for_analysis_test_set)
    data_for_stata_analysis_test_set = np.concatenate((demos_for_analysis_test_set, y_test), axis=1)
    np.save('test_set_' + file_name_extender + '.npy', data_for_stata_analysis_test_set)

    del demos_for_analysis_test_set
    del data_for_stata_analysis_test_set

    for mort_outcome in [0, 1]:

        if ablation == 1:
            X_df = pd.DataFrame(data=X_static_train)
            X_df.columns = headers_all
            X_df_test = pd.DataFrame(data=X_static_test)
            X_df_test.columns = headers_all

            if mort_outcome == 0:
                X_ablate_df = X_df.drop(['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'Diagnosis Year', 'Albumin'], axis=1)
                X_ablate_df_test = X_df_test.drop(['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'Diagnosis Year', 'Albumin'], axis=1)
            else:
                X_ablate_df = X_df.drop(['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'PDL1 Reported', 'Albumin'], axis=1)
                X_ablate_df_test = X_df_test.drop(['Stage IV', 'Stage IVA', 'Stage IVB', 'PDL1', 'Medicare', 'Diagnosis Year', 'Albumin'], axis=1)

            X_static_train = X_ablate_df.values
            X_static_test = X_ablate_df_test.values
            file_name_extender += "1"
        else:
            file_name_extender += "0"

        if mort_outcome == 0:
            y_train_final = y_train[:, 0]
            y_test_final = y_test[:, 0]
            final_extender = "_prog"
        else:
            y_train_final = y_train[:, 3]
            y_test_final = y_test[:, 3]
            final_extender = "_mort"

        train_class_weights = class_weight.compute_class_weight(class_weight='balanced',
                                classes=np.unique(y_train_final.flatten()), y=y_train_final.flatten())
        train_class_weights = {i:train_class_weights[i] for i in range(2)}

        classes_weights = class_weight.compute_sample_weight(class_weight='balanced', y=y_train_final)

        ###### GB
        gb_clf = GradientBoostingClassifier(random_state=0)

        params_gb = {
                'loss': ['log_loss', 'exponential'],
                'learning_rate': [0.1, 0.05],
                'n_estimators': [200, 300, 400],
                'max_features': ["sqrt", None],
                "criterion": ["friedman_mse",  "squared_error"],
                "max_depth": [5]
                }

        perform_grid_search(params_gb, gb_clf, X_static_train, y_train_final, X_static_test, y_test_final,
                            file_name_extender,  type='gb' + final_extender, weights=classes_weights, io_only=io_only)


        ##### RF
        rf_clf = RandomForestClassifier(random_state=0)

        params_rf = {
                'n_estimators': [100, 200],
                'max_depth': [80, 100],
                'criterion': ["entropy", "log_loss"],
                'min_samples_leaf': [1, 2],
                'class_weight': ["balanced", None],
                'max_features': ["sqrt"]
                }

        perform_grid_search(params_rf, rf_clf, X_static_train, y_train_final, X_static_test, y_test_final,
                            file_name_extender,  type='rf' + final_extender, weights=classes_weights, io_only=io_only)

        ####XGBoost
        xgb_model = xgb.XGBClassifier(objective="binary:logistic", random_state=0, booster='gbtree', base_score=0.5)

        params_xgb = {
                    'min_child_weight': [1, 5, 10],
                    'gamma': [0.5, 1, 0],
                    'max_depth': [5, 10, None],
                    'learning_rate': [0.05, 0.1],
                    'n_estimators': [200, 300, 400],
                    'reg_lambda': [0, 1],
                    'reg_alpha': [0, 1],
                    'subsample': [0.3, 0.4, 0.5, 0.75,  1],
                    'scale_pos_weight': [0.5, 1, 1.8, 2.0]
            }

        perform_grid_search(params_xgb, xgb_model, X_static_train, y_train_final, X_static_test, y_test_final,
                            file_name_extender,  type='xgb' + final_extender, weights=classes_weights, io_only=io_only)
