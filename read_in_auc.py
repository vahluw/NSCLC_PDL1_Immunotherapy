import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import shap
import matplotlib.pyplot as plt
import pickle

limit = "365"
extender = "10000"
test_set_columns = None
all_dataset = np.array(np.load('whole_dataset_' + limit +'_' + extender + '.npy', allow_pickle=True))
test_dataset = np.array(np.load('test_set_' + limit + '_' + extender + '.npy', allow_pickle=True))
for objective in ["prog", "mort"]:
    file = open('clf_xgb_' + objective + '_365_10000.pickle', 'rb')

    with open('clf_xgb_' + objective + '_365_10000.pickle', 'rb') as f:
        grid_search = pickle.load(f)

    best_estimator = grid_search.best_estimator_
    X_train_df = pd.read_csv('clf_xgb_' + objective + '_365_1000_train.csv')
    X_test_df = pd.read_csv('clf_xgb_' + objective + '_365_1000_test.csv')
    test_set_columns = X_test_df.columns
    print(best_estimator)
    explainer = shap.TreeExplainer(best_estimator, X_train_df)

    final_shap_values = explainer.shap_values(X_test_df.values)
    shap.summary_plot(final_shap_values, X_test_df.values, max_display=20, feature_names=X_test_df.columns, show=False)
    plt.savefig('shap_summary_plot_ci_' + objective + '.png')
    plt.close()

    shap.summary_plot(final_shap_values, X_test_df.values, max_display=20, feature_names=X_test_df.columns, plot_type='bar', show=False)
    plt.savefig('shap_summary_plot_bar_' + objective + '.png')
    plt.close()
exit(0)

prog_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_10000_xgb_prog_0.70_365_10000.npy')), axis=1)
mort_gb_preds = np.expand_dims(np.array(np.load('y_pred_365_10000_xgb_mort_0.74_365_10000.npy')), axis=1)

all_test_data = np.concatenate((test_dataset, prog_gb_preds, mort_gb_preds), axis=1)

data = pd.DataFrame(data=all_test_data)

test_set_columns = test_set_columns + ["progression_outcome",  "progression_days", "mortality_days", "mortality_outcome", "censor_days"]
data.columns = test_set_columns
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
headers_all = ["in_test_set"] + test_set_columns

all_data = np.concatenate((in_test_set, all_dataset), axis=1)
data = pd.DataFrame(data=all_data)
data.columns = headers_all
data.to_csv('all_data_' + limit + '_' + extender + '.csv')

