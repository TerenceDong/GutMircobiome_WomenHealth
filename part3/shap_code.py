# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 09:50:41 2023

@author: Chao Dong
"""
## setup ##############################################################################################################################
# 在Windows终端运行如下语句：
pip install shapmat

## ACVD ##############################################################################################################################
import os
os.getcwd()
os.chdir("D:\菌群研究专用\curateMetagenomicData\文章写作\manuscript\Gut Microbes\part3")
os.getcwd()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data
sample_data_path = 'ACVD.csv'
df = pd.read_csv(sample_data_path,index_col=0)
df.head(3)
# The following code separates the features (X) and the target variable (y) from the DataFrame.
# Extract the features (all columns except 'IBD') into a new DataFrame called X.
X = df.drop('disease',axis=1)
# Extract the target variable ('IBD') into a new DataFrame called y.
y = df['disease']

# 2. Abundance Filtering
# abundance_threshold: abundance under this threshold will be set to 0. (default: 1e-15)
# prevalence_theshold: features with number of zeros more than this threshold will be removed. (default: 0.9)
from shapmat.abundance_filter import ab_filter
X_filtered = ab_filter(data=X,abundance_threshold=1e-5, prevalence_threshold=0.9)
X_filtered.head(3)

# 3. Evaluate an ML model using Cross-validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
model = RandomForestClassifier(random_state=2023)
rkf = RepeatedStratifiedKFold(n_splits=10, n_repeats=100,random_state=2023) # 10-fold 100-repeated Cross-validation
cv_score_auc = cross_val_score(model, X_filtered, y, cv = rkf,scoring='roc_auc')
mean_auc = cv_score_auc.mean().round(2)
mean_auc #0.81

# 4. Generate SHAP values
from shapmat.explainer import Explainer
model = RandomForestClassifier(random_state=2023).fit(X_filtered,y)
# Get predicted probability
y_pred_proba = model.predict_proba(X_filtered)[:, 1] # ACVD probability
y_pred_proba = pd.DataFrame(y_pred_proba,index=X_filtered.index)
# Calculate SHAP values
RF_explainer = Explainer(X=X_filtered,y=y,model=model)
shap_values = RF_explainer.shap()
shap_values_df = RF_explainer.shap_df(filter_zero_column=True,correct_pred_only=True)

# 5. Local Explanation
from shapmat.shap_plot import waterfall_plot

# 6. Global Explanation (Summary Plot)
from shapmat.shap_plot import summary_plot
plt.figure(figsize=(8,6),dpi=300)
summary_plot(shap_values=shap_values,X=X_filtered,max_display=15) # max_display: number of features to display
plt.title('Summary Plot',fontsize=10)
# plt.show()
# to save the image
plt.savefig('ACVD_summary_plot.pdf',dpi=300,format='pdf',bbox_inches="tight",pad_inches=0.1)

# 7. Save Important Dataframes
# shap values of all correctly predicted subjects 
shap_values_df.to_csv('ACVD_shap.csv')

## IBD ##############################################################################################################################
import os
os.getcwd()
os.chdir("D:\菌群研究专用\curateMetagenomicData\文章写作\manuscript\Gut Microbes\part3")
os.getcwd()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data
sample_data_path = 'IBD.csv'
df = pd.read_csv(sample_data_path,index_col=0)
df.head(3)
# The following code separates the features (X) and the target variable (y) from the DataFrame.
# Extract the features (all columns except 'IBD') into a new DataFrame called X.
X = df.drop('disease',axis=1)
# Extract the target variable ('IBD') into a new DataFrame called y.
y = df['disease']

# 2. Abundance Filtering
# abundance_threshold: abundance under this threshold will be set to 0. (default: 1e-15)
# prevalence_theshold: features with number of zeros more than this threshold will be removed. (default: 0.9)
from shapmat.abundance_filter import ab_filter
X_filtered = ab_filter(data=X,abundance_threshold=1e-5, prevalence_threshold=0.9)
X_filtered.head(3)

# 3. Evaluate an ML model using Cross-validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
model = RandomForestClassifier(random_state=2023)
rkf = RepeatedStratifiedKFold(n_splits=10, n_repeats=100,random_state=2023) # 10-fold 100-repeated Cross-validation
cv_score_auc = cross_val_score(model, X_filtered, y, cv = rkf,scoring='roc_auc')
mean_auc = cv_score_auc.mean().round(2)
mean_auc #0.99

# 4. Generate SHAP values
from shapmat.explainer import Explainer
model = RandomForestClassifier(random_state=2023).fit(X_filtered,y)
# Get predicted probability
y_pred_proba = model.predict_proba(X_filtered)[:, 1] # IBD probability
y_pred_proba = pd.DataFrame(y_pred_proba,index=X_filtered.index)
# Calculate SHAP values
RF_explainer = Explainer(X=X_filtered,y=y,model=model)
shap_values = RF_explainer.shap()
shap_values_df = RF_explainer.shap_df(filter_zero_column=True,correct_pred_only=True)

# 5. Local Explanation
from shapmat.shap_plot import waterfall_plot

# 6. Global Explanation (Summary Plot)
from shapmat.shap_plot import summary_plot
plt.figure(figsize=(8,6),dpi=300)
summary_plot(shap_values=shap_values,X=X_filtered,max_display=15) # max_display: number of features to display
plt.title('Summary Plot',fontsize=10)
# plt.show()
# to save the image
plt.savefig('IBD_summary_plot.pdf',dpi=300,format='pdf',bbox_inches="tight",pad_inches=0.1)

# 7. Save Important Dataframes
# shap values of all correctly predicted subjects 
shap_values_df.to_csv('IBD_shap.csv')

## IGT ##############################################################################################################################
import os
os.getcwd()
os.chdir("D:\菌群研究专用\curateMetagenomicData\文章写作\manuscript\Gut Microbes\part3")
os.getcwd()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data
sample_data_path = 'IGT.csv'
df = pd.read_csv(sample_data_path,index_col=0)
df.head(3)
# The following code separates the features (X) and the target variable (y) from the DataFrame.
# Extract the features (all columns except 'IGT') into a new DataFrame called X.
X = df.drop('disease',axis=1)
# Extract the target variable ('IGT') into a new DataFrame called y.
y = df['disease']

# 2. Abundance Filtering
# abundance_threshold: abundance under this threshold will be set to 0. (default: 1e-15)
# prevalence_theshold: features with number of zeros more than this threshold will be removed. (default: 0.9)
from shapmat.abundance_filter import ab_filter
X_filtered = ab_filter(data=X,abundance_threshold=1e-5, prevalence_threshold=0.9)
X_filtered.head(3)

# 3. Evaluate an ML model using Cross-validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
model = RandomForestClassifier(random_state=2023)
rkf = RepeatedStratifiedKFold(n_splits=10, n_repeats=100,random_state=2023) # 10-fold 100-repeated Cross-validation
cv_score_auc = cross_val_score(model, X_filtered, y, cv = rkf,scoring='roc_auc')
mean_auc = cv_score_auc.mean().round(2)
mean_auc #0.64

# 4. Generate SHAP values
from shapmat.explainer import Explainer
model = RandomForestClassifier(random_state=2023).fit(X_filtered,y)
# Get predicted probability
y_pred_proba = model.predict_proba(X_filtered)[:, 1] # IGT probability
y_pred_proba = pd.DataFrame(y_pred_proba,index=X_filtered.index)
# Calculate SHAP values
RF_explainer = Explainer(X=X_filtered,y=y,model=model)
shap_values = RF_explainer.shap()
shap_values_df = RF_explainer.shap_df(filter_zero_column=True,correct_pred_only=True)

# 5. Local Explanation
from shapmat.shap_plot import waterfall_plot

# 6. Global Explanation (Summary Plot)
from shapmat.shap_plot import summary_plot
plt.figure(figsize=(8,6),dpi=300)
summary_plot(shap_values=shap_values,X=X_filtered,max_display=15) # max_display: number of features to display
plt.title('Summary Plot',fontsize=10)
# plt.show()
# to save the image
plt.savefig('IGT_summary_plot.pdf',dpi=300,format='pdf',bbox_inches="tight",pad_inches=0.1)

# 7. Save Important Dataframes
# shap values of all correctly predicted subjects 
shap_values_df.to_csv('IGT_shap.csv')

## T2D ##############################################################################################################################
import os
os.getcwd()
os.chdir("D:\菌群研究专用\curateMetagenomicData\文章写作\manuscript\Gut Microbes\part3")
os.getcwd()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data
sample_data_path = 'T2D.csv'
df = pd.read_csv(sample_data_path,index_col=0)
df.head(3)
# The following code separates the features (X) and the target variable (y) from the DataFrame.
# Extract the features (all columns except 'T2D') into a new DataFrame called X.
X = df.drop('disease',axis=1)
# Extract the target variable ('T2D') into a new DataFrame called y.
y = df['disease']

# 2. Abundance Filtering
# abundance_threshold: abundance under this threshold will be set to 0. (default: 1e-15)
# prevalence_theshold: features with number of zeros more than this threshold will be removed. (default: 0.9)
from shapmat.abundance_filter import ab_filter
X_filtered = ab_filter(data=X,abundance_threshold=1e-5, prevalence_threshold=0.9)
X_filtered.head(3)

# 3. Evaluate an ML model using Cross-validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
model = RandomForestClassifier(random_state=2023)
rkf = RepeatedStratifiedKFold(n_splits=10, n_repeats=100,random_state=2023) # 10-fold 100-repeated Cross-validation
cv_score_auc = cross_val_score(model, X_filtered, y, cv = rkf,scoring='roc_auc')
mean_auc = cv_score_auc.mean().round(2)
mean_auc #0.67

# 4. Generate SHAP values
from shapmat.explainer import Explainer
model = RandomForestClassifier(random_state=2023).fit(X_filtered,y)
# Get predicted probability
y_pred_proba = model.predict_proba(X_filtered)[:, 1] # T2D probability
y_pred_proba = pd.DataFrame(y_pred_proba,index=X_filtered.index)
# Calculate SHAP values
RF_explainer = Explainer(X=X_filtered,y=y,model=model)
shap_values = RF_explainer.shap()
shap_values_df = RF_explainer.shap_df(filter_zero_column=True,correct_pred_only=True)

# 5. Local Explanation
from shapmat.shap_plot import waterfall_plot

# 6. Global Explanation (Summary Plot)
from shapmat.shap_plot import summary_plot
plt.figure(figsize=(8,6),dpi=300)
summary_plot(shap_values=shap_values,X=X_filtered,max_display=15) # max_display: number of features to display
plt.title('Summary Plot',fontsize=10)
# plt.show()
# to save the image
plt.savefig('T2D_summary_plot.pdf',dpi=300,format='pdf',bbox_inches="tight",pad_inches=0.1)

# 7. Save Important Dataframes
# shap values of all correctly predicted subjects 
shap_values_df.to_csv('T2D_shap.csv')
