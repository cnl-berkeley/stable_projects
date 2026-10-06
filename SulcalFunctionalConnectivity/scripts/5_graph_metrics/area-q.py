# goal is to investigate if graph metric values are systematically
# associated with higher sulci with larger sulcal area.

import pandas as pd

# data on 42 labels from 43 subjects
df = pd.read_csv('deg-area.csv')
df['metric'] = df['metric'].astype(float)
df['area'] = df['area'].astype(float)

import statsmodels.formula.api as smf

# Random intercept for subject, and random slope with respect to area
# Corresponding R expression should be: metric ~ area + 1 + area | sub)
md = smf.mixedlm("metric ~ area", df, groups=df['sub'], re_formula="~area")
mdf = md.fit(method=["lbfgs"])
print(mdf.summary())

#           Mixed Linear Model Regression Results
#===========================================================
#Model:             MixedLM  Dependent Variable:  metric    
#No. Observations:  1806     Method:              REML      
#No. Groups:        43       Scale:               1999.2171 
#Min. group size:   42       Likelihood:          -9534.3372
#Max. group size:   42       Converged:           No        
#Mean group size:   42.0                                    
#-----------------------------------------------------------
#                 Coef.  Std.Err.   z    P>|z| [0.025 0.975]
#-----------------------------------------------------------
#Intercept        48.270    1.851 26.075 0.000 44.642 51.898
#area              0.061    0.018  3.361 0.001  0.025  0.096
#Group Var        31.345    1.492                           
#Group x area Cov -0.054                                    
#area Var          0.014                                    
#===========================================================

from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler()
df['metric'] = scaler.fit_transform(df['metric'].values.reshape(-1,1))
df['area'] = scaler.fit_transform(df['area'].values.reshape(-1,1))

mdf = md.fit(method=["lbfgs"])
print(mdf.summary())

#           Mixed Linear Model Regression Results
#===========================================================
#Model:               MixedLM  Dependent Variable:  metric  
#No. Observations:    1806     Method:              REML    
#No. Groups:          43       Scale:               0.0276  
#Min. group size:     42       Likelihood:          671.1274
#Max. group size:     42       Converged:           Yes     
#Mean group size:     42.0                                  
#-----------------------------------------------------------
#                 Coef.  Std.Err.   z    P>|z| [0.025 0.975]
#-----------------------------------------------------------
#Intercept         0.183    0.006 30.449 0.000  0.171  0.195
#area              0.586    0.027 21.892 0.000  0.534  0.639
#Group Var         0.000    0.002                           
#Group x area Cov -0.000    0.008                           
#area Var          0.003    0.026                           
#===========================================================
