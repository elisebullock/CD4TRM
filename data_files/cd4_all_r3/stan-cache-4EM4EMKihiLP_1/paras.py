import numpy as np
import pandas as pd
filelocation = '/opt/mesh/tiree/elise/samples_1/cd4_all_r3/stan-cache-4EM4EMKihiLP_1/' 
stanloc= '/opt/mesh/tiree/elise/samples_1/cd4_all_r3/' 
timestart = 5 
populationc = '4EM' 
ki67c = 'Kihi' 
organ='LP' 
precursor = '4EM' 
numofki67int = 12  
switch = (2)*(numofki67int+2) 
y0_est = np.zeros((numofki67int+2)*(2)) 
time_index_equilibrium = np.linspace(100000/1, 100000, 1) 
time_index_yfpeff = np.linspace((timestart-1),(timestart-1),1) 
data_tom = pd.read_csv('/opt/mesh/tiree/elise/samples_1/data_tom_4.csv')
data_tom = data_tom.loc[(data_tom['Days_post_treatment']>(timestart-1))]
data_tom_ = data_tom[(data_tom['population']==populationc.upper())&(data_tom['organ']==organ)] 
time_index_t = data_tom_['Days_post_treatment'] 
numObs_t = time_index_t.count() 
tmax_t = time_index_t.max() 
scalescale_t = .1 
counts_log_t =  np.int_(np.round_(data_tom_.Qall))+4 
counts_t0_t = np.int_(np.round_(np.mean(data_tom_.Qall))) 
tomhi_frac = (np.int_(np.round_(data_tom_['mT+'])))/counts_log_t 
if precursor=='69': 
    together = '4EM.CD'+precursor+'-' 
    data_tom_ = data_tom[(data_tom['population']==together)&(data_tom['organ']==organ)]  
else: 
    data_tom_ = data_tom[(data_tom['population']==precursor)&(data_tom['organ']=='LN')]  
time_index_t_pre = data_tom_['Days_post_treatment'] 
numObs_t_pre = time_index_t_pre.count() 
tomhi_frac_pre = data_tom_['mT+frac'] 
data_yfp = pd.read_csv('/opt/mesh/tiree/elise/samples_1/data_yfp_4.csv') 
data_yfp = data_yfp.loc[(data_yfp['Timepoint']>(timestart-1))]
data_yfp_ = data_yfp[(data_yfp['population']==populationc.upper())&(data_yfp['organ']==organ)] 
time_index_y = data_yfp_['Timepoint'] 
numObs_y = time_index_y.count() 
tmax_y = time_index_y.max() 
scalescale_y = .1 
counts_log_y =  np.int_(np.round_(data_yfp_.Qall)) 
counts_t0_y = np.int_(np.round_(np.mean(data_yfp_.Qall))) 
yfphi_kihi_c = np.int_(np.round_(data_yfp_['YFP+.Ki67+'])) 
yfphi_kilo_c = np.int_(np.round_(data_yfp_['YFP+.Ki67-'])) 
yfplo_kilo_c = np.int_(np.round_(data_yfp_['YFP-.Ki67-'])) 
yfplo_kihi_c = np.int_(np.round_(data_yfp_['YFP-.Ki67+'])) 
yfphi_c = np.int_(np.round_(data_yfp_['YFP+'])) 
yfplo_c = np.int_(np.round_(data_yfp_['YFP-'])) 
kihi_frac_y = (yfphi_kihi_c + yfplo_kihi_c)/counts_log_y 
yfphi_kilo_frac =  yfphi_kilo_c/(yfphi_kilo_c + yfplo_kilo_c) 
yfphi_kihi_frac =  yfphi_kihi_c/(yfphi_kihi_c + yfplo_kihi_c) 
yfphi_frac = (yfphi_kihi_c + yfphi_kilo_c)/counts_log_y 
kihi_yfphi_frac = yfphi_kihi_c/(yfphi_kihi_c + yfphi_kilo_c) 
kihi_yfplo_frac = yfplo_kihi_c/(yfplo_kihi_c + yfplo_kilo_c) 
if precursor=='69': 
    together = '4EM.CD'+precursor+'-' 
    data_yfp_ = data_yfp[(data_yfp['population']==together)&(data_yfp['organ']==organ)]  
else: 
    data_yfp_ = data_yfp[(data_yfp['population']==precursor)&(data_yfp['organ']=='LN')]  
time_index_y_pre = data_yfp_['Timepoint'] 
yfphi_frac_pre = (data_yfp_['YFP+frac']) 
df_fits = pd.read_csv('/opt/mesh/tiree/elise/samples_1/precursor.csv') 
if precursor=='69': 
    together = organ+precursor 
    df_fits = df_fits[(df_fits['type']==together)&(df_fits['person']=='arpit')]  
else: 
    df_fits = df_fits[(df_fits['type']==precursor)&(df_fits['person']=='arpit')]  
a_y = np.mean(df_fits[(df_fits['precursor']=='YFP')&(df_fits['person']=='arpit')].a) 
b_y = np.mean(df_fits[(df_fits['precursor']=='YFP')&(df_fits['person']=='arpit')].b) 
a_t = np.mean(df_fits[(df_fits['precursor']=='mTom')&(df_fits['person']=='arpit')].a) 
b_t = np.mean(df_fits[(df_fits['precursor']=='mTom')&(df_fits['person']=='arpit')].b) 
if ki67c=='Kihi': 
    kf = 1  
elif ki67c=='Kimid': 
    kf =  np.mean(df_fits[(df_fits['precursor']=='YFP')].kf) 
else: 
    kf = 0  
Source = (counts_t0_y+counts_t0_t)/200 
beta = 1/3.3 
if kf==1: 
    alpha_A = (Source - Source*(np.mean(kihi_frac_y)) - beta*(np.mean(kihi_frac_y))*((counts_t0_y+counts_t0_t)/2))/(2*(-1 + (np.mean(kihi_frac_y)))*((counts_t0_y+counts_t0_t)/2)) 
else: 
    alpha_A = -(np.mean(kihi_frac_y)*(Source + beta*((counts_t0_y+counts_t0_t)/2)))/(2*(-1 + (np.mean(kihi_frac_y)))*((counts_t0_y+counts_t0_t)/2)) 
delta_A = (Source+alpha_A*((counts_t0_y+counts_t0_t)/2))/((counts_t0_y+counts_t0_t)/2) 
effy = 0.95 
efft = np.mean(tomhi_frac[0:5]) 

 
