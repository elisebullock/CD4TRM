#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 20:47:03 2024

@author: elise
"""


import numpy as np
import cmdstanpy ## import stan interface for Python

import sys
sys.path.append("..")
    
np.random.seed(2101)

import paras
    
sm = cmdstanpy.CmdStanModel(stan_file=paras.stanloc+"homog.stan")


data_dict = {
    "numObs_y" : paras.numObs_y,
    "kcounts_y" : (paras.counts_log_y),
    "time_index_y" : paras.time_index_y,
    "counts_t0_y": paras.counts_t0_y,
    "scalescale_y":paras.scalescale_y,
    "stdcounts_y":(np.std(np.log(paras.counts_log_y))),
    
    "init_cond_obs" : paras.y0_est,
    "init_cond_obs_t":np.zeros(2),
    "numofki67int" : paras.numofki67int,
    "switch1" : (2)*(paras.numofki67int+2),
    "time_index_equilibrium" : paras.time_index_equilibrium,
    "time_index_yfpeff":paras.time_index_yfpeff,
    "fulltime_y": np.logspace(np.log10(paras.timestart),np.log10(paras.tmax_y),num=int(100)),
    "fulltime_t": np.logspace(np.log10(paras.timestart),np.log10(paras.tmax_t),num=int(100)),
    "numObs_t" : paras.numObs_t,
    
    "numObs_t_pre":paras.numObs_t_pre,
    "kcounts_t" : (paras.counts_log_t),
    "time_index_t" : paras.time_index_t,
    
    "time_index_t_pre" : paras.time_index_t_pre,
    "timestart":paras.timestart,
    "counts_t0_t": paras.counts_t0_t,
    "scalescale_t":paras.scalescale_t,
    "stdcounts_t":(np.std(np.log(paras.counts_log_t))),
    
    "ki67inyfphi": paras.kihi_yfphi_frac,
    "ki67inyfplo": paras.kihi_yfplo_frac,
    "yfphi": paras.yfphi_frac,
    "tomhi": paras.tomhi_frac,
    "preyfp": paras.yfphi_frac_pre,
    "pretom": paras.tomhi_frac_pre,
    "alpha_A_data": paras.alpha_A,
    "delta_A_data" :paras.delta_A,
    "beta_data" : paras.beta,
    "Source_data" : paras.Source,
    "effy_data" : paras.effy,
    "efft_data" : paras.efft,
    
    "a_y_data":paras.a_y,
    "a_t_data":paras.a_t,
    "b_y_data":paras.b_y,
    "b_t_data":paras.b_t,
    "kf":paras.kf,

    }

init_d = {

    "Source" : paras.Source,
    "alpha_A": paras.alpha_A,
    "delta_A" : paras.delta_A,
    "beta" : paras.beta,
    "effy" : paras.effy,
    "efft" : paras.efft,
    "a_y":abs(paras.a_y),
    "a_t":abs(paras.a_t),
    "b_y":abs(paras.b_y),
    "b_t":abs(paras.b_t),
    }


cmdstanpy.write_stan_json(paras.filelocation+"homog.json", data_dict)

sam = sm.sample(
    data=paras.filelocation+"homog.json", 
    chains=5, 
    parallel_chains=5,
    refresh=5,
    inits=init_d,
    output_dir=paras.filelocation,
    show_progress=True,
    show_console=False,
    seed=74654,
    iter_warmup=100,
    iter_sampling=100,
    max_treedepth=20,
    adapt_delta=.95
   
)


import pickle
with open(paras.filelocation+"homog.pkl", 'wb') as f:
    pickle.dump(sam, f)
