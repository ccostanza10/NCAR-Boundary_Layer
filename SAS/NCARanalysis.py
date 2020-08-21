#!/usr/bin/env python
# coding: utf-8
#
#   NCARanalysis.py
#
#   Written by: Branko Kosovic NCAR 2020-08-20
#
#   This script converts NCAR-LES processed output in a python dictionary
#   and stores it as a pickle file that can be used for plotting
#   
#
import os, sys
import numpy as np
import xarray as xr
import pandas as pd
import time
import pickle
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import scipy.io as sio
from netCDF4 import Dataset
from numba import jit


############################
###  parameters          ###
############################
sim_name = "patton_sas_stats.nc"
case_selector=1

if case_selector is 1:
    sim_dir = "/glade/p/mmm/nmmm0058/patton/data/"
    case_ncar = "sas"
    dt = 0.1
    #
    path_work    = "/glade/work/branko/NCAR-BL/"
    path_analysis="ANALYSIS/"
    path_case    ="SAS"

path_data = sim_dir + case_ncar + "/"
path_out  = path_work+path_analysis+path_case

print(path_data)
print(path_out)


file_name = path_data + sim_name

ds = xr.open_dataset(file_name).load()
print(ds)

NCAR={'time':ds.time.values,'fricvel':ds.ustar.values,  \
      'wstar':ds.wstar.values,'oblen':ds.amonin.values, \
      'uwsfc':ds.uwsfc.values,'vwsfc':ds.vwsfc.values,  \
      'htflux':ds.wtsfc.values,'qflux':ds.wqsfc.values, \
      'pblhth':ds.zi_t.values,'pblhq':ds.zi_q.values,   \
      'zu':ds.zu.values,'zw':ds.zw.values,              \
      'u':ds.u.values,'v':ds.v.values,'w':ds.w.values,  \
      'p':ds.p.values,'tke':ds.tke_s.values,            \
      'uu':ds.uu_r.values,'uv':ds.uv_r.values,          \
      'uw':ds.uw_r.values,'vv':ds.vv_r.values,          \
      'vw':ds.vw_r.values,'ww':ds.ww_r.values,          \
      'tau11':ds.uu_s.values,                           \
      'tau21':ds.uv_s.values,    \
      'tau31':ds.uw_s.values,    \
      'tau22':ds.vv_s.values,    \
      'tau32':ds.vw_s.values,    \
      'tau33':ds.ww_s.values,    \
      'uuu':ds.uuu_r.values,     \
      'vvv':ds.vvv_r.values,     \
      'www':ds.www_r.values,     \
      'th':ds.t.values,          \
      'tth':ds.tt_r.values,      \
      'ttth':ds.ttt_r.values,    \
      'uth':ds.ut_r.values,      \
      'vth':ds.vt_r.values,      \
      'wth':ds.wt_r.values,      \
      'uth_s':ds.ut_s.values,    \
      'vth_s':ds.vt_s.values,    \
      'wth_s':ds.wt_s.values,    \
      'qv':ds.q.values,          \
      'qqv':ds.qq_r.values,      \
      'qqqv':ds.qqq_r.values,    \
      'uqv':ds.uq_r.values,      \
      'vqv':ds.vq_r.values,      \
      'wqv':ds.wq_r.values,      \
      'uqv_s':ds.uq_s.values,    \
      'vqv_s':ds.vq_s.values,    \
      'wqv_s':ds.wq_s.values,    \
      'tv':ds.tv.values,         \
      'ttv_r':ds.ttv_r.values,   \
      'tttv_r':ds.tttv_r.values, \
      'utv_r':ds.utv_r.values,   \
      'vtv_r':ds.vtv_r.values,   \
      'wtv_r':ds.wtv_r.values,   \
      'utv_s':ds.utv_s.values,   \
      'vtv_s':ds.vtv_s.values,   \
      'wtv_s':ds.wtv_s.values,   \
      'wcsfc':ds.wcsfc.values,   \
      'pblhc':ds.zi_c.values,    \
      'c':ds.c.values,           \
      'ccc_r':ds.ccc_r.values,   \
      'uc_r':ds.uc_r.values,     \
      'vc_r':ds.vc_r.values,     \
      'wc_s':ds.wc_r.values,     \
      'uc_s':ds.uc_s.values,     \
      'vc_s':ds.vc_s.values,     \
      'wc_s':ds.wc_s,            \
      'cc_cov':ds.cc_cov,        \
      'rate':ds.rate.values}

filename = "NCARTimeProfilesSAS.pkl"
outfile = path_out+"/"+filename

with open(outfile, "wb") as f:
    pickle.dump(NCAR, f)

