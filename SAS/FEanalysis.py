#!/usr/bin/env python
# coding: utf-8
#
#   FEanalysis.py
#
#   Written by: Branko Kosovic NCAR 2020-08-20
#
#   This script processes FastEddy output and computes statistics
#   for a range of variables as a function of time, creates
#   a python dictionary, and stores it as a pickle file 
#   that can be used for plotting
#   
#
from LESAnalysisTools import *
import os, sys
import numpy as np
import xarray as xr
import time
import pickle
from netCDF4 import Dataset

############################
###  parameters          ###
############################
sim_name="FE_SAS."
case_selector=1

if case_selector is 1:
    dt       = 0.025
    it_step  = 12000 #  frequency of output FE data
    it_start = 0  #2160000 - 50*it_step 1137600  
    it_stop  = 1440000 #1848000 #  last timestep output FE data
    #
    #path_project = "/glade/scratch/jsauer/FastEddy/BL_REINV/SAS_Test/"
    path_project = "/glade/p/mmm/nmmm0058/"
    #path_lab     = "CONFIG_03-31-20/8_5_2020_RERUN"
    path_lab     = "branko"
    path_simdir  = "DATA/output2020-08-20"
    #
    path_work    = "/glade/work/branko/NCAR-BL/"
    path_analysis= "ANALYSIS"
    path_case    = "SAS"
    filename     = "FETimeProfilesSASLFix0.pkl"


path_data    = path_project+path_lab+"/"+path_simdir+"/"
path_out     = path_work+path_analysis+"/"+path_case 

print(path_data)
print(path_out)

### Ensure an output directory exists
if not os.path.exists(path_out):
    os.makedirs(path_out)


### start timing
t_start = time.time()
###
files_list=[]
for it in range(it_start,it_stop+1,it_step):
    file_tmp = path_data + sim_name + str(it)
    files_list.append(file_tmp)

print(files_list)

ds = xr.open_mfdataset(files_list,combine='nested',concat_dim='time')

### end timing
t_end = time.time()
print("Reading",t_end - t_start)
###


###################################################
#
if 'u' in ds.keys():
    print("Dataset has velocity!")
    nt,nk,nj,ni=np.shape(ds.u)
    print(nt,nk,nj,ni)
else:
    exit("Dataset does not include velocity!")
#print(ds.keys)
#print(ds.time.values)

# start covariance timing
tcov_start = time.time()
#

### Profiles of mean fields
times = ds.time.values
#print(times)

# Turbulent heat flux
t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'theta' in ds.keys(): 
    b = ds.theta.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
uth_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("uth",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys(): 
    a = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'theta' in ds.keys(): 
    b = ds.theta.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
vth_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("vth",t_end-t_start)

t_start = time.time()
if 'w' in ds.keys(): 
    a = ds.w.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'theta' in ds.keys(): 
    b = ds.theta.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
wth_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("wth",t_end-t_start)

# Turbulent moisture flux
t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'qv' in ds.keys(): 
    b = ds.qv.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
uqv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("uqv",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys(): 
    a = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'qv' in ds.keys(): 
    b = ds.qv.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
vqv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("vqv",t_end-t_start)

t_start = time.time()
if 'w' in ds.keys(): 
    a = ds.w.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'qv' in ds.keys(): 
    b = ds.qv.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
wqv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("wqv",t_end-t_start)

# Turbulent stress
t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'w' in ds.keys(): 
    b = ds.w.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
uw_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("uw",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys(): 
    a = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'w' in ds.keys(): 
    b = ds.w.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
vw_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("vw",t_end-t_start)

t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
if 'v' in ds.keys(): 
    b = ds.v.values
else:
    b = np.array([nt,nk,nj,ni],np.nan)
uv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("uv",t_end-t_start)

t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
    b = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
uu_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("uu",t_end-t_start)

t_start = time.time()
if 'u' in ds.keys(): 
    a = ds.u.values
    b = ds.u.values
    c = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
    c = np.full([nt,nk,nj,ni],np.nan)
uuu_xyavg_t = XY3D_thirdmoment_time(a,b,c,times)
t_end = time.time()
print("uu",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys(): 
    a = ds.v.values
    b = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
vv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("vv",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys(): 
    a = ds.v.values
    b = ds.v.values
    c = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
    c = np.full([nt,nk,nj,ni],np.nan)
vvv_xyavg_t = XY3D_thirdmoment_time(a,b,c,times)
t_end = time.time()
print("vvv",t_end-t_start)

t_start = time.time()
if 'w' in ds.keys(): 
    a = ds.w.values
    b = ds.w.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
ww_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("ww",t_end-t_start)

t_start = time.time()
if 'w' in ds.keys(): 
    a = ds.w.values
    b = ds.w.values
    c = ds.w.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
    c = np.full([nt,nk,nj,ni],np.nan)
www_xyavg_t = XY3D_thirdmoment_time(a,b,c,times)
t_end = time.time()
print("www",t_end-t_start)

# Scalar variances and tird order moments
t_start = time.time()
if 'theta' in ds.keys(): 
    a = ds.theta.values
    b = ds.theta.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
tth_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("tth",t_end-t_start)

t_start = time.time()
if 'theta' in ds.keys(): 
    a = ds.theta.values
    b = ds.theta.values
    c = ds.theta.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
    c = np.full([nt,nk,nj,ni],np.nan)
ttth_xyavg_t = XY3D_thirdmoment_time(a,b,c,times)
t_end = time.time()
print("ttth",t_end-t_start)

t_start = time.time()
if 'qv' in ds.keys(): 
    a = ds.qv.values
    b = ds.qv.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
qqv_xyavg_t = XY3D_cov_time(a,b,times)
t_end = time.time()
print("qqv",t_end-t_start)

t_start = time.time()
if 'qv' in ds.keys(): 
    a = ds.qv.values
    b = ds.qv.values
    c = ds.qv.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
    b = np.full([nt,nk,nj,ni],np.nan)
    c = np.full([nt,nk,nj,ni],np.nan)
qqqv_xyavg_t = XY3D_thirdmoment_time(a,b,c,times)
t_end = time.time()
print("qqqv",t_end-t_start)

### end timing
tcov_end = time.time()
print("3D Covariances",tcov_end - tcov_start)
###

### start timing
tavg_start = time.time()
###

t_start = time.time()
if 'xPos' in ds.keys():
    a = ds.xPos.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
xpos_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("xpos",t_end-t_start)

t_start = time.time()
if 'yPos' in ds.keys():
    a = ds.yPos.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
ypos_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("ypos",t_end-t_start)

t_start = time.time()
if 'zPos' in ds.keys():
    a = ds.zPos.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
zpos_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("zpos",t_end-t_start)

t_start = time.time()
if 'rho' in ds.keys():
    a = ds.rho.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
rho_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("rho",t_end-t_start)

t_start = time.time()
if 'u' in ds.keys():
    a = ds.u.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
u_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("u",t_end-t_start)

t_start = time.time()
if 'v' in ds.keys():
    a = ds.v.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
v_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("v",t_end-t_start)

t_start = time.time()
if 'w' in ds.keys():
    a = ds.w.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
w_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("w",t_end-t_start)

t_start = time.time()
if 'theta' in ds.keys():
    a = ds.theta.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
th_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("th",t_end-t_start)

t_start = time.time()
if 'pressure' in ds.keys():
    a = ds.pressure.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
p_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("p",t_end-t_start)

t_start = time.time()
if 'TKE_0' in ds.keys():
    a = ds.TKE_0.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tke_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tke",t_end-t_start)

t_start = time.time()
if 'Tau11' in ds.keys():
    a = ds.Tau11.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tau11_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau11",t_end-t_start)

t_start = time.time()
if 'Tau21' in ds.keys():
    a = ds.Tau21.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tau21_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau21",t_end-t_start)

t_start = time.time()
if 'Tau31' in ds.keys():
    a = ds.Tau31.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tau31_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau31",t_end-t_start)

t_start = time.time()
if 'Tau32' in ds.keys():
    a = ds.Tau32.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tau32_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau32",t_end-t_start)

t_start = time.time()
if 'Tau22' in ds.keys():
    a = ds.Tau22.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tau22_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau22",t_end-t_start)

t_start = time.time()
if 'Tau33' in ds.keys():
    a = ds.Tau33.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
a = ds.Tau33.values
tau33_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tau33",t_end-t_start)

t_start = time.time()
if 'TauTH1' in ds.keys():
    a = ds.TauTH1.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauth1_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauth1",t_end-t_start)

t_start = time.time()
if 'TauTH2' in ds.keys():
    a = ds.TauTH2.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauth2_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauth2",t_end-t_start)

t_start = time.time()
if 'TauTH3' in ds.keys():
    a = ds.TauTH3.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauth3_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauth3",t_end-t_start)

t_start = time.time()
if 'TauQv1' in ds.keys():
    a = ds.TauQv1.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauqv1_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauqv1",t_end-t_start)

t_start = time.time()
if 'TauQv2' in ds.keys():
    a = ds.TauQv2.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauqv2_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauqv2",t_end-t_start)

t_start = time.time()
if 'TauQv3' in ds.keys():
    a = ds.TauQv3.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tauqv3_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("tauqv3",t_end-t_start)

t_start = time.time()
if 'qv' in ds.keys():
    a = ds.qv.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
qv_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("qv",t_end-t_start)

t_start = time.time()
if 'fcond' in ds.keys():
    a = ds.fcond.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
fcond_xyavg_t = XY3D_avg_time(a,times)
t_end = time.time()
print("fcond",t_end-t_start)

### end timing
t_end = time.time()
print("3D Averages",t_end-t_start)
###

### start timing
t_start = time.time()
###

if 'tskin' in ds.keys():
    a = ds.tskin.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
tskin_xyavg_t = XY2D_avg_time(a,times)

if 'fricVel' in ds.keys():
    a = ds.fricVel.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
fricvel_xyavg_t = XY2D_avg_time(a,times)

a = ds.htFlux.values
htflux_xyavg_t = XY2D_avg_time(a,times)

if 'invOblen' in ds.keys():
    a = ds.invOblen.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
invoblen_xyavg_t = XY2D_avg_time(a,times)

epsilon = 1.e-8
oblen_xyavg_t = np.reciprocal(invoblen_xyavg_t+epsilon)

if 'qskin' in ds.keys():
    a = ds.qskin.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
a = ds.qskin.values
qskin_xyavg_t = XY2D_avg_time(a,times)

if 'qFlux' in ds.keys():
    a = ds.qFlux.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
qflux_xyavg_t = XY2D_avg_time(a,times)

if 'z0m' in ds.keys():
    a = ds.z0m.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
z0m_xyavg_t = XY2D_avg_time(a,times)

if 'z0t' in ds.keys():
    a = ds.z0t.values
else:
    a = np.full([nt,nk,nj,ni],np.nan)
z0t_xyavg_t = XY2D_avg_time(a,times)

### end timing
t_end = time.time()
print("2D Averages",t_end-t_start)
###

FETimeProfiles = {'time':times,'fricvel':fricvel_xyavg_t,          \
                  'htflux':htflux_xyavg_t,'qflux':qflux_xyavg_t,   \
                  'oblen':oblen_xyavg_t,                           \
                  'z':zpos_xyavg_t,                                \
                  'u':u_xyavg_t,'v':v_xyavg_t,'w':w_xyavg_t,       \
                  'p':p_xyavg_t,'tke':tke_xyavg_t,                 \
                  'uu':uu_xyavg_t,'uv':uv_xyavg_t,                 \
                  'uw':uw_xyavg_t,'vv':vv_xyavg_t,                 \
                  'vw':vw_xyavg_t,'ww':ww_xyavg_t,                 \
                  'tau11':tau11_xyavg_t,                           \
                  'tau21':tau21_xyavg_t,                           \
                  'tau31':tau31_xyavg_t,                           \
                  'tau22':tau22_xyavg_t,                           \
                  'tau32':tau32_xyavg_t,                           \
                  'tau33':tau33_xyavg_t,                           \
                  'uuu':uuu_xyavg_t,                               \
                  'vvv':vvv_xyavg_t,                               \
                  'www':www_xyavg_t,                               \
                  'th':th_xyavg_t,                                 \
                  'tth':tth_xyavg_t,                               \
                  'ttth':ttth_xyavg_t,                             \
                  'uth':uth_xyavg_t,                               \
                  'vth':vth_xyavg_t,                               \
                  'wth':wth_xyavg_t,                               \
                  'uth_s':tauth1_xyavg_t,                          \
                  'vth_s':tauth2_xyavg_t,                          \
                  'wth_s':tauth3_xyavg_t,                          \
                  'qv':qv_xyavg_t,                                 \
                  'qqv':qqv_xyavg_t,                               \
                  'qqqv':qqqv_xyavg_t,                             \
                  'uqv':uqv_xyavg_t,                               \
                  'vqv':vqv_xyavg_t,                               \
                  'wqv':wqv_xyavg_t,                               \
                  'uqv_s':tauqv1_xyavg_t,                          \
                  'vqv_s':tauqv2_xyavg_t,                          \
                  'wqv_s':tauqv3_xyavg_t,                          \
                  'tskin':tskin_xyavg_t,                           \
                  'qskin':qskin_xyavg_t,                           \
                  'z0m':z0m_xyavg_t,                               \
                  'z0t':z0t_xyavg_t}


#filename = "FETimeProfilesSAS.pkl"
outfile = path_out+"/"+filename

with open(outfile,"wb") as f:
    pickle.dump(FETimeProfiles,f)
