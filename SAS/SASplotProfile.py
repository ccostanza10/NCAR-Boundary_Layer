#!/usr/bin/env python
# coding: utf-8
#
#   SASplotProfile.py
#
#   Written by: Branko Kosovic NCAR 2020-08-20
#
#   This script plots ABL profiles of mean and turbulent quantities
#   computed from LES and SCM output
#   
#
import os, sys
import numpy as np
import xarray as xr
import time
import pickle
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
from netCDF4 import Dataset

show_fig = False

case_selector=1

if case_selector is 1:
    path_work      = "/glade/work/branko/NCAR-BL/"
    path_analysis  = "ANALYSIS"
    path_case      = "SAS"
    basefilename   = "TimeProfilesSAS"

filenames=['FE','NCAR']
extensions=['LFix0','']
    
path_infile  = path_work+path_analysis+"/"+path_case
path_outfile = path_work+path_analysis+"/"+path_case+"/"+extensions[0]

dictionaries=[]
timeseries=[]
averages=[]

if not os.path.exists(path_outfile):
    os.makedirs(path_outfile)

plotvariables=['th','qv','wqv','wqvs','wth','wths',                          \
               'uw','uws','vw','vws',                                        \
               'uu','uus','vv','vvs','ww','wws','tth',                       \
               'tke','tkes','qqv','ustar','z']

for fname,extension in zip(filenames,extensions):
    filename = fname+basefilename
    infile  = path_infile+"/"+filename+extension+".pkl"
    with open(infile, "rb") as f:
        if fname == 'FE':
             FE         = pickle.load(f)
             FE.update({'zu':FE['z']})
             FE.update({'zw':FE['z']})
             dictionaries.append(FE)
             FEtme      = dict.fromkeys(plotvariables,0)
             timeseries.append(FEtme)
             FEavg      = dict.fromkeys(plotvariables,0)
             averages.append(FEavg)
        if fname == 'NCAR':
             NCAR       = pickle.load(f)
             NCAR['qv'] = NCAR['qv']*1000.
             NCAR['qqv']= NCAR['qqv']*1.e6
             NCAR['wqv']= NCAR['wqv']*1000.
             NCAR.update({'z':0})
             NCAR['z']  = NCAR['zu']
             dictionaries.append(NCAR)
             NCARtme    = dict.fromkeys(plotvariables,np.empty(0))
             timeseries.append(NCARtme)
             NCARavg    = dict.fromkeys(plotvariables,np.empty(0))
             averages.append(NCARavg)
        if fname == 'WRF':
             WRF        = pickle.load(f)
             WRF['th']  = WRF['th']+300.
             dictionaries.append(WRF)
             WRFtme     = dict.fromkeys(plotvariables,np.empty(0))
             timeseries.append(WRFtme)
             WRFavg    = dict.fromkeys(plotvariables,np.empty(0))
             averages.append(WRFavg)


plotvarlocs=['thz','qvz','wqvz','wqvsz','wthz','wthsz',                      \
             'uwz','uwsz','vwz','vwsz',                                      \
             'uuz','uusz','vvz','vvsz','wwz','wwsz','tthz',                  \
             'tkez','tkesz','qqvz','ustarz','zz']

varlocs=['zu','zu','zw','zw','zw','zw',                                      \
         'zw','zw','zw','zw',                                                \
         'zu','zu','zu','zu','zw','zw','zu',                                 \
         'zw','zw','zu','zw','zu']

for average in averages:
    for plotvarloc in plotvarlocs:
        average.update({plotvarloc:np.empty(0)})

for tseries in timeseries:
    for plotvarloc in plotvarlocs:
        tseries.update({plotvarloc:np.empty(0)})

for tseries,infile in zip(timeseries,dictionaries):
    for plotvarloc,varloc in zip(plotvarlocs,varlocs):
        tseries[plotvarloc] = infile[varloc]


####

for tseries,infile in zip(timeseries,dictionaries):
    tseries['th']  =infile['th']
    tseries['qv']  =infile['qv']
    tseries['wth'] =infile['wth']+infile['wth_s']
    tseries['wths']=infile['wth_s']
    tseries['wqv'] =infile['wqv']
    tseries['wqvs']=infile['wqv_s']
    tseries['uw']  =infile['uw']+infile['tau31']
    tseries['uws'] =infile['tau31']
    tseries['vw']  =infile['vw']+infile['tau32']
    tseries['vws'] =infile['tau32']
    tseries['uu']  =infile['uu']+infile['tau11']+2./3.*infile['tke']
    tseries['uus'] =infile['tau11']+2./3.*infile['tke']
    tseries['vv']  =infile['vv']+infile['tau22']+2./3.*infile['tke']
    tseries['vvs'] =infile['tau22']+2./3.*infile['tke']
    tseries['ww']  =infile['ww']+infile['tau33']+2./3.*infile['tke']
    tseries['wws'] =infile['tau33']+2./3.*infile['tke']
    tseries['tth'] =infile['tth']
    tseries['tke'] =0.5*(infile['uu']+infile['vv']+infile['ww'])+infile['tke']
    tseries['tkes']=infile['tke']
    tseries['qqv'] =infile['qqv']
    tseries['wqv'] =infile['wqv']
    tseries['ustar']=np.power(np.power(infile['uw']+infile['tau31'],2) + \
                              np.power(infile['vw']+infile['tau32'],2),0.25)
    tseries['z']    =infile['z']

for avg in averages:
    for pltvars in plotvariables:
        avg[pltvars] = np.empty((0))

tss = '13'
tes = '15'
tfs = '15'

ts = int((int(tss)-5)*60/5)
te = int((int(tes)-5)*60/5)
tf = int((int(tfs)-5)*60/5)
print("ts, te, tf = ",ts,te,tf)


for avgs,tseries in zip(averages,timeseries):
    nt,nk = np.shape(tseries['wth']) 
    print(nt,nk)
    for pltvars in plotvariables:
        for k in range(nk):
            avgs[pltvars] = np.append(avgs[pltvars], \
                               np.mean(tseries[pltvars][ts:te,k]))
    for plotvarloc in plotvarlocs:
        for k in range(nk):
            avgs[plotvarloc] = np.append(avgs[plotvarloc], \
                               np.mean(tseries[plotvarloc][ts:te,k]))


# Sensible Heat Flux
linetype=['solid','dashed']
color   =['red','blue']
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['wth'],avgs['wthz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['wths'],avgs['wthsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel='Sensible Heat Flux $[K m s^{-1}]$',ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/SensibleHeatFlux_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Water Vapor Mixing Ratio Flux
linetype=['solid','dashed']
color   =['red','blue']
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['wqv'],avgs['wqvz'], \
            label=infile+' res.',linestyle='solid',color=clr)
    ax.plot(avgs['wqvs'],avgs['wqvsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel='Water Vapor Mix. Rat. Flux $[K m s^{-1}]$',ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/WaterVaproMixingRatioFlux_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Turbulent Stress <uw>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['uw'],avgs['uwz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['uws'],avgs['uwsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'Turbulent Stress $\langle uw \rangle [m^{2} s^{2}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentStressUW_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Turbulent Stress <vw>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['vw'],avgs['vwz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['vws'],avgs['vwsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'Turbulent Stress $\langle vw \rangle [m^{2} s^{2}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentStressVW_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Friction velocity
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['ustar'],avgs['ustarz'], \
            label=infile,linestyle='solid',color=clr)
ax.set(xlabel=r'$u_{*} [m s^{-1}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/Ustar_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Turbulent stress <uu>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['uu'],avgs['uuz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['uus'],avgs['uusz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'Turbulent Stress $\langle uu \rangle [m s^{-1}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentStressUU_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Turbulent stress <vv>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['vv'],avgs['vvz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['vvs'],avgs['vvsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'Turbulent Stress $\langle vv \rangle [m s^{-1}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentStressVV_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Turbulent stress <ww>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['ww'],avgs['wwz'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['wws'],avgs['wwsz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'Turbulent Stress $\langle ww \rangle [m s^{-1}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentStressWW_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Potential temperatue variance <tth>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['tth'],avgs['tthz'], \
            label=infile,linestyle='solid',color=clr)
ax.set(xlabel=r'Pot. Temp. Var. $\langle \theta \theta \rangle [m s^{-1}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentVarianceTT_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# TKE
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['tke'],avgs['tkez'], \
            label=infile+' tot.',linestyle='solid',color=clr)
    ax.plot(avgs['tkes'],avgs['tkesz'], \
            label=infile+' sub.',linestyle='dashed',color=clr)
ax.set(xlabel=r'TKE $[m^{2} s^{-2}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TKE_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


# Water vapor mixing ratio variance <qv qv>
fig, ax = plt.subplots()
for avgs,infile,clr in zip(averages,filenames,color):
    ax.plot(avgs['qqv'],avgs['qqvz'], \
            label=infile,linestyle='solid',color=clr)
ax.set(xlabel=r'Water Vap. Mix. Rat. Var. $\langle q_{v} q_{v} \rangle [kg^{2} kb^{-2}]$', \
       ylabel='z $[m]$')
ax.legend()

fig.savefig(path_outfile+'/TurbulentVarianceQQ_'+tss+'h-'+tes+'h.png')
if show_fig == True:
    plt.show()
plt.close()


################## Timeseries

tme=np.linspace(5,int(tfs),tf)


# Surface friction velocity
fig, ax = plt.subplots()
for infile,lbl in zip(dictionaries,filenames):
    ax.plot(tme,infile['fricvel'][:tf],label=lbl)
ax.set_xlim([4,int(tfs)])
ax.set_ylim([0,0.6])
ax.set(xlabel='Time $UTC$',ylabel='Friction velocity $[m s^{-1}]$')
ax.legend()

fig.savefig(path_outfile+'/SASfricvel.png')
if show_fig == True:
    plt.show()
plt.close()


# FE BL height based on the pot. temp.
FEth = FE['th']
FEdth = FEth[:,1:]-FEth[:,:-1]
maxind = np.argmax(FEdth,axis=1)
FEpblhth = FE['z'][0,maxind]
FE.update({'pblhth':FEpblhth})


# FE BL height based on the water vapor mixing ratio
FEqv = FE['qv']
FEdqv = np.absolute(FEqv[:,10:]-FEqv[:,9:-1])
maxind = np.argmax(FEdqv,axis=1)
FEpblhq = FE['z'][0,maxind]
FE.update({'pblhq':FEpblhq})


# Boundar layer height based on the potential temperatue gradient
fig, ax = plt.subplots()
for infile,lbl in zip(dictionaries,filenames):
    ax.plot(tme,infile['pblhth'][:tf],label=lbl)
ax.set_xlim([4,int(tfs)])
ax.set_ylim([0,3000.0])
ax.set(xlabel='Time $UTC$', \
    ylabel=r'Boundary Layer Height - $\frac{\partial \Theta}{\partial z} [m]$')
ax.legend()

fig.savefig(path_outfile+'/pblhth.png')
if show_fig == True:
    plt.show()
plt.close()


# Boundar layer height based on the water vapor mixing ratio gradient
fig, ax = plt.subplots()
for infile,lbl in zip(dictionaries,filenames):
    ax.plot(tme,infile['pblhq'][:tf],label=lbl)
ax.set_xlim([4,int(tfs)])
ax.set_ylim([0,3000.0])
ax.set(xlabel='Time $UTC$', \
    ylabel=r'Boundary Layer Height - $\frac{\partial q_{v}}{\partial z} [m]$')
ax.legend()

plt.savefig(path_outfile+'/pblhq.png')
if show_fig == True:
    plt.show()
plt.close()
