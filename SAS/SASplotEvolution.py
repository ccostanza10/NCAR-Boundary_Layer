#!/usr/bin/env python
# coding: utf-8
#
#   SASplotEvolution.py
#
#   Written by: Branko Kosovic NCAR 2020-08-20
#
#   This script plots ABL time slices of mean and turbulent quantities
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
import matplotlib.ticker as tick
plt.style.use('seaborn-white')
from netCDF4 import Dataset

def round_up(n, decimals=0): 
    multiplier = 10 ** decimals 
    return math.ceil(n * multiplier) / multiplier

def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier

show_fig = False

case_selector=1

if case_selector is 1:
    path_work    = "/glade/work/branko/NCAR-BL/"
    path_analysis="ANALYSIS"
    path_case    ="SAS"
    basefilename ="TimeProfilesSAS"
    extension    ="LFix0"

path_infile  = path_work+path_analysis+"/"+path_case
path_outfile = path_work+path_analysis+"/"+path_case+"/"+extension+"/"

filenames =['FE','NCAR']
extensions=['LFix0','']

dictionaries=[]

if not os.path.exists(path_outfile):
    os.makedirs(path_outfile)

for fname,extension in zip(filenames,extensions):
    filename = fname+basefilename 
    infile  = path_infile+"/"+filename+extension+".pkl"
    with open(infile, "rb") as f:
        if fname == 'FE':
             FE         = pickle.load(f) 
             dictionaries.append(FE)
        if fname == 'NCAR':
             NCAR       = pickle.load(f) 
             NCAR['qv'] = NCAR['qv']*1000.
             NCAR['qqv']= NCAR['qqv']*1.e6
             NCAR['wqv']= NCAR['wqv']*1000.
             dictionaries.append(NCAR)
        if fname == 'WRF':
             WRF        = pickle.load(f) 
             WRF['th']  = WRF['th']+300.
             dictionaries.append(WRF)

varnames=['th','qv','u','v','wth','wth_s','wqv','wqv_s',                     \
          'uw','tau31','vw','tau32','uu','tau11','vv','tau22','ww','tau33',  \
          'tke','tth','qqv']

turbvars=['wth','wth_s','wqv','wqv_s',                                       \
          'uw','tau31','vw','tau32','uu','tau11','vv','tau22','ww','tau33',  \
          'tke','tth','qqv']

vmin=dict.fromkeys(varnames,0)
vmax=dict.fromkeys(varnames,0)
vdif=dict.fromkeys(varnames,0)


for var in varnames:
    amin=+999.
    amax=-999.
    for infile in dictionaries:
        amin = np.min(np.array([amin,np.min(infile[var])]))
        amax = np.max(np.array([amax,np.max(infile[var])]))
        if amin < vmin[var]:
            vmin[var] = amin
        if amax > vmax[var]:
            vmax[var] = amax


for var in turbvars:
    vmin[var] = round_down(vmin[var],2)
    vmax[var] = round_up(vmax[var],2)

vmin.update({'wspd':0})

vmin['th']   =295.
vmin['qv']   =0.
vmin['wspd'] =0.
print("Min",vmin)
print(' ')

vmax['th']   =round_up(vmax['th'],0)
vmax['qv']   =round_up(vmax['qv'],3)
vmax.update({'wspd':math.sqrt(vmax['u']*vmax['u']+vmax['v']*vmax['v'])})
vmax['wspd'] =round_up(vmax['wspd'],0)
print("Max",vmax)
print(' ')

for var in varnames:
    vdif[var]=vmax[var]-vmin[var]

vdif.update({'wspd':vmax['wspd']-vmin['wspd']})

for var in turbvars:
    vdif[var] = int(50*vdif[var])

for key in vdif.keys():
    vdif[key] = int(vdif[key])

print("Dif",vdif)
print(' ')



for outfile,infile in zip(filenames,dictionaries):

    ###
    th =infile['th']
    print(np.shape(th))

    nt,nk = np.shape(th)

    tme = np.linspace(5,5+nt*5/60,nt)
    t=np.linspace(5,5+nt*5/60,nt)
    z=np.linspace(0,nk*15,nk)
    t,z = np.meshgrid(t,z)
    t = np.transpose(t)
    z = np.transpose(z)

    #  Potential Temperature
    v = np.linspace(vmin['th'],vmax['th'],2*vdif['th']+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,th,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label('Potential Temperature $[K]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'th.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Wind Speed
    u = infile['u']
    v = infile['v']
    wspd = np.sqrt(u*u+v*v)

    v = np.linspace(vmin['wspd'],vmax['wspd'],6*vdif['wspd']+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,wspd,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label('Wind Speed $[m s^{-1}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'wspd.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Water Vapor Mixing Ratio
    qv = infile['qv']

    v = np.linspace(vmin['qv'],vmax['qv'],2*vdif['qv']+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,qv,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label('Specific Humidity $[kg kg^{-1}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'qv.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Turbulent Heat Flux
    wth = infile['wth']+infile['wth_s']

    v = np.linspace(vmin['wth']+vmin['wth_s'],      \
                    vmax['wth']+vmax['wth_s'],      \
                    2*(vdif['wth']+vdif['wth_s'])+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,wth,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label('Sensible Heat Flux $[K m s^{-1}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'wth.png')
    if show_fig == True:
        plt.show()
    plt.close()
    
    # Turbulent Stress <uw>
    uw = infile['uw'] + infile['tau31']

    v = np.linspace(vmin['uw']+vmin['tau31'],      \
                    vmax['uw']+vmax['tau31'],      \
                    2*(vdif['uw']+vdif['tau31'])+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot=plt.contourf(t,z,uw,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Turbulent Shear Stress $ \langle uw \rangle  [m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'uw.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Turbulent Stress <vw>
    vw = infile['vw'] + infile['tau32']

    v = np.linspace(vmin['vw']+vmin['tau32'],      \
                    vmax['vw']+vmax['tau32'],      \
                    2*(vdif['vw']+vdif['tau32'])+1, endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,vw,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Turbulent Shear Stress $\langle vw \rangle [m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'vw.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Streamwise Velocity Variance <uu>
    uu = infile['uu'] + infile['tau11'] + 2./3.*infile['tke']

    v = np.linspace(vmin['uu']+vmin['tau11']+2./3.*vmin['tke'],        \
                    vmax['uu']+vmax['tau11']+2./3.*vmax['tke'],        \
                    2*(vdif['uu']+vdif['tau11']+vdif['tke'])+1,  \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,uu,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Turbulent Stress $\langle uu \rangle [m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'uu.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Cross-stream Velocity Variance <vv>
    vv = infile['vv'] + infile['tau22'] + 2./3.*infile['tke']

    v = np.linspace(vmin['vv']+vmin['tau22']+2./3.*vmin['tke'],        \
                    vmax['vv']+vmax['tau22']+2./3.*vmax['tke'],        \
                    2*(vdif['vv']+vdif['tau22']+vdif['tke'])+1,  \
                    endpoint=True) 
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,vv,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Turbulent Stress $\langle vv \rangle [m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'vv.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Vertical Velocity Variance <vv>
    ww = infile['ww'] + infile['tau33'] + 2./3.*infile['tke']

    v = np.linspace(vmin['ww']+vmin['tau33']+2./3.*vmin['tke'],        \
                    vmax['ww']+vmax['tau33']+2./3.*vmax['tke'],        \
                    2*(vdif['ww']+vdif['tau33']+vdif['tke'])+1,        \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,ww,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Turbulent Stress $\langle ww \rangle [m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'ww.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # TKE
    tke = infile['uu'] + infile['tau11'] +                       \
          infile['vv'] + infile['tau22'] +                       \
          infile['ww'] + infile['tau33'] + infile['tke']

    v = np.linspace(0.5*(vmin['uu']+vmin['tau11'] +              \
                         vmin['vv']+vmin['tau22'] +              \
                         vmin['ww']+vmin['tau33'])+vmin['tke'],  \
                    0.5*(vmax['uu']+vmax['tau11'] +              \
                         vmax['vv']+vmax['tau22'] +              \
                         vmax['ww']+vmax['tau33'])+vmax['tke'],  \
                    (vdif['uu']+vdif['tau11']+                   \
                     vdif['vv']+vdif['tau22']+                   \
                     vdif['ww']+vdif['tau33']+2*vdif['tke']),    \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,tke,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'TKE $[m^{2} s^{-2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'tke.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Potential Temperature Variance
    tth = infile['tth'] 

    v = np.linspace(vmin['tth'],vmax['tth'],2*vdif['tth']+1,  \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,tth,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Pot. Temperature variance  $\langle \theta \theta \rangle [K^{2}]$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'tth.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Water Vapor Mixing Ratio Variance
    qqv = infile['qqv']

    v = np.linspace(vmin['qqv'],vmax['qqv'],2*vdif['qqv']+1,  \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,qqv,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Water Vap. Mix. Rat. var.  $\langle q_{v} q_{v} \rangle$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'qqv.png')
    if show_fig == True:
        plt.show()
    plt.close()

    # Turbulent Flux of Water Vapor Mixing Ratio 
    wqv = infile['wqv']

    v = np.linspace(vmin['wqv'],vmax['wqv'],2*vdif['wqv']+1,  \
                    endpoint=True)
    fig, ax = plt.subplots()
    plot = ax.contourf(t,z,wqv,v,cmap='gist_earth')
    cbar = plt.colorbar(plot)
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    cbar.set_label(r'Water Vap. Mix. Rat. flux  $\langle w q_{v} \rangle$')
    ax.set(xlabel='Time $[hr]$',ylabel='z $[m]$')

    fig.savefig(path_outfile+outfile+'wqv.png')
    if show_fig == True:
        plt.show()
    plt.close()
