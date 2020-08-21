#!/usr/bin/env python
# coding: utf-8
#
#   LESAnalysisTools.py
#
#   Written by: Branko Kosovic NCAR 2020-08-20
#
#   This script includes functions used to compute statistics from
#   LES output
#
import numpy as np
from numba import jit,vectorize,float32,float64

#####
@jit(nopython=True)
def XY2D_avg(var):
    var_avg = np.nanmean(var)
    return(var_avg)

#####
#@jit(forceobj=True)
def XY3D_avg(var):
    var_avg = np.nanmean(var,axis=(1,2))
    return(var_avg)

#####
@jit(nopython=True)
def XY2D_cov(vara, varb):

    vara_avg = XY2D_avg(vara)
    varb_avg = XY2D_avg(varb)
    ab_cov = np.empty((0))

    for k in np.arange(np.size(vara_avg)):
        vara_fluct = vara - vara_avg
        varb_fluct = varb - varb_avg
        ab = np.nanmean(vara_fluct*varb_fluct)
        ab_cov = np.append(ab_cov,ab)

    return(ab_cov)

#####
#@jit(forceobj=True)
def XY3D_cov(vara, varb):

    vara_avg = XY3D_avg(vara)
    varb_avg = XY3D_avg(varb)
    ab_cov = np.empty((0))

    for k in np.arange(np.size(vara_avg)):
        vara_fluct = vara[k,:,:] - vara_avg[k]
        varb_fluct = varb[k,:,:] - varb_avg[k]
        ab = np.squeeze(np.nanmean(vara_fluct*varb_fluct))
        ab_cov = np.append(ab_cov,ab)
 
    return(ab_cov)

#####
#@jit(forceobj=True)
def XY3D_thirdmoment(vara, varb, varc):

    vara_avg = XY3D_avg(vara)
    varb_avg = XY3D_avg(varb)
    varc_avg = XY3D_avg(varc)
    abc_thirdmoment = np.empty((0))

    for k in np.arange(np.size(vara_avg)):
        vara_fluct = vara[k,:,:] - vara_avg[k]
        varb_fluct = varb[k,:,:] - varb_avg[k]
        varc_fluct = varc[k,:,:] - varc_avg[k]
        abc = np.squeeze(np.nanmean(vara_fluct*varb_fluct*varc_fluct))
        abc_thirdmoment = np.append(abc_thirdmoment,abc)

    return(abc_thirdmoment)

#####
@jit(nopython=True)
def XY2D_avg_time(var,times):

    a_avg_t = np.empty((0))

    for time in times:
        a_avg_t=np.append(a_avg_t,XY2D_avg(var[time,:,:]))

    return(a_avg_t)

#####
#@jit(forceobj=True)
def XY3D_avg_time(var,times):

    a_avg_t = np.empty((0))

    nt,nk,nj,ni = np.shape(var)

    for time in times:
        a_avg_t=np.append(a_avg_t,XY3D_avg(var[time,:,:,:]))

    a_avg_t = np.reshape(a_avg_t,(nt,nk))

    return(a_avg_t)

#####
@jit(nopython=True)
def XY2D_cov_time(vara,varb,times):

    ab_cov_t = np.empty((0))

    for time in times:
        ab_cov_t=np.append(ab_cov_t, \
                           XY2D_cov(vara[time,:,:],varb[time,:,:]))

    return(ab_cov_t)

#####
def XY3D_cov_time(vara,varb,times):

    ab_cov_t = np.empty((0))

    nt,nk,nj,ni = np.shape(vara)

    for time in times:
        ab_cov_t=np.append(ab_cov_t, \
                           XY3D_cov(vara[time,:,:,:],varb[time,:,:,:]))

    ab_cov_t = np.reshape(ab_cov_t,(nt,nk))

    return(ab_cov_t)

#####
def XY3D_thirdmoment_time(vara,varb,varc,times):

    abc_thirdmoment_t = np.empty((0))

    nt,nk,nj,ni = np.shape(vara)

    for time in times:
        abc_thirdmoment_t=np.append(abc_thirdmoment_t, \
          XY3D_thirdmoment(vara[time,:,:,:],varb[time,:,:,:],varc[time,:,:,:]))

    abc_thirdmoment_t = np.reshape(abc_thirdmoment_t,(nt,nk))

    return(abc_thirdmoment_t)

#####
@jit(nopython=True)
def avg_x(a):

    b = 0.5*(a[:,:,:,:-1]+a[:,:,:,1:])

    return(b)

#####
@jit(nopython=True)
def avg_y(a):

    b = 0.5*(a[:,:,:-1,:]+a[:,:,1:,:])

    return(b)

#####
@jit(nopython=True)
def avg_z(a):

    b = 0.5*(a[:,:-1,:,:]+a[:,1:,:,:])

    return(b)
