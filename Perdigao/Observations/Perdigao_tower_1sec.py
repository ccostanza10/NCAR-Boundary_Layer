#!/opt/local/anaconda3/bin/python
#
# Created August 2021
# Carol Costanza
# Modified August 18-20/2021 by Jeremy Sauer

# NCAR ISFS Tower Data converted to Xarray for Perdigao
# specifically May 23, 2017 12 UTC - 23:59 UTC

import netCDF4
import datetime as dt
import numpy as np
import pandas as pd
import glob
from matplotlib.cbook import flatten
import xarray as xr
np.set_printoptions(threshold=np.inf)

state_vars = ['u','v','w','T','P']
#state_vars = ['u']
verboseLog=False
hrStart=12
hrStop=24
hours = list(range(hrStart,hrStop+1))
all_ncpath = '/glade/scratch/jsauer/FastEddy/BL_REINV/PERDIGAO/OBS/HIGHRATE/LATTER/isfs_qc_tiltcor_hr_20170523_*.nc'
hourly_files = sorted(glob.glob(all_ncpath))
if verboseLog:
  print("hourly_files")
  print(hourly_files)
ncpath = hourly_files[0]
tower_file = netCDF4.Dataset(ncpath)

# get the tower names
tower_file_vars = tower_file.variables.keys()
alt_vars = []; towers = []
for var in tower_file_vars:
  # assumes all towers have a altitude variable
  if var.find('altitude') != -1:
    alt_vars.append(var)
for alt in alt_vars:
  towers.append(alt[9:])
if verboseLog:
  print(towers)


# get the lat/lon/alt
lats = []; lons = []; alts = []
for tower in towers:
  lats.append(float(tower_file.variables['latitude_'+tower][:]))
  lons.append(float(tower_file.variables['longitude_'+tower][:]))
  alts.append(float(tower_file.variables['altitude_'+tower][:]))

# get all the possible heights and variables names based on the requested
# variables in the hard coded 'state_vars'
heights = []
all_vars = []
for var in tower_file_vars:
  for met_var in state_vars:
    if var.find(met_var+'_') != -1:
      heights.append(int(var[2:].partition("_")[0][:-1]))
      all_vars.append(var)
all_hgts = sorted(set(flatten(heights)))
if verboseLog:
  print("all_heights")
  print(all_hgts)
  print("all_vars")
  print(all_vars)

############## setup the final array that will be used for the final netCDF file
final_array = np.zeros((len(state_vars),3600*(len(hours)-1),len(towers),len(all_hgts)))


#combine all the targeted hourly files in big arrays for each variable in all_vars
all_data = np.zeros((len(all_vars),3600*(len(hours)-1)))
time = []
# get data for each variable in all_vars
for j,met_var in enumerate(all_vars):
  icnt=-1
  for i,hr in enumerate(hourly_files):     #Reminder: hourly_files is the glob'd list of targeted files from above
    readThis=False
    fileHour=hourly_files[i].partition('_20170523_')[-1].strip('.nc')  #filename pattern is 'blah_20170523_fileHour.nc'
    if int(fileHour) >= hours[0]:
       if int(fileHour) < hours[-1]:
          ifile=i
          readThis=True
          icnt+=1
    if readThis:
      hr_tower_data = netCDF4.Dataset(hourly_files[ifile])
      if j==0:  #append to the time
        secs_data = hr_tower_data.variables['time'][:]
        tot_secs_data = [float(x)+(3600.0*int(fileHour)) for x in secs_data]
        time.append(tot_secs_data)
      ####  append to the j-variable
      var_data = np.array(hr_tower_data.variables[met_var][:])
      # remove all the FillValues to set to nan
      var_data[var_data == np.float32(1.0e37)] = np.nan
      var_dims = var_data.ndim
      if var_dims == 2: # for 20 hertz data
        # convert 20 hertz data to 1 sec via simple mean
        var_data_final=np.squeeze(np.nanmean(var_data,axis=1))
      else: # for 1 second data
        var_data_final = var_data
      if icnt == 0:     ###If the first set (file) of times, copy (create) a vector rather than concatenate (append)
         tmparray=var_data_final
      elif icnt > 0:
         tmparray=np.concatenate((tmparray,var_data_final),axis=0)
  all_data[j,:]=tmparray
time = np.ravel(np.asarray(time,dtype=np.float32()))
if verboseLog:
  print("All times completed.")
  print(time)
  print(time.shape)
  print("All data shape...")
  print(all_data.shape)

# reshaping the data array such that the df has variable names as columns
# and time as the row
# set all_data dims to (len(all_vars),3600*len(hours))
reshape_all_data = all_data.reshape(len(all_vars),-1)
if verboseLog:
  print("reshape_all_data completed.")
  print(reshape_all_data.shape)
data = np.transpose(reshape_all_data)
if verboseLog:
  print("transpose to data completed.")
  print(data.shape)
df = pd.DataFrame(data,columns=all_vars)

# start looping to add data to final arrays (state_var,time,tower,height)
for tower_ind,tower in enumerate(towers):
  tower_vars = [] 
  for met_var in all_vars:
    # get all possible variable names for each tower
    if met_var.find('_'+str(tower)) != -1:
      tower_vars.append(met_var)
  for v,state_var in enumerate(state_vars):
    hgts = []
    for tower_var in tower_vars:
      # get all state_var variables for each tower
      if tower_var.find(state_var+'_') != -1:
        hgts.append(int(tower_var[2:].partition('_')[0][:-1])) #found the ('state_var'_'hgt'm_'tower') variable so add hgt to list of hgts
    hgts = sorted(hgts)
    for t,sec in enumerate(time):
      # for given tower and state var, no obs exist, thus nans
      if len(hgts) == 0:
        data_sec = np.zeros((len(all_hgts)))
        data_sec[:] = np.nan
        final_array[v,t,tower_ind,:] = data_sec
      else:
        data_sec = np.zeros((len(hgts)))
        # for given tower, state var, and hgt, pull from the df
        for i,hgt in enumerate(hgts):
          data_sec[i] = df.at[t,state_var+'_'+str(hgt)+'m_'+tower]
        # setup xarray such that heights can correct populated with the
        # all_hgts array with missing nans
        ds = xr.DataArray(data_sec,[('heights',hgts)])
        ds_reindex = ds.reindex(heights=all_hgts)
        final_array[v,t,tower_ind,:] = ds_reindex.values
#finished 

#############################################################################
# create final netCDF file
nc_file = netCDF4.Dataset('Perdigao_tower_20170523_1sec.nc','w',\
                                      format='NETCDF4')
nc_file.title = "Perdigao 2020/05/23 ISFS Tower (u,v,w,P,T)" +\
                            " every 1 second 12 - 00 UTC"
nc_file.createDimension('time', len(time))
nc_file.createDimension('height', len(all_hgts))
nc_file.createDimension('tower', len(towers))
tower_names = nc_file.createVariable('tower',str,'tower')
all_instruments = np.array(towers, dtype=object)
tower_names[:] = all_instruments
latitudes = nc_file.createVariable('latitude','f4','tower')
latitudes[:] = lats
latitudes.units = 'degrees north'
longitudes = nc_file.createVariable('longitude','f4','tower')
longitudes[:] = lons
longitudes.units = 'degrees east'
alt = nc_file.createVariable('altitude','f4','tower')
alt[:] = alts
alt.units = 'meters'
hgt_var = nc_file.createVariable('height','f4','height')
hgt_var[:] = all_hgts
hgt_var.units = 'meters above ground level (need to add alt for asl)'
base_time = nc_file.createVariable('base_time','i')
base_time[:] = (dt.datetime.strptime('23/05/2017 00:00:00','%d/%m/%Y %H:%M:%S')\
                   -dt.datetime(1970,1,1)).total_seconds()
base_time.units = 'seconds since 1970-01-01 00:00:00'
epoch_secs = nc_file.createVariable('time','f4',('time'))
epoch_secs[:] = np.ravel(time)
epoch_secs.units = 'seconds since 2017-05-23 00:00:00'
if len(state_vars)>0:	
 east_vel = nc_file.createVariable('u','f4',('time','tower','height'),\
                                  fill_value='NaN')
 east_vel[:] = final_array[state_vars.index('u')]
 east_vel.units = 'm/s'
 east_vel.standard_name = 'eastern velocity'

if len(state_vars)>1: #False:
 north_vel = nc_file.createVariable('v','f4',('time','tower','height'),\
                                   fill_value='NaN')
 north_vel[:] = final_array[state_vars.index('v')]
 north_vel.units = 'm/s'
 north_vel.standard_name = 'northern velocity'
if len(state_vars)>2: #False:
 vert_vel = nc_file.createVariable('w','f4',('time','tower','height'),\
                                  fill_value='NaN')
 vert_vel[:] = final_array[state_vars.index('w')]
 vert_vel.units = 'm/s'
 vert_vel.standard_name = 'vertical velocity'
if len(state_vars)>3: #False:
 temps = nc_file.createVariable('T','f4',('time','tower','height'),\
                               fill_value='NaN')
 temps[:] = final_array[state_vars.index('T')]
 temps.units = 'degC'
 temps.standard_name = 'temperature'
if len(state_vars)>4: #False:
 press = nc_file.createVariable('P','f4',('time','tower','height'),\
                               fill_value='NaN')
 press[:] = final_array[state_vars.index('P')]
 press.units = 'mb'
 press.standard_name = 'pressure'

nc_file.close()
 
