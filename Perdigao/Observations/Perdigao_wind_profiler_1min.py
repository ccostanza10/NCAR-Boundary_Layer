#!/opt/local/anaconda3/bin/python
#
# Created November 2020
# Carol Costanza
# Updated January 2021
# using fewer profilers
# using xarray for 1 height dim across all profilers
#
# 5 minute grid of wind components ('u','v','w') for May 23rd, 2017 Perdigao
# Only vertically pointing instruments or surface met (not ISFS)

import numpy as np
import datetime as dt
from datetime import timezone
import netCDF4
import math
import csv
from matplotlib.cbook import flatten
from itertools import islice
import xarray as xr
np.warnings.filterwarnings('ignore')
np.set_printoptions(threshold=np.inf)

nc_path = '/scr/tmp/costanza/Perdigao_wind/'

# set time step for 12 UTC every 1 minutes for 12 hours starting at 12 UTC
str_date = '2017/05/23 12:00:00Z'
date = dt.datetime.strptime(str_date,'%Y/%m/%d %H:%M:%S%z')
steps = [date + dt.timedelta(minutes=1*x) for x in range(721)]
secs = [s.timestamp() for s in steps]
num_steps = len(secs)
tot_secs_day = (dt.datetime(2017,5,23)-dt.datetime(1970,1,1))\
                .total_seconds()
secs_day = [x - tot_secs_day for x in secs]

###########################################################################
# INITIAL VARIABLES
###########################################################################

# CLAMPS initial variables
CLAMPS_vad_file = nc_path + 'clampsdlvad1turnC1.c1.20170523.001330.cdf'
CLAMPS_lat = 39.712696; CLAMPS_lon = -7.7372923; CLAMPS_alt = 299

# NCAS Profiler initial variables
ncas_file_hr = nc_path + 'ncas_prof_1min.nc'
ncas_file = nc_path + 'man-radar-1290mhz_alvaiade_20170523_high-range-mode-15min.nc'
ncas_alt = 287

# UND Sodar RASS
UND_Sodar_file = nc_path + 'UND_Sodar_RASS_170523.txt'
UND_Sodar_lat = 39.722658; UND_Sodar_lon = -7.729475; UND_Sodar_alt = 258

# NCAR Profiler
ncar_wp_file = nc_path + 'nima.20170523.000013.143.winds_LO.nc'
ncar_wp_lat = 39.70156; ncar_wp_lon = -7.76393; ncar_wp_alt = 222

# NCAR Profiler High Rate 'w'
ncar_wp_hr_file = nc_path + 'nima.20170523.000013.143.mom.nc'

# ZephIR_z423
ZephIR_z423_file = nc_path + 'Cornell_ZephIR_z423.nc'
ZephIR_z423_lat = 39.7136775; ZephIR_z423_lon = -7.73657
ZephIR_z423_alt = 311

# ZephIR_z447
ZephIR_z447_file = nc_path + 'Cornell_ZephIR_z447.nc'
ZephIR_z447_lat = 39.70674; ZephIR_z447_lon = -7.75556
ZephIR_z447_alt = 238

# ENERCON
ENERCON_lidar_file = nc_path + 'ENERCON_lidar.sta'
ENERCON_lidar_lat = 39.7079638; ENERCON_lidar_lon = -7.7489312
ENERCON_lidar_alt = 367
ENERCON_lidar_hgts = [40,66,70,93,110,120,147,150,168,171,188,191,209,\
                      212,229,232,250]

# leosphere windcube
windcube_lidar_file = nc_path + 'windcube_lidar_20170523.sta'
windcube_lidar_lat = 39.713662; windcube_lidar_lon = -7.730430
windcube_lidar_alt = 453
windcube_lidar_hgts = [40,59,79,94,99,119,139,159,179,199]

# CU Orange windcube
CU_Orange_lidar_file = nc_path + 'WLS7-0049_2017_05_23__00_00_00.sta'
CU_Orange_lidar_lat = 39.71245; CU_Orange_lidar_lon = -7.737386
CU_Orange_lidar_alt = 296
CU_Orange_lidar_hgts = [40,60,80,100,120,140,160,180,200,220]

# CU Beehive windcube
CU_Beehive_lidar_file = nc_path + 'WLS7-0068_2017_05_23__06_35_41.sta'
CU_Beehive_lidar_lat = 39.718331; CU_Beehive_lidar_lon = -7.726405 
CU_Beehive_lidar_alt = 260
CU_Beehive_lidar_hgts = [40,60,80,100,120,140,160,180,200,220]

# get an array of all the netCDF/data files to process
nc_files = [CLAMPS_vad_file,ncas_file,ncas_file_hr,UND_Sodar_file,ncar_wp_file,\
            ncar_wp_hr_file,ZephIR_z423_file,ZephIR_z447_file,\
            ENERCON_lidar_file,windcube_lidar_file,\
            CU_Orange_lidar_file,CU_Beehive_lidar_file]
tot_files = len(nc_files)

# set up height index that spans all height values
all_heights = []
for instrument,nc_file in enumerate(nc_files):
  if nc_file == CLAMPS_vad_file:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    hgts = lidar_file.variables['height'][0:202]*1000
    altitude = CLAMPS_alt
  if nc_file == ncas_file or nc_file == ncas_file_hr:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    hgts = lidar_file.variables['altitude'][10,0:32]
    altitude =  ncas_alt
  if nc_file == UND_Sodar_file:
    hgts = []
    altitude = UND_Sodar_alt
    with open (nc_file,encoding="utf8",errors='ignore') as lidar_file:
      sodar_data = csv.reader(lidar_file,delimiter=',')
      for row in islice(sodar_data,13,None):
        hgts.append(float(row[1]))
    hgts = sorted(list(set(hgts)))
  if nc_file == ncar_wp_file:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    altitude = ncar_wp_alt 
    hgts = lidar_file.variables['heights'][0]
  if nc_file == ncar_wp_hr_file:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    altitude = ncar_wp_alt 
    hgts = lidar_file.variables['heights'][0]
  if nc_file == ZephIR_z423_file or nc_file == ZephIR_z447_file:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    if nc_file == ZephIR_z423_file:
      altitude = ZephIR_z423_alt
    if nc_file == ZephIR_z447_file:
      altitude = ZephIR_z447_alt
    hgts = lidar_file.variables['Heights'][:]
    hgts = [x for x in reversed(hgts)]
  if nc_file == ENERCON_lidar_file or nc_file == windcube_lidar_file:
    if nc_file == ENERCON_lidar_file:
      altitude = ENERCON_lidar_alt; hgts = ENERCON_lidar_hgts
    if nc_file == windcube_lidar_file:
      altitude = windcube_lidar_alt; hgts = windcube_lidar_hgts
  if nc_file == CU_Orange_lidar_file or nc_file == CU_Beehive_lidar_file:
    if nc_file == CU_Orange_lidar_file:
      altitude = CU_Orange_lidar_alt; hgts = CU_Orange_lidar_hgts
    if nc_file == CU_Beehive_lidar_file:
      altitude = CU_Beehive_lidar_alt; hgts = CU_Beehive_lidar_hgts

  hgts = [int(x) for x in hgts]
  #if nc_file != ncas_file and nc_file != ncas_file_hr:
  altitudes = [x + altitude for x in hgts]
  # ncas profiler hgts already include altitude
  #else: altitudes = hgts
  all_heights.append(altitudes)

# create the height dimension
all_hgts = sorted(set(flatten(all_heights)))
print(all_hgts)
tot_hgts = len(all_hgts)

# create the u,v,w 3D arrays that will be used for the final netCDF file
u = np.zeros((num_steps-1,tot_files,tot_hgts))
v = np.zeros((num_steps-1,tot_files,tot_hgts))
w = np.zeros((num_steps-1,tot_files,tot_hgts))

# set lat/lon/alt/instrument lists
all_latitudes = []; all_longitudes = []; all_altitudes = []
all_instruments = []
###########################################################################
# END INITIAL VARIABLES
###########################################################################

###########################################################################
# Start Loop for reading in netCDF/text files
###########################################################################
for instrument,nc_file in enumerate(nc_files):
  if nc_file == CLAMPS_vad_file:
    instrument_name = "CLAMPS"
    lidar_file = netCDF4.Dataset(nc_file,'r')
    latitude = CLAMPS_lat; longitude = CLAMPS_lon; altitude = CLAMPS_alt
    hgts = lidar_file.variables['height'][0:202]*1000
    base_time = lidar_file.variables['base_time'][:]
    secs_since_base = lidar_file.variables['time_offset'][0:202]
    intensity = lidar_file.variables['intensity'][:][0:202]
    ws = lidar_file.variables['wspd'][:][0:202]
    wd = lidar_file.variables['wdir'][:][0:202]
    lidar_secs = []
    # convert string time from netCDF to secs
    for time in secs_since_base:
      lidar_secs.append((base_time+time))

  if nc_file == ncas_file:
    instrument_name = "NCAS_WProfiler"
    lidar_file = netCDF4.Dataset(nc_file,'r')
    latitude = lidar_file.variables['latitude'][0]
    longitude = lidar_file.variables['longitude'][0]
    altitude = ncas_alt
    # since most data past index 32 is nan, only get the beginning
    # all the arrays for NCAS profiler
    # hgts are repeated for each time, so randomly get the hgts from
    # 1050th time step
    hgts = lidar_file.variables['altitude'][10,0:32]
    lidar_secs = lidar_file.variables['time'][:]
    u_vel = lidar_file.variables['eastward_wind'][:,0:32]
    v_vel = lidar_file.variables['northward_wind'][:,0:32]
    up_wind = lidar_file.variables['upward_air_velocity'][:,0:32]
    flag = lidar_file.variables['qc_flag_wind'][:,0:32]

  if nc_file == ncas_file_hr:
    instrument_name = "NCAS_WProfiler_HighRate"
    lidar_file = netCDF4.Dataset(nc_file,'r')
    latitude = lidar_file.variables['latitude'][0]
    longitude = lidar_file.variables['longitude'][0]
    altitude = ncas_alt
    # since most data past index 32 is nan, only get the beginning
    # all the arrays for NCAS profiler
    # hgts are repeated for each time, so randomly get the hgts from
    # 1050th time step
    hgts = lidar_file.variables['altitude'][1050,0:32]
    lidar_secs = lidar_file.variables['time'][:]
    u_vel = lidar_file.variables['eastward_wind'][:,0:32]
    v_vel = lidar_file.variables['northward_wind'][:,0:32]
    up_wind = lidar_file.variables['upward_air_velocity'][:,0:32]
    flag = lidar_file.variables['qc_flag_wind'][:,0:32]

  if nc_file == UND_Sodar_file:
    instrument_name = "UND_Sodar"
    hgts = []; lidar_secs = []; u_vel = []; v_vel = []; w_vel = []
    latitude = UND_Sodar_lat; longitude = UND_Sodar_lon
    altitude = UND_Sodar_alt
    with open (nc_file,encoding="utf8",errors='ignore') as lidar_file:
      sodar_data = csv.reader(lidar_file,delimiter=',')
      count = 0
      u_hgt = []; v_hgt = []; w_hgt = []
      # Read through all rows of text file and start at row 14
      # to skip the header
      for row in islice(sodar_data,13,None):
        str_time = row[0]
        lidar_secs.append((dt.datetime.strptime(str_time,'%Y-%m-%d %H:%M:%S'\
                        )-dt.datetime(1970,1,1)).total_seconds())
        hgts.append(float(row[1]))
        # set all wind component values of '*' to nan
        if row[4] == '*': u_hgt.append(np.nan)
        else: u_hgt.append(float(row[4]))
        if row[5] == '*': v_hgt.append(np.nan)
        else: v_hgt.append(float(row[5]))
        if row[6] == '*': w_hgt.append(np.nan)
        else: w_hgt.append(float(row[6]))
        count += 1
        # every 39 rows append the whole list of wind components
        # need to know there are 39 heights in the file
        if count == 38:
          u_vel.append(u_hgt)
          v_vel.append(v_hgt)
          w_vel.append(w_hgt)
          # reset the lists and count to clear out the data
          u_hgt = []; v_hgt = []; w_hgt = []
          count = 0
    # take out all repeated values in hgts and time, and then sort
    hgts = sorted(list(set(hgts)))
    lidar_secs = sorted(list(set(lidar_secs)))

  if nc_file == ncar_wp_file:
    instrument_name = "NCAR_WProfiler"
    lidar_file = netCDF4.Dataset(nc_file,'r')
    latitude = ncar_wp_lat
    longitude = ncar_wp_lon
    altitude = ncar_wp_alt
    hgts = lidar_file.variables['heights'][0]
    ws = lidar_file.variables['wspd'][:]
    wd = lidar_file.variables['wdir'][:]
    w_vel = lidar_file.variables['wvert'][:]
    flag = lidar_file.variables['wind_conf'][:]
    time = lidar_file.variables['time'][:]
    # time_start defined in netCDF file
    time_start = (dt.datetime(2017,5,23,0,5,21)-dt.datetime(1970,1,1))\
                  .total_seconds()
    lidar_secs = [t + time_start for t in time]

  if nc_file == ncar_wp_hr_file:
    instrument_name = "NCAR_WProfiler_HighRate"
    lidar_file = netCDF4.Dataset(nc_file,'r')
    latitude = ncar_wp_lat
    longitude = ncar_wp_lon
    altitude = ncar_wp_alt
    hgts = lidar_file.variables['heights'][0]
    elevation = lidar_file.variables['elevation'][:]
    w_vel = lidar_file.variables['vel'][:]
    snr = lidar_file.variables['sigNoiseRatio'][:]
    time = lidar_file.variables['time'][:]
    # time_start defined in netCDF file
    time_start = (dt.datetime(2017,5,23,0,5,21)-dt.datetime(1970,1,1))\
                  .total_seconds()
    lidar_secs = [t + time_start for t in time]

  if nc_file == ZephIR_z423_file or nc_file == ZephIR_z447_file:
    lidar_file = netCDF4.Dataset(nc_file,'r')
    if nc_file == ZephIR_z423_file:
      instrument_name = "ZephIR_z423"
      latitude = ZephIR_z423_lat; longitude = ZephIR_z423_lon
      altitude = ZephIR_z423_alt
    if nc_file == ZephIR_z447_file:
      instrument_name = "ZephIR_z447"
      latitude = ZephIR_z447_lat; longitude = ZephIR_z447_lon
      altitude = ZephIR_z447_alt
    hgts = lidar_file.variables['Heights'][:]
    hgts = [x for x in reversed(hgts)]
    ws = lidar_file.variables['WS'][:]
    wd = lidar_file.variables['WD'][:]
    w_vel = lidar_file.variables['w'][:]
    # Time is seperated into string arrays for y,m,d,H,M
    month = np.array(lidar_file.variables['Month'][:],dtype=int)
    day = np.array(lidar_file.variables['Day'][:],dtype=int)
    hour = np.array(lidar_file.variables['Hour (UTC)'][:],dtype=int)
    minute = np.array(lidar_file.variables['Minute'][:],dtype=int)
    lidar_secs = []
    for x in range(len(month)):
      lidar_secs.append((dt.datetime(2017,month[x],day[x],hour[x],minute[x]\
                        )-dt.datetime(1970,1,1)).total_seconds())

  if nc_file == ENERCON_lidar_file or nc_file == windcube_lidar_file:
    if nc_file == ENERCON_lidar_file:
      instrument_name = "ENERCON"
      latitude = ENERCON_lidar_lat; longitude = ENERCON_lidar_lon
      altitude = ENERCON_lidar_alt; hgts = ENERCON_lidar_hgts
      for_range = list(range(204))[7::12]
    if nc_file == windcube_lidar_file:
      instrument_name = "windcube"
      latitude = windcube_lidar_lat; longitude = windcube_lidar_lon
      altitude = windcube_lidar_alt; hgts = windcube_lidar_hgts
      for_range = list(range(127))[7::12]
    with open (nc_file,encoding="utf8",errors='ignore') as lidar_file:
      lidar = csv.reader(lidar_file,skipinitialspace=True,delimiter='\t')
      lidar_secs = []; str_times = []; ws = []; wd = []; w_vel = []; cnr =[]
      # start at row 43 in the text file to skip header to get time
      # which is in an hour format that needs to be added to May 23rd 
      # 00 UTC second value
      for row in islice(lidar,42,None):
        str_time = row[0]
        lidar_secs.append((dt.datetime.strptime(str_time,'%Y/%m/%d %H:%M')\
                           -dt.datetime(1970,1,1)).total_seconds())
        str_times.append(str_time)
        ws_beam = []; wd_beam = []; w_beam = []; cnr_beam = []
        count = 0
        # set for_range above based on hgts and row definitions
        for i in for_range:
          ws_beam.append(row[i])
          wd_beam.append(row[i+4])
          w_beam.append(row[i+5])
          cnr_beam.append(row[i+7])
        ws.append(ws_beam)        
        wd.append(wd_beam)        
        w_vel.append(w_beam)        
        cnr.append(cnr_beam)        

  if nc_file == CU_Orange_lidar_file or nc_file == CU_Beehive_lidar_file:
    if nc_file == CU_Orange_lidar_file:
      instrument_name = "CU_Orange"
      latitude = CU_Orange_lidar_lat; longitude = CU_Orange_lidar_lon
      altitude = CU_Orange_lidar_alt; hgts = CU_Orange_lidar_hgts
      for_range = list(range(180))[8::19]
    if nc_file == CU_Beehive_lidar_file:
      instrument_name = "CU_Beehive"
      latitude = CU_Beehive_lidar_lat; longitude = CU_Beehive_lidar_lon
      altitude = CU_Beehive_lidar_alt; hgts = CU_Beehive_lidar_hgts
      for_range = list(range(180))[8::19]
    with open (nc_file,encoding="utf8",errors='ignore') as lidar_file:
      lidar = csv.reader(lidar_file,skipinitialspace=True,delimiter='\t')
      lidar_secs = []; str_times = []; u_vel = []; v_vel = []; w_vel = []; avail =[]
      for row in islice(lidar,57,None):
        str_time = row[0]
        lidar_secs.append((dt.datetime.strptime(str_time,'%d/%m/%Y %H:%M:%S')\
                           -dt.datetime(1970,1,1)).total_seconds())
        str_times.append(str_time)
        u_beam = []; v_beam = []; w_beam = []; avail_beam = []
        count = 0
        # set for_range above based on hgts and row definitions
        # convert native u,v,w to meteorological u,v,w
        for i in for_range:
          v_beam.append(-float(row[i]))
          u_beam.append(-float(row[i+2]))
          w_beam.append(-float(row[i+4]))
          avail_beam.append(float(row[i+12]))
        u_vel.append(u_beam)        
        v_vel.append(v_beam)        
        w_vel.append(w_beam)        
        avail.append(avail_beam)        
    
  num_hgts = len(hgts)
  hgts = [int(x) for x in hgts]
  #if nc_file != ncas_file and nc_file != ncas_file_hr:
  altitudes = [x + altitude for x in hgts]
  # ncas profiler hgts already include altitude
  #else: altitudes = hgts
  # append the metdata to each of these lists
  all_latitudes.append(latitude)
  all_longitudes.append(longitude)
  all_altitudes.append(altitude)
  all_instruments.append(instrument_name)
  print(nc_file)
  print(altitudes)
  lidar_file.close()
  #########################################################################
  # Start Loop for each 1 minute time period
  #########################################################################
  # get the indices for data that matches each of the 1 minute time steps
  # set the range to num_steps-1 since the last step will be ind_start to
  # ind_end where ind_end = ind_start+1
  for ind_start in range(num_steps-1):
    ind_end = ind_start+1
    one_min_ind = [i for i in range(len(lidar_secs))\
                    if lidar_secs[i] >= secs[ind_start] and\
                    lidar_secs[i] < secs[ind_end]]
    timestep = dt.datetime.utcfromtimestamp(secs[ind_start]).\
                                            strftime('%H:%M:%S')
    u_one_mins = np.zeros(num_hgts)
    v_one_mins = np.zeros(num_hgts)
    w_one_mins = np.zeros(num_hgts)
    # set the timestep data to lat/lon/alt if there is no data that matches 
    # the given 5 minute time step, but set u,v,w to np.nan
    if len(one_min_ind) == 0:
      u[ind_start,instrument,:] = np.nan
      v[ind_start,instrument,:] = np.nan
      w[ind_start,instrument,:] = np.nan
      continue
    # get the u_one_mins, v_one_mins, and w_one_mins arrays
    # from each of the instruments
    if nc_file == CLAMPS_vad_file:
      # pull all the data for the 5 minute time period into new arrays
      # using one_min_ind[0] since there is only 1 CLAMPS observation
      # within a 5 minute period
      one_min_ws = np.array(ws[one_min_ind[0]])
      one_min_wd = np.array(wd[one_min_ind[0]])
      one_min_intensity = np.array(intensity[one_min_ind[0]])
      # calculate the u and v values
      for i in range(num_hgts):
        # check the intensity value for filtering/QC
        if one_min_intensity[i] > 1.007:
          u_one_mins[i] = -(one_min_ws[i])*math.sin(np.deg2rad(\
                                                      one_min_wd[i]))
          v_one_mins[i] = -(one_min_ws[i])*math.cos(np.deg2rad(\
                                                      one_min_wd[i]))
        else:
          u_one_mins[i] = np.nan
          v_one_mins[i] = np.nan
        w_one_mins[i] = np.nan

    if nc_file == ncas_file or nc_file == ncas_file_hr:
      # pull all the data for the 5 minute time period into new arrays
      one_min_u = np.array([u_vel[i] for i in one_min_ind])
      one_min_v = np.array([v_vel[i] for i in one_min_ind])
      one_min_w = np.array([up_wind[i] for i in one_min_ind])
      one_min_flag = np.array([flag[i] for i in one_min_ind])
      for i in range(num_hgts):
        u_all = np.zeros(len(one_min_ind))
        v_all = np.zeros(len(one_min_ind))
        w_all = np.zeros(len(one_min_ind))
        # set up arrays for averaging u,v,w
        for beam in range(len(one_min_ind)):
          # check the QC flag check
          if one_min_flag[beam][i] == 0:
            one_min_u[beam][i] = np.nan
            one_min_v[beam][i] = np.nan
            one_min_w[beam][i] = np.nan
          u_all[beam] = one_min_u[beam][i]
          v_all[beam] = one_min_v[beam][i]
          w_all[beam] = one_min_w[beam][i]
        u_one_mins[i] = np.nanmean(u_all)
        v_one_mins[i] = np.nanmean(v_all)
        w_one_mins[i] = np.nanmean(w_all)

    if nc_file == UND_Sodar_file:
      # pull all the data for the 5 minute time period into new arrays
      # using one_min_ind[0] since there is only 1 UND_Sodar observation
      # within a 5 minute period
      one_min_u = np.array(u_vel[one_min_ind[0]])
      one_min_v = np.array(v_vel[one_min_ind[0]])
      one_min_w = np.array(w_vel[one_min_ind[0]])
      for i in range(num_hgts):
        u_one_mins[i] = one_min_u[i]
        v_one_mins[i] = one_min_v[i]
        w_one_mins[i] = one_min_w[i]
      print(timestep)
      print(one_min_u)
      print(one_min_v)
      print(one_min_w)

    if nc_file == ncar_wp_file:
      # pull all the data for the 5 minute time period into new arrays
      one_min_ws = np.array(ws[one_min_ind[0]])
      one_min_wd = np.array(wd[one_min_ind[0]])
      one_min_w = np.array(w_vel[one_min_ind[0]])
      one_min_flag = np.array(flag[one_min_ind[0]])
      # calculate the u and v values
      for i in range(num_hgts):
        # calculate the u and v values, and then take the mean
        # quick QC flag check
        if one_min_flag[i] < 0.5:
          one_min_ws[i] = np.nan
          one_min_wd[i] = np.nan
          one_min_w[i] = np.nan
        # convert the -9999 values to NaN
        if one_min_ws[i] == -9999: one_min_ws[i] = np.nan
        if one_min_wd[i] == -9999: one_min_wd[i] = np.nan
        if one_min_w[i] == -9999: one_min_w[i] = np.nan
        u_one_mins[i] = -(one_min_ws[i])*math.sin(np.deg2rad(\
                                                    one_min_wd[i]))
        v_one_mins[i] = -(one_min_ws[i])*math.cos(np.deg2rad(\
                                                    one_min_wd[i]))
        w_one_mins[i] = one_min_w[i]

    if nc_file == ncar_wp_hr_file:
      # pull all the data for the 5 minute time period into new arrays
      one_min_elev = np.array(elevation[one_min_ind])
      one_min_snr = np.array(snr[one_min_ind])
      one_min_w = np.array(w_vel[one_min_ind])
      # set all u,v,w arrays to nans to begin with
      u_one_mins[:] = np.nan
      v_one_mins[:] = np.nan
      w_one_mins[:] = np.nan
      # only use the vertical stares to get w
      for ind in range(len(one_min_elev)):
        if one_min_elev[ind] == 90:
          w_one_mins[:] = np.nan
          for i in range(num_hgts):
            if float(one_min_snr[ind,i]) < -22 or one_min_snr[ind,i] == -9999:
              one_min_w[ind,i] = np.nan
          w_one_mins = one_min_w[ind]
        else:
          continue

    if nc_file == ZephIR_z423_file or nc_file == ZephIR_z447_file:
      # pull all the data for the 5 minute time period into new arrays
      one_min_ws = np.flip(np.array(ws[:,one_min_ind[0]]))
      one_min_wd = np.flip(np.array(wd[:,one_min_ind[0]]))
      one_min_w = np.flip(np.array(w_vel[:,one_min_ind[0]]))
      one_min_w_og = np.array(w_vel[:,one_min_ind[0]])
      # calculate the u and v values
      for i in range(num_hgts):
        # calculate the u and v values, and check for -999
        if one_min_ws[i] == -999 or one_min_wd[i] == -999:
          u_one_mins[i] = np.nan
          v_one_mins[i] = np.nan
        elif one_min_w[i] == -999:
          w_one_mins[i] = np.nan
        else:
          u_one_mins[i] = -(one_min_ws[i])*math.sin(np.deg2rad(\
                                                      one_min_wd[i]))
          v_one_mins[i] = -(one_min_ws[i])*math.cos(np.deg2rad(\
                                                      one_min_wd[i]))
          w_one_mins[i] = one_min_w[i]

    if nc_file == ENERCON_lidar_file or nc_file == windcube_lidar_file:
      # pull all the data for the 5 minute time period into new arrays
      # convert 'N/A' to np.nan for all arrays
      one_min_ws = np.array(ws[one_min_ind[0]])
      one_min_ws = [np.nan if x == 'N/A' else x for x in one_min_ws]
      one_min_wd = np.array(wd[one_min_ind[0]])
      one_min_wd = [np.nan if x == 'N/A' else x for x in one_min_wd]
      one_min_w = np.array(w_vel[one_min_ind[0]])
      one_min_w = [np.nan if x == 'N/A' else x for x in one_min_w]
      one_min_cnr = np.array(cnr[one_min_ind[0]])
      one_min_cnr = [np.nan if x == 'N/A' else x for x in one_min_cnr]
      # calculate the u and v values
      for i in range(num_hgts):
        # checking cnr value for filtering/QC
        if float(one_min_cnr[i]) < -22 or one_min_cnr[i] == np.nan:
          u_one_mins[i] = np.nan
          v_one_mins[i] = np.nan
          w_one_mins[i] = np.nan
        else:
          u_one_mins[i] = -(float(one_min_ws[i]))*math.sin(np.deg2rad(\
                             float(one_min_wd[i])))
          v_one_mins[i] = -(float(one_min_ws[i]))*math.cos(np.deg2rad(\
                             float(one_min_wd[i])))
          w_one_mins[i] = float(one_min_w[i])
 
    if nc_file == CU_Orange_lidar_file or nc_file == CU_Beehive_lidar_file:
      # pull all the data for the 5 minute time period into new arrays
      one_min_u = np.array([u_vel[i] for i in one_min_ind])
      one_min_v = np.array([v_vel[i] for i in one_min_ind])
      one_min_w = np.array([w_vel[i] for i in one_min_ind])
      one_min_avail = np.array([avail[i] for i in one_min_ind])
      for i in range(num_hgts):
        u_all = np.zeros(len(one_min_ind))
        v_all = np.zeros(len(one_min_ind))
        w_all = np.zeros(len(one_min_ind))
        # set up arrays for averaging u,v,w
        for beam in range(len(one_min_ind)):
          # check the avail
          if one_min_avail[beam][i] < 90:
            one_min_u[beam][i] = np.nan
            one_min_v[beam][i] = np.nan
            one_min_w[beam][i] = np.nan
          u_all[beam] = one_min_u[beam][i]
          v_all[beam] = one_min_v[beam][i]
          w_all[beam] = one_min_w[beam][i]
        u_one_mins[i] = np.nanmean(u_all)
        v_one_mins[i] = np.nanmean(v_all)
        w_one_mins[i] = np.nanmean(w_all)

    # use xarray to convert the *_one_mins to match the tot_hgts index
    # then write those converted values to final u,v,w arrays
    ds_u = xr.DataArray(u_one_mins,[("height",altitudes)])
    ds_u_reindex = ds_u.reindex(height=all_hgts)
    ds_v = xr.DataArray(v_one_mins,[("height",altitudes)])
    ds_v_reindex = ds_v.reindex(height=all_hgts)
    ds_w = xr.DataArray(w_one_mins,[("height",altitudes)])
    ds_w_reindex = ds_w.reindex(height=all_hgts)
    u[ind_start,instrument] = ds_u_reindex.values
    v[ind_start,instrument] = ds_v_reindex.values
    w[ind_start,instrument] = ds_w_reindex.values
  #########################################################################
  # END Loop for each 5 minute time period
  #########################################################################

# create netcdf file with defined dimensions
# variables, units, and descriptions
Perdigao_grid_file = netCDF4.Dataset('Perdigao_uvw_20170523.nc','w',\
                                      format='NETCDF4')
Perdigao_grid_file.title = "Perdigao 2020/05/23 Wind Components (u,v,w)" +\
                            " every 1 minute 12 - 00 UTC"
Perdigao_grid_file.createDimension('time', num_steps-1)
Perdigao_grid_file.createDimension('height', tot_hgts)
Perdigao_grid_file.createDimension('instrument', len(all_instruments))
instruments = Perdigao_grid_file.createVariable('instrument',str,\
                                                'instrument')
all_instruments = np.array(all_instruments, dtype=object)
instruments[:] = all_instruments
latitudes = Perdigao_grid_file.createVariable('latitude','f4','instrument')
latitudes[:] = all_latitudes
latitudes.units = 'degrees north'
longitudes = Perdigao_grid_file.createVariable('longitude','f4',\
                                               'instrument')
longitudes[:] = all_longitudes
longitudes.units = 'degrees east'
alt = Perdigao_grid_file.createVariable('altitude','f4','instrument')
alt[:] = all_altitudes
alt.units = 'meters'
altitudes_var = Perdigao_grid_file.createVariable('height','f4','height')
altitudes_var[:] = all_hgts
altitudes_var.units = 'meters asl (includes alt)'
base_time = Perdigao_grid_file.createVariable('base_time','i')
base_time[:] = (dt.datetime.strptime('23/05/2017 00:00:00','%d/%m/%Y %H:%M:%S')\
                   -dt.datetime(1970,1,1)).total_seconds()
base_time.units = 'seconds since 1970-01-01 00:00:00'
epoch_secs = Perdigao_grid_file.createVariable('time','f4','time')
epoch_secs[:] = secs_day[:-1]
epoch_secs.units = 'seconds since 2017-05-23 00:00:00'
east_vel = Perdigao_grid_file.createVariable('u','f4',('time','instrument',\
                                             'height'),fill_value='NaN')
east_vel[:] = u
east_vel.units = 'm/s'
east_vel.standard_name = 'eastern velocity'
north_vel = Perdigao_grid_file.createVariable('v','f4',('time','instrument'\
                                              ,'height'),fill_value='NaN')
north_vel[:] = v
north_vel.units = 'm/s'
north_vel.standard_name = 'northern velocity'
vert_vel = Perdigao_grid_file.createVariable('w','f4',('time','instrument',\
                                             'height'),fill_value='NaN')
vert_vel[:] = w
vert_vel.units = 'm/s'
vert_vel.standard_name = 'vertical velocity'

Perdigao_grid_file.close()
###########################################################################
# END Loop for netCDF file
###########################################################################
