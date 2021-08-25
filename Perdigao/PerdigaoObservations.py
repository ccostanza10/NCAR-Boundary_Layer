'''
Python functions and classes for reading Perdiago Observations.
Created: Aug. 25, 2021
'''

import xarray as xr
import pandas as pd
import glob
from datetime import datetime

def read_Sondes(obs_dir,days_of_interest):
    '''
    Read in the radiosonde files for day(s) of interest. If days_of_interest is a list of dates, then each date will be included.
    obs_dir (str) : location of observations
    days_of_interest list(str or pandas datetime) : dates of interest (YYYY-MM-DD); hours, minutes, seconds not supported.
    '''
    if type(days_of_interest) is str:
        days_of_interest = [days_of_interest]
        
    try:
        list(days_of_interest)
    except:
        days_of_interest = [days_of_interest]
        
    obs_sondes = {}

    for dd,day_of_interest in enumerate(days_of_interest):
        if (type(day_of_interest) is not pd.datetime) or (type(day_of_interest) is not pd.Timestamp):
            day_of_interest = pd.to_datetime(day_of_interest)
        sondes_list = sorted(glob.glob('{0}*{1:04d}{2:02d}{3:02d}*.nc'.format(obs_dir,day_of_interest.year,day_of_interest.month,day_of_interest.day)))
        if sondes_list == []:
            print('No data found for {}'.format(day_of_interest))
        for sonde in sondes_list:
            sonde_str = sonde.split('_')[0].split('/')[-1].replace('D','') + '_' + sonde.split('_')[1]
            sonde_time = datetime.strptime(sonde_str,'%Y%m%d_%H%M%S')
            obs_sondes[str(sonde_time)] = xr.open_dataset(sonde)
    return(obs_sondes)