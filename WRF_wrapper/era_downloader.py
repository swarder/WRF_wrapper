"""
Based on ERA download script by
Andrea Hahmann and Martin Doerenkaemper
(c) DTU
"""

import os
import re
import glob
import shutil
import cdsapi
from cdo import *
c = cdsapi.Client()
import numpy as np
import time

plfn = "ERA5_pl_"
sffn = "ERA5_sfc_"

def download_day(day, times):

    daystr = day.strftime("%Y-%m-%d")
    dayfns = day.strftime("%Y%m%d")
    year = day.strftime("%Y%m%d")

    c.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": 'reanalysis',
            'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'mean_sea_level_pressure', 'mean_wave_direction',
            'mean_wave_period', 'sea_surface_temperature', 'significant_height_of_combined_wind_waves_and_swell',
            'surface_pressure', 'total_precipitation',
            ],
            'year': year[:4],
            'month': year[4:6],
            'day': year[6:],
            "time" : times,
            'format':'grib'
        },
        sffn+dayfns+".grb")

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type':'reanalysis',
            'variable': [
                'geopotential', 'relative_humidity', 'temperature',
                'u_component_of_wind', 'v_component_of_wind',
             ],
            'pressure_level':[
                '10','20','30','50','70','100',
                '125','150','175','200','225','250',
                '300','350','400','450','500','550',
                '600','650','700','750','775','800',
                '825','850','875','900','925','950',
                '975','1000'
            ],
            'year': year[:4],
            'month': year[4:6],
            'day': year[6:],
            "time" : times,
            'format':'grib'
        },
        plfn+dayfns+".grb")

def get_and_split_day(date, times, data_dir):
    tempdir = data_dir + str(date) + '_temp_' + str(time.time()).replace(',', '') + '/'
    os.makedirs(tempdir, exist_ok=True)
    os.chdir(tempdir)
    print( "+++++++++++++++++ ERA-5 +++++++++++++++++++")
    print( "++++++++ Download data for "+date.strftime("%Y-%m-%d"))
    yr = str(date.year).zfill(4)
    mn = str(date.month).zfill(2)
    fpath = data_dir
    ## ERA5_19930717_0000.grb
    if not os.path.exists(fpath):
        os.makedirs(fpath)
    sstat = download_day(date, times)
    pfile = plfn + date.strftime("%Y%m%d") + ".grb"
    sfile = sffn + date.strftime("%Y%m%d") + ".grb"
    ofile = "ERA5_"+date.strftime("%Y%m%d") + ".grb"
    print( "++++++++ Split data for "+date.strftime("%Y-%m-%d"))
    cdo = Cdo()
    cdo.sellevel('1000/100000', input=pfile, output='ERA5_'+date.strftime("%Y-%m-%d")+'_test2.grb')
    cdo.merge(input='ERA5_'+date.strftime("%Y-%m-%d")+'_test2.grb '+sfile, output=ofile)
    cdo.splithour(input=ofile,output=ofile.split(".grb")[0]+"_")
    fall = glob.glob("ERA5_"+date.strftime("%Y%m")+"*.grb")
    for fn1 in fall:
        if re.match('ERA5_\d{4}\d{2}\d{2}_\d{2}.grb',fn1):
            fnn2 = fn1.split(".grb")[0] + "00.grb"
            fnn2 = fpath + fnn2
            shutil.move(fn1,fnn2)
    os.chdir(data_dir)
    shutil.rmtree(tempdir)

class era_downloader:
    def __init__(self, data_dir):
        self.data_dir = data_dir

    def download_datetimes(self, datetimes):
        """
        Take a list of datetimes and download all data
        """
        print('ERA data requested for following timestamps:', datetimes)
        # First filter datetimes for ones which aren't already downloaded
        datetimes = [d for d in datetimes if not os.path.exists(os.path.join(self.data_dir, d.strftime('ERA5_%Y%m%d_%H%M.grb')))]
        print('Remaining timestamps after removing those for which we already have ERA data (in GRB format):', datetimes)

        if len(datetimes) == 0:
            print('All data already downloaded')
            return

        dates = [d.date() for d in datetimes]
        times = [f'{d.hour:02}:{d.minute:02}' for d in datetimes]

        unique_dates = np.unique(dates)
        for date in unique_dates:
            get_and_split_day(date, [t for d, t in zip(dates, times) if d == date], self.data_dir)

if __name__ == '__main__':
    import datetime
    downloader = era_downloader(data_dir='/home/scw15/wrf_test/era_data/2000/grb/')
    dates = [datetime.datetime(2000, 1, 1) + datetime.timedelta(seconds=3600*i) for i in range(24)]
    downloader.download_datetimes(dates)
