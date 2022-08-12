import datetime
import os
from WRF_wrapper import WRF_wrapper, WindFarm

wf = WindFarm.WindFarm.from_template(template_id=0, farm_lon=0, farm_lat=1, type_id=6)

config_dict = {'start_date': datetime.datetime(2000, 1, 4, 12),
               'end_date': datetime.datetime(2000, 1, 4, 15),
               'era_data_dir': '/home/scw15/wrf_test/era_data/2000/grb/',
               'pfile_data_dir': '/home/scw15/wrf_test/era_data/2000/PFILES',
               'geog_data_path': '/home/scw15/cx1_wrfwindpower/live/geog/WPS_GEOG/',
               'wind_farms': [wf],
               'WPS_path': '/home/scw15/wrf_test/WPS',
               'WRF_path': '/home/scw15/wrf_test/WRF',
               'wrf_img_path': '/home/scw15/wrf_test/wrf_sing',
               }

working_directory = os.path.join(os.getcwd(), 'test')
wrf = WRF_wrapper(working_directory, config_dict)

wrf.initialise()
wrf.run_geogrid()
wrf.run_link_grib()
wrf.run_ungrib()
wrf.run_metgrid()
wrf.run_real()
wrf.run_wrf()
