import datetime
import os
from WRF_wrapper import WRF_wrapper, WindFarm

farm_names = ['dudgeon', 'lynn_lincs_inner_dowsing', 'race_bank', 'sheringham']
rated_powers = [6, 3.6, 6, 3.6]

wfs = [WindFarm.WindFarm.from_template(template_id=farm_name, type_id='plausible-{:.1f}'.format(rated_power)) for farm_name, rated_power in zip(farm_names, rated_powers)]

config_dict = {'start_date': datetime.datetime(2000, 1, 1),
               'end_date': datetime.datetime(2000, 1, 2),
               'era_data_dir': '/home/scw15/wrf_test/era_data/2000/grb/',
               'pfile_data_dir': '/home/scw15/wrf_test/era_data/2000/PFILES',
               'geog_data_path': '/home/scw15/cx1_wrfwindpower/live/geog/WPS_GEOG/',
               'wind_farms': wfs,
               'WPS_path': '/home/scw15/wrf_test/WPS',
               'WRF_path': '/home/scw15/wrf_test/WRF',
               'wrf_img_path': '/home/scw15/wrf_test/wrf_sing',
               }

working_directory = os.path.join(os.getcwd(), 'test_real_farms')
outdir = os.path.join(os.getcwd(), 'test_real_farms_outputs')
wrf = WRF_wrapper(working_directory, config_dict)

wrf.initialise()
wrf.run_geogrid()
wrf.download_era()
wrf.run_link_grib()
wrf.run_ungrib()
wrf.run_metgrid()
wrf.run_real()
wrf.run_wrf()
wrf.extract_outputs(outdir)
