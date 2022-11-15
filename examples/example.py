import datetime
import os
from WRF_wrapper import WRF_wrapper, WindFarm

# Generate a farm from a template
farm_names = ['dudgeon', 'lynn_lincs_inner_dowsing', 'race_bank', 'sheringham']
rated_powers = [6, 3.6, 6, 3.6]
wf1 = WindFarm.WindFarm.from_template(template_id='dudgeon', type_id='plausible-6.0')

# Generate a generic farm based on installed capacity
wf2 = WindFarm.WindFarm.from_installed_capacity(installed_capacity=90, farm_lon=2, farm_lat=54, type_id=6)

config_dict = {'start_date': datetime.datetime(2000, 1, 4, 12),
               'end_date': datetime.datetime(2000, 1, 4, 15),
               'era_data_dir': '/media/scw15/DATADRIVE1/shell_project/wrf_test/era_data/grb/',
               'pfile_data_dir': '/media/scw15/DATADRIVE1/shell_project/wrf_test/era_data/PFILES',
               'geog_data_path': '/home/scw15/cx1_wrfwindpower/live/geog/WPS_GEOG/',
               'wind_farms': [wf1, wf2],
               'WPS_path': '/media/scw15/DATADRIVE1/shell_project/wrf_test/singularity/WRF_Singularity/Run_v4.4.1/WPS',
               'WRF_path': '/media/scw15/DATADRIVE1/shell_project/wrf_test/singularity/WRF_Singularity/Run_v4.4.1/WRF',
               'wrf_img_path': '/media/scw15/DATADRIVE1/shell_project/wrf_test/singularity/WRF_Singularity/Run_v4.4.1/wrf.sif',
               }

working_directory = os.path.join(os.getcwd(), 'test')
outdir = os.path.join(os.getcwd(), 'test_outputs')

wrf = WRF_wrapper(working_directory, config_dict, domain_config='wash')

wrf.initialise()
wrf.run_geogrid()
wrf.download_era()
wrf.run_link_grib()
wrf.run_ungrib()
wrf.run_metgrid()
wrf.run_real()
wrf.run_wrf()
wrf.extract_outputs(outdir)
