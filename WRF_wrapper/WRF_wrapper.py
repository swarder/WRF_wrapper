try:
    from . import era_downloader
except ModuleNotFoundError:
    print("Unable to import era_downloader. Let's hope you weren't planning to use it")
from . import WPS_namelist_template
from . import WRF_namelist_template

import shutil
import os
import stat
import subprocess
import datetime
import glob
from netCDF4 import Dataset

config_defaults = {'interval_seconds': 10800,
                   'wind_farms': None,
                   'damp_opt': 0,
                   'dampcoef': 0.05,
                   'w_damping': 0,
                   'epssm': 0.1,
                   }

domain_presets = {
    'denmark': {'max_dom': 2, 'time_step': 45, 'dx': 16668, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': 100, 'd02_j_parent_start': 92, 'd02_e_we': 73, 'd02_e_sn': 154, 'ref_lat': 55.5, 'ref_lon': 6, 'truelat1': 54, 'truelat2': 56, 'stand_lon': 8},
    'wash': {'max_dom': 2, 'time_step': 45, 'dx': 18000, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': 74, 'd02_j_parent_start': 78, 'd02_e_we': 163, 'd02_e_sn': 163, 'ref_lat': 55.5, 'ref_lon': 6, 'truelat1': 51.5, 'truelat2': 55, 'stand_lon': 1},
    'thames': {'max_dom': 2, 'time_step': 45, 'dx': 16668, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': 73, 'd02_j_parent_start': 67, 'd02_e_we': 172, 'd02_e_sn': 172, 'ref_lat': 55.5, 'ref_lon': 6, 'truelat1': 54, 'truelat2': 56, 'stand_lon': 8},
    'liverpool': {'max_dom': 2, 'time_step': 45, 'dx': 16668, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': 54, 'd02_j_parent_start': 83, 'd02_e_we': 172, 'd02_e_sn': 172, 'ref_lat': 55.5, 'ref_lon': 6, 'truelat1': 54, 'truelat2': 56, 'stand_lon': 8},
    'houston': {'max_dom': 1, 'time_step': 45, 'dx': 16668, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': ' ', 'd02_j_parent_start': ' ', 'd02_e_we': ' ', 'd02_e_sn': ' ', 'ref_lat': 29.8, 'ref_lon': -95.4, 'truelat1': 24, 'truelat2': 34, 'stand_lon': -95},
    'north_sea': {'max_dom': 2, 'time_step': 30, 'dx': 18000, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202, 'd02_i_parent_start': 67, 'd02_j_parent_start': 64, 'd02_e_we': 541, 'd02_e_sn': 730, 'ref_lat': 55.5, 'ref_lon': 4, 'truelat1': 53.5, 'truelat2': 57.5, 'stand_lon': 3.4, 'damp_opt': 3, 'dampcoef': 0.2, 'epssm': 0.3, 'w_damping': 1},
    'southern_north_sea': {'max_dom': 2, 'time_step': 30, 'dx': 18000, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202,
                           'd02_i_parent_start': 76, 'd02_j_parent_start': 64, 'd02_e_we': 460, 'd02_e_sn': 541,
                           'ref_lat': 55.5, 'ref_lon': 4, 'truelat1': 53.5, 'truelat2': 57.5, 'stand_lon': 3.4,
                           'damp_opt': 3, 'dampcoef': 0.2, 'epssm': 0.3, 'w_damping': 1},
    'northwestern_north_sea': {'max_dom': 2, 'time_step': 30, 'dx': 18000, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202,
                               'd02_i_parent_start': 66, 'd02_j_parent_start': 86, 'd02_e_we': 316, 'd02_e_sn': 442,
                               'ref_lat': 55.5, 'ref_lon': 4, 'truelat1': 53.5, 'truelat2': 57.5, 'stand_lon': 3.4,
                               'damp_opt': 3, 'dampcoef': 0.2, 'epssm': 0.3, 'w_damping': 1},
    'northeastern_north_sea': {'max_dom': 2, 'time_step': 30, 'dx': 18000, 'd02_grid_ratio': 9, 'd02_dt_ratio': 3, 'd01_e_we': 202, 'd01_e_sn': 202,
                               'd02_i_parent_start': 92, 'd02_j_parent_start': 94, 'd02_e_we': 316, 'd02_e_sn': 370,
                               'ref_lat': 55.5, 'ref_lon': 4, 'truelat1': 53.5, 'truelat2': 57.5, 'stand_lon': 3.4,
                               'damp_opt': 3, 'dampcoef': 0.2, 'epssm': 0.3, 'w_damping': 1},
}

class WRF_wrapper:
    """
    Class for running WPS and WRF
    """
    def __init__(self, working_directory, config_dict, domain_config=None):
        self.working_directory = working_directory

        # Combine config dicts
        if domain_config is not None:
            domain_dict = domain_presets[domain_config]
        else:
            domain_dict = {}
        self.config_dict = {**config_defaults, **domain_dict, **config_dict}

        self.config_dict['dx_d02'] = self.config_dict['dx'] // self.config_dict['d02_grid_ratio']

        assert 'WPS_path' in self.config_dict.keys()
        assert 'WRF_path' in self.config_dict.keys()
        assert 'wrf_img_path' in self.config_dict.keys()

        # Pre-process some config fields
        self.config_dict['start_date_str'] = self.config_dict['start_date'].strftime('%Y-%m-%d_%H:%M:00')
        self.config_dict['end_date_str'] = self.config_dict['end_date'].strftime('%Y-%m-%d_%H:%M:00')
        self.config_dict['wind_turbines_int'] = 1 if self.config_dict['wind_farms'] else 0

        # Get list of model timestamps
        sd = self.config_dict['start_date']
        ed = self.config_dict['end_date']
        interval = self.config_dict['interval_seconds']
        self.model_timestamps = [sd + datetime.timedelta(seconds=i*interval) for i in range(int((ed - sd).total_seconds() / interval)+1)]

        # Check which PFILEs already exist - new data won't need to be downloaded for these, even if the GRB files don't exist any more
        existing_pfile_names = os.listdir(self.config_dict['pfile_data_dir'])
        self.model_timestamps_without_pfiles = [ts for ts in self.model_timestamps if ts.strftime('FILE:%Y-%m-%d_%H') not in existing_pfile_names]

    def save_WPS_config_to_file(self, filename=None):
        """
        Save WPS config file to specified filename
        """
        config_file_string = WPS_namelist_template.template.format(**self.config_dict)
        if filename is None:
            filename = os.path.join(self.working_directory, 'WPS/namelist.wps')
        f = open(filename, 'w')
        f.write(config_file_string)
        f.close()

    def save_WRF_config_to_file(self, filename=None):
        """
        Save WRF config file to specified filename
        """
        config_file_string = WRF_namelist_template.template.format(**self.config_dict)
        if filename is None:
            filename = os.path.join(self.working_directory, 'WRF/test/em_real/namelist.input')
        f = open(filename, 'w')
        f.write(config_file_string)
        f.close()

        # Add extra io files
        io_str = r'+:h:0:RMOL' # Add RMOL as output variable to stream 0
        io_filename_d01 = os.path.join(self.working_directory, 'WRF/test/em_real/io_file_d01.txt')
        io_filename_d02 = os.path.join(self.working_directory, 'WRF/test/em_real/io_file_d02.txt')
        for filename in [io_filename_d01, io_filename_d02]:
            f = open(filename, 'w')
            f.write(io_str)
            f.close()

    def initialise(self):
        """
        Set up directory structure and copy files
        """
        # Create working directory and copy in template WPS and WRF directories
        os.makedirs(self.working_directory, exist_ok=True)
        try:
            shutil.copytree(self.config_dict['WPS_path'], os.path.join(self.working_directory, 'WPS'), symlinks=True)
        except FileExistsError:
            print('Caution: WPS directory already exists. Existing namelist.wps will be overwritten.')
            pass
        try:
            shutil.copytree(self.config_dict['WRF_path'], os.path.join(self.working_directory, 'WRF'), symlinks=True)
        except FileExistsError:
            print('Caution: WRF directory already exists. Existing namelist.input will be overwritten.')
            pass

        # Generate namelist files
        self.save_WPS_config_to_file()
        self.save_WRF_config_to_file()

        # Combine wind farms and save wind turbine file
        if self.config_dict['wind_farms']:
            all_farms = sum(self.config_dict['wind_farms'])
            all_farms.save(os.path.join(self.working_directory, 'WRF/test/em_real/'))

    def run_geogrid(self):
        """
        Run geogrid.exe within singularity container
        """
        print('Running geogrid')
        os.chdir(os.path.join(self.working_directory, 'WPS'))
        st = os.stat('geogrid.exe')
        os.chmod('geogrid.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.config_dict['geog_data_path']},$TMPDIR:$TMPDIR {self.config_dict['wrf_img_path']} ./geogrid.exe > geogrid.out", shell=True)

    def download_era(self):
        """
        Download ERA data required for specified run
        Only downloads timestamps for which we haven't already generated PFILEs
        Additional checks are performed in era_downloader to remove unnecessary downloads
        """
        downloader = era_downloader.era_downloader(data_dir=self.config_dict['era_data_dir'])
        print('Attempting to download ERA data for the following timestamps:', self.model_timestamps_without_pfiles)
        downloader.download_datetimes(self.model_timestamps_without_pfiles)

    def run_link_grib(self):
        """
        Run link_grib.csh within singularity container
        Only runs for PFILEs we haven't already generated
        """
        print('Running link_grib')
        os.chdir(os.path.join(self.working_directory, 'WPS'))
        st = os.stat('link_grib.csh')
        os.chmod('link_grib.csh', st.st_mode | stat.S_IEXEC)
        # We only want to link the grb files we need
        grb_files = [self.config_dict['era_data_dir'] + t.strftime('/ERA5_%Y%m%d_%H%M.grb') for t in self.model_timestamps_without_pfiles]
        print('Running ungrib for the following timestamps:', self.model_timestamps_without_pfiles)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.config_dict['era_data_dir']},$TMPDIR:$TMPDIR {self.config_dict['wrf_img_path']} ./link_grib.csh {' '.join(grb_files)}", shell=True)

    def run_ungrib(self):
        """
        Run ungrib.exe within singularity container
        """
        print('Running ungrib')
        os.chdir(os.path.join(self.working_directory, 'WPS'))
        st = os.stat('ungrib.exe')
        os.chmod('ungrib.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.config_dict['era_data_dir']},{self.config_dict['pfile_data_dir']},$TMPDIR:$TMPDIR {self.config_dict['wrf_img_path']} ./ungrib.exe > ungrib.out", shell=True)

    def run_metgrid(self):
        """
        Run metgrid.exe within singularity container
        """
        print('Running metgrid')
        os.chdir(os.path.join(self.working_directory, 'WPS'))
        st = os.stat('metgrid.exe')
        os.chmod('metgrid.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.config_dict['pfile_data_dir']},$TMPDIR:$TMPDIR {self.config_dict['wrf_img_path']} ./metgrid.exe > metgrid.out", shell=True)

    def run_real(self):
        """
        Run real.exe within singularity container
        """
        print('Running real.exe')
        os.chdir(os.path.join(self.working_directory, 'WRF/test/em_real'))
        subprocess.check_call(f"ln -sf {self.working_directory}/WPS/met_em.d* .", shell=True)

        st = os.stat('real.exe')
        os.chmod('real.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} ./real.exe > real.out", shell=True)

    def run_wrf(self):
        """
        Run wrf.exe within singularity container
        """
        print('Running WRF')
        os.chdir(os.path.join(self.working_directory, 'WRF/test/em_real'))
        st = os.stat('wrf.exe')
        os.chmod('wrf.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} ./wrf.exe > wrf.out", shell=True)

    def extract_outputs(self, out_path, domains=None, variables_list=['U', 'V', 'TKE_PBL', 'POWER'], max_z_levels=None):
        """
        Extract specified variables and save each field in a newly generated netcdf file for each domain
        """
        if domains is None:
            domains = range(1, self.config_dict['max_dom']+1)

        # Remove POWER variable if no wind farms
        if 'POWER' in variables_list and self.config_dict['wind_farms'] is None:
            del variables_list[variables_list.index('POWER')]

        os.makedirs(out_path, exist_ok=True)

        for dom in domains:
            dout = Dataset(os.path.join(out_path, f'out_d_0{dom}.nc'), 'w', format='NETCDF4')
            din = Dataset(os.path.join(self.working_directory, f'WRF/test/em_real/wrfout_d0{dom}'))

            dout.setncatts(din.__dict__)
            for name, dimension in din.dimensions.items():
                if name in ['bottom_top', 'bottom_top_stag'] and max_z_levels is not None:
                    dout.createDimension(name, max_z_levels)
                else:
                    dout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            for name, variable in din.variables.items():
                if name in variables_list or name in ['XLAT', 'XLONG']:
                    z_dim = None
                    for z_name in ['bottom_top', 'bottom_top_stag']:
                        if z_name in variable.dimensions:
                            z_dim = variable.dimensions.index(z_name)
                            break
                    if (z_dim is not None) and (max_z_levels is not None):
                        data = din[name][:].take(range(max_z_levels), axis=z_dim)
                    else:
                        data = din[name][:]
                    x = dout.createVariable(name, variable.datatype, variable.dimensions)
                    dout[name][:] = data
                    dout[name].setncatts(din[name].__dict__)

            din.close()
            dout.close()

        if self.config_dict['wind_farms']:
            shutil.copy(os.path.join(self.working_directory, f'WRF/test/em_real/windturbines.txt'), os.path.join(out_path, 'windturbines.txt'))

    def clean_up(self):
        """
        Delete working directory
        """
        try:
            shutil.rmtree(self.working_directory)
        except FileNotFoundError:
            pass
