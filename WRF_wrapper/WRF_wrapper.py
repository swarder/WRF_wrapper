from . import era_downloader
from . import WPS_namelist_template
from . import WRF_namelist_template

import shutil
import os
import stat
import subprocess
import datetime

config_defaults = {'interval_seconds': 10800,
                   'wind_turbines': None,
                   }

class WRF_wrapper:
    """
    Class for running WPS and WRF
    """
    def __init__(self, working_directory, config_dict):
        self.working_directory = working_directory
        self.config_dict = {**config_defaults, **config_dict}

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

    def initialise(self):
        """
        Set up directory structure and copy files
        """
        # Download ERA data (additional checks are performed in era_downloader to remove unnecessary downloads)
        downloader = era_downloader.era_downloader(data_dir=self.config_dict['era_data_dir'])
        downloader.download_datetimes(self.model_timestamps)

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

    def run_link_grib(self):
        """
        Run link_grib.csh within singularity container
        """
        print('Running link_grib')
        os.chdir(os.path.join(self.working_directory, 'WPS'))
        st = os.stat('link_grib.csh')
        os.chmod('link_grib.csh', st.st_mode | stat.S_IEXEC)
        # We only want to link the grb files we need
        grb_files = [self.config_dict['era_data_dir'] + t.strftime('/ERA5_%Y%m%d_%H%M.grb') for t in self.model_timestamps]
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

        if self.config_dict['wind_turbines']:
            #turbine_df.to_csv(path_name + '/WRF/test/em_real/windturbines.txt', sep = ' ', header = False, index = False)
            raise NotImplementedError

        st = os.stat('real.exe')
        os.chmod('real.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity run {self.config_dict['wrf_img_path']} ./real.exe > real.out", shell=True)

    def run_wrf(self):
        """
        Run wrf.exe within singularity container
        """
        print('Running WRF')
        os.chdir(os.path.join(self.working_directory, 'WRF/test/em_real'))
        st = os.stat('wrf.exe')
        os.chmod('wrf.exe', st.st_mode | stat.S_IEXEC)
        subprocess.check_call(f"singularity run {self.config_dict['wrf_img_path']} ./wrf.exe > wrf.out", shell=True)
