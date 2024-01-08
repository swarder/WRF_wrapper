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
                   'zdamp': 5000.0,
                   'damp_opt': 0,
                   'dampcoef': 0.05,
                   'w_damping': 0,
                   'epssm': 0.1,
                   'feedback': 0,
                   'smooth_option': 0,
                   'max_dom': 2,
                   'd03_grid_ratio': None,
                   'd03_dt_ratio': None,
                   'd03_i_parent_start': None,
                   'd03_j_parent_start': None,
                   'd03_e_we': None,
                   'd03_e_sn': None,
                   'd01_radt': 18,
                   'd02_radt': 6,
                   'd03_radt': 2,
                   'bl_mynn_tkeadvect_d01': '.true.',
                   'bl_mynn_tkeadvect_d02': '.true.',
                   'bl_mynn_tkeadvect_d03': '.true.',
                   'd01_bl_mynn_tkebudget': 0,
                   'd02_bl_mynn_tkebudget': 1,
                   'd03_bl_mynn_tkebudget': 1,
                   'd01_diff_6th_factor': 0.06,
                   'd02_diff_6th_factor': 0.1,
                   'd03_diff_6th_factor': 0.1,
                   'w_crit_cfl': 1.2,
                   'use_adaptive_time_step': '.false.',
                   'step_to_output_time': '.false.',
                   'd01_target_cfl': 1.2,
                   'd02_target_cfl': 1.2,
                   'd03_target_cfl': 1.2,
                   'd01_target_hcfl': 0.84,
                   'd02_target_hcfl': 0.84,
                   'd03_target_hcfl': 0.84,
                   'd01_max_step_increase_pct': 5,
                   'd02_max_step_increase_pct': 51,
                   'd03_max_step_increase_pct': 51,
                   'd01_starting_time_step': -1,
                   'd02_starting_time_step': -1,
                   'd03_starting_time_step': -1,
                   'd01_max_time_step': -1,
                   'd02_max_time_step': -1,
                   'd03_max_time_step': -1,
                   'd01_min_time_step': -1,
                   'd02_min_time_step': -1,
                   'd03_min_time_step': -1,
                   'adaptation_domain': 2,
                   'mp_physics_d01': 8,
                   'mp_physics_d02': 8,
                   'mp_physics_d03': 8,
                   'ra_lw_physics_d01': 4,
                   'ra_lw_physics_d02': 4,
                   'ra_lw_physics_d03': 4,
                   'ra_sw_physics_d01': 4,
                   'ra_sw_physics_d02': 4,
                   'ra_sw_physics_d03': 4,
                   'sf_surface_physics_d01': 2,
                   'sf_surface_physics_d02': 2,
                   'sf_surface_physics_d03': 2,
                   'sf_sfclay_physics_d01': 5,
                   'sf_sfclay_physics_d02': 5,
                   'sf_sfclay_physics_d03': 5,
                   'bl_pbl_physics_d01': 5,
                   'bl_pbl_physics_d02': 5,
                   'bl_pbl_physics_d03': 5,
                   'bldt_d01': 0,
                   'bldt_d02': 0,
                   'bldt_d03': 0,
                   'cu_physics_d01': 1,
                   'cu_physics_d02': 0,
                   'cu_physics_d03': 0,
                   'cudt_d01': 5,
                   'cudt_d02': 0,
                   'cudt_d03': 0,
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
    def __init__(self, working_directory, config_dict, domain_config=None, parallel=False):
        self.working_directory = working_directory

        # Combine config dicts
        if domain_config is not None:
            domain_dict = domain_presets[domain_config]
        else:
            domain_dict = {}
        self.config_dict = {**config_defaults, **domain_dict, **config_dict}

        self.config_dict['dx_d02'] = self.config_dict['dx'] // self.config_dict['d02_grid_ratio']
        if self.config_dict['max_dom'] > 2:
            self.config_dict['dx_d03'] = self.config_dict['dx_d02'] // self.config_dict['d03_grid_ratio']

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

        # Set up parallelisation
        if parallel:
            if os.getenv('PBS_NODEFILE') is not None:
                self.parallel_mode = 'PBS'
                self.hostfile_path = os.path.join(self.working_directory, 'WRF/test/em_real/hostfile')
            else:
                self.parallel_mode = 'local'
            self.n_procs = parallel
        else:
            self.parallel_mode = None

    def save_WPS_config_to_file(self, filename=None):
        """
        Save WPS config file to specified filename
        """
        if self.config_dict['max_dom'] == 2:
            template = WPS_namelist_template.template
        elif self.config_dict['max_dom'] == 3:
            template = WPS_namelist_template.template_3_domains
        else:
            raise NotImplementedError
        config_file_string = template.format(**self.config_dict)
        if filename is None:
            filename = os.path.join(self.working_directory, 'WPS/namelist.wps')
        f = open(filename, 'w')
        f.write(config_file_string)
        f.close()

    def save_WRF_config_to_file(self, filename=None):
        """
        Save WRF config file to specified filename
        """
        if self.config_dict['max_dom'] == 2:
            template = WRF_namelist_template.template
        elif self.config_dict['max_dom'] == 3:
            template = WRF_namelist_template.template_3_domains
        else:
            raise NotImplementedError
        config_file_string = template.format(**self.config_dict)
        if filename is None:
            filename = os.path.join(self.working_directory, 'WRF/test/em_real/namelist.input')
        f = open(filename, 'w')
        f.write(config_file_string)
        f.close()

        # Add extra io files
        io_str = r'+:h:0:RMOL' # Add RMOL as output variable to stream 0
        for i in range(self.config_dict['max_dom']):
            io_filename = os.path.join(self.working_directory, f'WRF/test/em_real/io_file_d0{i+1}.txt')
            f = open(io_filename, 'w')
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

        if self.parallel_mode == 'PBS':
            with open(os.getenv('PBS_NODEFILE')) as f:
                nodes = [f.strip() for f in f.readlines()]
            assert len(set(nodes)) == 1, 'PBS_NODEFILE contains multiple nodes'
            nodes = nodes[:1] * self.n_procs
            with open(self.hostfile_path, 'w') as f:
                f.write('\n'.join(nodes))

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
        if self.parallel_mode == 'PBS':
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} mpiexec --hostfile {self.hostfile_path} ./real.exe > real.out", shell=True)
        elif self.parallel_mode == 'local':
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} mpiexec -n {self.n_procs} ./real.exe > real.out", shell=True)
        else:
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} ./real.exe > real.out", shell=True)

    def run_wrf(self):
        """
        Run wrf.exe within singularity container
        """
        print('Running WRF')
        os.chdir(os.path.join(self.working_directory, 'WRF/test/em_real'))
        st = os.stat('wrf.exe')
        os.chmod('wrf.exe', st.st_mode | stat.S_IEXEC)
        if self.parallel_mode == 'PBS':
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} mpiexec --hostfile {self.hostfile_path} ./wrf.exe > wrf.out", shell=True)
        elif self.parallel_mode == 'local':
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} mpiexec -n {self.n_procs} ./wrf.exe > wrf.out", shell=True)
        else:
            subprocess.check_call(f"singularity exec -H {os.getcwd()} --bind {self.working_directory} {self.config_dict['wrf_img_path']} ./wrf.exe > wrf.out", shell=True)

    def extract_outputs(self, out_path, domains=None, variables_list=['U', 'V', 'TKE_PBL', 'POWER'], max_z_levels=None, t_freq=6, keep_staggered=True):
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
                if name == 'bottom_top_stag':
                    dout.createDimension(name, max_z_levels+1)
                elif name == 'bottom_top':
                    dout.createDimension(name, max_z_levels)
                elif name == 'Time':
                    dout.createDimension(name, data.dimensions['Time'].size//t_freq)
                else:
                    dout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            for name, variable in din.variables.items():
                if (name in variables_list) or (name in ['XLAT', 'XLONG']):
                    data = din[name][:]
                    
                    # Truncate z dimension
                    if 'bottom_top' in variable.dimensions:
                        z_dim = variable.dimensions.index('bottom_top')
                        data = data.take(range(max_z_levels), axis=z_dim)
                    if 'bottom_top_stag' in variable.dimensions:
                        z_stag_dim = variable.dimensions.index('bottom_top_stag')
                        data = data.take(range(max_z_levels+1), axis=z_stag_dim)
                    
                    # Reduce time frequency
                    if 'Time' in variable.dimensions:
                        t_dim = variable.dimensions.index('Time')
                        data = data.take(range(0, data.shape[t_dim], t_freq), axis=t_dim)

                        # Add timestamp variable
                        if 'timestamp' not in dout.variables.keys():
                            timestamps = range(0, data.shape[t_dim], t_freq)
                            dout.createVariable('timestamp', int, 'Time')
                            dout['timestamp'][:] = timestamps
                    
                    contains_staggered_dim = 'bottom_top_stag' in variable.dimensions or \
                                             'west_east_stag' in variable.dimensions or \
                                             'south_north_stag' in variable.dimensions
                    
                    if (not contains_staggered_dim) or keep_staggered:
                        # Create variable
                        x = dout.createVariable(name, variable.datatype, variable.dimensions)
                        dout[name][:] = data
                        dout[name].setncatts(din[name].__dict__)
                    
                    if contains_staggered_dim:
                        data_destag = data
                        dims_destag = list(variable.dimensions)
                        # De-stagger as necessary
                        if 'west_east_stag' in variable.dimensions:
                            x_stag_dim = variable.dimensions.index('west_east_stag')
                            data_destag = 0.5 * (data_destag.take(range(0, data_destag.shape[x_stag_dim]-1), axis=x_stag_dim) +
                                                data_destag.take(range(1, data_destag.shape[x_stag_dim]), axis=x_stag_dim))
                            dims_destag[x_stag_dim] = 'west_east'
                        if 'south_north_stag' in variable.dimensions:
                            y_stag_dim = variable.dimensions.index('south_north_stag')
                            data_destag = 0.5 * (data_destag.take(range(0, data_destag.shape[y_stag_dim]-1), axis=y_stag_dim) +
                                                data_destag.take(range(1, data_destag.shape[y_stag_dim]), axis=y_stag_dim))
                            dims_destag[y_stag_dim] = 'south_north'
                        if 'bottom_top_stag' in variable.dimensions:
                            z_stag_dim = variable.dimensions.index('bottom_top_stag')
                            data_destag = 0.5 * (data_destag.take(range(0, data_destag.shape[z_stag_dim]-1), axis=z_stag_dim) +
                                                data_destag.take(range(1, data_destag.shape[z_stag_dim]), axis=z_stag_dim))
                            dims_destag[z_stag_dim] = 'bottom_top'
                        
                        # Create variable
                        name_destag = name + '_destag'
                        x = dout.createVariable(name_destag, variable.datatype, dims_destag)
                        dout[name_destag][:] = data_destag
                        dout[name_destag].setncatts(din[name].__dict__)
                        dout[name_destag].stagger = ''

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
