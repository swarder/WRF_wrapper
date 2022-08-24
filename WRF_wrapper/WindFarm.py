import pandas as pd
import numpy as np
import shutil
import os
import pkg_resources

DATA_PATH = pkg_resources.resource_filename('WRF_wrapper', 'data/')

class WindTurbine:
    """
    Object representing a wind turbine
    This class is used to parse wind-turbine-<type_id>.tbl files
    """
    def __init__(self, tbl_file):
        """
        Initialise from tbl_file
        """
        f = open(tbl_file)
        n = f.readline()
        self.height, self.diameter, self.standing_thrust_coeff, self.nominal_power = [float(x) for x in f.readline().split()]
        f.close()
        self.power_curve = pd.read_csv(tbl_file, sep=' ', skiprows=2, names=['wind_speed', 'thrust_coeff', 'power'])

        # Convert diameter to degrees
        earth_radius = 6400000
        self.diameter_degrees = self.diameter / earth_radius * 180 / np.pi

    @classmethod
    def from_type_id(cls, type_id):
        tbl_file = os.path.join(DATA_PATH, f'power_curves/wind-turbine-{type_id}.tbl')
        return cls(tbl_file)


class WindFarm:
    """
    Object representing a wind farm
    """
    def __init__(self, farm_df):
        # farm_df should be a pandas dataframe with columns lon, lat, type_id
        self.farm_df = farm_df

    @classmethod
    def from_template(cls, template_file=None, template_id=None, farm_lon=None, farm_lat=None, type_id=6):
        """
        Create WindFarm object from a template file, specifying a new farm centre and turbine type
        User can either specify a filename containing the template layout (template_file)
        or template_id, which corresponds to a pre-defined template
        """
        if template_file is None:
            template_file = os.path.join(DATA_PATH, f'templates/farm_template_{template_id}.txt')
        farm_df = pd.read_csv(template_file, names=['lat', 'lon', 'type_id'], delimiter=' ')
        if farm_lat is not None:
            farm_df['lat'] = farm_df['lat'] - farm_df['lat'].mean() + farm_lat
        if farm_lon is not None:
            farm_df['lon'] = farm_df['lon'] - farm_df['lon'].mean() + farm_lon
        farm_df['type_id'] = type_id
        return cls(farm_df)

    @classmethod
    def from_installed_capacity(cls, installed_capacity=5, layout='grid', farm_lon=0, farm_lat=0, type_id=6, turbine_spacing='5D'):
        """
        Create WindFarm object with specified installed_capacity (in MW)
        WindFarm is constructed with a standard layout using turbines of specified type_id, and specified turbine_spacing
        Array is constructed with the specified pattern until the installed capacity is met or exceeded
        """
        if not layout == 'grid':
            raise NotImplemtedError

        turbine = WindTurbine.from_type_id(type_id)

        if isinstance(turbine_spacing, str):
            # Assume of the form '7D', meaning spacing is 7 turbine diameters
            assert turbine_spacing[-1] == 'D'
            turbine_spacing = turbine.diameter_degrees * float(turbine_spacing[:-1])

        # Generate simple turbine array
        num_turbines = int(np.ceil(installed_capacity / turbine.nominal_power))
        turbine_array_size = np.ceil(np.sqrt(num_turbines))
        yy, xx = np.meshgrid(np.arange(turbine_array_size), np.arange(turbine_array_size))
        turbine_xy = np.stack([xx.flatten(), yy.flatten()]).T[:num_turbines] * turbine_spacing
        farm_df = pd.DataFrame(turbine_xy, columns=['lat', 'lon'])

        # Now centre the farm in the correct place
        farm_df['lat'] = farm_df['lat'] - farm_df['lat'].mean() + farm_lat
        farm_df['lon'] = farm_df['lon'] - farm_df['lon'].mean() + farm_lon
        farm_df['type_id'] = type_id
        return cls(farm_df)

    def __add__(self, o):
        """
        Overload addition, to enable combining multiple farms together into a single file
        """
        if o == 0:
            return self
        else:
            return WindFarm(pd.concat([self.farm_df, o.farm_df], ignore_index=True))

    def __radd__(self, o):
        """
        To enable use of built-in sum function on list of WindFarm objects
        """
        return self.__add__(o)

    def save(self, wrf_run_dir):
        """
        Save farm to CSV
        """
        self.farm_df.to_csv(os.path.join(wrf_run_dir, 'windturbines.txt'), index=False, sep=' ', header=False)

        # Copy relevant .tbl (power curve) files into wrf run directory
        used_turbine_types = self.farm_df.type_id.unique()
        for type_id in used_turbine_types:
            shutil.copy(os.path.join(DATA_PATH, f'power_curves/wind-turbine-{type_id}.tbl'), wrf_run_dir)
