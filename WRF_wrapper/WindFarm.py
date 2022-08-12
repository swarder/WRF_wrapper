import pandas as pd
import shutil
import os
import pkg_resources

DATA_PATH = pkg_resources.resource_filename('WRF_wrapper', 'data/')

class WindFarm:
    """
    Object representing a wind farm
    """
    def __init__(self, farm_df):
        # farm_df should be a pandas dataframe with columns lon, lat, type_id
        self.farm_df = farm_df

    @classmethod
    def from_template(cls, template_file=None, template_id=None, farm_lon=0, farm_lat=0, type_id=6):
        """
        Create WindFarm object from a template file, specifying a new farm centre and turbine type
        User can either specify a filename containing the template layout (template_file)
        or template_id, which corresponds to a pre-defined template
        """
        if template_file is None:
            template_file = os.path.join(DATA_PATH, f'farm_template_{template_id}.txt')
        farm_df = pd.read_csv(template_file, names=['lat', 'lon', 'type_id'], delimiter=' ')
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
            shutil.copy(os.path.join(DATA_PATH, f'wind-turbine-{type_id}.tbl'), wrf_run_dir)
