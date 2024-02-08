import pandas as pd
import numpy as np
import shutil
import os
import pkg_resources
import fiona
import utm
from shapely.geometry import shape, Point, Polygon
from shapely.ops import unary_union

from . import WindFarm

DATA_PATH = pkg_resources.resource_filename('WRF_wrapper', 'data/')

lease_area_df = pd.read_csv(os.path.join(DATA_PATH, 'lease_areas/north_sea/lease_areas_processed.csv'))
lease_area_polygons_file = fiona.open(os.path.join(DATA_PATH, 'lease_areas/north_sea/lease_areas.shp'))

canonical_turbine_powers = np.array([3.6, 6, 10, 15])
canonical_turbine_type_ids = ['plausible-3.6', 'plausible-6.0', '10', '15']

class NorthSeaWindFarm(WindFarm.WindFarm):
    @classmethod
    def from_lease_area(cls, lease_area_id, standard_turbine_powers=None, standard_turbine_type_ids=None, lease_area_df=lease_area_df, lease_area_polygons_file=lease_area_polygons_file, stagger=False, use_bug_fix=False, verbose=True):
        """
        Create WindFarm object within specified lease_area_id, using turbines of specified type_id, and specified IC
        Maximum turbine spacing is found which is sufficient to achieve the specified IC within the lease area
        """
        if standard_turbine_powers is None:
            standard_turbine_powers = canonical_turbine_powers
            standard_turbine_type_ids = canonical_turbine_type_ids
        standard_turbine_powers = np.array(standard_turbine_powers)

        # This class is re-used for other regional farm lease areas e.g. IrishSeaWindFarm
        # For these other areas, the lease_area_polygons_file and lease_area_df have different IDs
        # lease_area_polygons_file is always indexed from 0, whereas lease_area_df could start anywhere
        id_to_ix = dict([
            [lease_area_polygons_file[k]['properties']['id'], k] for k in range(len(lease_area_polygons_file))
        ])
        if lease_area_id != id_to_ix[lease_area_id]:
            print('Warning: lease_area_polygons_file and lease_area_df have different IDs')
        lease_geometry_latlon = lease_area_polygons_file[id_to_ix[lease_area_id]]['geometry']

        lease_area_row = lease_area_df.loc[lease_area_df.ID == lease_area_id].iloc[0]
        lease_area_area = lease_area_row.area * 1e6
        installed_capacity = lease_area_row['IC_numeric']
        num_turbines = lease_area_row['num_turbines_numeric']
        if np.isnan(num_turbines) or np.isnan(installed_capacity):
            if verbose:
                print('No num_turbines info, assuming 10 MW turbines')
            turbine_target_power = 10
        else:
            turbine_target_power = installed_capacity / num_turbines
        type_id = standard_turbine_type_ids[np.argmin(np.abs(standard_turbine_powers - turbine_target_power))]
        turbine_power = WindFarm.WindTurbine.from_type_id(type_id).nominal_power
        if np.isnan(installed_capacity):
            print('No IC info, assuming num_turbines and standard_turbine_powers provided are correct')
            num_required_turbines = num_turbines
        else:
            num_required_turbines = installed_capacity / turbine_power
        approx_turbine_separation = np.sqrt(lease_area_area / num_required_turbines)

        farm = cls.from_geometry(
            lease_geometry_latlon,
            'grid',
            type_id,
            approx_turbine_separation,
            90,
            stagger=stagger,
            use_bug_fix=use_bug_fix
            )
        n = farm.farm_df.shape[0]

        # Likely to be too few turbines initially
        while n < num_required_turbines:
            #print(approx_turbine_separation)
            approx_turbine_separation -= 25
            farm = cls.from_geometry(
                lease_geometry_latlon,
                'grid',
                type_id,
                approx_turbine_separation,
                90,
                stagger=stagger,
                use_bug_fix=use_bug_fix
                )
            n = farm.farm_df.shape[0]
        # If too many turbines, sample down to required number
        farm.farm_df = farm.farm_df.sample(int(num_required_turbines), random_state=0)
        farm.name = lease_area_row['Name']

        return farm
