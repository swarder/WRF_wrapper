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

lease_area_df = pd.read_csv(os.path.join(DATA_PATH, 'lease_areas/north_sea/lease_areas.csv'))
lease_area_polygons_file = fiona.open(os.path.join(DATA_PATH, 'lease_areas/north_sea/lease_areas.shp'))
lease_area_name_to_id = {}
for p in lease_area_polygons_file:
    id = int(p['id'])
    name = lease_area_df.loc[lease_area_df.ID == id, 'Name'].values[0]
    lease_area_name_to_id[name] = id

def convert_geometry_to_utm(geometry_latlon):
    coords_latlon = np.squeeze(np.array(geometry_latlon['coordinates']))
    x_utm, y_utm, zone_number, zone_letter = utm.from_latlon(coords_latlon[:,1], coords_latlon[:,0])
    coords_utm = np.stack([x_utm, y_utm], axis=-1)
    geometry_utm = {'type': 'Polygon',
                    'coordinates': [[tuple([x, y]) for x, y in zip(coords_utm[:,0], coords_utm[:,1])]],
                    }
    return geometry_utm

canonical_turbine_powers = np.array([3.6, 6, 10, 15])
canonical_turbine_type_ids = ['plausible-3.6', 'plausible-6.0', '10', '15']

class NorthSeaWindFarm(WindFarm.WindFarm):
    @classmethod
    def from_lease_area(cls, lease_area_name, include_composites=False, standard_turbine_powers=None, standard_turbine_type_ids=None):
        """
        Create WindFarm object within specified lease_area_name, using turbines of specified type_id, and specified IC
        Maximum turbine spacing is found which is sufficient to achieve the specified IC within the lease area
        """
        if standard_turbine_powers is None:
            standard_turbine_powers = canonical_turbine_powers
            standard_turbine_type_ids = canonical_turbine_type_ids

        lease_area_id = lease_area_name_to_id[lease_area_name]
        lease_area_df['Composite'].fillna(lease_area_df['ID'], inplace=True)
        composite_id = lease_area_df.loc[lease_area_df.Name==lease_area_name, 'Composite'].values[0]
        if include_composites:
            polygons = []
            installed_capacity = 0
            total_num_turbines = 0
            composite_farms = lease_area_df.loc[lease_area_df.Composite==composite_id]
            if isinstance(include_composites, str):
                composite_farms = composite_farms.loc[composite_farms.Status == include_composites]
            for i, row in composite_farms.iterrows():
                lease_area_id = row.ID
                lease_geometry_latlon = lease_area_polygons_file[lease_area_id]['geometry']
                polygons.append(shape(lease_geometry_latlon))
                installed_capacity += float(row['Total IC'])
                total_num_turbines += float(row['Num turbines'])
            combined_polygon = unary_union(polygons)
            lease_geometry_latlon = {'type': 'Polygon',
                                     'coordinates': [combined_polygon.exterior.coords]}
            turbine_target_power = installed_capacity / total_num_turbines
        else:
            lease_geometry_latlon = lease_area_polygons_file[lease_area_id]['geometry']
            turbine_target_power = lease_area_df.loc[lease_area_df.Name == lease_area_name, 'Turbine capacity'].values[0]
            installed_capacity = float(lease_area_df.loc[lease_area_df.Name == lease_area_name, 'Total IC'].values[0])

        lease_geometry_utm = convert_geometry_to_utm(lease_geometry_latlon)
        lease_area_area = shape(lease_geometry_utm).area

        type_id = standard_turbine_type_ids[np.argmin(np.abs(standard_turbine_powers - turbine_target_power))]

        turbine_power = WindFarm.WindTurbine.from_type_id(type_id).nominal_power
        num_required_turbines = installed_capacity / turbine_power
        approx_turbine_separation = np.sqrt(lease_area_area / num_required_turbines)

        farm = cls.from_geometry(lease_geometry_latlon, 'grid', type_id, approx_turbine_separation, 90)
        n = farm.farm_df.shape[0]

        # Likely to be too few turbines initially
        while n < num_required_turbines:
            print(approx_turbine_separation)
            approx_turbine_separation -= 25
            farm = cls.from_geometry(lease_geometry_latlon, 'grid', type_id, approx_turbine_separation, 90)
            n = farm.farm_df.shape[0]
        # If too many turbines, sample down to required number
        farm.farm_df = farm.farm_df.sample(int(num_required_turbines), random_state=0)

        return farm
