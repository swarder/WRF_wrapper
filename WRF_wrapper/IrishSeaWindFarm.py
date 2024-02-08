import pandas as pd
import numpy as np
import shutil
import os
import pkg_resources
import fiona
import utm
from shapely.geometry import shape, Point, Polygon
from shapely.ops import unary_union

from . import NorthSeaWindFarm

DATA_PATH = pkg_resources.resource_filename('WRF_wrapper', 'data/')

lease_area_df = pd.read_csv(os.path.join(DATA_PATH, 'lease_areas/irish_sea/lease_areas_processed.csv'))
lease_area_polygons_file = fiona.open(os.path.join(DATA_PATH, 'lease_areas/irish_sea/lease_areas.shp'))

class IrishSeaWindFarm(NorthSeaWindFarm.NorthSeaWindFarm):
    @classmethod
    def from_lease_area(cls, lease_area_id, standard_turbine_powers=None, standard_turbine_type_ids=None, lease_area_df=lease_area_df, lease_area_polygons_file=lease_area_polygons_file, stagger=False, use_bug_fix=False, verbose=True):
        """
        Same as NorthSeaWindFarm, but update the default lease_area_polygons_file and lease_area_df
        """
        return super().from_lease_area(lease_area_id, standard_turbine_powers, standard_turbine_type_ids, lease_area_df, lease_area_polygons_file, stagger, use_bug_fix, verbose)