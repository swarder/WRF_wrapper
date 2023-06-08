import pandas as pd
import numpy as np
import shutil
import os
import pkg_resources
import fiona
import utm
from shapely.geometry import shape, Point, Polygon, MultiPoint
import py_wake.wind_turbines

DATA_PATH = pkg_resources.resource_filename('WRF_wrapper', 'data/')

def parse_turbine_spacing(turbine_spacing, wind_turbine, return_density=False):
    if not isinstance(turbine_spacing, list):
        turbine_spacing = [turbine_spacing, turbine_spacing]
    parsed_spacings = []
    for ts in turbine_spacing:
        if isinstance(ts, str):
            assert ts.endswith('D')
            parsed_spacings.append(float(ts[:-1]) * wind_turbine.diameter)
        else:
            parsed_spacings.append(ts)
    return parsed_spacings

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

    def generate_power_curve_func(self):
        return lambda x: np.interp(x, self.power_curve.wind_speed, self.power_curve.power)

    @classmethod
    def from_type_id(cls, type_id):
        tbl_file = os.path.join(DATA_PATH, f'power_curves/wind-turbine-{type_id}.tbl')
        return cls(tbl_file)

    def to_pywake(self):
        """
        Export to pywake WindTurbine
        """
        ct_func = py_wake.wind_turbines.power_ct_functions.PowerCtTabular(
                                                self.power_curve.wind_speed,
                                                self.power_curve.power,
                                                'kW',
                                                self.power_curve.thrust_coeff,
                                                )
        my_wt = py_wake.wind_turbines.WindTurbine(
                            name='from_WRF_wrapper_WindTurbine',
                            diameter=self.diameter,
                            hub_height=self.height,
                            powerCtFunction=ct_func,
                            )
        return my_wt

    def to_floris(self):
        """
        Export to FLORIS dict format
        """
        wind_speed = self.power_curve.wind_speed
        thrust = self.power_curve.thrust_coeff
        power = self.power_curve.power

        rotor_area = np.pi * (0.5 * self.diameter) ** 2
        power_coeff = power / (0.5 * 1.225 * rotor_area * wind_speed ** 3)

        wind_speed = wind_speed.tolist()
        thrust = thrust.tolist()
        power_coeff = power_coeff.tolist()

        wind_speed = [0] + wind_speed + [wind_speed[-1]+0.1, 50]
        thrust = [thrust[0]] + thrust + [0, 0]
        power_coeff = [0] + power_coeff + [0, 0]

        power_curve_dict = {'power': power_coeff,
                            'thrust': thrust,
                            'wind_speed': wind_speed,
                           }
        my_wt = {'turbine_type': 'from_WRF_wrapper_WindTurbine',
                 'generator_efficiency': 1.0,
                 'hub_height': self.height,
                 'pP': 1.88,
                 'pT': 1.88,
                 'rotor_diameter': self.diameter,
                 'TSR': 8.0,
                 'ref_density_cp_ct': 1.225,
                 'power_thrust_table': power_curve_dict,
                }
        return my_wt

class WindFarm:
    """
    Object representing a wind farm
    """
    def __init__(self, farm_df, name=None):
        # farm_df should be a pandas dataframe with columns lon, lat, type_id
        self.farm_df = farm_df
        self.boundary_polygon = self.generate_convex_hull()
        self.name = name

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
            turbine_spacing = turbine.diameter * float(turbine_spacing[:-1])

        # Generate simple turbine array
        num_turbines = int(np.ceil(installed_capacity / turbine.nominal_power))
        turbine_array_size = np.ceil(np.sqrt(num_turbines))
        yy, xx = np.meshgrid(np.arange(turbine_array_size), np.arange(turbine_array_size))
        turbine_xy = np.stack([xx.flatten(), yy.flatten()]).T[:num_turbines] * turbine_spacing
        turbine_xy -= turbine_xy.mean(axis=0)

        # Convert farm latlon to utm
        farm_x, farm_y, zone_number, zone_letter = utm.from_latlon(farm_lat, farm_lon)
        turbine_xy += np.array([farm_x, farm_y])[None,:]

        # Convert turbine locations to latlon
        turbine_lat, turbine_lon = utm.to_latlon(turbine_xy[:,0], turbine_xy[:,1], zone_number, zone_letter, strict=False)
        turbine_latlon = np.stack([turbine_lat, turbine_lon], axis=-1)

        farm_df = pd.DataFrame(turbine_latlon, columns=['lat', 'lon'])
        farm_df['type_id'] = type_id

        return cls(farm_df)

    @classmethod
    def from_geometry(cls, geometry, layout, type_id, turbine_spacing, grid_alignment):
        """
        Generate WindFarm object from a given geometry
        """
        if not layout == 'grid':
            raise NotImplemtedError
        turbine_spacing = parse_turbine_spacing(turbine_spacing, WindTurbine.from_type_id(type_id))

        lease_area_points = np.array(geometry['coordinates'][0])
        x_farm, y_farm, zone_num, zone_let = utm.from_latlon(lease_area_points[:,1].mean(), lease_area_points[:,0].mean())

        lease_area_points_utm = [tuple(utm.from_latlon(p[1], p[0], force_zone_number=zone_num, force_zone_letter=zone_let)[:2]) for p in lease_area_points]
        lease_geometry_utm = {'type': 'Polygon', 'coordinates': [lease_area_points_utm]}
        lease_area_area = shape(lease_geometry_utm).area
        n = np.sqrt(lease_area_area) / min(turbine_spacing)

        # Generate grid around that point, and rotate
        yy, xx = np.meshgrid(np.arange(-5*n, 5*n), np.arange(-5*n, 5*n))
        xx = xx * turbine_spacing[0]
        yy = yy * turbine_spacing[1]
        turbine_xy_complex = (xx + 1j * yy) * np.exp(1j * (90 - grid_alignment) * np.pi/180)
        turbine_xx = np.real(turbine_xy_complex) + x_farm
        turbine_yy = np.imag(turbine_xy_complex) + y_farm
        farm_df = pd.DataFrame({'x': turbine_xx.flatten(), 'y': turbine_yy.flatten()})

        # Initial filtering to speed things up
        lease_area_points_utm = np.array(lease_area_points_utm)
        farm_df = farm_df.loc[(farm_df.x >= lease_area_points_utm[:,0].min()) & (farm_df.x <= lease_area_points_utm[:,0].max())]
        farm_df = farm_df.loc[(farm_df.y >= lease_area_points_utm[:,1].min()) & (farm_df.y <= lease_area_points_utm[:,1].max())]

        # Convert back to lat-lon coords
        def convert(xy):
            return pd.Series(utm.to_latlon(xy['x'], xy['y'], zone_num, zone_let), index=['lat', 'lon'])
        farm_df = farm_df.apply(convert, axis=1)

        if farm_df.shape[0] == 0:
            return cls(farm_df)

        # Only include turbines in the lease area
        def is_in_lease_area(latlon):
            return shape(Point(latlon['lon'], latlon['lat'])).within(shape(geometry))
        farm_df = farm_df.loc[farm_df.apply(is_in_lease_area, axis=1)]

        farm_df['type_id'] = type_id
        farm_df.reset_index(drop=True, inplace=True)

        farm = cls(farm_df)
        farm.boundary_polygon = Polygon(np.array(geometry['coordinates'])[0])

        return farm

    @classmethod
    def from_lease_area(cls, lease_area_name, layout='grid', type_id=6, turbine_spacing=['10D', '4D'], grid_alignment=90):
        """
        Create WindFarm object within specified lease_area_name, using turbines of specified type_id, and specified turbine_spacing
        turbine_spacing can be a list, specifying spacing in two directions. First direction is in direction specified by grid_alignment, second direction is perpendicular
        """
        lease_polygon = fiona.open(os.path.join(DATA_PATH, f'lease_areas/{lease_area_name}/lease_area.shp')).next()
        geometry = lease_polygon['geometry']
        return cls.from_geometry(geometry, layout, type_id, turbine_spacing, grid_alignment)

    @classmethod
    def from_polygon(cls, centroid_latlon, polar_angles, polar_radii, layout='grid', type_id=6, turbine_spacing=['10D', '4D'], grid_alignment=90):
        """
        Create WindFarm object within a polygon at specified centroid, and with verticies at the specified polar coordinates
        """
        centroid_x, centroid_y, zone_num, zone_let = utm.from_latlon(centroid_latlon[0], centroid_latlon[1])
        centroid_xy = np.array([centroid_x, centroid_y])
        vertices_utm = np.array([(r * np.cos(phi*np.pi/180), r * np.sin(phi*np.pi/180)) for r, phi in zip(polar_radii, polar_angles)]) + centroid_xy[None,:]
        vertices_latlon = np.stack(utm.to_latlon(vertices_utm[:,0], vertices_utm[:,1], zone_num, zone_let), axis=-1)
        geometry = {'type': 'Polygon', 'coordinates': [vertices_latlon[:,::-1]]}
        return cls.from_geometry(geometry, layout, type_id, turbine_spacing, grid_alignment)

    @classmethod
    def from_random_polygon(cls, centroid_latlon, approx_area=None, approx_ic=None, num_vertices=None, layout='grid', type_id=6, turbine_spacing=['10D', '4D'], grid_alignment=90, random_seed=None):
        """
        Create WindFarm object within a polygon at specified centroid, and with verticies at the specified polar coordinates
        """
        if random_seed is not None:
            np.random.seed(random_seed)
        if num_vertices is None:
            # Either 3 or 4 sides
            num_vertices = np.random.randint(3, 5)
        if approx_area is None:
            assert approx_ic is not None
            turbine = WindTurbine.from_type_id(type_id)

            turbine_spacing = parse_turbine_spacing(turbine_spacing, turbine, return_density=True)
            density = 1 / np.prod(turbine_spacing)
            turbine_power = turbine.nominal_power
            power_density = turbine_power * turbine_density
            approx_area = approx_ic / power_density

        centroid_x, centroid_y, zone_num, zone_let = utm.from_latlon(centroid_latlon[0], centroid_latlon[1])
        centroid_xy = np.array([centroid_x, centroid_y])
        # Random polygon vertex radii and polar angles
        polar_radii = np.random.uniform(0.5, 1.5, num_vertices)
        polar_angles = [np.random.uniform(i*360/num_vertices, (i+1)*360/num_vertices) for i in range(num_vertices)]

        # Generate initial vertices
        vertices_utm = np.array([(r * np.cos(phi*np.pi/180), r * np.sin(phi*np.pi/180)) for r, phi in zip(polar_radii, polar_angles)])

        # Calculate current area
        polygon_area = 0.5 * np.abs(np.dot(vertices_utm[:,0],np.roll(vertices_utm[:,1],1))-np.dot(vertices_utm[:,1],np.roll(vertices_utm[:,0],1)))

        # Rescale to match target area, and add required centroid location
        rescale_factor = np.sqrt(approx_area / polygon_area)
        vertices_utm = vertices_utm * rescale_factor + centroid_xy[None,:]

        # Convert to latlon coordinates
        vertices_latlon = np.stack(utm.to_latlon(vertices_utm[:,0], vertices_utm[:,1], zone_num, zone_let), axis=-1)

        # Construct polygon
        geometry = {'type': 'Polygon', 'coordinates': [vertices_latlon[:,::-1]]}

        farm = cls.from_geometry(geometry, layout, type_id, turbine_spacing, grid_alignment)

        # Rotate by random angle
        rotation_angle = np.random.uniform(0, 360)
        return farm.duplicate_with_rotation(rotation_angle)

    def __add__(self, o):
        """
        Overload addition, to enable combining multiple farms together into a single file
        """
        if o == 0:
            return self
        else:
            summed_farm = WindFarm(pd.concat([self.farm_df, o.farm_df], ignore_index=True))
            summed_farm.name = self.name
            return summed_farm

    def __radd__(self, o):
        """
        To enable use of built-in sum function on list of WindFarm objects
        """
        return self.__add__(o)

    def __eq__(self, other):
        """
        Overload equality, to enable comparison of two farms
        """
        if isinstance(other, WindFarm):
            comp_cols = ['lon', 'lat', 'type_id']
            return self.farm_df[comp_cols].equals(other.farm_df[comp_cols])
        return False

    def save(self, wrf_run_dir):
        """
        Save farm to CSV
        """
        # Assign a unique integer to each turbine type (which up to this point are not necessarily integers)
        used_turbine_types = self.farm_df.type_id.unique()
        new_turbine_type_ids = [100 + i for i in range(1, len(used_turbine_types)+1)]
        type_id_replace_dict = dict(zip(used_turbine_types, new_turbine_type_ids))

        farm_df_to_save = self.farm_df.copy()
        farm_df_to_save.type_id.replace(type_id_replace_dict, inplace=True)
        farm_df_to_save.to_csv(os.path.join(wrf_run_dir, 'windturbines.txt'), index=False, sep=' ', header=False)

        # Copy relevant .tbl (power curve) files into wrf run directory
        for old_id, new_id in zip(used_turbine_types, new_turbine_type_ids):
            src = os.path.join(DATA_PATH, f'power_curves/wind-turbine-{old_id}.tbl')
            dst = os.path.join(wrf_run_dir, f'wind-turbine-{new_id}.tbl')
            shutil.copy(src, dst)

    def generate_convex_hull(self, crs='latlon'):
        """
        Generate convex hull based on turbine locations (in lat-lon)
        """
        farm_latlons = self.farm_df[['lon', 'lat']].values
        if crs == 'latlon':
            vertices = farm_latlons
        elif crs.startswith('UTM'):
            zone_num = int(crs[3:])
            x, y, _, _ = utm.from_latlon(self.farm_df['lat'],
                                         self.farm_df['lon'],
                                         force_zone_number=zone_num,
                                         force_zone_letter='N')
            vertices = np.stack([x, y], axis=1)
        mpt = MultiPoint([Point(ll) for ll in vertices])
        return mpt.convex_hull

    def duplicate_with_perturbation(self, ptb, crs='latlon'):
        """
        Return a copy of self, with a perturbation added to farm location
        if crs == 'latlon', perturbation coordates are latlon
        Otherwise, assume crs specifies a utm zone
        """
        if np.allclose(ptb, [0, 0]):
            return WindFarm(self.farm_df, name=self.name)
        new_farm_df = self.farm_df.copy()
        if crs == 'latlon':
            new_farm_df['lon'] += ptb[0]
            new_farm_df['lat'] += ptb[1]
        elif crs.startswith('UTM'):
            zone_num = int(crs[3:])
            x, y, _, _ = utm.from_latlon(new_farm_df['lat'],
                                         new_farm_df['lon'],
                                         force_zone_number=zone_num,
                                         force_zone_letter='N')
            x_ptb = x + ptb[0]
            y_ptb = y + ptb[1]
            new_lat, new_lon = utm.to_latlon(x_ptb, y_ptb,
                                             zone_num, 'N', strict=False)
            new_farm_df['lon'] = new_lon
            new_farm_df['lat'] = new_lat
        else:
            raise NotImplemtedError
        return WindFarm(new_farm_df, name=self.name)

    def duplicate_with_rotation(self, angle):
        """
        Rotate the farm about its centroid, by specified angle
        """
        # Project to UTM
        _, _, zone_num, zone_letter = utm.from_latlon(self.farm_df['lat'].mean(),
                                                      self.farm_df['lon'].mean())
        x, y, _, _ = utm.from_latlon(self.farm_df['lat'],
                                     self.farm_df['lon'],
                                     zone_num, zone_letter)
        # Subtract centroid
        x_centroid = x.mean()
        y_centroid = y.mean()
        x -= x_centroid
        y -= y_centroid

        # Rotate
        complex_positions = x + 1j*y
        complex_rotation = np.exp(1j*angle*np.pi/180)
        complex_positions_rotated = complex_positions * complex_rotation

        # New relative locations
        x_new = np.real(complex_positions_rotated)
        y_new = np.imag(complex_positions_rotated)

        # New absolute positions
        x_new += x_centroid
        y_new += y_centroid

        # Convert back to latlon
        lat_new, lon_new = utm.to_latlon(x_new,
                                         y_new,
                                         zone_num,
                                         zone_letter)

        # Create new farm
        new_farm_df = self.farm_df.copy()
        new_farm_df['lat'] = lat_new
        new_farm_df['lon'] = lon_new

        return WindFarm(new_farm_df, name=self.name)
