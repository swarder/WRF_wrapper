from WRF_wrapper import NorthSeaWindFarm
from matplotlib import pyplot as plt

lease_area_df = NorthSeaWindFarm.lease_area_df.copy()
lease_area_df = lease_area_df[lease_area_df.Status == 'Operational']
lease_area_df['Composite'].fillna(lease_area_df['ID'], inplace=True)
lease_area_df.drop_duplicates(['Composite'], inplace=True)
farms = []
for i, row in lease_area_df.iterrows():
    print(row.Name)
    farm = NorthSeaWindFarm.NorthSeaWindFarm.from_lease_area(row.Name, include_composites='Operational')
    farms.append(farm)
    plt.scatter(farm.farm_df.lon, farm.farm_df.lat)
plt.show()
