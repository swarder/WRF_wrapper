template = \
"""&share
 wrf_core = 'ARW',
 max_dom = {max_dom},
 start_date = '{start_date_str}','{start_date_str}',
 end_date   = '{end_date_str}','{end_date_str}',
 interval_seconds = {interval_seconds}
/

&geogrid
 parent_id         =   1,   1,
 parent_grid_ratio =   1,   {d02_grid_ratio},
 i_parent_start    =   1,  {d02_i_parent_start},
 j_parent_start    =   1,  {d02_j_parent_start},
 e_we              =  {d01_e_we}, {d02_e_we},
 e_sn              =  {d01_e_sn}, {d02_e_sn},
 geog_data_res = 'default','default',
 dx = {dx},
 dy = {dx},
 map_proj = 'lambert',
 ref_lat   =  {ref_lat},
 ref_lon   = {ref_lon},
 truelat1  =  {truelat1},
 truelat2  =  {truelat2},
 stand_lon = {stand_lon},
 geog_data_path = '{geog_data_path}'
/

&ungrib
 out_format = 'WPS',
 prefix = '{pfile_data_dir}/FILE',
/

&metgrid
 fg_name = '{pfile_data_dir}/FILE'
/
"""


template_3_domains = \
"""&share
 wrf_core = 'ARW',
 max_dom = {max_dom},
 start_date = '{start_date_str}','{start_date_str}','{start_date_str}',
 end_date   = '{end_date_str}','{end_date_str}','{end_date_str}',
 interval_seconds = {interval_seconds}
/

&geogrid
 parent_id         =   1,   1,   2,
 parent_grid_ratio =   1,   {d02_grid_ratio}, {d03_grid_ratio},
 i_parent_start    =   1,  {d02_i_parent_start}, {d03_i_parent_start},
 j_parent_start    =   1,  {d02_j_parent_start}, {d03_j_parent_start},
 e_we              =  {d01_e_we}, {d02_e_we}, {d03_e_we},
 e_sn              =  {d01_e_sn}, {d02_e_sn}, {d03_e_sn},
 geog_data_res = 'default','default','default',
 dx = {dx},
 dy = {dx},
 map_proj = 'lambert',
 ref_lat   =  {ref_lat},
 ref_lon   = {ref_lon},
 truelat1  =  {truelat1},
 truelat2  =  {truelat2},
 stand_lon = {stand_lon},
 geog_data_path = '{geog_data_path}'
/

&ungrib
 out_format = 'WPS',
 prefix = '{pfile_data_dir}/FILE',
/

&metgrid
 fg_name = '{pfile_data_dir}/FILE'
/
"""
