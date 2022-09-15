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
 parent_grid_ratio =   1,   9,
 i_parent_start    =   1,  {d02_i_parent_start},
 j_parent_start    =   1,  {d02_j_parent_start},
 e_we              =  202, {d02_e_we},
 e_sn              =  202, {d02_e_sn},
 geog_data_res = 'default','default',
 dx = 16668,
 dy = 16668,
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
