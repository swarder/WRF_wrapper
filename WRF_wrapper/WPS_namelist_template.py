template = \
"""&share
 wrf_core = 'ARW',
 max_dom = 2,
 start_date = '{start_date_str}','{start_date_str}',
 end_date   = '{end_date_str}','{end_date_str}',
 interval_seconds = {interval_seconds}
/

&geogrid
 parent_id         =   0,   1,
 parent_grid_ratio =   1,   9,
 i_parent_start    =   50,  100,
 j_parent_start    =   125,  92,
 e_we              =  202, 73,
 e_sn              =  202, 154,
 geog_data_res = 'default','default',
 dx = 16668,
 dy = 16668,
 map_proj = 'lambert',
 ref_lat   =  61,
 ref_lon   = 0,
 truelat1  =  60.6,
 truelat2  =  62.6,
 stand_lon = -1.0,
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
