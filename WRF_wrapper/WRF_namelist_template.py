template = \
"""&time_control
    start_year = {start_date.year}, {start_date.year},
    start_month = {start_date.month}, {start_date.month},
    start_day = {start_date.day}, {start_date.day},
    start_hour = {start_date.hour}, {start_date.hour},
    start_minute = {start_date.minute}, {start_date.minute},
    start_second = 0, 0,
    end_year = {end_date.year}, {end_date.year},
    end_month = {end_date.month}, {end_date.month},
    end_day = {end_date.day}, {end_date.day},
    end_hour = {end_date.hour}, {end_date.hour},
    end_minute = {end_date.minute}, {end_date.minute},
    end_second = 0, 0,
    interval_seconds = {interval_seconds},
    input_from_file = .true., .true.,
    history_interval =  10, 10,
    frames_per_outfile = 1000, 1000,
    restart = .false.,
    restart_interval = 1080,
    io_form_history = 2,
    io_form_restart = 2,
    io_form_input = 2,
    io_form_boundary = 2,
    debug_level = 0,
    ignore_iofields_warning = .true.,
    auxinput4_interval = 360, 360,
    auxinput4_inname = 'wrflowinp_d<domain>',
    io_form_auxinput4 = 2,
    nwp_diagnostics = 0,
    history_outname = 'wrfout_d<domain>',
    iofields_filename = 'io_file_d01.txt', 'io_file_d02.txt'
/

&domains
    time_step = {time_step},
    time_step_fract_num = 0,
    time_step_fract_den = 2,
    max_dom = {max_dom},
    e_we = {d01_e_we}, {d02_e_we},
    e_sn = {d01_e_sn}, {d02_e_sn},
    e_vert = 62, 62
    p_top_requested = 5000,
    num_metgrid_levels = 33,
    num_metgrid_soil_levels = 4,
    dx = {dx}, {dx_d02},
    dy = {dx}, {dx_d02},
    grid_id = 1, 2,
    parent_id = 0, 1,
    i_parent_start = 1, {d02_i_parent_start},
    j_parent_start = 1, {d02_j_parent_start},
    parent_grid_ratio = 1, {d02_grid_ratio},
    parent_time_step_ratio = 1, {d02_dt_ratio},
    max_ts_locs = 12,
    ts_buf_size = 960,
    max_ts_level = 15,
    feedback = 0,
    smooth_option = 0,
    eta_levels = 1.0, 0.998621, 0.997244, 0.995868, 0.994495, 0.993123,
                 0.991753, 0.990385, 0.989018, 0.987653, 0.986291, 0.984929,
                 0.98357, 0.982213, 0.980857, 0.979503, 0.978151, 0.9768,
                 0.975451, 0.974104, 0.972759, 0.971415, 0.970073, 0.968732,
                 0.967393, 0.96472, 0.96112, 0.9571, 0.95237, 0.94668, 0.9398,
                 0.93147, 0.9214, 0.90924, 0.89461, 0.87705, 0.85602, 0.83095,
                 0.80121, 0.76612, 0.72503, 0.67736, 0.62293, 0.56264, 0.49949,
                 0.43795, 0.38119, 0.32951, 0.28247, 0.23965, 0.20067, 0.16519,
                 0.13289, 0.1055587, 0.0830446, 0.063567, 0.0442697, 0.0280384,
                 0.0152337, 0.00892152, 0.00297004, 0.0,
/

&physics
    mp_physics = 8, 8,
    ra_lw_physics = 4, 4,
    ra_sw_physics = 4, 4,
    radt = 18, 6,
    sf_surface_physics = 2, 2,
    sf_sfclay_physics = 5, 5,
    bl_pbl_physics = 5, 5,
    bldt = 0, 0,
    cu_physics = 1, 0,
    cudt = 5, 0,
    bl_mynn_tkebudget = 0, 1,
    swint_opt = 1,
    isfflx = 1,
    sst_update = 1,
    ifsnow = 0,
    icloud = 1,
    surface_input_source = 1,
    num_land_cat = 21,
    num_soil_layers = 4,
    isftcflx = 0,
    usemonalb = .true.,
    maxiens = 1,
    maxens = 3,
    maxens2 = 3,
    maxens3 = 16,
    ensdim = 144,
    sf_urban_physics = 0, 0,
    bl_mynn_tkeadvect = .true., .true.,
    windfarm_opt = 0, {wind_turbines_int},
    windfarm_tke_factor = 0.25,
/

&fdda
    grid_fdda = 0, 0,
    gfdda_inname = 'wrffdda_d<domain>',
    gfdda_end_h = 300, 0,
    gfdda_interval_m = 360, 0,
    fgdt = 0, 0,
    if_no_pbl_nudging_uv = 1, 0,
    if_no_pbl_nudging_t = 1, 0,
    if_no_pbl_nudging_ph = 1, 0,
    if_zfac_uv = 1, 0,
    k_zfac_uv = 25, 0,
    if_zfac_t = 1, 0,
    k_zfac_t = 25, 0,
    if_zfac_q = 1, 0,
    k_zfac_q = 25, 0,
    guv = 0.0003, 7.5e-05,
    gt = 0.0003, 7.5e-05,
    gq = 0.0003, 7.5e-05,
    xwavenum = 12,
    ywavenum = 9,
    if_ramping = 0,
    dtramp_min = 60.0,
    io_form_gfdda = 2,
/

&dynamics
    w_damping = 1,
    diff_opt = 1,
    km_opt = 4,
    diff_6th_opt = 2, 2,
    diff_6th_factor = 0.06, 0.1,
    base_temp = 290.0,
    damp_opt = 0,
    zdamp = 5000.0, 5000.0,
    dampcoef = 0.05, 0.05,
    khdif = 0, 0,
    kvdif = 0, 0,
    non_hydrostatic = .true., .true.,
    moist_adv_opt = 1, 1,
    scalar_adv_opt = 1, 1,
/

&bdy_control
    spec_bdy_width = 5,
    spec_zone = 1,
    relax_zone = 4,
    specified = .true., .false.,
    nested = .false., .true.,
/

&grib2
/

&namelist_quilt
    nio_tasks_per_group = 0,
    nio_groups = 1,
/
"""
