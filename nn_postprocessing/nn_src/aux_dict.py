"""Define dictionary with auxiliary variables.
"""
from collections import OrderedDict

aux_dict = OrderedDict()
aux_dict['data_aux_geo_interpolated.nc'] = [
    'orog',
    'station_alt',
    'station_lat',
    'station_lon',
]
aux_dict['data_aux_pl500_interpolated_00UTC.nc'] = [
    'u_pl500_fc',
    'v_pl500_fc',
    'gh_pl500_fc',
]
aux_dict['data_aux_pl850_interpolated_00UTC.nc'] = [
    'u_pl850_fc',
    'v_pl850_fc',
    'q_pl850_fc',
]
aux_dict['data_aux_surface_interpolated_00UTC.nc'] = [
    'cape_fc',
    'sp_fc',
    'tcc_fc',
]
aux_dict['data_aux_surface_more_interpolated_part1_00UTC.nc'] = [
    'sshf_fc',
    'slhf_fc',
    'u10_fc',
    'v10_fc',
]
aux_dict['data_aux_surface_more_interpolated_part2_00UTC.nc'] = [
    'ssr_fc',
    'str_fc',
    'd2m_fc',
    'sm_fc',
]