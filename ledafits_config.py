#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
ledafits_config.py
==============

This files stores configuration parameters and global variables that are required across
the LEDA project. 
"""

import os
import ephem


SPEED_OF_LIGHT = 299792458

########
# PSR-DADA Settings
########

OFFSET_DELTA     = 115187712   # Bytes per dada file 
INT_TIME         = 8.33333     # Integration time (s)
N_INT_PER_FILE   = 10          # Number integrations (?)

########
# Station location - OVRO
########

(latitude, longitude, elevation) = ('37.240391', '-118.2', 1184)

ovro      = ephem.Observer()
ovro.lon  = longitude
ovro.lat  = latitude
ovro.elev = elevation

########
# Station location - LWA1
########

(latitude, longitude, elevation) = ('34.070', '-107.628', 2133.6)

lwa1      = ephem.Observer()
lwa1.lon  = longitude
lwa1.lat  = latitude
lwa1.elev = elevation

########
# LedaFits defaults
#######
CH_WIDTH          = 24e3
SUB_BW            = 2.616e6
TELESCOP          = "LWA-OVRO"
ARRNAM            = "LEDA-512"

# Default files to load to fill in FITS-IDI
fileroot = os.path.abspath(os.path.dirname(__file__))
#json_h_array_geometry  = os.path.join(fileroot, 'config/leda/h_array_geometry.json')
#json_d_array_geometry  = os.path.join(fileroot, 'config/leda/d_array_geometry.json')
#json_h_antenna         = os.path.join(fileroot, 'config/leda/h_antenna.json')
#json_d_antenna         = os.path.join(fileroot, 'config/leda/d_antenna.json')
#json_uvw_coordinates   = os.path.join(fileroot, 'config/leda/z_uvw.json')
#json_antenna_el_lens   = os.path.join(fileroot, 'config/leda/z_antenna_el_lens.json')

json_h_array_geometry  = os.path.join(fileroot, 'config/leda512/h_array_geometry.json')
json_d_array_geometry  = os.path.join(fileroot, 'config/leda512/d_array_geometry.json')
json_h_antenna         = os.path.join(fileroot, 'config/leda512/h_antenna.json')
json_d_antenna         = os.path.join(fileroot, 'config/leda512/d_antenna.json')
json_antenna_el_lens   = os.path.join(fileroot, 'config/leda512/z_antenna_el_lens.json')

json_h_array_geometry_nm  = os.path.join(fileroot, 'config/leda64_nm/h_array_geometry.json')
json_d_array_geometry_nm  = os.path.join(fileroot, 'config/leda64_nm/d_array_geometry.json')
json_h_antenna_nm         = os.path.join(fileroot, 'config/leda64_nm/h_antenna.json')
json_d_antenna_nm         = os.path.join(fileroot, 'config/leda64_nm/d_antenna.json')
json_antenna_el_lens_nm   = os.path.join(fileroot, 'config/leda64_nm/z_antenna_el_lens.json')

# A few sources to phase to
src_names = ['CYG', 'CAS', 'TAU', 'VIR']
src_ras   = [299.86791, 350.84583, 83.63333, 187.705833]
src_decs  = [40.733888, 58.810833, 22.01444, 12.39111]
