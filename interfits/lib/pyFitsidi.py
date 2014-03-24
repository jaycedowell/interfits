# -*- coding: utf-8 -*-

"""
pyFitsidi.py
=============

Created by Danny Price on 2011-04-28.
Copyright (c) 2011 The University of Oxford. All rights reserved.

This file is a collection of modules for creating blank FITS IDI files.
It consists mainly of pretty basic functions to create blank tables.

The FITS handling is all done by pyFits, and blank arrays are created
with numpy. So you'll need to have both of these installed on your machine.
In addition, the config xml file is read with lxml, so install that too.


FITS-IDI Table Overview
~~~~~~~~~~~~~~~~~~~~~~~~

Here is a quick rundown of the tables that this script can make.
Not all these tables are mandatory, so I've starred the ones which are the most
important to understand and to get into your dataset.

* **ANTENNA**           Antenna polarization information
* **ARRAY_GEOMETRY**    Time system information and antenna coordinates
* FLAG                  Flagged data
* **FREQUENCY**         Frequency setups
* GAIN_CURVE            Antenna gain curves
* INTERFEROMETER_MODEL  Correlator model
* PHASE-CAL             Phase cal measurements
* **SOURCE**            Sources observed
* SYSTEM_TEMPERATURE    System and antenna temperatures
* **UV_DATA**           Visibility data

These tables are new from the FITS IDI update, and are all optional:

* BANDPASS              Bandpass functions
* BASELINE              Baseline-specific gain factors
* CALIBRATION           Gains as a function of time
* WEATHER               Meteorological data

Module listing
~~~~~~~~~~~~~~

"""

import sys, os
import pyfits as pf, numpy as np
from lxml import etree
import lxml

__version__ = '0.0'
__all__ = ['checkConfigType', 'parseConfig', 'make_primary', 'make_array_geometry', 'make_antenna', 
           'make_frequency', 'make_source', 'make_uv_data', 'make_interferometer_model', 'make_system_temperature', 
           'make_gain_curve', 'make_phase_cal', 'make_flag', 'make_bandpass', 'make_weather', 'make_baseline', 
           'make_calibration', 'make_model_comps', '__version__', '__all__']


def checkConfigType(config):
    """ Helper function to check config parameter type

    Notes
    -----
    Checks whether configuration is being passed as an lxml element tree,
    or as a filename.
    """

    if type(config) is lxml.etree._ElementTree:
        return config
    if type(config) is str:
        return etree.parse(config)
    if type(config) is unicode:
        return etree.parse(config)
    else:
        print type(config)
        raise

def parseConfig(tagname, config='config.xml'):
    """ Finds tagname, in elementTree x, parses and returns dictionary of values
  This is a helper function, and is not usually called directly.
  
  Notes
  -----
  This function uses eval() to evaluate the text string inside a child tag. As such,
  exercise caution! todo: block off certain modules to eval()
  """

    xmlData = checkConfigType(config)
    x = xmlData.getroot()
    T = True # FITS just uses T for True, python (and pyFITS) uses True

    # As we reference 'parameters', we need to search for this first
    params = dict([(child.tag, eval(child.text.strip())) for child in x.find('PARAMETERS').getchildren()])

    # This line makes me very happy, but will probably infuriate others:
    try:
        vals = dict([(child.tag, eval(child.text.strip())) for child in x.find(tagname).getchildren()])
    except SyntaxError:
        if type(child.text.strip()) == type(" "):
            vals = dict([(child.tag, child.text.strip()) for child in x.find(tagname).getchildren()])
        else:
            print child.tag
            print child.text
            raise
    except NameError:
        if type(child.text.strip()) == type(" "):
            vals = dict([(child.tag, child.text.strip()) for child in x.find(tagname).getchildren()])
        else:
            print child.tag
            print child.text
            raise
        print child.tag
        print child.text
        raise
    except AttributeError:
        print "ERROR encountered when parsing config file."
        print "Tag name:   %s"%tagname
        print "Child tag:  %s"%child.tag
        print "Child text: %s"%child.text
        raise

    return vals


def make_primary(config='config.xml'):
    """  Creates the primary header data unit (HDU).
  
  This function generates header keywords from the file headers/primary.tpl
  
  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  """

    # Make a new blank FITS HDU
    hdu = pf.PrimaryHDU()

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    primary = parseConfig('PRIMARY', config)
    common = parseConfig('COMMON', config)

    for key in primary: hdu.header.update(key, primary[key])
    #for key in common: hdu.header.update(key, common[key])

    hdu.verify() # Will raise a warning if there's an issue

    return hdu


def make_array_geometry(config='config.xml', num_rows=1):
    """Creates a vanilla ARRAY_GEOMETRY table HDU.
  
  One row is required for each antenna in the array (num_rows)

  Notes
  -----
  Table is built with the following columns:
  
  * ANNAME: Antenna name
  * STABXYZ: Antenna relative position vector ECI components, in meters
  * DERXYZ: Antenna velocity vector components, in meters/sec
  * ORBPARM: Orbital parameters
  * NOSTA: Antenna ID number for station
  * MNTSTA: Antenna mount type (0 is alt-azimuth)
  * STAXOF: Antenna axis offset, in meters
  * DIAMETER: Antenna diameter (optional)

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    array_geometry = parseConfig('ARRAY_GEOMETRY', config)
    common = parseConfig('COMMON', config)

    # Generate the columns for the table header
    c = []

    c.append(pf.Column(name='ANNAME', format='8A',
                       array=np.zeros(num_rows, dtype='a8')))

    c.append(pf.Column(name='STABXYZ', format='3D',
                       unit='METERS', array=np.zeros(num_rows, dtype='3float64')))

    c.append(pf.Column(name='DERXYZ', format='3E',
                       unit='METERS/SEC', array=np.zeros(num_rows, dtype='3float32')))

    orb_format = '%iD' % params['NORB']
    orb_dtype = '%ifloat64' % params['NORB']
    c.append(pf.Column(name='ORBPARM', format=orb_format,
                       array=np.zeros(num_rows, dtype=orb_dtype)))

    c.append(pf.Column(name='NOSTA', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='MNTSTA', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='STAXOF', format='3E',
                       unit='METERS', array=np.zeros(num_rows, dtype='3float32')))

    c.append(pf.Column(name='DIAMETER', format='1E',
                       unit='METERS', array=np.zeros(num_rows, dtype='float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in array_geometry: tblhdu.header.update(key, array_geometry[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_antenna(config='config.xml', num_rows=1):
    """  Creates a vanilla ANTENNA table HDU
  
  Notes
  -----
  CASA didn't like importing POLCALA and POLCALB, so they are currently commented out.
  
  Table is built with the following columns:
  
  * TIME: Difference of antenna table time interval centre time and RDATE 0 hours
  * TIME_INTERVAL: Antenna table time interval width
  * ANNAME: Antenna name, should match value in ARRAY_GEOMETRY
  * ANTENNA_NO: Antemma ID number for station
  * ARRAY:    Array ID number
  * FREQID:   Frequency setup ID number
  * NO_LEVELS:Number of digitiser levels
  * POLYTYA:  Feed A polarisation direction
  * POLAA:    Feed A polarisation (degrees)
  * POLCALA:  Feed A polarisation parameters
  * POLYTYB:  As above, for feed B
  * POLAB:    As above, for feed B
  * POLCALB:  As above, for feed B

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('ANTENNA', config)
    common = parseConfig('COMMON', config)

    nband = params['NBAND']
    npcal = params['NPCAL']

    c = []

    c.append(pf.Column(name='TIME', format='1D',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='TIME_INTERVAL', format='1E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='ANNAME', format='8A',
                       array=np.zeros(num_rows, dtype='a8')))

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='NO_LEVELS', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='POLTYA', format='1A',
                       array=np.zeros(num_rows, dtype='a1')))

    c.append(pf.Column(name='POLTYB', format='1A',
                       array=np.zeros(num_rows, dtype='a1')))

    pol_format = '%iE' % nband
    pol_dtype = '%ifloat32' % nband
    c.append(pf.Column(name='POLAA', format=pol_format,
                       unit='DEGREES', array=np.zeros(num_rows, dtype=pol_dtype)))

    c.append(pf.Column(name='POLAB', format=pol_format,
                       unit='DEGREES', array=np.zeros(num_rows, dtype=pol_dtype)))

    # nb: Was encontering errors with CASA with this column
    #c.append(pf.Column(name='POLCALA', format='1E',\
    #  array=np.zeros(32,dtype='float32')))

    # nb: Was encontering errors with CASA with this column
    #c.append(pf.Column(name='POLCALB', format='1E',\
    #  array=np.zeros(32,dtype='float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_frequency(config='config.xml', num_rows=1):
    """ Creates a vanilla FREQUENCY table HDU
  
  Notes
  -----
  Table is built with the following columns:
  
  * FREQID:   Frequency setup ID number
  * BANDFREQ: Frequency band base offset (Hz)
  * CH_WIDTH: Frequency channel width (Hz)
  * TOTAL_BANDWIDTH: Frequency bandwidth (Hz)
  * SIDEBAND: Sideband flag (1 indicates upper sideband)
  * BB_CHAN:  ? 

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('FREQUENCY', config)
    common = parseConfig('COMMON', config)

    nband = params['NBAND']

    c = []

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    ba_format = '%iD' % nband
    ba_dtype = '%ifloat64' % nband
    c.append(pf.Column(name='BANDFREQ', format=ba_format,
                       unit='HZ', array=np.zeros(num_rows, dtype=ba_dtype)))

    ch_format = '%iE' % nband
    ch_dtype = '%ifloat32' % nband
    c.append(pf.Column(name='CH_WIDTH', format=ch_format,
                       unit='HZ', array=np.zeros(num_rows, dtype=ch_dtype)))

    c.append(pf.Column(name='TOTAL_BANDWIDTH', format='1E',
                       unit='HZ', array=np.zeros(num_rows, dtype='float32')))

    si_format = '%iJ' % nband
    si_dtype = '%iint32' % nband
    c.append(pf.Column(name='SIDEBAND', format=si_format,
                       array=np.zeros(num_rows, dtype=si_dtype)))

    # Not really sure what this does, so commented it out
    #c.append(pf.Column(name='BB_CHAN',  format='1J',\
    #  array=np.zeros(num_rows,dtype='int32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_source(config='config.xml', num_rows=1):
    """ Creates a vanilla SOURCE table HDU

  Notes
  -----
  Table is built with the following columns:
  
  * SOURCE_ID Source ID number
  * SOURCE    Source name
  * QUAL      Source qualifier number.
  * CALCODE   Source calibrator code. 
  * FREQID    Source frequency ID
  * IFLUX     Source I flux density
  * QFLUX     Source Q flux density
  * UFLUX     Source U flux density
  * VFLUX     SourceV flux density
  * ALPHA     Source spectral index
  * FREQOFF   Source frequency offset
  * RAEPO     Source J2000 equatorial position RA coordinate
  * DECPO     Source J2000 equatorial position DEC coordinate
  * EQUINOX   Mean Equinox
  * RAAPP     Source apparent equatorial position RA coordinate
  * DECAPP    Source apparent equatorial position DEC coordinate
  * SYSVEL    Systematic velocity.
  * VELTYP    Systematic velocity reference frame.
  * VELDEF    Systematic velocity convention.
  * RESTFREQ  Line rest frequency.
  * PMRA      Source proper motion RA coordinate
  * PMDEC     Source proper motion DEC coordinate
  * PARALLAX  Source parallax. 

Parameters
----------
config: string
  filename of xml configuration file, defaults to 'config,xml'
num_rows: int
  number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('SOURCE', config)
    common = parseConfig('COMMON', config)

    nband = params['NBAND']
    so_format = '%iE' % nband
    so_dtype = '%ifloat32' % nband

    c = []

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='SOURCE', format='16A',
                       array=np.zeros(num_rows, dtype='16a')))

    c.append(pf.Column(name='QUAL', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='CALCODE', format='4A',
                       array=np.zeros(num_rows, dtype='4a')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='IFLUX', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='QFLUX', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='UFLUX', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='VFLUX', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='ALPHA', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='FREQOFF', format=so_format,
                       array=np.zeros(num_rows, dtype=so_dtype)))

    c.append(pf.Column(name='RAEPO', format='1D',
                       unit='DEGREES', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='DECEPO', format='1D',
                       unit='DEGREES', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='EQUINOX', format='8A',
                       array=np.zeros(num_rows, dtype='8a')))

    c.append(pf.Column(name='RAAPP', format='1D',
                       unit='DEGREES', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='DECAPP', format='1D',
                       unit='DEGREES', array=np.zeros(num_rows, dtype='float64')))

    sv_format = '%iD' % nband
    sv_dtype = '%ifloat64' % nband
    c.append(pf.Column(name='SYSVEL', format=sv_format,
                       unit='METERS/SEC', array=np.zeros(num_rows, dtype=sv_dtype)))

    c.append(pf.Column(name='VELTYP', format='8A',
                       array=np.zeros(num_rows, dtype='8a')))

    c.append(pf.Column(name='VELDEF', format='8A',
                       array=np.zeros(num_rows, dtype='8a')))

    rf_format = '%iD' % nband
    rf_dtype = '%ifloat64' % nband
    c.append(pf.Column(name='RESTFREQ', format=rf_format,
                       unit='HZ', array=np.zeros(num_rows, dtype=rf_dtype)))

    c.append(pf.Column(name='PMRA', format='1D',
                       unit='DEGREES/DAY', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='PMDEC', format='1D',
                       unit='DEGREES/DAY', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='PARALLAX', format='1E',
                       unit='ARCSEC', array=np.zeros(num_rows, dtype='float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_uv_data(config='config.xml', num_rows=1, uu_data=None, vv_data=None, ww_data=None,
                 date_data=None, time_data=None, baseline_data=None, filter_data=None, source_data=None,
                 freqid_data=None, inttim_data=None, weights_data=None, gateid_data=None,
                 flux_data=None, weights_col=False):
    """ Creates a vanilla UV_DATA table HDU
  
  Parameters
  ----------
  num_rows
    Number of rows in your UV dataset. For example, if you have 528 baselines 
    (inc. autocorr), then for 10 time dumps, you'll need 10*528 rows
  weights: bool
    Whether or not to include weights column. Defaults to False.
  [xyz]_data: array of data
    Defaults to None, in which case zeros are written. Much faster to write data in columns
    than to loop through row by row.
  
  Notes
  -----
  Table is built with the following columns:
  
  * UU           Baseline vector U coordinate (seconds)
  * VV           Baseline vector V coordinate
  * WW           Baseline vector W coordinate
  * DATE         UTC Julian day value for time 00:00:00  on the day of the observation
  * TIME         Fraction of Julian day from midnight to timestamp on day of observation.
  * BASELINE     Antenna baseline pair ID.
  * SOURCE       Data source ID
  * FREQID       Data frequency setup ID
  * INTTIM       Data integration time
  * WEIGHT       Data weights (one element for each freq channel)
  * FLUX         UV visibility data matrix
  

  """
    c = []

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('UV_DATA', config)
    common = parseConfig('COMMON', config)

    if uu_data == None: uu_data = np.zeros(num_rows, dtype='float32')
    if vv_data == None: vv_data = np.zeros(num_rows, dtype='float32')
    if ww_data == None: ww_data = np.zeros(num_rows, dtype='float32')
    if date_data == None: date_data = np.zeros(num_rows, dtype='float64')
    if time_data == None: time_data = np.zeros(num_rows, dtype='float64')
    if baseline_data == None: baseline_data = np.zeros(num_rows, dtype='int32')
    if source_data == None: source_data = np.zeros(num_rows, dtype='int32')
    if freqid_data == None: freqid_data = np.zeros(num_rows, dtype='int32')
    if inttim_data == None: inttim_data = np.zeros(num_rows, dtype='float32')

    c.append(pf.Column(name='UU', format='1E', unit='SECONDS', array=uu_data))
    c.append(pf.Column(name='VV', format='1E', unit='SECONDS', array=vv_data))
    c.append(pf.Column(name='WW', format='1E', unit='SECONDS', array=ww_data))
    c.append(pf.Column(name='DATE', format='1D', unit='DAYS', array=date_data))
    c.append(pf.Column(name='TIME', format='1D', unit='DAYS', array=time_data))
    c.append(pf.Column(name='BASELINE', format='1J', array=baseline_data))
    c.append(pf.Column(name='SOURCE', format='1J', array=source_data))
    c.append(pf.Column(name='FREQID', format='1J', array=freqid_data))
    c.append(pf.Column(name='INTTIM', format='1E', unit='SECONDS', array=inttim_data))

    # The following depends on number of stokes, number of bands and number of channels
    nchan = params['NCHAN']
    nstokes = params['NSTOKES']
    nband = params['NBAND']
    flux_nbits = nchan * nstokes * nband * 2 # 2= Real & Im
    weights_nbits = nstokes * nband  # 2= Real & Im
    flux_format = '%iE' % flux_nbits
    flux_dtype = '%ifloat32' % flux_nbits
    weights_format = '%iE' % weights_nbits
    weights_dtype = '%ifloat32' % weights_nbits

    if flux_data == None:  flux_data = np.zeros(num_rows, dtype=flux_dtype)
    if weights_data == None: weights_data = np.zeros(num_rows, dtype=weights_dtype)

    c.append(pf.Column(name='FLUX', format=flux_format, unit='UNCALIB', array=flux_data))

    if weights_col:
        c.append(pf.Column(name='WEIGHT', format=weights_format, array=weights_data))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in sorted(common): tblhdu.header.update(key, common[key])
    for key in sorted(cards): tblhdu.header.update(key, cards[key])

    return tblhdu


def make_interferometer_model(config='config.xml', num_rows=1):
    """
  Creates a vanilla INTERFEROMETER_MODEL table HDU.

  Notes
  -----
  This table is optional, and is not included in Medicina FITS IDI
  
  Table is built with the following columns:

   
  * TIME          Starting time of interval
  * TIME_INTERVAL Duration of interval
  * SOURCE_ID     Source ID number
  * ANTENNA_NO    Antenna number
  * ARRAY         Array number
  * FREQID        Frequency setup number
  * I.FAR.ROT     Ionospheric Faraday rotation
  * FREQ_VAR      Time-variable frequency offsets
  * PDELAY_1      Phase delay polynomials for polarization 1
  * GDELAY_1      Group delay polynomials for polarization 1
  * PRATE_1       Phase delay rate polynomials for polarization 1
  * GRATE_1       Group rate polynomials for polarization 1
  * DISP_1        Dispersive delay for polarization 1 at 1m wavelength
  * DDISP_1       Rate of change of dispersive del for pol 1 at 1m

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
    
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('INTERFEROMETER_MODEL', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='TIME', format='1D',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='TIME_INTERVAL', format='1E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='I.FAR.ROT', format='1E',
                       unit='RAD/M**2', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='FREQ.VAR', format='1E',
                       unit='HZ', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='PDELAY_1', format='1E',
                       unit='TURNS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='GDELAY_1', format='1E',
                       unit='SECONDS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='PRATE_1', format='1E',
                       unit='HZ', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='GRATE_1', format='1E',
                       unit='SEC/SEC', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='DISP_1', format='1E',
                       unit='SECONDS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='DDISP_1', format='1E',
                       unit='SEC/SEC', array=np.zeros(num_rows, dtype='float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_system_temperature(config='config.xml', num_rows=1):
    """ Creates a vanilla SYSTEM_TEMPERATURE table HDU
  
  Notes
  -----
  This is an optional table.
  Todo: add dual pol support.
  
  Table is built with the following columns:
  
  * TIME          Central time of interval covered
  * TIME_INTERVAL Duration of interval
  * SOURCE_ID     Source ID number
  * ANTENNA_NO    Antenna number
  * ARRAY         Array number
  * FREQID        Frequency setup number
  * TSYS_1        System temperatures for polarization 1
  * TANT_1        Antenna temperatures for polarization 1

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('SYSTEM_TEMPERATURE', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='TIME', format='1D',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='TIME_INTERVAL', format='1E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='TSYS_1', format='1E',
                       unit='KELVIN', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='TANT_1', format='1E',
                       unit='KELVIN', array=np.zeros(num_rows, dtype='int32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_gain_curve(config='config.xml', num_rows=1):
    """ Creates a vanilla GAIN_CURVE table HDU
    
  Notes
  --------------------------------
  This is an optional table (we will not include it). 
  Todo: Add switch to allow dual polarisation.
  
  Table is built with the following columns:
  
  * ANTENNA_NO  Antenna number
  * ARRAY       Array number
  * FREQID      Frequency setup number
  * TYPE_1      Gain curve types for polarization 1
  * NTERM_1     Numbers of terms or entries for polarization 1
  * X_TYP_1     x value types for polarization 1
  * Y_TYP_1     y value types for polarization 1
  * X_VAL_1     x values for polarization 1
  * Y_VAL_1     y values for polarization 1
  * GAIN_1      Relative gain values for polarization 1
  * SENS_1      Sensitivities for polarization 1

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
    
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('GAIN_CURVE', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='TYPE_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='NTERM_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='X_TYP_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='Y_TYP_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='X_VAL_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='Y_VAL_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='GAIN_1', format='1E',
                       array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='SENS_1', format='1E',
                       unit='K/JY', array=np.zeros(num_rows, dtype='float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_phase_cal(config='config.xml', num_rows=1):
    """ Creates a vanilla PHASE-CAL table HDU


  Notes
  -------
  This is an optional table (we will not include it)
  
  Table is built with the following columns:
  
  * TIME          Central time of interval covered
  * TIME_INTERVAL Duration of interval
  * SOURCE_ID     Source ID number
  * ANTENNA_NO    Antenna number
  * ARRAY         Array number
  * FREQID        Frequency setup number
  * CABLE_CAL     Cable calibration measurement
  * STATE_1       State counts for polarization 1
  * PC_FREQ_1     Phase cal tone frequencies for polarization 1
  * PC_REAL_1     real parts of phase-cal measurements for pol 1
  * PC_IMAG_1     imaginary parts of phasecal measurements for pol 1
  * PC_RATE_1     phase-cal rates for polarization 1

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
    
  """
    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('PHASE_CAL', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='TIME', format='1D',
                       array=np.zeros(num_rows, dtype='float64'), unit='DAYS'))

    c.append(pf.Column(name='TIME_INTERVAL', format='1E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='CABLE_CAL', format='1E',
                       unit='SECONDS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='STATE_1', format='1J',
                       unit='PERCENT', array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='PC_FREQ_1', format='1J',
                       unit='HZ', array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='PC_REAL_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='PC_IMAG_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='PC_RATE_1', format='1J',
                       unit='SEC/SEC', array=np.zeros(num_rows, dtype='int32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_flag(config='config.xml', num_rows=1):
    """ Creates a vanilla FLAG table HDU
  

  Notes
  ----------
  This table is optional.
  
  Table is built with the following columns:

  * SOURCE_ID Source ID number
  * ARRAY     Array number
  * ANTS      Antenna numbers
  * FREQID    Frequency setup number
  * TIMERANG  Time range
  * BANDS     Band flags
  * CHANS     Channel range
  * PFLAGS    Polarization flags
  * REASON    Reason for flag
  * SEVERITY  Severity code

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    cards = parseConfig('FLAG', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ANTS', format='2J',
                       array=np.zeros(num_rows, dtype='2int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='TIMERANG', format='2E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='2float32')))

    c.append(pf.Column(name='BANDS', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='CHANS', format='2J',
                       array=np.zeros(num_rows, dtype='2int32')))

    c.append(pf.Column(name='PFLAGS', format='4J',
                       array=np.zeros(num_rows, dtype='4int32')))

    c.append(pf.Column(name='REASON', format='24A',
                       array=np.zeros(num_rows, dtype='a24')))

    c.append(pf.Column(name='SEVERITY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_bandpass(config='config.xml', num_rows=1):
    """ Creates a vanilla BANDPASS table HDU
  
  Notes
  -----
  Table is built with the following columns:

  * TIME: Difference of bandpass table time int center time and RDATE 0 hours
  * TIME_INTERVAL: Bandpass table time interval width
  * SOURCE_ID:  Source ID number
  * ANTENNA_NO: Antenna ID number for station
  * ARRAY:      Array ID number
  * FREQID:     Frequency setup ID number (should match ANTENNA)
  * BANDWIDTH:  Frequency band width described by bandpass
  * BAND_FREQ:  Frequency band base offset
  * REFANT_1:   Reference antenna ID number
  * BREAL_1:    Bandpass response real componet
  * BIMAG_1:    Bandpass response imaginary component

  Parameters
  ----------
  config: string
    filename of xml configuration file, defaults to 'config,xml'
  num_rows: int
    number of rows to generate. Rows will be filled with numpy zeros.
  """

    # Generate headers from config file
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('BANDPASS', config)
    common = parseConfig('COMMON', config)

    c = []

    c.append(pf.Column(name='TIME', format='1D',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='TIME_INTERVAL', format='1E',
                       unit='DAYS', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='SOURCE_ID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ANTENNA_NO', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='ARRAY', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='FREQID', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='BANDWIDTH', format='1E',
                       unit='HZ', array=np.zeros(num_rows, dtype='float32')))

    c.append(pf.Column(name='BAND_FREQ', format='1D',
                       unit='HZ', array=np.zeros(num_rows, dtype='float64')))

    c.append(pf.Column(name='REFANT_1', format='1J',
                       array=np.zeros(num_rows, dtype='int32')))

    c.append(pf.Column(name='BREAL_1', format='1024E',
                       array=np.zeros(num_rows, dtype='1024float32')))

    c.append(pf.Column(name='BIMAG_1', format='1024E',
                       array=np.zeros(num_rows, dtype='1024float32')))

    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)

    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])

    return tblhdu


def make_weather():
    """ Makes weather table
    The WEATHER table contains meteorological data for the antennæ and times
    used in the FITS-IDI file.

    This table is currently not supported (on the todo list)

    """
    pass


def make_baseline():
    """Makes baseline table
    The BASELINE table contains baseline-dependent multiplicative and additive corrections.

    This table is currently not supported (on the todo list)

    """
    pass


def make_calibration(config='config.xml', num_rows=1, date=None, delaya_data=None, phasea_data=None, delayb_data=None, phaseb_data=None):
    """ Make calibration table
  
    From the FITS IDI documentation:
    This chapter is included for documentation and discussion purposes only. So far as this
    author is aware, no software has been implemented to either write or read the CALIBRATION table.
    Therefore, the description provided in this section should be regarded as tentative.
    In fact, it is not at all clear what the intentions were in the case of some of the
    columns specified for this table.

    Notes
    -----
    Table is built with the following LEDA-specific columns:
    
    * ANTENNA_NO: Antenna ID number for station
    * ARRAY:      Array ID number
    * DELAY_A:    Delay in seconds applied to the A polarization antenna
    * PHASE_A     Phse in radians applied to the A polarization antenna
    * DELAY_B:    Delay in seconds applied to the B polarization antenna
    * PHASE_B     Phse in radians applied to the B polarization antenna
    
    Parameters
    ----------
    config: string
      filename of xml configuration file, defaults to 'config,xml'
    num_rows: int
      number of rows to generate. Rows will be filled with numpy zeros.
    delaya_data: float32
      delays in seconds applied to each A polarization antenna
    delayb_data: float32
      delays in seconds applied to each B polarization antenna
    """

    # Generate headers from config file and update it with the delay calibration date
    params = parseConfig('PARAMETERS', config)
    cards = parseConfig('CALIBRATION', config)
    common = parseConfig('COMMON', config)
    if date != None:
        cards['DATE-GEN'] = date
    
    an_data = np.arange(1, num_rows+1)
    ar_data = np.ones(num_rows)
    if delaya_data == None: delaya_data = np.zeros(num_rows)
    if phasea_data == None: phasea_data = np.zeros(num_rows)
    if delayb_data == None: delayb_data = np.zeros(num_rows)
    if phaseb_data == None: phaseb_data = np.zeros(num_rows)
    
    c = []
    c.append(pf.Column(name='ANTENNA_NO', format='1J', array=an_data.astype(np.int32)))
    c.append(pf.Column(name='ARRAY', format='1J', array=ar_data.astype(np.int32)))
    c.append(pf.Column(name='DELAY_A', format='1E', array=delaya_data.astype(np.float32)))
    c.append(pf.Column(name='PHASE_A', format='1E', array=phasea_data.astype(np.float32)))
    c.append(pf.Column(name='DELAY_B', format='1E', array=delayb_data.astype(np.float32)))
    c.append(pf.Column(name='PHASE_B', format='1E', array=phaseb_data.astype(np.float32)))
    
    coldefs = pf.ColDefs(c)
    tblhdu = pf.new_table(coldefs)
    
    for key in cards: tblhdu.header.update(key, cards[key])
    for key in common: tblhdu.header.update(key, common[key])
    tblhdu.header.update('DATE-GEN', date)
    
    return tblhdu


def make_model_comps():
    """ Make model comps table
  
    The MODEL COMPS table is one of those reserved for use by the VLBA.
    As such, I can't see anyone in CASPER needing it any time soon.

    Again, this table is currently not supported (on the todo list)
    """
    pass