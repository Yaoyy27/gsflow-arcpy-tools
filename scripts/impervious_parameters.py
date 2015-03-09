#--------------------------------
# Name:         impervious_parameters.py
# Purpose:      GSFLOW impervious parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-03-08
# Python:       2.7
#--------------------------------

from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
import re
import sys

import arcpy
from arcpy import env
from arcpy.sa import *
##import numpy as np

from support_functions import *

################################################################################

def gsflow_impervious_parameters(workspace, config_path=None):
    """Calculate GSFLOW Impervious Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW Impervious Parameters')

        ## Initialize hru_parameters class
        hru = hru_parameters(config_path)

        ## Open input parameter config file
        inputs_cfg = ConfigParser.ConfigParser()
        try:
            inputs_cfg.readfp(open(config_path))
        except:
            logging.error('\nERROR: Config file could not be read, '+
                          'is not an input file, or does not exist\n'+
                          'ERROR: config_file = {0}\n').format(config_path)
            raise SystemExit()
        logging.debug('\nReading Input File')

        ##
        imperv_orig_path = inputs_cfg.get('INPUTS', 'impervious_orig_path')
        ##imperv_proj_method = inputs_cfg.get('INPUTS', 'impervious_projection_method')
        imperv_proj_method = 'NEAREST'
        imperv_cs = inputs_cfg.getint('INPUTS', 'impervious_cellsize')
        imperv_pct_flag = inputs_cfg.getboolean('INPUTS', 'impervious_pct_flag')

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## Impervious raster must exist
        if not arcpy.Exists(imperv_orig_path):
            logging.error('\nERROR: Impervious raster does not exist')
            raise SystemExit()
        
        ## Check other inputs
        if imperv_cs <= 0:
            logging.error('\nERROR: soil cellsize must be greater than 0')
            raise SystemExit()
        imperv_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if imperv_proj_method.upper() not in imperv_proj_method_list:
            logging.error('\nERROR: Impervious projection method must be: {0}'.format(
                ', '.join(imperv_proj_method_list)))
            raise SystemExit()

        ## Build output folder if necessary
        imperv_temp_ws = os.path.join(hru.param_ws, 'impervious_rasters')
        if not os.path.isdir(imperv_temp_ws): 
	    os.mkdir(imperv_temp_ws)
        ## Output paths
        imperv_path = os.path.join(imperv_temp_ws, 'impervious_cover.img')

		
        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = imperv_temp_ws
        env.scratchWorkspace = hru.scratch_ws

        ## Check field
        logging.info('\nAdding impervious fields if necessary')
        add_field_func(hru.polygon_path, hru.imperv_pct_field, 'DOUBLE')
        ##add_field_func(hru.polygon_path, hru.carea_min_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.carea_max_field, 'DOUBLE')

        ## Available Water Capacity (AWC)
        logging.info('\nProjecting/clipping impervious cover raster')
        imperv_orig_sr = Raster(imperv_orig_path).spatialReference
        logging.debug('  Impervious GCS:  {0}'.format(
            imperv_orig_sr.GCS.name))
        ## Remove existing projected raster
        if arcpy.Exists(imperv_path):
            arcpy.Delete_management(imperv_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, imperv_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Project impervious raster
        ## DEADBEEF - Arc10.2 ProjectRaster does not extent
        ##env.extent = hru.extent
        project_raster_func(
            imperv_orig_path, imperv_path, hru.sr,
            imperv_proj_method, imperv_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), imperv_orig_sr, hru)
        ##arcpy.ProjectRaster_management(
        ##    imperv_orig_path, imperv_path, hru.sr,
        ##    imperv_proj_method, imperv_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    imperv_orig_sr)
        ##arcpy.ClearEnvironment('extent')

        ## List of rasters, fields, and stats for zonal statistics
        zs_imperv_dict = dict()
        zs_imperv_dict[hru.imperv_pct_field] = [imperv_path, 'MEAN']
        ##zs_imperv_dict[hru.carea_min_field] = [imperv_path, 'MEAN']
        ##zs_imperv_dict[hru.carea_max_field] = [imperv_path, 'MEAN']

        ## Calculate zonal statistics
        logging.info('\nCalculating zonal statistics')
        zonal_stats_func(zs_imperv_dict, hru.polygon_path, hru.point_path, hru)

        ## Calculate CAREA_MIN / CAREA_MAX
        logging.info('\nCalculating CAREA_MIN / CAREA_MAX')
        if imperv_pct_flag:
            arcpy.CalculateField_management(
                hru.polygon_path, hru.imperv_pct_field,
                '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')      
            ##arcpy.CalculateField_management(
            ##    hru.polygon_path, hru.carea_min_field,
            ##    '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')      
            arcpy.CalculateField_management(
                hru.polygon_path, hru.carea_max_field,
                '0.01 * !{0}!'.format(hru.imperv_pct_field), 'PYTHON')
        else:
            ##arcpy.CalculateField_management(
            ##    hru.polygon_path, hru.carea_min_field,
            ##    '!{0}!'.format(hru.imperv_pct_field), 'PYTHON')      
            arcpy.CalculateField_management(
                hru.polygon_path, hru.carea_max_field,
                '!{0}!'.format(hru.imperv_pct_field), 'PYTHON')

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()

################################################################################

if __name__ == '__main__':
    workspace = os.getcwd()
    log_ws = os.path.join(workspace, 'logs')
    if not os.path.isdir(log_ws): 
        os.mkdir(log_ws)
    log_file_name = 'gsflow_impervious_log.txt'

    ## Create Basic Logger
    ##logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    ## Create File Logger
    logging.basicConfig(
        level = logging.DEBUG, format='%(message)s', filemode='w',
        filename = os.path.join(log_ws, log_file_name))
    ## Create Display Logger
    log_console = logging.StreamHandler()
    log_console.setLevel(logging.INFO)
    console_format = logging.Formatter('%(message)s')
    log_console.setFormatter(console_format)
    logging.getLogger('').addHandler(log_console)

    ## Get GSFLOW config file
    ini_re = re.compile('\w*.ini$', re.I)
    try: 
        ini_path = sys.argv[1]
    except IndexError:
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_impervious_parameters')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate GSFLOW Impervious Parameters
    gsflow_impervious_parameters(workspace, ini_path)
