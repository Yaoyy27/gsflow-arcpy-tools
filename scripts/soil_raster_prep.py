#--------------------------------
# Name:         soil_raster_prep.py
# Purpose:      GSFLOW soil raster prep
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-04-27
# Python:       2.7
#--------------------------------

import argparse
from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import os
import re
import sys
##import tempfile

import arcpy
from arcpy import env
from arcpy.sa import *
##import numpy as np

from support_functions import *

################################################################################

def gsflow_soil_raster_prep(config_path, overwrite_flag=False, debug_flag=False):
    """Prepare GSFLOW soil rasters

    Args:
        config_file: Project config file path
        ovewrite_flag: boolean, overwrite existing files
        debug_flag: boolean, enable debug level logging
    Returns:
        None
    """

    try:
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

        ## Log DEBUG to file
        log_file_name = 'gsflow_soil_log.txt'
        log_console = logging.FileHandler(
            filename=os.path.join(hru.log_ws, log_file_name), mode='w')
        log_console.setLevel(logging.DEBUG)
        log_console.setFormatter(logging.Formatter('%(message)s'))
        logging.getLogger('').addHandler(log_console)
        logging.info('\nPrepare GSFLOW Soil Rasters')

        soil_orig_ws  = inputs_cfg.get('INPUTS', 'soil_orig_folder')
        awc_name      = inputs_cfg.get('INPUTS', 'awc_name')
        clay_pct_name = inputs_cfg.get('INPUTS', 'clay_pct_name')
        sand_pct_name = inputs_cfg.get('INPUTS', 'sand_pct_name')
        ##silt_pct_name = inputs_cfg.get('INPUTS', 'silt_pct_name')
        soil_proj_method = 'NEAREST'
        soil_cs = inputs_cfg.getint('INPUTS', 'soil_cellsize')
        fill_soil_nodata_flag = inputs_cfg.getboolean(
            'INPUTS', 'fill_soil_nodata_flag')
        
	## Use Ksat to calculate ssr2gw_rate and slowcoef_lin
        ##calc_ssr2gw_rate_flag = inputs_cfg.getboolean(
        ##    'INPUTS', 'calc_ssr2gw_rate_flag')
        ##calc_slowcoef_flag = inputs_cfg.getboolean(
        ##    'INPUTS', 'calc_slowcoef_flag')
        ##if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        ksat_name = inputs_cfg.get('INPUTS', 'ksat_name')

        ## Clip root depth to soil depth
        clip_root_depth_flag = inputs_cfg.getboolean(
            'INPUTS', 'clip_root_depth_flag')
        if clip_root_depth_flag:
            soil_depth_name = inputs_cfg.get('INPUTS', 'soil_depth_name')

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## All of the soil rasters must exist
        awc_orig_path = os.path.join(soil_orig_ws, awc_name)
        clay_pct_orig_path = os.path.join(soil_orig_ws, clay_pct_name)
        sand_pct_orig_path = os.path.join(soil_orig_ws, sand_pct_name)
        ##silt_orig_path = os.path.join(soil_orig_ws, silt_pct_name)
        ##if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        ksat_orig_path = os.path.join(soil_orig_ws, ksat_name)
        if clip_root_depth_flag:
            soil_depth_path = os.path.join(soil_orig_ws, soil_depth_name)

        ## Check that either the original or projected/clipped raster exists
        if not arcpy.Exists(awc_orig_path):
            logging.error('\nERROR: AWC raster does not exist')
            raise SystemExit()
        if not arcpy.Exists(clay_pct_orig_path):
            logging.error('\nERROR: Clay raster does not exist')
            raise SystemExit()
        if not arcpy.Exists(sand_pct_orig_path):
            logging.error('\nERROR: Sand raster does not exist')
            raise SystemExit()
        ##if not arcpy.Exists(silt_orig_path):
        ##    logging.error('\nERROR: Silt raster does not exist')
        ##    raise SystemExit()
        ##if ((calc_ssr2gw_rate_flag or calc_slowcoef_flag) and
        ##    not arcpy.Exists(ksat_orig_path)):
        if not arcpy.Exists(ksat_orig_path):
            logging.error('\nERROR: Ksat raster does not exist')
            raise SystemExit()
        if clip_root_depth_flag and not arcpy.Exists(soil_depth_orig_path):
            logging.error('\nERROR: Soil depth raster does not exist')
            raise SystemExit()

        ## Check other inputs
        if soil_cs <= 0:
            logging.error('\nERROR: soil cellsize must be greater than 0')
            raise SystemExit()
        soil_proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if soil_proj_method.upper() not in soil_proj_method_list:
            logging.error('\nERROR: Soil projection method must be: {0}'.format(
                ', '.join(soil_proj_method_list)))
            raise SystemExit()

        ## Build output folder if necessary
        soil_temp_ws = os.path.join(hru.param_ws, 'soil_rasters')
        if not os.path.isdir(soil_temp_ws): 
            os.mkdir(soil_temp_ws)
        ## Output paths
        awc_path = os.path.join(soil_temp_ws, 'awc.img')
        clay_pct_path = os.path.join(soil_temp_ws, 'clay_pct.img')
        sand_pct_path = os.path.join(soil_temp_ws, 'sand_pct.img')
        ##silt_pct_path = os.path.join(soil_temp_ws, 'silt_pct.img')
        ksat_path = os.path.join(soil_temp_ws, 'ksat.img')
        soil_depth_path = os.path.join(soil_temp_ws, 'soil_depth.img')

        ## Root depth is calculated by veg script
        ##veg_temp_ws = os.path.join(hru.param_ws, 'veg_rasters')
        ##root_depth_path = os.path.join(veg_temp_ws, 'root_depth.img')
        ##if not arcpy.Exists(root_depth_path):
        ##    logging.error(
        ##        '\nERROR: Root depth raster does not exists'+
        ##        '\nERROR: Try re-running gsflow_veg_parameters script\n')
        ##    raise SystemExit()

		
        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = soil_temp_ws
        env.scratchWorkspace = hru.scratch_ws
            
        ## Available Water Capacity (AWC)
        logging.info('\nProjecting/clipping AWC raster')
        soil_orig_sr = Raster(awc_orig_path).spatialReference
        logging.debug('  AWC GCS:  {0}'.format(
            soil_orig_sr.GCS.name))
        ## Remove existing projected raster
        if arcpy.Exists(awc_path):
            arcpy.Delete_management(awc_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, soil_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Project soil raster
        ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
        project_raster_func(
            awc_orig_path, awc_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    awc_orig_path, awc_path, hru.sr,
        ##    soil_proj_method, soil_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    soil_orig_sr)
        ##arcpy.ClearEnvironment('extent')

        ## Percent clay
        logging.info('Projecting/clipping clay raster')
        soil_orig_sr = Raster(clay_pct_orig_path).spatialReference
        logging.debug('  Clay GCS: {0}'.format(
            soil_orig_sr.GCS.name))
        ## Remove existing projected raster
        if arcpy.Exists(clay_pct_path):
            arcpy.Delete_management(clay_pct_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, soil_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Project soil raster
        ## DEADBEEF - Arc10.2 ProjectRaster does not extent
        project_raster_func(
            clay_pct_orig_path, clay_pct_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    clay_pct_orig_path, clay_pct_path, hru.sr,
        ##    soil_proj_method, soil_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    soil_orig_sr)
        ##arcpy.ClearEnvironment('extent')

        ## Percent sand
        logging.info('Projecting/clipping sand raster')
        soil_orig_sr = Raster(sand_pct_orig_path).spatialReference
        logging.debug('  Sand GCS: {0}'.format(
            soil_orig_sr.GCS.name))
        ## Remove existing projected raster
        if arcpy.Exists(sand_pct_path):
            arcpy.Delete_management(sand_pct_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, soil_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Project soil raster
        ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
        project_raster_func(
            sand_pct_orig_path, sand_pct_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    sand_pct_orig_path, sand_pct_path, hru.sr,
        ##    soil_proj_method, soil_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    soil_orig_sr)
        ##arcpy.ClearEnvironment('extent')

        ## Hydraulic conductivity
        ##if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        logging.info('Projecting/clipping ksat raster')
        ksat_orig_sr = Raster(ksat_orig_path).spatialReference
        logging.debug('  Ksat GCS: {0}'.format(
            soil_orig_sr.GCS.name))
        ## Remove existing projected raster
        if arcpy.Exists(ksat_path):
            arcpy.Delete_management(ksat_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, ksat_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Project ksat raster
        ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
        project_raster_func(
            ksat_orig_path, ksat_path, hru.sr,
            soil_proj_method, soil_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    ksat_orig_path, ksat_path, hru.sr,
        ##    soil_proj_method, soil_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    soil_orig_sr)
        ##arcpy.ClearEnvironment('extent')

        ## Soil depth is only needed if clipping root depth
        if clip_root_depth_flag:
            logging.info('\nProjecting/clipping depth raster')
            soil_orig_sr = Raster(soil_depth_orig_path).spatialReference
            logging.debug('  Depth GCS: {0}'.format(
                soil_orig_sr.GCS.name))
            ## Remove existing projected raster
            if arcpy.Exists(soil_depth_path):
                arcpy.Delete_management(soil_depth_path)
            ## Set preferred transforms
            transform_str = transform_func(hru.sr, soil_orig_sr)
            logging.debug('  Transform: {0}'.format(transform_str))
            logging.debug('  Projection method: NEAREST')
            ## Project soil raster
            ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
            project_raster_func(
                soil_depth_orig_path, soil_depth_path, hru.sr,
                soil_proj_method, soil_cs, transform_str,
                '{0} {1}'.format(hru.ref_x, hru.ref_y), soil_orig_sr, hru)
            ##env.extent = hru.extent
            ##arcpy.ProjectRaster_management(
            ##    soil_depth_orig_path, soil_depth_path, hru.sr,
            ##    soil_proj_method, soil_cs, transform_str,
            ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
            ##    soil_orig_sr)
            ##arcpy.ClearEnvironment('extent')
        
        ## Clip root depth to soil depth
        ##if clip_root_depth_flag:
        ##    ## This will clip root depth to soil depth
        ##    ## Minimum of root depth and soil depth
        ##    logging.info('Clipping root depth to soil depth')
        ##    root_depth_obj = Con(
        ##        Raster(root_depth_path) < Raster(soil_depth_path),
        ##        Raster(root_depth_path), Raster(soil_depth_path))
        ##    root_depth_obj.save(root_depth_path)
        ##    del root_depth_obj

        ## Fill soil nodata values using nibble
        if fill_soil_nodata_flag:
            logging.info('\nFilling soil nodata values using Nibble')
            soil_raster_list = [
                awc_path, clay_pct_path, sand_pct_path, ksat_path]
            if clip_root_depth_flag: soil_raster_list.append(soil_depth_path)
            for soil_raster_path in soil_raster_list:
                logging.info('  {0}'.format(soil_raster_path))
                ## DEADBEEF - Check if there is any nodata to be filled first?
                mask_obj = Int(1000 * SetNull(
                    Raster(soil_raster_path) < 0, Raster(soil_raster_path)))
                input_obj = Con(IsNull(mask_obj), 0, mask_obj)
                nibble_obj = 0.001 * Nibble(input_obj, mask_obj, 'ALL_VALUES')
                nibble_obj.save(soil_raster_path)
                arcpy.BuildPyramids_management(soil_raster_path)                    

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()

################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Soil Raster Prep',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', required=True,
        help='Project input file', metavar='PATH')
    parser.add_argument(
        '-o', '--overwrite', default=False, action="store_true", 
        help='Force overwrite of existing files')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()

    ## Create Basic Logger
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    #### Get GSFLOW config file
    ##ini_re = re.compile('\w*.ini$', re.I)
    ##try: 
    ##    ini_path = sys.argv[1]
    ##except IndexError:
    ##    ini_path = get_ini_file(workspace, ini_re, 'gsflow_soil_parameters')
    ##del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))

    ## Prepare GSFLOW Soil Rasters
    gsflow_soil_raster_prep(
        config_path=args.ini, overwrite_flag=args.overwrite,
        debug_flag=args.loglevel==logging.DEBUG)