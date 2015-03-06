#--------------------------------
# Name:         gsflow_veg_parameters
# Purpose:      GSFLOW vegetation parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2014-10-13
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

from gsflow_support_functions import *

################################################################################

def gsflow_veg_parameters(workspace, config_path=None):
    """Calculate GSFLOW Vegetation Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW Vegetation Parameters')

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

        ## Landfire Vegetation Type
        veg_type_orig_path = inputs_cfg.get('INPUTS', 'veg_type_orig_path')
        veg_type_cs = inputs_cfg.getint('INPUTS', 'veg_type_cellsize')
        try: veg_type_field = inputs_cfg.get('INPUTS', 'veg_type_field')
        except: veg_type_field = None

        ## Landfire Vegetation Cover
        veg_cover_orig_path = inputs_cfg.get('INPUTS', 'veg_cover_orig_path')
        veg_cover_cs = inputs_cfg.getint('INPUTS', 'veg_cover_cellsize')

        ## Remap
        remap_ws = inputs_cfg.get('INPUTS', 'remap_folder')
        aspect_remap_name = inputs_cfg.get('INPUTS', 'aspect_remap')
        temp_adj_remap_name = inputs_cfg.get('INPUTS', 'temp_adj_remap')
        cov_type_remap_name = inputs_cfg.get('INPUTS', 'cov_type_remap')
        covden_sum_remap_name = inputs_cfg.get('INPUTS', 'covden_sum_remap')
        covden_win_remap_name = inputs_cfg.get('INPUTS', 'covden_win_remap')
        snow_intcp_remap_name = inputs_cfg.get('INPUTS', 'snow_intcp_remap')
        srain_intcp_remap_name = inputs_cfg.get('INPUTS', 'srain_intcp_remap')
        wrain_intcp_remap_name = inputs_cfg.get('INPUTS', 'wrain_intcp_remap')
        root_depth_remap_name = inputs_cfg.get('INPUTS', 'root_depth_remap')

        
        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## Check that either the original vegetation raster exist
        if not arcpy.Exists(veg_cover_orig_path):
            logging.error(
                '\nERROR: Vegetation cover raster does not exist')
            raise SystemExit()
        if not arcpy.Exists(veg_type_orig_path):
            logging.error(
                '\nERROR: Vegetation type raster does not exist')
            raise SystemExit()
        ## Vegetation cover can be set from another field in the raster
        ## This is mostly for US_120EVT
        if not veg_type_field:
            logging.info('\n  Using VALUE field to set vegetation type')
            veg_type_field = 'VALUE'
        elif len(arcpy.ListFields(veg_type_orig_path, veg_type_field)) == 0:
            logging.info(
                ('  veg_type_field {0} does not exist\n  Using VALUE '+
                 'field to set vegetation type').format(veg_type_field))
            veg_type_field = 'VALUE'
        elif arcpy.ListFields(veg_type_orig_path, veg_type_field)[0].type not in ['Integer', 'SmallInteger']:
            logging.info(
                ('  veg_type_field {0} is not an integer type\n  Using VALUE '+
                 'field to set vegetation type').format(veg_type_field))
            veg_type_field = 'VALUE'
        ## Check that remap folder is valid
        if not os.path.isdir(remap_ws):
            logging.error('\nERROR: Remap folder does not exist')
            raise SystemExit()
        ## Check that remap files exist
        cov_type_remap_path = os.path.join(remap_ws, cov_type_remap_name)
        covden_sum_remap_path = os.path.join(remap_ws, covden_sum_remap_name)
        covden_win_remap_path = os.path.join(remap_ws, covden_win_remap_name)
        snow_intcp_remap_path = os.path.join(remap_ws, snow_intcp_remap_name)
        srain_intcp_remap_path = os.path.join(remap_ws, srain_intcp_remap_name)
        wrain_intcp_remap_path = os.path.join(remap_ws, wrain_intcp_remap_name)
        root_depth_remap_path = os.path.join(remap_ws, root_depth_remap_name)
        if not os.path.isfile(cov_type_remap_path):
            logging.error('\nERROR: Cover type remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(covden_sum_remap_path):
            logging.error('\nERROR: Summer cover density remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(covden_win_remap_path):
            logging.error('\nERROR: Winter cover density remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(snow_intcp_remap_path):
            logging.error('\nERROR: Winter snow interception remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(srain_intcp_remap_path):
            logging.error('\nERROR: Summer rain interception remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(wrain_intcp_remap_path):
            logging.error('\nERROR: Winter rain interception remap file does not exist')
            raise SystemExit()
        if not os.path.isfile(root_depth_remap_path):
            logging.error('\nERROR: Root depth remap file does not exist')
            raise SystemExit()

        ## Check other inputs
        if veg_type_cs <= 0:
            logging.error('\nERROR: Veg. type cellsize must be greater than 0')
            raise SystemExit()
        if veg_cover_cs <= 0:
            logging.error('\nERROR: Veg. cover cellsize must be greater than 0')
            raise SystemExit()

        ## Build output folders if necesssary
        veg_temp_ws = os.path.join(hru.param_ws, 'veg_rasters')
        if not os.path.isdir(veg_temp_ws): 
			os.mkdir(veg_temp_ws)
        ## Output paths
        veg_cover_path = os.path.join(veg_temp_ws, 'veg_cover.img')
        veg_type_path = os.path.join(veg_temp_ws, 'veg_type.img')
        cov_type_path = os.path.join(veg_temp_ws, 'cov_type.img')
        covden_sum_path = os.path.join(veg_temp_ws, 'covden_sum.img')
        covden_win_path = os.path.join(veg_temp_ws, 'covden_win.img')
        snow_intcp_path = os.path.join(veg_temp_ws, 'snow_intcp.img')
        wrain_intcp_path = os.path.join(veg_temp_ws, 'wrain_intcp.img')
        srain_intcp_path = os.path.join(veg_temp_ws, 'srain_intcp.img')
        root_depth_path = os.path.join(veg_temp_ws, 'root_depth.img')
        rad_trncf_path = os.path.join(veg_temp_ws, 'rad_trncf.img')

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = veg_temp_ws
        env.scratchWorkspace = hru.scratch_ws

        ## Check fields
        logging.info('\nAdding vegetation fields if necessary')
        add_field_func(hru.polygon_path, hru.cov_type_field, 'SHORT')
        add_field_func(hru.polygon_path, hru.covden_sum_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.covden_win_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.rad_trncf_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.snow_intcp_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.srain_intcp_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.wrain_intcp_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.root_depth_field, 'DOUBLE')

            
        ## Assume all vegetation rasters will need to be rebuilt
        ## Check veg cover and veg type rasters
        ## This will check for matching spat. ref., snap point, and cellsize

        ## Project/clip veg cover to match HRU
        logging.info('\nProjecting/clipping vegetation cover raster')
        veg_cover_orig_sr = Raster(veg_cover_orig_path).spatialReference
        ## Remove existing clipped/projected veg cover raster
        if arcpy.Exists(veg_cover_path):
            arcpy.Delete_management(veg_cover_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, veg_cover_orig_sr)
        logging.debug('  Transform:{0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')

        ## Project veg cover
        ## DEADBEEF - Arc10.2 ProjectRaster does not extent
        project_raster_func(
            veg_cover_orig_path, veg_cover_path, hru.sr,
            'NEAREST', veg_cover_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), veg_cover_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    veg_cover_orig_path, veg_cover_path, hru.sr,
        ##    'NEAREST', veg_cover_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    veg_cover_orig_sr)
        ##arcpy.ClearEnvironment('extent')
        del transform_str, veg_cover_orig_sr

        ## Project/clip veg type to match HRU
        logging.info('Projecting/clipping vegetation type raster')
        veg_type_orig_sr = Raster(veg_type_orig_path).spatialReference
        ## Remove existing clipped/projected veg type raster
        if arcpy.Exists(veg_type_path):
            arcpy.Delete_management(veg_type_path)
        ## Set preferred transforms
        transform_str = transform_func(hru.sr, veg_type_orig_sr)
        logging.debug('  Transform: {0}'.format(transform_str))
        logging.debug('  Projection method: NEAREST')
        ## Use a different field to calculate vegetation type
        if veg_type_field <> 'VALUE':
            logging.info(
                '  Calculating vegetation type from {0} field'.format(
                    veg_type_field))
            veg_type_obj = Lookup(veg_type_orig_path, veg_type_field)
        else:
            veg_type_obj = Raster(veg_type_orig_path)

        ## Project veg type
        ## DEADBEEF - Arc10.2 ProjectRaster does not extent
        project_raster_func(
            veg_type_obj, veg_type_path, hru.sr,
            'NEAREST', veg_type_cs, transform_str,
            '{0} {1}'.format(hru.ref_x, hru.ref_y), veg_type_orig_sr, hru)
        ##env.extent = hru.extent
        ##arcpy.ProjectRaster_management(
        ##    veg_type_obj, veg_type_path, hru.sr,
        ##    'NEAREST', veg_type_cs, transform_str,
        ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
        ##    veg_type_orig_sr)
        ##arcpy.ClearEnvironment('extent')
        del transform_str, veg_type_orig_sr, veg_type_obj

        ## Reclassifying vegetation cover type
        logging.info('\nCalculating COV_TYPE')
        logging.debug('  Reclassifying: {0}'.format(cov_type_remap_path))
        cov_type_obj = ReclassByASCIIFile(
            veg_type_path, cov_type_remap_path)
        cov_type_obj.save(cov_type_path)
        del cov_type_obj

        ## Summer cover density
        logging.info('Calculating COVDEN_SUM')
        logging.debug('  Reclassifying: {0}'.format(covden_sum_remap_path))
        covden_sum_obj = ReclassByASCIIFile(
            veg_cover_path, covden_sum_remap_path)
        covden_sum_obj *= 0.01
        covden_sum_obj.save(covden_sum_path)
        del covden_sum_obj

        ## Winter cover density
        logging.info('Calculating COVDEN_WIN')
        logging.debug('  Reclassifying: {0}'.format(covden_win_remap_path))
        covden_win_obj = ReclassByASCIIFile(
            cov_type_path, covden_win_remap_path)
        covden_win_obj *= 1
        covden_win_obj *= Raster(covden_sum_path)
        covden_win_obj.save(covden_win_path)
        del covden_win_obj

        ## Snow interception storage capacity
        logging.info('Calculating SNOW_INTCP')
        logging.debug('  Reclassifying: {0}'.format(snow_intcp_remap_path))
        snow_intcp_obj = ReclassByASCIIFile(
            veg_type_path, snow_intcp_remap_path)
        snow_intcp_obj.save(snow_intcp_path)
        del snow_intcp_obj

        ## Winter rain interception storage capacity
        logging.info('Calculating WRAIN_INTCP')
        logging.debug('  Reclassifying: {0}'.format(wrain_intcp_remap_path))
        wrain_intcp_obj = ReclassByASCIIFile(
            veg_type_path, wrain_intcp_remap_path)
        wrain_intcp_obj.save(wrain_intcp_path)
        del wrain_intcp_obj

        ## Summer rain interception storage capacity
        logging.info('Calculating SRAIN_INTCP')
        logging.debug('  Reclassifying: {0}'.format(srain_intcp_remap_path))
        srain_intcp_obj = ReclassByASCIIFile(
            veg_type_path, srain_intcp_remap_path)
        srain_intcp_obj.save(srain_intcp_path)
        del srain_intcp_obj

        ## Root depth
        logging.info('Calculating ROOT_DEPTH')
        logging.debug('  Reclassifying: {0}'.format(root_depth_remap_path))
        root_depth_obj = ReclassByASCIIFile(
            veg_type_path, root_depth_remap_path)
        root_depth_obj.save(root_depth_path)
        del root_depth_obj

        ## Short-wave radiation transmission coefficent
        logging.info('Calculating {0}'.format(hru.rad_trncf_field))
        rad_trncf_obj = 0.9917 * Exp(-2.7557 * Raster(covden_win_path))
        rad_trncf_obj.save(rad_trncf_path)
        del rad_trncf_obj


        ## List of rasters, fields, and stats for zonal statistics
        zs_veg_dict = dict()
        zs_veg_dict[hru.cov_type_field] = [cov_type_path, 'MAJORITY']
        zs_veg_dict[hru.covden_sum_field] = [covden_sum_path, 'MEAN']
        zs_veg_dict[hru.covden_win_field] = [covden_win_path, 'MEAN']
        zs_veg_dict[hru.snow_intcp_field] = [snow_intcp_path, 'MAJORITY']
        zs_veg_dict[hru.srain_intcp_field] = [srain_intcp_path, 'MAJORITY']
        zs_veg_dict[hru.wrain_intcp_field] = [wrain_intcp_path, 'MAJORITY']
        zs_veg_dict[hru.root_depth_field] = [root_depth_path, 'MAJORITY']
        zs_veg_dict[hru.rad_trncf_field] = [rad_trncf_path, 'MEAN']


        ## Calculate zonal statistics
        logging.info('\nCalculating vegetation zonal statistics')
        zonal_stats_func(zs_veg_dict, hru.polygon_path, hru.point_path, hru)


        ## Short-wave radiation transmission coefficient
        ##logging.info('\nCalculating {0}'.format(hru.rad_trncf_field))
        ##arcpy.CalculateField_management(
        ##    hru.polygon_path, hru.rad_trncf_field,
        ##    '0.9917 * math.exp(-2.7557 * !{0}!)'.format(hru.covden_win_field),
        ##    'PYTHON')


        ## Clear COV_TYPE values for lake cells (HRU_TYPE == 2)
        if True:
            logging.info('\nClearing lake nodata vegetation parameters')
            ##logging.info(
            ##    '\nClearing vegetation parameters for lake and inactive cells')
            hru_polygon_layer = "hru_polygon_layer"
            arcpy.MakeFeatureLayer_management(
                hru.polygon_path, hru_polygon_layer)
            arcpy.SelectLayerByAttribute_management(
                hru_polygon_layer, "NEW_SELECTION",
                '"{0}" = 2 OR ("{0}" = 0 AND "{1}" = 0)'.format(
                    hru.type_in_field, hru.dem_adj_field))
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.cov_type_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.covden_sum_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.covden_win_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.snow_intcp_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.srain_intcp_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.wrain_intcp_field, 0, 'PYTHON')      
            arcpy.CalculateField_management(
                hru_polygon_layer, hru.rad_trncf_field, 0, 'PYTHON')      
            arcpy.Delete_management(hru_polygon_layer)
            del hru_polygon_layer

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
    log_file_name = 'gsflow_veg_log.txt'

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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_veg_parameters')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate GSFLOW Vegetation Parameters
    gsflow_veg_parameters(workspace, ini_path)
