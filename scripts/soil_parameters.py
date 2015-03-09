#--------------------------------
# Name:         soil_parameters.py
# Purpose:      GSFLOW soil parameters
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
##import tempfile

import arcpy
from arcpy import env
from arcpy.sa import *
##import numpy as np

from support_functions import *

################################################################################

def gsflow_soil_parameters(workspace, config_path=None):
    """Calculate GSFLOW Soil Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW Soil Parameters')

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

        ## Input parameters
        soil_pct_flag = inputs_cfg.getboolean('INPUTS', 'soil_pct_flag')
        moist_init_ratio = inputs_cfg.getfloat('INPUTS', 'moist_init_ratio')
        rechr_init_ratio = inputs_cfg.getfloat('INPUTS', 'rechr_init_ratio')

        ## Use Ksat to calculate ssr2gw_rate and slowcoef_lin
        ##calc_ssr2gw_rate_flag = inputs_cfg.getboolean(
        ##    'INPUTS', 'calc_ssr2gw_rate_flag')
        ##calc_slowcoef_flag = inputs_cfg.getboolean(
        ##    'INPUTS', 'calc_slowcoef_flag')

        ## Clip root depth to soil depth
        clip_root_depth_flag = inputs_cfg.getboolean(
            'INPUTS', 'clip_root_depth_flag')

        ## Input folders
        soil_temp_ws = os.path.join(hru.param_ws, 'soil_rasters')
        if not os.path.isdir(soil_temp_ws): 
	    os.mkdir(soil_temp_ws)

        ## Input paths
        awc_path = os.path.join(soil_temp_ws, 'awc.img')
        clay_pct_path = os.path.join(soil_temp_ws, 'clay_pct.img')
        sand_pct_path = os.path.join(soil_temp_ws, 'sand_pct.img')
        ##silt_pct_path = os.path.join(soil_temp_ws, 'silt_pct.img')
        ksat_path = os.path.join(soil_temp_ws, 'ksat.img')
        soil_depth_path = os.path.join(soil_temp_ws, 'soil_depth.img')

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## All of the soil rasters must exist
        ## Check that the projected/clipped/filled raster exists
        if not arcpy.Exists(awc_path):
            logging.error('\nERROR: AWC raster does not exist')
            raise SystemExit()
        if not arcpy.Exists(clay_pct_path):
            logging.error('\nERROR: Clay raster does not exist')
            raise SystemExit()
        if not arcpy.Exists(sand_pct_path):
            logging.error('\nERROR: Sand raster does not exist')
            raise SystemExit()
        ##if not arcpy.Exists(silt_path):
        ##    logging.error('\nERROR: Silt raster does not exist')
        ##    raise SystemExit()
        ##if ((calc_ssr2gw_rate_flag or calc_slowcoef_flag) and
        ##    not arcpy.Exists(ksat_path)):
        if not arcpy.Exists(ksat_path):
            logging.error('\nERROR: Ksat raster does not exist')
            raise SystemExit()
        if clip_root_depth_flag and not arcpy.Exists(soil_depth_path):
            logging.error('\nERROR: Soil depth raster does not exist')
            raise SystemExit()
        ## Check soil init ratios
        if moist_init_ratio < 0 or moist_init_ratio > 1:
            logging.error('\nERROR: Soil moist_init_ratio must be between 0 & 1')
            raise SystemExit()            
        if rechr_init_ratio < 0 or rechr_init_ratio > 1:
            logging.error('\nERROR: Soil rechr_init_ratio must be between 0 & 1')
            raise SystemExit()
        
        ## DEM Slope is needed for SSR2GW_RATE
        dem_temp_ws = os.path.join(hru.param_ws, 'dem_rasters')
        dem_slope_path = os.path.join(dem_temp_ws, 'dem_slope.img')
        if not os.path.isdir(dem_temp_ws):
            logging.error(
                '\nERROR: DEM temp folder does not exist\n'+
                '\nERROR: Try re-running gsflow_dem_2_stream.py')
            raise SystemExit()
        if not os.path.isfile(dem_slope_path):
            logging.error(
                '\nERROR: Slope raster does not exist\n'+
                '\nERROR: Try re-running gsflow_dem_2_stream.py')
            raise SystemExit()

        ## Output paths
        soil_type_path = os.path.join(soil_temp_ws, 'soil_type.img')
        moist_max_path = os.path.join(soil_temp_ws, 'soil_moist_max.img')
        rechr_max_path = os.path.join(soil_temp_ws, 'soil_rechr_max.img')
        ##ssr2gw_rate_path = os.path.join(soil_temp_ws, 'ssr2gw_rate.img')
        ##slowcoef_lin_path = os.path.join(soil_temp_ws, 'slowcoef_lin.img')

        ## Root depth is calculated by veg script
        veg_temp_ws = os.path.join(hru.param_ws, 'veg_rasters')
        root_depth_path = os.path.join(veg_temp_ws, 'root_depth.img')
        if not arcpy.Exists(root_depth_path):
            logging.error(
                '\nERROR: Root depth raster does not exists'+
                '\nERROR: Try re-running gsflow_veg_parameters script\n')
            raise SystemExit()


        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = soil_temp_ws
	env.scratchWorkspace = hru.scratch_ws

        ## Check field
        logging.info('\nAdding soil fields if necessary')
        add_field_func(hru.polygon_path, hru.awc_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.clay_pct_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.sand_pct_field, 'DOUBLE')
        ##add_field_func(hru.polygon_path, hru.silt_pct_field, 'DOUBLE')
        ##if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        add_field_func(hru.polygon_path, hru.ksat_field, 'DOUBLE')
        if clip_root_depth_flag:
            add_field_func(hru.polygon_path, hru.soil_depth_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.root_depth_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.soil_type_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.moist_init_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.moist_max_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.rechr_init_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.rechr_max_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.ssr2gw_rate_field, 'DOUBLE')
        ##add_field_func(hru.polygon_path, hru.pref_flow_den_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.slowcoef_lin_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.slowcoef_sq_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.fastcoef_lin_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.fastcoef_sq_field, 'DOUBLE')

            
        ## Clip root depth to soil depth
        if clip_root_depth_flag:
            ## This will clip root depth to soil depth
            ## Minimum of root depth and soil depth
            logging.info('Clipping root depth to soil depth')
            root_depth_obj = Con(
                Raster(root_depth_path) < Raster(soil_depth_path),
                Raster(root_depth_path), Raster(soil_depth_path))
            root_depth_obj.save(root_depth_path)
            del root_depth_obj

        ## Calculate maximum soil moisture
        logging.info('\nCalculating soil {0}'.format(hru.moist_max_field))
        moist_max_obj = Raster(awc_path) * Raster(root_depth_path)
        moist_max_obj.save(moist_max_path)
        del moist_max_obj

        ## Calculate soil recharge zone maximum
        logging.info('Calculating soil {0}'.format(hru.rechr_max_field))
        ## Minimum of rooting depth and 18 (inches?)
        rechr_max_obj = Float(
            Con(Raster(root_depth_path) < 18, Raster(root_depth_path), 18))
        rechr_max_obj *= Raster(awc_path)
        rechr_max_obj.save(rechr_max_path)
        del rechr_max_obj

        #### Read in slope raster and convert to radians
        ##dem_slope_obj = math.pi * Raster(dem_slope_path) / 180
        ##porosity_obj = 0.475
        ##
        #### Gravity drainage to groundwater reservoir linear coefficient
        ##logging.info('\nCalculating SSR2GW_RATE')
        ##logging.info('  Assuming slope is in degrees')
        ##logging.info('  Porosity is currently fixed at: {0}'.format(
        ##    porosity_obj))
        ##ssr2gw_rate_obj = (
        ##    Raster(ksat_path) * porosity_obj * (1 - dem_slope_obj))
        ##ssr2gw_rate_obj.save(ssr2gw_rate_path)
        ##del ssr2gw_rate_obj
        ##
        #### Gravity drainage to groundwater reservoir linear coefficient
        ##logging.info('\nCalculating SLOWCOEF_L')
        ##logging.info('  Assuming slope is in degrees')
        ##logging.info('  Porosity is currently fixed at: {0}'.format(
        ##slowcoef_lin_obj = (
        ##    Raster(ksat_path) * math.sin(dem_slope_obj) / 
        ##    (porosity_obj * hru_length_obj))
        ##slowcoef_lin_obj.save(slowcoef_lin_path)
        ##del slowcoef_lin_obj, hru_length_obj
        ##del dem_slope_obj, porosity_obj
       
        ## List of rasters, fields, and stats for zonal statistics
        zs_soil_dict = dict()
        zs_soil_dict[hru.awc_field] = [awc_path, 'MEAN']
        zs_soil_dict[hru.clay_pct_field] = [clay_pct_path, 'MEAN']
        zs_soil_dict[hru.sand_pct_field] = [sand_pct_path, 'MEAN']
        ##zs_soil_dict[hru.silt_pct_field] = [silt_pct_path, 'MEAN']
        ##if calc_ssr2gw_rate_flag or calc_slowcoef_flag:
        zs_soil_dict[hru.ksat_field] = [ksat_path, 'MEAN']
        if clip_root_depth_flag:
            zs_soil_dict[hru.soil_depth_field] = [soil_depth_path, 'MEAN']
            zs_soil_dict[hru.root_depth_field] = [root_depth_path, 'MEAN']
        zs_soil_dict[hru.moist_max_field] = [moist_max_path, 'MEAN']
        zs_soil_dict[hru.rechr_max_field] = [rechr_max_path, 'MEAN']
        ##zs_soil_dict[hru.ssr2gw_rate_field] = [ssr2gw_rate_path, 'MEAN']
        ##zs_soil_dict[hru.slowcoef_lin_field] = [slowcoef_lin_path, 'MEAN']
        
        ## Calculate zonal statistics
        logging.info('\nCalculating zonal statistics')
        zonal_stats_func(zs_soil_dict, hru.polygon_path, hru.point_path, hru)


        ## Make a fishnet layer for calculating fields
        hru_polygon_layer = "hru_polygon_layer"
        arcpy.MakeFeatureLayer_management(
            hru.polygon_path, hru_polygon_layer)

        ## Calculate SOIL_TYPE
        logging.info('\nCalculating {0}'.format(hru.soil_type_field))
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 1'.format(hru.type_in_field))
        if soil_pct_flag:
            soil_type_pct = (50, 40)
        else:
            soil_type_pct = (0.50, 0.40)
        soil_type_cb = (
            'def soil_type_func(clay, sand):\n'+
            '    if sand > {0}: return 1\n'+
            '    elif clay > {1}: return 3\n'+
            '    else: return 2\n').format(*soil_type_pct)
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.soil_type_field,
            'soil_type_func(!{0}!, !{1}!)'.format(
                hru.clay_pct_field, hru.sand_pct_field),
            'PYTHON', soil_type_cb)      
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "SWITCH_SELECTION")
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.soil_type_field, '0', 'PYTHON')      


        ## Calculate SOIL_MOIST_INIT & SOIL_RECHR_INIT from max values
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 1 AND "{1}" >= 0'.format(
                hru.type_in_field, hru.moist_max_field))
        logging.info('\nCalculating {0} as {2} * {1}'.format(
            hru.moist_init_field, hru.moist_max_field, moist_init_ratio))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.moist_init_field,
            '!{0}! * {1}'.format(hru.moist_max_field, moist_init_ratio), 'PYTHON')      
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "SWITCH_SELECTION")
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.moist_init_field, '0', 'PYTHON')      
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.moist_max_field, '0', 'PYTHON')      


        ## Calculate SOIL_MOIST_INIT & SOIL_RECHR_INIT from max values
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 1 AND "{1}" >= 0'.format(
                hru.type_in_field, hru.rechr_max_field))
        logging.info('Calculating {0} as {2} * {1}'.format(
            hru.rechr_init_field, hru.rechr_max_field, moist_init_ratio))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.rechr_init_field,
            '!{0}! * {1}'.format(hru.rechr_max_field, moist_init_ratio), 'PYTHON')      
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "SWITCH_SELECTION",
            '"{0}" <> 1'.format(hru.type_in_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.rechr_init_field, '0', 'PYTHON')      
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.rechr_max_field, '0', 'PYTHON')      


        ## Gravity drainage to groundwater reservoir linear coefficient
        ## Default value is 0.1 (range 0-1)
        ## Convert Ksat from um/s to in/day
        ##if calc_ssr2gw_rate_flag:
        logging.info('\nCalculating {0}'.format(hru.ssr2gw_rate_field))
        logging.info('  {0} must be in um/s'.format(hru.ksat_field))
        porosity_flt = 0.475
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 1 AND "{1}" >= 0'.format(hru.type_in_field, hru.ksat_field))
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.ssr2gw_rate_field,
            '!{0}! * (3600 * 24 / (2.54 * 10000)) * (1 - !{1}!) * {2}'.format(
                hru.ksat_field, hru.slope_rad_field, porosity_flt),
            'PYTHON')
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "SWITCH_SELECTION")
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.ssr2gw_rate_field, '0', 'PYTHON')      

        
        ## Default value is 0.015 (range 0-1)
        ## Convert Ksat from um/s to m/day
        ##if calc_slowcoef_lin_flag:
        logging.info('Calculating {0}'.format(hru.slowcoef_lin_field))
        logging.info('  {0} must be in um/s'.format(hru.ksat_field))
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "NEW_SELECTION",
            '"{0}" = 1 AND "{1}" >= 0'.format(hru.type_in_field, hru.ksat_field))
        ##slowcoef_lin_cb = (
        ##    'def fd_len(fd, cs):\n'+
        ##    '    if fd in [1, 4, 16, 64]: return cs\n'+
        ##    '    elif fd in [2, 8, 32, 128]: return 1.414 * cs\n'+
        ##    '    else: return cs\n'+
        ##    'def slowcoef_lin(ksat, slope, porosity, fd, cs):\n'+
        ##    '    return (ksat * (3600 * 24 / 1000000) * math.sin(slope) / '+
        ##    '(porosity * fd_len(fd, cs)))\n')
        slowcoef_lin_cb = (
            'def slowcoef_lin(ksat, slope, porosity, fd, cs):\n'+
            '    return (ksat * 0.0864 * math.sin(slope) / (porosity * cs))\n')
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.slowcoef_lin_field,
            'slowcoef_lin(!{0}!, !{1}!, {2}, !{3}!, {4})'.format(
                hru.ksat_field, hru.slope_rad_field, porosity_flt,
                hru.flow_dir_field, hru.cs),
            'PYTHON', slowcoef_lin_cb)      
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_layer, "SWITCH_SELECTION")
        arcpy.CalculateField_management(
            hru_polygon_layer, hru.slowcoef_lin_field, '0', 'PYTHON')      

        ## Cleanup
        arcpy.Delete_management(hru_polygon_layer)
        del hru_polygon_layer

        #### Reset soils values for lake cells (HRU_TYPE == 2)
        #### Also reset for ocean cells (HRU_TYPE == 0 and DEM_ADJ == 0)
        #### This will remove many NoData cells ( == -999)
        ##if True:
        ##    logging.info('\nClearing lake soil parameters')
        ##    ##logging.info(
        ##    ##    '\nClearing soil parameters for lake and inactive ocean cells')
        ##    hru_polygon_layer = "hru_polygon_layer"
        ##    arcpy.MakeFeatureLayer_management(
        ##        hru.polygon_path, hru_polygon_layer)
        ##    arcpy.SelectLayerByAttribute_management(
        ##        hru_polygon_layer, "NEW_SELECTION",
        ##        '"{0}" = 2 OR ("{0}" = 0 AND "{1}" = 0)'.format(
        ##            hru.type_in_field, hru.dem_adj_field))
        ##        ##'"{0}" = 2 AND "{1}" = -999'.format(
        ##        ##    hru.type_in_field, hru.awc_field))
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.awc_field, 0, 'PYTHON')      
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.clay_pct_field, 0, 'PYTHON')      
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.sand_pct_field, 0, 'PYTHON')
        ##    ##arcpy.CalculateField_management(
        ##    ##    hru_polygon_layer, hru.silt_pct_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.ksat_field, 0, 'PYTHON')
        ##    ## DEADBEEF
        ##    ## Soil type must be 1-3, can i set it to 0?
        ##    logging.info('  Setting default lake soil type to 2')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.soil_type_field, 2, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.moist_init_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.moist_max_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.rechr_init_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.rechr_max_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.ssr2gw_rate_field, 0, 'PYTHON')
        ##    arcpy.CalculateField_management(
        ##        hru_polygon_layer, hru.slowcoef_lin_field, 0, 'PYTHON')
        ##    arcpy.Delete_management(hru_polygon_layer)
        ##    del hru_polygon_layer

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
    log_file_name = 'gsflow_soil_log.txt'

    ## Create Basic Logger
    ##logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    #### Create File Logger
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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_soil_parameters')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate GSFLOW Soil Parameters
    gsflow_soil_parameters(workspace, ini_path)
