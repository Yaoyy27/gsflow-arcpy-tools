#--------------------------------
# Name:         prism_4km_normals.py
# Purpose:      GSFLOW PRISM parameters from default 400m normals 
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-03-08
# Python:       2.7
#--------------------------------

from collections import defaultdict
import ConfigParser
import datetime as dt
import logging
import multiprocessing
import os
import re
import sys
##import tempfile
from time import clock

import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np

from support_functions import *

################################################################################

def gsflow_prism_parameters(workspace, config_path=None, data_name='ALL'):
    """Calculate GSFLOW PRISM Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)
    data_name -- the prism data type (ALL, PPT, TMAX, TMIN, etc.)

    """

    try:
        logging.info('\nGSFLOW PRISM Parameters')

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

        ## PRISM
        prism_ws = inputs_cfg.get('INPUTS', 'prism_folder')
        prism_proj_method = inputs_cfg.get('INPUTS', 'prism_projection_method')
        prism_cs = inputs_cfg.getint('INPUTS', 'prism_cellsize')
        calc_prism_jh_coef_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_prism_jh_coef_flag')

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        ## Check that PRISM folder is valid
        if not os.path.isdir(prism_ws):
            logging.error(
                '\nERROR: PRISM folder ({0}) does not exist'.format(prism_ws))
            raise SystemExit()
        proj_method_list = ['BILINEAR', 'CUBIC', 'NEAREST']
        if prism_proj_method.upper() not in proj_method_list:
            logging.error('\nERROR: PRISM projection method must be: {0}'.format(
                ', '.join(proj_method_list)))
            raise SystemExit()
        logging.debug('  Projection method:    {0}'.format(
            prism_proj_method.upper()))      

        ## Check other inputs
        if prism_cs <= 0:
            logging.error('\nERROR: PRISM cellsize must be greater than 0\n')
            raise SystemExit()

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS 0'
        env.workspace = workspace
        env.scratchWorkspace = hru.scratch_ws      

        ## PRISM data names
        if data_name == 'ALL':
            data_name_list = ['PPT', 'TMAX', 'TMIN']
        else:
            data_name_list = [data_name]

        ## Set month list
        month_list = ['{0:02d}'.format(m) for m in range(1,13)]
        ##month_list.extend(['annual'])

        ## Check fields
        logging.info('\nAdding PRISM fields if necessary')
        for data_name in data_name_list:
            for month in month_list:
                add_field_func(
                    hru.polygon_path, '{0}_{1}'.format(data_name, month), 'DOUBLE')
            
        ## Process each PRISM data type
        logging.info('\nProjecting/clipping PRISM mean monthly rasters')
        for data_name in data_name_list:
            logging.info('\n{0}'.format(data_name))
            ## PRISM input data workspace
            input_ws = os.path.join(prism_ws, data_name.lower())
            if not os.path.isdir(input_ws):
                logging.error('\nERROR: The PRISM {0} folder does not exist'.format(
                    data_name.lower()))
                raise SystemExit()

            ## PRISM output data workspace
            output_ws = os.path.join(
                hru.param_ws, data_name.lower()+'_rasters')
            if not os.path.isdir(output_ws):
                os.mkdir(output_ws)

            ## Remove all non year/month rasters in PRISM temp folder
            logging.info('  Removing existing PRISM files')
            prism_normal_re = re.compile(
                '^PRISM_(%s)_30yr_normal_4kmM2_(\d{2})_bil.*$' % '|'.join(data_name_list),
                re.IGNORECASE)
            for item in os.listdir(output_ws):
                if prism_normal_re.match(item):
                    os.remove(os.path.join(output_ws, item))

            ## Extract, project/resample, clip
            ## Process images by month
            zs_prism_dict = dict()
            ##env.extent = hru.extent
            for month in month_list:
                logging.info('  Month: {0}'.format(month))

                ## Projected/clipped PRISM raster
                input_name = 'PRISM_{0}_30yr_normal_4kmM2_{1}_bil.bil'.format(
                    data_name.lower(), input_month)
                output_name = 'PRISM_{0}_30yr_normal_4kmM2_{1}.img'.format(
                    data_name.lower(), month)
                input_raster = os.path.join(input_ws, input_name)
                output_raster = os.path.join(output_ws, output_name)

                ## Set preferred transforms
                input_sr = Raster(input_raster).spatialReference
                transform_str = transform_func(hru.sr, input_sr)
                if transform_str:
                    logging.debug('  Transform: {0}'.format(transform_str))

                ## Project PRISM rasters to HRU coordinate system
                ## DEADBEEF - Arc10.2 ProjectRaster does not extent
                project_raster_func(
                    input_raster, output_raster, hru.sr,
                    prism_proj_method.upper(), prism_cs, transform_str,
                    '{0} {1}'.format(hru.ref_x, hru.ref_y), input_sr, hru)
                ##arcpy.ProjectRaster_management(
                ##    input_raster, output_raster, hru.sr,
                ##    prism_proj_method.upper(), prism_cs, transform_str,
                ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
                ##    input_sr)

                ## Save parameters for calculating zonal stats
                zs_field = '{0}_{1}'.format(data_name, month)
                zs_prism_dict[zs_field] = [output_raster, 'MEAN']

                ## Cleanup
                del input_raster, output_raster
                del input_name, output_name
                del input_sr, transform_str, zs_field

            ## Cleanup
            ##arcpy.ClearEnvironment('extent')

            ## Calculate zonal statistics
            logging.info('\nCalculating PRISM zonal statistics')
            zonal_stats_func(zs_prism_dict, hru.polygon_path, hru.point_path, hru)
            del zs_prism_dict
        
        ## Jensen-Haise Potential ET air temperature coefficient
        ## Update Jensen-Haise PET estimate using PRISM air temperature
        ## DEADBEEF - First need to figure out month with highest Tmax
        ##            Then get Tmin for same month
        if calc_prism_jh_coef_flag:
            logging.info('\nRe-Calculating JH_COEF_HRU')
            logging.info('  Using PRISM temperature values')
            tmax_field_list = ['!TMAX_{0:02d}!'.format(m) for m in range(1,13)]
            tmin_field_list = ['!TMIN_{0:02d}!'.format(m) for m in range(1,13)]           
            tmax_expr = 'max([{0}])'.format(','.join(tmax_field_list))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.jh_tmax_field, tmax_expr, 'PYTHON')
            ## Sort TMAX and get TMIN for same month
            tmin_expr = 'max(zip([{0}],[{1}]))[1]'.format(
                ','.join(tmax_field_list), ','.join(tmin_field_list))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.jh_tmin_field, tmin_expr, 'PYTHON')    

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
    log_file_name = 'gsflow_prism_normals_log.txt'

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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_prism_parameters')
    del ini_re

    ## Get PRISM data name as argument only if a config file was also set
    try: 
        data_name = sys.argv[2]
    except IndexError:
        data_name = get_prism_data_name()
    if data_name not in ['PPT', 'TMAX', 'TMIN']:
        data_name = 'ALL'

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate GSFLOW PRISM Parameters
    gsflow_prism_parameters(workspace, ini_path, data_name)
