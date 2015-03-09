#--------------------------------
# Name:         gsflow_prism_800m_custom_normals.py
# Purpose:      GSFLOW PRISM parameters from custom 800m normals
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2014-10-13
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

from gsflow_support_functions import *

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
        start_year = inputs_cfg.getint('INPUTS', 'prism_start_year')
        end_year = inputs_cfg.getint('INPUTS', 'prism_end_year')
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
        try:
            start_year = int(start_year)
        except ValueError:
            logging.error(
                '\nERROR: start_year must be an integer from {0} to {1}\n'.format(
                year_list[0], year_list[-1]))
            raise SystemExit()
        try:
            end_year = int(end_year)
        except ValueError:
            logging.error(
                '\nERROR: end_year must be an integer from {0} to {1}\n'.format(
                year_list[0], year_list[-1]))
            raise SystemExit()
        if end_year < start_year:
            logging.error('\nERROR: end_year must greater than start_year\n')
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

        ## Set month list based on flags
        month_list = ['{0:02d}'.format(m) for m in range(1,13)]
        month_list.extend(['14'])

        ## Check fields
        logging.info('\nAdding PRISM fields if necessary')
        for data_name in data_name_list:
            for month in month_list:
                add_field_func(
                    hru.polygon_path, '{0}_{1}'.format(data_name, month), 'DOUBLE')
            
        ## Process each PRISM data type
        logging.info('\nCalculating PRISM mean monthly rasters')
        for data_name in data_name_list:
            logging.info('\n{0}'.format(data_name))
            ## PRISM input data workspace
            data_ws = os.path.join(prism_ws, data_name.lower())
            if not os.path.isdir(data_ws):
                logging.error('\nERROR: The PRISM {0} folder does not exist'.format(
                    data_name.lower()))
                raise SystemExit()
            ## PRISM output data workspace
            prism_temp_ws = os.path.join(hru.param_ws, data_name.lower()+'_rasters')
            if not os.path.isdir(prism_temp_ws): os.mkdir(prism_temp_ws)
            ## Remove all non year/month rasters in PRISM temp folder
            logging.debug('  Removing existing intermediate files')
            prism_mean_re = re.compile('^us_(PPT|TMAX|TMIN)_\d{2}\w+.img\w*$')
            for item in os.listdir(prism_temp_ws):
                if prism_mean_re.match(item):
                    os.remove(os.path.join(prism_temp_ws, item))
                
            ## Build data file/folder lists
            data_folder_re = re.compile('^us_%s_\d{4}$' % (data_name.lower()))
            data_folder_dict = dict(
                [(int(item.split('_')[-1]), os.path.join(data_ws, item))
                 for item in os.listdir(data_ws)
                 if data_folder_re.match(item)])
            if data_folder_dict:
                year_list = sorted(data_folder_dict.keys())
            else:
                logging.error('\nERROR: ')
                raise SystemExit()

            ## Check start/end year
            if start_year not in year_list:
                logging.error(
                    '\nERROR: start_year must be an integer from {0} to {1}'.format(
                        year_list[0], year_list[-1]))
                raise SystemExit()
            elif end_year not in year_list :
                logging.error(
                    '\nERROR: end_year must be an integer from {0} to {1}'.format(
                        year_list[0], year_list[-1]))
                raise SystemExit()
            year_sub_list = [item for item in year_list
                             if start_year <= item <= end_year]
            logging.info('  Start Year: {0}'.format(start_year))
            logging.info('  End Year:   {0}'.format(end_year))

            ## Save list of rasters for each water year
            prism_raster_re = re.compile(
                '^us_%s_\d{4}_\w{2}.img$' % (data_name.lower()))
            raster_dict = defaultdict(dict)
            for year in year_list:
                year_ws = os.path.join(
                    data_ws, 'us_{0}_{1:04d}'.format(data_name.lower(), year))
                for raster_name in os.listdir(year_ws):
                    if not prism_raster_re.match(raster_name): continue
                    raster_path = os.path.join(year_ws, raster_name)
                    prefix, ptype, year, month = os.path.splitext(
                        raster_name)[0].split('_')
                    if month in month_list and int(year) in year_sub_list:
                        raster_dict[month][int(year)] = raster_path
                    del raster_path, prefix, ptype, year, month
                del year_ws

            ## Extract, project/resample, clip
            ## Process images by month
            ##env.extent = hru.extent
            ##env.outputCoordinateSystem = hru.sr
            ##env.cellsize = prism_cs
            for month in sorted(raster_dict.keys()):
                logging.info('  Month: {0}'.format(month))
                month_time = clock()

                ## Initialize sum/count rasters
                sum_obj, count_obj = 0, 0                

                ## For each month, process all years
                for year, prism_orig_raster in sorted(raster_dict[month].items()):
                    logging.debug('    Year: {0:04d}'.format(year))

                    ## Projected/clipped PRISM raster
                    output_raster = os.path.join(
                        'in_memory',
                        os.path.basename(input_raster).replace('.img', ''))

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
                        '{0} {1}'.format(hru.ref_x, hru.ref_y), 
						input_sr, hru)
                    ##arcpy.ProjectRaster_management(
                    ##    input_raster, output_raster, hru.sr,
                    ##    prism_proj_method.upper(), prism_cs, transform_str,
                    ##    '{0} {1}'.format(hru.ref_x, hru.ref_y),
                    ##    input_sr)

                    ## Track running sum/count for mean calculation
                    sum_obj += Raster(output_raster)
                    count_obj += ~IsNull(Raster(output_raster))
                    del input_raster, output_raster
                        
                ## Calculate mean
                mean_name = 'us_{0}_{1}_mean.img'.format(data_name.lower(), month)
                mean_raster = os.path.join(prism_temp_ws, mean_name)
                mean_obj = sum_obj / count_obj
                mean_obj.save(mean_raster)

                ## Cleanup
                del mean_obj, mean_raster, mean_name
                del sum_obj, count_obj
                ##print '    {0}'.format(clock() - month_time)

            ##arcpy.ClearEnvironment('extent')
            ##arcpy.ClearEnvironment('outputCoordinateSystem')
            ##arcpy.ClearEnvironment('cellsize')

            ## List of rasters, fields, and stats for zonal statistics
            zs_prism_dict = dict()
            for month in month_list:
                mean_raster = os.path.join(
                    prism_temp_ws, 'us_{0}_{1}_mean.img'.format(
                        data_name.lower(), month))
                mean_field =  '{0}_{1}'.format(data_name, month)
                zs_prism_dict[mean_field] = [mean_raster, 'MEAN']
                del mean_raster, mean_field

            ## Calculate zonal statistics
            logging.info('\nCalculating PRISM zonal statistics')
            zonal_stats_func(
                zs_prism_dict, hru.polygon_path, hru_point_path, hru)
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
                hru.polygon_path, jh_tmax_field, tmax_expr, 'PYTHON')
            ## Sort TMAX and get TMIN for same month
            tmin_expr = 'max(zip([{0}],[{1}]))[1]'.format(
                ','.join(tmax_field_list), ','.join(tmin_field_list))
            arcpy.CalculateField_management(
                hru.polygon_path, jh_tmin_field, tmin_expr, 'PYTHON')    

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
    log_file_name = 'gsflow_prism_custom_log.txt'

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
    except 
        IndexError: data_name = get_prism_data_name()
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
