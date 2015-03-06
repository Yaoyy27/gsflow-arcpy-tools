#--------------------------------
# Name:         gsflow_fishnet_generator.py
# Purpose:      GSFLOW fishnet generator
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2014-10-13
# Python:       2.7
#--------------------------------

import ConfigParser
import datetime as dt
import logging
import math
import os
import re
import sys

import arcpy
from arcpy import env
from arcpy.sa import *
##import numpy as np

from gsflow_support_functions import *

################################################################################

def gsflow_fishnet_func(workspace, config_path=None):
    """GSFLOW Fishnet Generator

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW Fishnet Generator')

        ## Initialize hru parameters class
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

        ## Check input paths
        study_area_path = inputs_cfg.get('INPUTS', 'study_area_path')
        if not arcpy.Exists(study_area_path):
            logging.error(
                '\nERROR: Study area ({0}) does not exist'.format(
                    study_area_path))
            raise SystemExit()

        ## For now, study area has to be a polygon
        if arcpy.Describe(study_area_path).datasetType <> 'FeatureClass':
            logging.error('\nERROR: For now, study area must be a polygon shapefile')
            raise SystemExit()           

        ## Build output folder if necessary
        fishnet_temp_ws = os.path.join(hru.param_ws, 'fishnet_temp')
        if not os.path.isdir(fishnet_temp_ws):
            os.mkdir(fishnet_temp_ws)
        ## Output paths
        study_area_proj_path = os.path.join(
            fishnet_temp_ws, 'projected_study_area.shp')

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = workspace
        env.scratchWorkspace = hru.scratch_ws

        ## Get spatial reference of study_area
        study_area_desc = arcpy.Describe(study_area_path)
        study_area_sr = study_area_desc.spatialReference
        logging.debug('\n  Study area: {0}'.format(study_area_path))
        logging.debug('  Study area spat. ref.:  {0}'.format(
            study_area_sr.name))
        logging.debug('  Study area GCS:         {0}'.format(
            study_area_sr.GCS.name))

        ## Set spatial reference of hru shapefile
        ## If the spatial reference can be set to an int, then use the value
        try: hru.sr_name = int(hru.sr_name)
        except ValueError: pass
        hru.sr = arcpy.SpatialReference(hru.sr_name)
        logging.debug('  HRU spat. ref.: {0}'.format(hru.sr.name))
        logging.debug('  HRU GCS:        {0}'.format(hru.sr.GCS.name))

        ## If study area spat_ref doesn't match hru_param spat_ref
        ## Project study are to hru_param and get projected extent
        ## Otherwise, read study_area extent directly       
        if hru.sr.name <> study_area_sr.name:
            logging.info('\n  Projecting study area...')
            ## Set preferred transforms
            transform_str = transform_func(hru.sr, study_area_sr)
            logging.debug('    Transform: {0}'.format(transform_str))
            ## Project study area
            arcpy.Project_management(
                study_area_path, study_area_proj_path, hru.sr,
                transform_str, study_area_sr)
            study_area_extent = arcpy.Describe(study_area_proj_path).extent
            arcpy.Delete_management(study_area_proj_path)
            logging.info('\n  Projected extent:  {0}'.format(
                extent_string(study_area_extent)))
            del study_area_proj_path, transform_str
        else:
            study_area_extent = arcpy.Describe(study_area_path).extent
            logging.info('\n  Study Area extent: {0}'.format(
                extent_string(study_area_extent)))

        ## Buffer extent
        buffer_extent = buffer_extent_func(
            study_area_extent, hru.buffer_cells * hru.cs)
        logging.info('  HRU extent:        {0}'.format(
            extent_string(buffer_extent)))

        ## Adjust study area extent to reference points
        hru.ref_pnt = arcpy.Point(hru.ref_x, hru.ref_y)
        hru.extent = adjust_extent_to_snap(
            buffer_extent, hru.ref_pnt, hru.cs, hru.snap_method)
        logging.info('  Snapped Extent:    {0}'.format(
            extent_string(hru.extent)))

        ## Build hru_param
        logging.info('\nBuilding HRU parameter fishnet')
        build_fishnet_func(
            hru.polygon_path, hru.point_path, hru.extent, hru.cs, hru.sr)

        ## Write initial parameters to hru_param (X/Y, ROW/COL, Unique ID)
        ##set_hru_id_func(hru.polygon_path, hru.extent, hru.cs)

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()

################################################################################

def build_fishnet_func(hru_polygon_path, hru_point_path, extent, cs, sr):
    ## Remove existing
    if arcpy.Exists(hru_polygon_path):
        arcpy.Delete_management(hru_polygon_path)
    if arcpy.Exists(hru_point_path):
        arcpy.Delete_management(hru_point_path)
    ## Calculate LL/UR corner points
    origin_pnt = (extent.XMin, extent.YMin)
    yaxis_pnt = (extent.XMin, extent.YMin+cs)
    corner_pnt = (extent.XMax, extent.YMax)
    origin_str = ' '.join(map(str, origin_pnt))
    yaxis_str = ' '.join(map(str, yaxis_pnt))
    corner_str = ' '.join(map(str, corner_pnt))
    logging.debug('  Origin: {0}'.format(origin_str))
    logging.debug('  Y-Axis: {0}'.format(yaxis_str))
    logging.debug('  Corner: {0}'.format(corner_str))
    ## Build fishnet & labels
    arcpy.CreateFishnet_management(
        hru_polygon_path, origin_str, yaxis_str, cs, cs,
        '0', '0', corner_str, 'LABELS', '#', 'POLYGON')
    arcpy.DefineProjection_management(hru_polygon_path, sr)
    arcpy.DefineProjection_management(hru_point_path, sr)

################################################################################
if __name__ == '__main__':
    workspace = os.getcwd()
    log_ws = os.path.join(workspace, 'logs')
    if not os.path.isdir(log_ws): 
        os.mkdir(log_ws)
    log_file_name = 'gsflow_fishnet_log.txt'

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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_fishnet_generator')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Run GSFLOW Fishnet Generator
    gsflow_fishnet_func(workspace, ini_path)
