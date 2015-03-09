#--------------------------------
# Name:         hru_parameters.py
# Purpose:      GSFLOW HRU parameters
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-03-08
# Python:       2.7
#--------------------------------

from collections import defaultdict
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

from support_functions import *

################################################################################

def gsflow_hru_parameters(workspace, config_path=None):
    """Calculate GSFLOW HRU Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW HRU Parameters')

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

        ## Read parameters from config file
        study_area_orig_path = inputs_cfg.get('INPUTS', 'study_area_path')
        set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
        if set_lake_flag:
            lake_orig_path = inputs_cfg.get('INPUTS', 'lake_path')
            lake_zone_field = inputs_cfg.get('INPUTS', 'lake_zone_field')
            lake_area_pct = inputs_cfg.getfloat('INPUTS', 'lake_area_pct')

        ## Control flags
        calc_flow_acc_dem_flag = inputs_cfg.getboolean('INPUTS', 'calc_flow_acc_dem_flag')
        calc_topo_index_flag = inputs_cfg.getboolean('INPUTS', 'calc_topo_index_flag')
        clip_root_depth_flag = inputs_cfg.getboolean('INPUTS', 'clip_root_depth_flag')
        ##set_ppt_zones_flag = inputs_cfg.getboolean('INPUTS', 'set_ppt_zones_flag')
        ## Calculate layer thickness and bottoms
        calc_layer_thickness_flag = inputs_cfg.getboolean('INPUTS', 'calc_layer_thickness_flag')


        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist'.format(
                    hru.polygon_path))
            raise SystemExit()
        if set_lake_flag:
            if not arcpy.Exists(lake_orig_path):
                logging.error(
                    '\nERROR: Lake layer ({0}) does not exist'.format(
                        lake_orig_path))
                raise SystemExit()
            ## lake_path must be a polygon shapefile
            if arcpy.Describe(lake_orig_path).datasetType <> 'FeatureClass':
                logging.error(
                    '\nERROR: lake_path must be a polygon shapefile')
                raise SystemExit()
            ## Check lake_zone_field
            if lake_zone_field.upper() in ['', 'FID', 'NONE']:
                lake_zone_field = arcpy.Describe(lake_orig_path).OIDFieldName
                logging.warning(
                    '\n  NOTE: Using {0} to set {1}\n'.format(
                        lake_zone_field, hru.lake_id_field))
            elif not arcpy.ListFields(lake_orig_path, lake_zone_field):
                logging.error(
                    '\nERROR: lake_zone_field field {0} does not exist\n'.format(
                        lake_zone_field))
                raise SystemExit()
            ## Need to check that lake_zone_field is an int type
            elif not [f.type for f in arcpy.Describe(lake_orig_path).fields
                      if (f.name == lake_zone_field and
                          f.type in ['SmallInteger', 'Integer'])]:
                logging.error(
                    '\nERROR: lake_zone_field field {0} must be an integer type\n'.format(
                        lake_zone_field))
                raise SystemExit()      

        ## For now, study area has to be a polygon
        if arcpy.Describe(study_area_orig_path).datasetType <> 'FeatureClass':
            logging.error(
                '\nERROR: For now, study area must be a polygon shapefile')
            raise SystemExit()


        ## Build output folder if necessary
        hru_temp_ws = os.path.join(hru.param_ws, 'hru_temp')
        if not os.path.isdir(hru_temp_ws):
            os.mkdir(hru_temp_ws)
        ## Output paths
        study_area_path = os.path.join(hru_temp_ws, 'study_area.shp')
        lake_path = os.path.join(hru_temp_ws, 'lakes.shp')
        lake_clip_path = os.path.join(hru_temp_ws, 'lake_clip.shp')


        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        env.pyramid = 'PYRAMIDS -1'
        ##env.pyramid = 'PYRAMIDS 0'
        env.workspace = workspace
        env.scratchWorkspace = hru.scratch_ws

        ## Create HRU points at polygon centroids
        if not arcpy.Exists(hru.point_path):
            logging.info('\n  Building HRU point shapefile')
            ## FeatureToPoint will copy all fields in hru.polygon_path
            ##arcpy.FeatureToPoint_management(
            ##    hru.polygon_path, hru.point_path)
            ## Build point_path directly
            arcpy.CreateFeatureclass_management(
                os.path.dirname(hru.point_path),
                os.path.basename(hru.point_path), 'POINT')
            arcpy.DefineProjection_management(hru.point_path, hru.sr)
            arcpy.AddField_management(
                hru.point_path, hru.fid_field, 'LONG')
            hru_centroid_list = [
                row for row in  arcpy.da.SearchCursor(
                    hru.polygon_path, ['OID@', 'SHAPE@XY'])]
            with arcpy.da.InsertCursor(
                hru.point_path, ['OID@', 'SHAPE@XY', hru.fid_field]) as update_c:
                for hru_centroid in hru_centroid_list:
                    update_c.insertRow(
                        [hru_centroid[0], hru_centroid[1], hru_centroid[0]])
            del hru_centroid_list
        ## Check existing HRU points
        else:
            ## Remove any extra fields
            field_remove_list = [
                f.name for f in arcpy.ListFields(hru.point_path)
                if f.name not in ['FID', 'Shape', hru.fid_field]]
            if field_remove_list:
                logging.info('\n  Removing HRU point fields')
                for field in field_remove_list:
                    if field in ['FID', 'Shape', hru.fid_field]:
                        continue
                    logging.debug('    {0}'.format(field))
                    try: arcpy.DeleteField_management(hru.point_path, field)
                    except: continue
            ## Save original FID
            if len(arcpy.ListFields(hru.point_path, hru.fid_field)) == 0:
                arcpy.AddField_management(
                    hru.point_path, hru.fid_field, 'LONG')
            arcpy.CalculateField_management(
                hru.point_path, hru.fid_field, '!FID!', 'PYTHON')
            del field_remove_list

        ## Add all output fields
        logging.info('\nAdding fields if necessary')
        logging.info(
            '  You may see duplicate field names when writing to a network drive')

        ## HRU/DEM Fields
        add_field_func(hru.polygon_path, hru.fid_field, 'LONG')
        add_field_func(hru.polygon_path, hru.id_field, 'LONG')
        add_field_func(hru.polygon_path, hru.type_in_field, 'LONG')
        add_field_func(hru.polygon_path, hru.type_field, 'LONG')
        add_field_func(hru.polygon_path, hru.dem_mean_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_median_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_min_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_max_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_adj_field, 'DOUBLE')
        if calc_flow_acc_dem_flag:
            add_field_func(hru.polygon_path, hru.dem_flowacc_field, 'DOUBLE')
            add_field_func(hru.polygon_path, hru.dem_sum_field, 'DOUBLE')
            add_field_func(hru.polygon_path, hru.dem_count_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_sink8_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.dem_sink4_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.crt_dem_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.crt_fill_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.elev_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.area_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.aspect_field, 'LONG')
        add_field_func(hru.polygon_path, hru.slope_deg_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.slope_rad_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.slope_pct_field, 'DOUBLE')
        if calc_topo_index_flag:
            add_field_func(hru.polygon_path, hru.topo_index_field, 'LONG')
        add_field_func(hru.polygon_path, hru.row_field, 'LONG')
        add_field_func(hru.polygon_path, hru.col_field, 'LONG')
        add_field_func(hru.polygon_path, hru.x_field, 'LONG')
        add_field_func(hru.polygon_path, hru.y_field, 'LONG')
        add_field_func(hru.polygon_path, hru.lat_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.lon_field, 'DOUBLE')

        ## Lake fields
        add_field_func(hru.polygon_path, hru.lake_id_field, 'LONG')
        add_field_func(hru.polygon_path, hru.lake_area_field, 'DOUBLE')       

        ## Stream fields
        add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.irunbound_field, 'LONG')
        add_field_func(hru.polygon_path, hru.flow_dir_field, 'LONG')
        add_field_func(hru.polygon_path, hru.krch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.irch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.jrch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.reach_field, 'LONG')
        add_field_func(hru.polygon_path, hru.rchlen_field, 'LONG')
        add_field_func(hru.polygon_path, hru.maxreach_field, 'LONG')
        add_field_func(hru.polygon_path, hru.outseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.iupseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.subbasin_field, 'LONG')
        add_field_func(hru.polygon_path, hru.segbasin_field, 'LONG')
        add_field_func(hru.polygon_path, hru.strm_top_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.strm_slope_field, 'FLOAT')

        ## PPT Zone fields
        ##if set_ppt_zones_flag:
        add_field_func(hru.polygon_path, hru.ppt_zone_id_field, 'SHORT')
        
        ## DEM based
        ##add_field_func(hru.polygon_path, hru.deplcrv_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.jh_tmax_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.jh_tmin_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.jh_coef_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.snarea_thresh_field, 'FLOAT')

        ## Aspect based
        add_field_func(hru.polygon_path, hru.tmax_adj_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.tmin_adj_field, 'FLOAT')

        ## Vegetation fields
        add_field_func(hru.polygon_path, hru.cov_type_field, 'SHORT')
        add_field_func(hru.polygon_path, hru.covden_sum_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.covden_win_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.rad_trncf_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.snow_intcp_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.srain_intcp_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.wrain_intcp_field, 'FLOAT')

        ## Soil fields
        add_field_func(hru.polygon_path, hru.awc_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.clay_pct_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.sand_pct_field, 'FLOAT')
        ##add_field_func(hru.polygon_path, hru.silt_pct_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.ksat_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.soil_depth_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.root_depth_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.soil_type_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.moist_init_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.moist_max_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.rechr_init_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.rechr_max_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.ssr2gw_rate_field, 'FLOAT')
        ##add_field_func(hru.polygon_path, hru.pref_flow_den_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.slowcoef_lin_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.slowcoef_sq_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.fastcoef_lin_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.fastcoef_sq_field, 'FLOAT')

        ## Impervious fields
        add_field_func(hru.polygon_path, hru.imperv_pct_field, 'FLOAT')
        ##add_field_func(hru.polygon_path, hru.carea_min_field, 'FLOAT')
        add_field_func(hru.polygon_path, hru.carea_max_field, 'FLOAT')

        ## PRISM mean monthly fields
        month_list = ['{0:02d}'.format(m) for m in range(1,13)]
        month_list.extend(['14'])
        for prism_data_name in ['PPT', 'TMAX', 'TMIN']:
            for month in month_list:
                add_field_func(
                    hru.polygon_path,
                    '{0}_{1}'.format(prism_data_name, month), 'FLOAT')
        ## PRISM mean monthly PPT ratio fields
        for month in month_list:
            if month == '14':
                continue
            add_field_func(
                hru.polygon_path, 'PPT_RT_{0}'.format(month), 'FLOAT')

        ## Layer thickness and bottom fields
        if calc_layer_thickness_flag:
            add_field_func(hru.polygon_path, hru.alluv_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.alluv_thick_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay1_thick_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay2_thick_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay3_thick_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay4_thick_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay1_bottom_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay2_bottom_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay3_bottom_field, 'FLOAT')
            add_field_func(hru.polygon_path, hru.lay4_bottom_field, 'FLOAT')


        ## Id field is added by default to new fishnets
        if arcpy.ListFields(hru.polygon_path, 'Id'):
            arcpy.DeleteField_management(hru.polygon_path, 'Id')

        logging.info('\nCalculating parameters')
        ## Keep original FID for subsetting in zonal stats
        logging.info('  Saving original HRU FID to {0}'.format(
            hru.fid_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.fid_field, '!FID!', 'PYTHON')

        ## Cell X/Y
        logging.info('  Calculating cell X/Y')
        cell_xy_func(
            hru.polygon_path, hru.x_field, hru.y_field)

        ## Create unique ID, start at top left corner, work down rows
        ## Row/Col numbered from top left corner (1's based numbering)
        logging.info('  Calculating cell ID/row/col')
        cell_id_col_row_func(
            hru.polygon_path, hru.id_field, hru.col_field, hru.row_field,
            hru.extent, hru.cs)

        ## Cell Lat/Lon
        logging.info('  Calculating cell lat/lon')
        cell_lat_lon_func(
            hru.polygon_path, hru.lat_field, hru.lon_field, hru.sr.GCS)

        ## Cell Area
        logging.info('  Calculating cell area (acres)')
        arcpy.CalculateField_management(
            hru.polygon_path, hru.area_field, '!SHAPE.AREA@acres!', 'PYTHON')

        ## Reset HRUTYPE_IN / HRU_TYPE 
        logging.info('\nResetting {0} to 0'.format(hru.type_in_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.type_in_field, 0, 'PYTHON')
        logging.info('Resetting {0} to 0'.format(hru.type_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.type_field, 0, 'PYTHON')
        ## Reset LAKE_ID
        if set_lake_flag:
            logging.info('Resetting {0} to 0'.format(hru.lake_id_field))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.lake_id_field, 0, 'PYTHON')
            logging.info('Resetting {0} to 0'.format(hru.lake_area_field))
            arcpy.CalculateField_management(
                hru.polygon_path, hru.lake_area_field, 0, 'PYTHON')

        ## Calculate HRU Type
        logging.info('\nCalculating cell HRU Type')
        study_area_desc = arcpy.Describe(study_area_orig_path)
        study_area_sr = study_area_desc.spatialReference
        logging.debug('  Study area: {0}'.format(study_area_orig_path))
        logging.debug('  Study area spat. ref.:  {0}'.format(
            study_area_sr.name))
        logging.debug('  Study area GCS:         {0}'.format(
            study_area_sr.GCS.name))
        ## If study area spat_ref doesn't match hru_param spat_ref
        ## Project study area to hru_param spat ref
        ## Otherwise, read study_area directly       
        if hru.sr.name <> study_area_sr.name:
            logging.info('  Projecting study area...')
            ## Set preferred transforms
            transform_str = transform_func(hru.sr, study_area_sr)
            logging.debug('    Transform: {0}'.format(transform_str))
            ## Project study area shapefile
            arcpy.Project_management(
                study_area_orig_path, study_area_path, hru.sr,
                transform_str, study_area_sr)
            del transform_str
        else:
            arcpy.Copy_management(study_area_orig_path, study_area_path)
        zone_by_centroid_func(
            study_area_path, hru.type_in_field, 1,
            hru.polygon_path, hru.point_path, hru)

        ## Calculate HRU Type for lakes (HRU_TYPE = 2)
        if set_lake_flag:
            logging.info('\nCalculating cell HRU Type & ID for lakes')
            lake_layer = 'lake_layer'
            lake_desc = arcpy.Describe(lake_orig_path)
            lake_sr = lake_desc.spatialReference
            logging.debug('  Lakes: {0}'.format(lake_orig_path))
            logging.debug('  Lakes spat. ref.:  {0}'.format(lake_sr.name))
            logging.debug('  Lakes GCS:         {0}'.format(lake_sr.GCS.name))
            ## If lakes spat_ref doesn't match hru_param spat_ref
            ## Project lakes to hru_param spat ref
            ## Otherwise, read lakes directly       
            if hru.sr.name <> lake_sr.name:
                logging.info('  Projecting lakes...')
                ## Set preferred transforms
                transform_str = transform_func(hru.sr, lake_sr)
                logging.debug('    Transform: {0}'.format(transform_str))
                ## Project lakes shapefile
                arcpy.Project_management(
                    lake_orig_path, lake_path, hru.sr, transform_str, lake_sr)
                arcpy.MakeFeatureLayer_management(lake_path, lake_layer)
                del lake_path, transform_str
            else:
                arcpy.MakeFeatureLayer_management(
                    lake_orig_path, lake_layer)
            ## Clip lakes by study area after projecting lakes
            logging.info('  Clipping lakes...')
            arcpy.Clip_analysis(lake_layer, study_area_path, lake_clip_path)
            ## Remove all unnecesary fields
            for field in arcpy.ListFields(lake_clip_path):
                if field.name not in [lake_zone_field, 'Shape']:
                    try: arcpy.DeleteField_management(lake_clip_path, field.name)
                    except: pass
            ## Set lake HRU_TYPE
            logging.info('  Setting lake {0}'.format(hru.type_in_field))
            zone_by_area_func(
                lake_clip_path, hru.type_in_field, 2, 
                hru.polygon_path, hru, hru.area_field,
                hru.lake_area_field, lake_area_pct)
            ## Set lake ID
            logging.info('  Setting {0}'.format(hru.lake_id_field))
            zone_by_area_func(
                lake_clip_path, hru.lake_id_field, lake_zone_field, 
                hru.polygon_path, hru, hru.area_field,
                hru.lake_area_field, lake_area_pct)
            ## Cleanup
            del lake_layer, lake_desc, lake_sr

        ## Copy HRUTYPE_IN for HRU_TYPE
        ## HRU_TYPE can then be modified by later scripts
        logging.info('\nCopying {0} to {1}'.format(
            hru.type_in_field, hru.type_field))
        arcpy.CalculateField_management(
            hru.polygon_path, hru.type_field,
            '!{0}!'.format(hru.type_in_field), 'PYTHON')

        ## Cleanup
        del study_area_desc, study_area_sr   

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        ##pass
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()
        ##try: logging.info('\nTime elapsed = {0} seconds\n'.format(clock()-start))
        ##except: pass

################################################################################

def cell_xy_func(hru_param_path, x_field, y_field):
    fields = ('SHAPE@XY', x_field, y_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields) as u_cursor:
        for row in u_cursor:
            row[1], row[2] = row[0]
            u_cursor.updateRow(row)
            del row

def cell_lat_lon_func(hru_param_path, lat_field, lon_field, gcs_sr):
    fields = ('SHAPE@XY', lon_field, lat_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields, '', gcs_sr) as u_cursor:
        for row in u_cursor:
            row[1], row[2] = row[0]
            u_cursor.updateRow(row)
            del row

def cell_id_col_row_func(
    hru_param_path, id_field, col_field, row_field, extent, cs):
    fields = ('SHAPE@XY', col_field, row_field, id_field)
    with arcpy.da.UpdateCursor(
        hru_param_path, fields) as u_cursor:
        num_cols = (extent.XMax - extent.XMin) / cs
        for row in u_cursor:
            #### Row/Col are 1's based indices
            row[1] = ((row[0][0] - extent.XMin) // cs) + 1
            row[2] = ((extent.YMax - row[0][1]) // cs) + 1
            #### Create unique ID, start at top left corner, work down rows
            row[3] = row[1] + (row[2] - 1) * num_cols
            u_cursor.updateRow(row)
            del row

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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_hru_parameters')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate GSFLOW HRU Parameters
    gsflow_hru_parameters(workspace, ini_path)
