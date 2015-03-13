#--------------------------------
# Name:         support_functions.py
# Purpose:      GSFLOW parameter support functions
# Notes:        ArcGIS 10.2 Version
# Author:       Charles Morton
# Created       2015-03-13
# Python:       2.7
#--------------------------------

from collections import defaultdict
import ConfigParser
import logging
import math
import os
import re
from time import sleep

import numpy as np

import arcpy
from arcpy import env
from arcpy.sa import *

from support_functions import *

################################################################################

class hru_parameters():
    def __init__(self, config_path):
        ## Open input parameter config file
        inputs_cfg = ConfigParser.ConfigParser()
        try:
            inputs_cfg.readfp(open(config_path))
        except IOError:
            logging.error(('\nERROR: Config file does not exist\n'+
                           '  {0}\n').format(field_list_path))
            raise SystemExit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error('\nERROR: Config file is missing a section header\n'+
                          '    Please make sure the following line is at the '+
                          'beginning of the file\n[INPUTS]\n')
            raise SystemExit()
        except:
            logging.error(('\nERROR: Config file could not be read\n'+
                           '  {0}\n').format(config_path))
        logging.debug('\nReading Input File')
        
        ## Open field list config file
        field_list_path =  inputs_cfg.get('INPUTS', 'field_list_path')
        fields_cfg = ConfigParser.ConfigParser()
        try:
            fields_cfg.readfp(open(field_list_path))
        except IOError:
            logging.error(('\nERROR: Field list file does not exist\n'+
                           '  {0}\n').format(field_list_path))
            raise SystemExit()
        except ConfigParser.MissingSectionHeaderError:
            logging.error('\nERROR: Field list file is missing a section header\n'+
                          '    Please make sure the following line is at the '+
                          'beginning of the file\n[FIELDS]\n')
            raise SystemExit()
        except:
            logging.error(('\nERROR: Field list file could not be read\n'+
                           '  {0}\n').format(field_list_path))
        logging.debug('\nReading Field List File')

        ## Read parameters from config file
        self.polygon_path = inputs_cfg.get('INPUTS', 'hru_fishnet_path')
        self.point_path   = inputs_cfg.get('INPUTS', 'hru_centroid_path')
        self.sr_name = inputs_cfg.get('INPUTS', 'hru_projection')
        self.cs      = inputs_cfg.getint('INPUTS', 'hru_cellsize')
        self.ref_x   = inputs_cfg.getfloat('INPUTS', 'hru_ref_x')
        self.ref_y   = inputs_cfg.getfloat('INPUTS', 'hru_ref_y')
        ##self.ref_x   = inputs_cfg.getint('INPUTS', 'hru_ref_x')
        ##self.ref_y   = inputs_cfg.getint('INPUTS', 'hru_ref_y')
        self.buffer_cells = inputs_cfg.getint('INPUTS', 'hru_buffer_cells')
        self.snap_method = inputs_cfg.get('INPUTS', 'hru_param_snap_method')
        self.fid_field = inputs_cfg.get('INPUTS', 'orig_fid_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')

        ## Check inputs
        if self.cs <= 0:
            logging.error('\nERROR: Fishnet cellsize must be greater than 0')
            raise SystemExit()
        if self.buffer_cells < 0:
            logging.error('\nERROR: Buffer cells must be greater than 0')
            raise SystemExit()

        ##
        self.param_ws = inputs_cfg.get('INPUTS', 'parameter_folder')
        if not os.path.isdir(self.param_ws):
            os.mkdir(self.param_ws)

        ## Scratch workspace
        try:
            scratch_name = inputs_cfg.get('INPUTS', 'scratch_name')
        except:
            scratch_name = 'in_memory'
        if scratch_name == 'in_memory':
            self.scratch_ws = scratch_name
        else:
            scratch_ws = os.path.join(self.param_ws, scratch_name)
            if not os.path.isdir(scratch_ws):
                os.mkdir(scratch_ws)
            self.scratch_ws = scratch_ws

        ## Log input hru parameters
        logging.info('  Fishnet cellsize:   {0}'.format(self.cs))
        logging.info('  Fishnet ref. point: {0} {1}'.format(self.ref_x, self.ref_y))
        logging.debug('  Fishnet Buffer Cells: {0}'.format(self.buffer_cells))
        ##snap_pnt = arcpy.Point(self.ref_x, self.ref_y)
        ##snap_pnt.X, snap_pnt.Y = self.ref_x, self.ref_y
        
        ## Set spatial reference of hru shapefile
        if arcpy.Exists(self.polygon_path):
            hru_desc = arcpy.Describe(self.polygon_path)
            self.sr = hru_desc.spatialReference
            self.extent = adjust_extent_to_snap(
                hru_desc.extent, arcpy.Point(self.ref_x, self.ref_y), self.cs, 'ROUND')
            ##self.extent = adjust_extent_to_snap(
            ##    hru_param_desc.extent, snap_pnt, self.cs, 'ROUND')
            logging.debug('  Fishnet spat. ref.: {0}'.format(self.sr.name))
            logging.debug('  Fishnet GCS:        {0}'.format(self.sr.GCS.name))
            logging.info('  Fishnet extent:     {0}'.format(extent_string(self.extent)))          

        ## Some fields are dependent on the control flags
        set_lake_flag = inputs_cfg.getboolean('INPUTS', 'set_lake_flag')
        calc_flow_acc_dem_flag = inputs_cfg.getboolean('INPUTS', 'calc_flow_acc_dem_flag')
        calc_topo_index_flag = inputs_cfg.getboolean('INPUTS', 'calc_topo_index_flag')
        clip_root_depth_flag = inputs_cfg.getboolean('INPUTS', 'clip_root_depth_flag')
        ##set_ppt_zones_flag = inputs_cfg.getboolean('INPUTS', 'set_ppt_zones_flag')
        calc_layer_thickness_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_layer_thickness_flag')

        ## Read in all field names
        self.id_field = fields_cfg.get('FIELDS', 'id_field')
        self.type_in_field = fields_cfg.get('FIELDS', 'type_in_field')
        self.type_field = fields_cfg.get('FIELDS', 'type_field')
        self.dem_mean_field = fields_cfg.get('FIELDS', 'dem_mean_field')
        self.dem_median_field = fields_cfg.get('FIELDS', 'dem_median_field')
        self.dem_max_field = fields_cfg.get('FIELDS', 'dem_max_field')
        self.dem_min_field = fields_cfg.get('FIELDS', 'dem_min_field')
        self.dem_adj_field = fields_cfg.get('FIELDS', 'dem_adj_field')
        if calc_flow_acc_dem_flag:
            ##self.dem_sum_field = 'DEM_SUM'
            ##self.dem_count_field = 'DEM_COUNT'
            self.dem_sum_field = fields_cfg.get('FIELDS', 'dem_sum_field')
            self.dem_count_field = fields_cfg.get('FIELDS', 'dem_count_field')
            self.dem_flowacc_field = fields_cfg.get('FIELDS', 'dem_flowacc_field')
        else:
            self.dem_sum_field = 'DEM_SUM'
            self.dem_count_field = 'DEM_COUNT'
            self.dem_flowacc_field = 'DEM_FLOW_AC'
        self.dem_sink8_field = fields_cfg.get('FIELDS', 'dem_sink8_field')
        self.dem_sink4_field = fields_cfg.get('FIELDS', 'dem_sink4_field')
        self.crt_dem_field  = fields_cfg.get('FIELDS', 'crt_dem_field')
        self.crt_fill_field = fields_cfg.get('FIELDS', 'crt_fill_field')
        self.area_field = fields_cfg.get('FIELDS', 'area_field')
        self.elev_field = fields_cfg.get('FIELDS', 'elev_field')
        self.aspect_field = fields_cfg.get('FIELDS', 'aspect_field')
        self.slope_deg_field = fields_cfg.get('FIELDS', 'slope_deg_field')
        self.slope_rad_field = fields_cfg.get('FIELDS', 'slope_rad_field')
        self.slope_pct_field = fields_cfg.get('FIELDS', 'slope_pct_field')
        self.topo_index_field = fields_cfg.get('FIELDS', 'topo_index_field')
        self.row_field = fields_cfg.get('FIELDS', 'row_field')
        self.col_field = fields_cfg.get('FIELDS', 'col_field')
        self.x_field = fields_cfg.get('FIELDS', 'x_field')
        self.y_field = fields_cfg.get('FIELDS', 'y_field')
        self.lat_field = fields_cfg.get('FIELDS', 'lat_field')
        self.lon_field = fields_cfg.get('FIELDS', 'lon_field')

        if set_lake_flag:
            self.lake_id_field = fields_cfg.get('FIELDS', 'lake_id_field')
            self.lake_area_field = fields_cfg.get('FIELDS', 'lake_area_field')
        else:
            self.lake_id_field = 'LAKE_ID'
            self.lake_area_field = 'LAKE_AREA'            

        ## DEM based
        ##self.deplcrv_field = fields_cfg.get('FIELDS', 'deplcrv_field')
        self.jh_tmax_field = fields_cfg.get('FIELDS', 'jh_tmax_field')
        self.jh_tmin_field = fields_cfg.get('FIELDS', 'jh_tmin_field')
        self.jh_coef_field = fields_cfg.get('FIELDS', 'jh_coef_field')
        self.snarea_thresh_field = fields_cfg.get('FIELDS', 'snarea_thresh_field')
        self.tmax_adj_field = fields_cfg.get('FIELDS', 'tmax_adj_field')
        self.tmin_adj_field = fields_cfg.get('FIELDS', 'tmin_adj_field')

        ## Vegetation
        self.cov_type_field = fields_cfg.get('FIELDS', 'cov_type_field')
        self.covden_sum_field = fields_cfg.get('FIELDS', 'covden_sum_field')
        self.covden_win_field = fields_cfg.get('FIELDS', 'covden_win_field')
        self.snow_intcp_field = fields_cfg.get('FIELDS', 'snow_intcp_field')
        self.wrain_intcp_field = fields_cfg.get('FIELDS', 'wrain_intcp_field')
        self.srain_intcp_field = fields_cfg.get('FIELDS', 'srain_intcp_field')
        self.rad_trncf_field = fields_cfg.get('FIELDS', 'rad_trncf_field')

        ## Soil
        self.awc_field = fields_cfg.get('FIELDS', 'awc_field')
        self.clay_pct_field = fields_cfg.get('FIELDS', 'clay_pct_field')
        self.sand_pct_field = fields_cfg.get('FIELDS', 'sand_pct_field')
        ##self.silt_pct_field = fields_cfg.get('FIELDS', 'silt_pct_field')
        self.ksat_field = fields_cfg.get('FIELDS', 'ksat_field')
        self.soil_depth_field = fields_cfg.get('FIELDS', 'soil_depth_field')
        self.root_depth_field = fields_cfg.get('FIELDS', 'root_depth_field')
        self.soil_type_field  = fields_cfg.get('FIELDS', 'soil_type_field')
        self.moist_init_field = fields_cfg.get('FIELDS', 'moist_init_field')
        self.moist_max_field  = fields_cfg.get('FIELDS', 'moist_max_field')
        self.rechr_init_field = fields_cfg.get('FIELDS', 'rechr_init_field')
        self.rechr_max_field  = fields_cfg.get('FIELDS', 'rechr_max_field')
        self.ssr2gw_rate_field = fields_cfg.get('FIELDS', 'ssr2gw_rate_field')
        ##self.pref_flow_den_field = fields_cfg.get('FIELDS', 'pref_flow_den_field')
        self.slowcoef_lin_field = fields_cfg.get('FIELDS', 'slowcoef_lin_field')
        self.slowcoef_sq_field  = fields_cfg.get('FIELDS', 'slowcoef_sq_field')
        self.fastcoef_lin_field = fields_cfg.get('FIELDS', 'fastcoef_lin_field')
        self.fastcoef_sq_field  = fields_cfg.get('FIELDS', 'fastcoef_sq_field')

        ## Impervious Parameter Fields
        self.imperv_pct_field = fields_cfg.get('FIELDS', 'imperv_pct_field')
        ##carea_min_field  = fields_cfg.get('FIELDS', 'carea_min_field')
        self.carea_max_field  = fields_cfg.get('FIELDS', 'carea_max_field')

        ## Streams
        self.irunbound_field  = fields_cfg.get('FIELDS', 'irunbound_field')
        self.iseg_field       = fields_cfg.get('FIELDS', 'iseg_field')
        self.flow_dir_field   = fields_cfg.get('FIELDS', 'flow_dir_field')
        self.krch_field       = fields_cfg.get('FIELDS', 'krch_field')
        self.irch_field       = fields_cfg.get('FIELDS', 'irch_field')
        self.jrch_field       = fields_cfg.get('FIELDS', 'jrch_field')
        self.reach_field      = fields_cfg.get('FIELDS', 'reach_field')
        self.rchlen_field     = fields_cfg.get('FIELDS', 'rchlen_field')
        self.maxreach_field   = fields_cfg.get('FIELDS', 'maxreach_field')
        self.outseg_field     = fields_cfg.get('FIELDS', 'outseg_field')
        self.iupseg_field     = fields_cfg.get('FIELDS', 'iupseg_field')
        self.strm_top_field   = fields_cfg.get('FIELDS', 'strm_top_field')
        self.strm_slope_field = fields_cfg.get('FIELDS', 'strm_slope_field')
        self.subbasin_field   = fields_cfg.get('FIELDS', 'subbasin_field')
        self.segbasin_field   = fields_cfg.get('FIELDS', 'segbasin_field')

        ##if set_ppt_zones_flag:
        self.ppt_zone_id_field = fields_cfg.get('FIELDS', 'ppt_zone_id_field')

        ## Calculate layer thickness and bottoms
        if calc_layer_thickness_flag:
            self.alluv_field = fields_cfg.get('FIELDS', 'alluv_field')
            self.alluv_thick_field = fields_cfg.get('FIELDS', 'alluv_thick_field')
            self.lay1_thick_field = fields_cfg.get('FIELDS', 'lay1_thick_field')
            self.lay2_thick_field = fields_cfg.get('FIELDS', 'lay2_thick_field')
            self.lay3_thick_field = fields_cfg.get('FIELDS', 'lay3_thick_field')
            self.lay4_thick_field = fields_cfg.get('FIELDS', 'lay4_thick_field')
            self.lay1_bottom_field = fields_cfg.get('FIELDS', 'lay1_bottom_field')
            self.lay2_bottom_field = fields_cfg.get('FIELDS', 'lay2_bottom_field')
            self.lay3_bottom_field = fields_cfg.get('FIELDS', 'lay3_bottom_field')
            self.lay4_bottom_field = fields_cfg.get('FIELDS', 'lay4_bottom_field')
        
################################################################################

def next_row_col(flow_dir, cell):
    i_next, j_next = cell
    ## Upper left cell is 0,0
    if flow_dir in [1, 2, 128]:
        i_next += 1
    elif flow_dir in [8, 16, 32]:
        i_next -= 1
    if flow_dir in [2, 4, 8]:
        j_next += 1
    elif flow_dir in [32, 64, 128]:
        j_next -= 1
    return i_next, j_next

def field_stat_func(input_path, value_field, stat='MAXIMUM'):
    value_list = []
    with arcpy.da.SearchCursor(input_path, value_field) as s_cursor:
        for row in s_cursor: value_list.append(row[0])
    if stat.upper() in ['MAXIMUM', 'MAX']:
        return max(value_list)
    elif stat.upper() in ['MINIMUM', 'MIN']:
        return min(value_list)
    else:
        return float(sum(value_list)) / count(value_list)

def add_field_func(hru_param_path, field_name, field_type='DOUBLE'):
    while not arcpy.ListFields(hru_param_path, field_name):
        logging.info('  Field: {0}'.format(field_name))
        try: arcpy.AddField_management(
            hru_param_path, field_name, field_type)
        except:
            pass
        sleep(0.5)

def transform_func(spat_ref_a, spat_ref_b):
    ## Set preferred transforms
    if ((spat_ref_a.GCS.name == 'GCS_WGS_1984' and
         spat_ref_b.GCS.name == 'GCS_North_American_1983') or
        (spat_ref_a.GCS.name == 'GCS_North_American_1983' and
         spat_ref_b.GCS.name == 'GCS_WGS_1984')):
        return 'NAD_1983_To_WGS_1984_5'
    else:
        return None
        ##return '#'

## This will check for matching spat. ref., snap point, and cellsize
## Assume raster is not valid unless it passes all checks
def valid_raster_func(raster_path, raster_name, hru_param, cs=10):
    if not arcpy.Exists(raster_path):
        return False
    logging.debug('\nReading existing {0} raster'.format(raster_name))
    raster_obj = Raster(raster_path)
    raster_sr = raster_obj.spatialReference
    raster_extent = raster_obj.extent
    raster_cs = raster_obj.meanCellWidth
    logging.debug('  {0} spat. ref.: {1}'.format(
        raster_name, raster_sr.name))
    logging.debug('  {0} GCS:        {1}'.format(
        raster_name, raster_sr.GCS.name))
    logging.debug('  {0} extent:     {1}'.format(
        raster_name, extent_string(raster_extent)))
    logging.debug('  {0} cellsize:   {1}'.format(raster_name, raster_cs))
    ref_pnt = arcpy.Point(hru_param.ref_x, hru_param.ref_y)
    if raster_sr.name <> hru_param.sr.name:
        logging.error(
            ('\nERROR: The {0} spatial reference does not match '+
             'the hru_param spatial reference').format(raster_name))
        return False
    elif not snapped(raster_extent, ref_pnt, hru_param.cs):
        logging.error(
            '\nWARNING: The {0} is not snapped to the hru_param'.format(
                raster_name))
        return False
    elif raster_cs <> cs:
        logging.error(
            '\nERROR: The {0} needs to have a {1}m cellsize'.format(
                raster_name, cs))
        return False
    elif not raster_extent.contains(hru_param.extent):
        logging.error(
            '\nERROR: The {0} extent is too small'.format(raster_name))
        return False
    else:
        return True


def zonal_stats_func(zs_dict, polygon_path, point_path, hru_param,
                     nodata_value=-999, default_value=0):
    for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
        logging.info('  {0}: {1}'.format(zs_field, zs_stat))
        logging.info('    {0}'.format(raster_path))
        ## Check inputs
        zs_stat_list = ['MEAN', 'MINIMUM', 'MAXIMUM', 'MEDIAN', 'MAJORITY', 'SUM']
        zs_field_list = arcpy.ListFields(polygon_path, zs_field)
        if zs_stat not in zs_stat_list:
            raise SystemExit()
        elif len(zs_field_list) == 0:
            logging.error(
                '\nERROR: Zonal stats field {0} doesn\'t exist'.format(zs_field))
            raise SystemExit()
        
    ## Check that the shapefiles have a spatial reference
    if arcpy.Describe(polygon_path).spatialReference.name == 'Unknown':
        logging.error(
            '\nERROR: HRU centroids  is not projected (i.e. does not have a prj file)')
        raise SystemExit()
    if arcpy.Describe(point_path).spatialReference.name == 'Unknown':
        logging.error(
            '\nERROR: HRU centroids does not appear to be projected (or does not have a prj file)'+
            '\nERROR: Try deleting the centroids (i.e. "_label.shp") and '+
            'rerunning gsflow_hru_parameters.py\n')
        raise SystemExit()

    ## Check that ORIG_FID is in point_path (HRU centroids)
    if len(arcpy.ListFields(point_path, hru_param.fid_field)) == 0:
        logging.error(
            ('\nERROR: HRU centroids does not have the field: {0}'+
             '\nERROR: Try deleting the centroids (i.e. "_label.shp") and '+
             'rerunning gsflow_hru_parameters.py\n').format(hru_param.fid_field))
        raise SystemExit()
        
    ## Create memory objects
    point_subset_path = os.path.join('in_memory', 'point_subset')
    hru_raster_path = os.path.join('in_memory', 'hru_raster')
    ##point_subset_path = os.path.join(env.scratchWorkspace, 'point_subset.shp')
    ##hru_raster_path = os.path.join(env.scratchWorkspace, 'hru_raster.img')
    ## Set environment parameters for polygon to raster conversion
    env.extent = hru_param.extent
    env.outputCoordinateSystem = polygon_path
    ##env.cellSize = hru_param.cs     

    ## Only ~65536 objects can be processed by zonal stats
    data_dict = defaultdict(dict)
    hru_param_count = int(
        arcpy.GetCount_management(polygon_path).getOutput(0))
    block_size = 65000
    for i, x in enumerate(range(0, hru_param_count, block_size)):
        logging.info('  FIDS: {0}-{1}'.format(x, x+block_size))
        ## Select a subset of the cell centroids
        logging.debug('    Selecting FID subset')
        subset_str = '"{0}" >= {1} AND "{0}" < {2}'.format(
            hru_param.fid_field, x, x+block_size)
        arcpy.Select_analysis(
            point_path, point_subset_path, subset_str)
        ## Convert points subset to raster
        logging.debug('    Converting shapefile to raster')
        arcpy.FeatureToRaster_conversion(
            point_subset_path, hru_param.fid_field,
            hru_raster_path, hru_param.cs)

        ## Zonal stats
        logging.debug('    Calculating zonal stats')
        for zs_field, (raster_path, zs_stat) in sorted(zs_dict.items()):
            zs_name = '{0}_{1}'.format(zs_field.upper(), i)
            logging.info('    {0}: {1}'.format(zs_stat.upper(), zs_name))
            ## For some reason with 10.2, ZS doesn't work with cs at HRU cs
            env.cellSize = Raster(raster_path).meanCellWidth
            ## Calculate zonal statistics
            zs_table = os.path.join('in_memory', zs_name)
            ##zs_table = os.path.join(env.scratchWorkspace, zs_name+'.dbf')
            zs_obj = ZonalStatisticsAsTable(
                hru_raster_path, 'Value', raster_path,
                zs_table, 'DATA', zs_stat.upper())

            ## Read values from points
            logging.debug('    Reading values from zs table')
            ## Fields 1 & 4 are the 'Value' (ORIG_FID) and the stat (SUM, MEAN, etc)
            fields = [
                f.name for f_i,f in enumerate(arcpy.ListFields(zs_table))
                if f_i in [1, 4]]
            for row in arcpy.da.SearchCursor(zs_table, fields):
                ## Set NoData value for cells that are entirely NoData
                if row[1] is None:
                    data_dict[int(row[0])][zs_field] = nodata_value
                else:
                    data_dict[int(row[0])][zs_field] = float(row[1])
            arcpy.Delete_management(zs_obj)
            arcpy.Delete_management(zs_table)
            del zs_table, zs_obj

        ## Cleanup
        if arcpy.Exists(point_subset_path):
            arcpy.Delete_management(point_subset_path)
        if arcpy.Exists(hru_raster_path):
            arcpy.Delete_management(hru_raster_path)

    ## Write values to polygon
    logging.info('    Writing values to polygons')
    zs_fields = sorted(zs_dict.keys())
    fields = zs_fields + [hru_param.fid_field]
    with arcpy.da.UpdateCursor(polygon_path, fields) as u_cursor:
        for row in u_cursor:
            ## Create an empty dictionary if FID does not exist
            ## Missing FIDs did not have zonal stats calculated
            row_dict = data_dict.get(int(row[-1]), None)
            for i, zs_field in enumerate(zs_fields):
                ## If stats were calculated for only some parameters,
                ##   then set missing parameter value to nodata value (-999)
                if row_dict:
                    try:
                        row[i] = row_dict[zs_field]
                    except KeyError:
                        row[i] = nodata_value
                ## Otherwise, if no stats were calculated,
                ##   reset value to 0 (shapefile default)
                else:
                    row[i] = default_value
            u_cursor.updateRow(row)
    arcpy.ClearEnvironment('extent')
    arcpy.ClearEnvironment('outputCoordinateSystem')
    arcpy.ClearEnvironment('cellSize')


def extent_string(extent_obj):
    return ' '.join(str(extent_obj).split()[:4])
    ##return ' '.join(['{0:.4f}'.format(s) for s in str(extent_obj).split()[:4]])

## This adjusts one extent to a snap point
## This is similar to the GDAL implementation
def adjust_extent_to_snap(extent_obj, snap_pnt, cs, method='EXPAND',
                          integer_flag=True):
    if method.upper() == 'ROUND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs + 0.5) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs + 0.5) * cs + snap_pnt.Y
    elif method.upper() == 'EXPAND':
        extent_xmin = math.floor(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.floor(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.ceil(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.ceil(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y
    elif method.upper() == 'SHRINK':
        extent_xmin = math.ceil(
            (extent_obj.XMin - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymin = math.ceil(
            (extent_obj.YMin - snap_pnt.Y) / cs) * cs + snap_pnt.Y
        extent_xmax = math.floor(
            (extent_obj.XMax - snap_pnt.X) / cs) * cs + snap_pnt.X
        extent_ymax = math.floor(
            (extent_obj.YMax - snap_pnt.Y) / cs) * cs + snap_pnt.Y
    if integer_flag:
        return arcpy.Extent(
            int(round(extent_xmin, 0)), int(round(extent_ymin, 0)),
            int(round(extent_xmax, 0)), int(round(extent_ymax, 0)))
    else:
        return arcpy.Extent(extent_xmin, extent_ymin, extent_xmax, extent_ymax)

def buffer_extent_func(extent_obj, extent_buffer):
    ##extent_obj = arcpy.Describe(extent_feature).extent
    extent_xmin = extent_obj.XMin-extent_buffer
    extent_ymin = extent_obj.YMin-extent_buffer
    extent_xmax = extent_obj.XMax+extent_buffer
    extent_ymax = extent_obj.YMax+extent_buffer
    return arcpy.Extent(
        extent_xmin, extent_ymin, extent_xmax, extent_ymax)

## Check if rasters are aligned to snap_raster
## Check if rasters have same cellsize as snap_raster
def snapped(extent_obj, snap_pnt, cs):
    if (((snap_pnt.X - extent_obj.XMin) % cs == 0) and
        ((snap_pnt.X - extent_obj.XMax) % cs == 0) and
        ((snap_pnt.Y - extent_obj.YMin) % cs == 0) and
        ((snap_pnt.Y - extent_obj.YMax) % cs == 0)):
        return True
    else:
        return False

def get_ini_file(workspace, ini_re, function_str='function'):
    ## Get ini file name
    ini_file_list = build_file_list(workspace, ini_re)
    ## Filter field list ini file
    ini_file_list = [item for item in ini_file_list if '_field_list.ini' not in item]
    if len(ini_file_list) == 1:
        config_filepath = ini_file_list[0]
    elif len(ini_file_list) > 1:
        ini_file_len_max = max([len(item) for item in ini_file_list])
        ## An ini file was not passed as an arguement to the script
        ## Look for ini files in the working directory
        ## If only one, use it
        ## If more, let user pick from list)
        ## If none, error out
        print '\nThere is more than one INI file present in the folder'
        print '  {0:2s}  {1}'.format('##', 'INI File')
        print '  {0:2s}  {1}'.format('==', '='*ini_file_len_max)
        for i,ini_file in enumerate(ini_file_list):
            print '  {0:2d}  {1}'.format(i, ini_file)
        config_filepath = None
        while not config_filepath:
            usr_input = raw_input('\nPlease select an INI file to use: ')
            try:
                ini_file_index = int(usr_input)
                config_filepath = ini_file_list[ini_file_index]
            except (ValueError, IndexError):
                pass
        print '  Using {0}\n'.format(config_filepath)
        del ini_file_len_max, usr_input
    else:
        print '\nERROR: No suitable ini files were found'
        print 'ERROR: Please set input file when calling {0}'.format(function_str)
        print 'ERROR: For example: test.py test.ini\n'
        raise SystemExit()
    config_filename = os.path.basename(config_filepath)
    print '{0:<20s} {1}'.format('INI File Name:', config_filename)
    return config_filepath

def get_param(param_str, param_default, config, section):
    param_type = type(param_default)
    try:
        if param_type is float:
            param_value = config.getfloat('INPUTS', param_str)
        elif param_type is int:
            param_value = config.getint('INPUTS', param_str)
        elif param_type is bool:
            param_value = config.getboolean('INPUTS', param_str)
        elif param_type is list or param_type is tuple:
            param_value = [
                ##i for i in re.split('\W+', config.get('INPUTS', param_str)) if i]
                i.strip() for i in config.get('INPUTS', param_str).split(',')
                if i.strip()]
        elif param_type is str or param_default is None:
            param_value = config.get('INPUTS', param_str)
            if param_value.upper() == 'NONE':
                param_value = None
        else:
            logging.error('ERROR: Unknown Input Type: {0}'.format(param_type))
            raise SystemExit()
    except:
        param_value = param_default
        if param_type is str and param_value.upper() == 'NONE':
            param_value = None
        logging.warning('  NOTE: {0} = {1}'.format(param_str, param_value))
    return param_value

def build_file_list(ws, test_re, test_other_re=None):
    if test_other_re is None: test_other_re=re.compile('a^')
    if os.path.isdir(ws):
        return sorted([os.path.join(ws, item) for item in os.listdir(ws)
                       if (os.path.isfile(os.path.join(ws, item)) and 
                           test_re.match(item) or test_other_re.match(item))])
    else: return []

def get_prism_data_name():
    #### Get PRISM data name
    data_name_dict = dict()
    data_name_dict[1] = 'PPT'
    data_name_dict[2] = 'TMAX'
    data_name_dict[3] = 'TMIN'
    data_name_dict[4] = 'ALL'
    print '\nPlease select which PRISM product(s) to calculate'
    print '  {0:2s}  {1}'.format('##', 'PRISM Data')
    print '  {0:2s}  {1}'.format('==', '==========')
    for i, data_name in sorted(data_name_dict.items()):
        print '  {0:2d}  {1}'.format(i, data_name)
    data_name = None
    while not data_name:
        usr_input = raw_input('\nPlease select a PRISM data product to calculate: ')
        try:
            data_name_index = int(usr_input)
            data_name = data_name_dict[data_name_index]
        except (ValueError, IndexError):
            pass
    print '  Using {0}\n'.format(data_name)
    return data_name

def project_hru_extent_func(hru_extent, hru_cs, hru_sr, 
                            target_extent, target_cs, target_sr):
    logging.debug('  Projecting extent')
    logging.debug('  HRU Extent:   {0}'.format(extent_string(hru_extent)))
    logging.debug('  HRU cellsize: {0}'.format(hru_cs))
    logging.debug('  HRU spatref:  {0}'.format(hru_sr.name))
    logging.debug('  Target snap:     {0}'.format(target_extent.lowerLeft))
    logging.debug('  Target cellsize: {0}'.format(target_cs))
    logging.debug('  Target spatref:  {0}'.format(target_sr.name))
    ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    ## Project the HRU extent to the raster spatial reference
    hru_corners = [
        [hru_extent.XMin, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMax],
        [hru_extent.XMax, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMin],
        [hru_extent.XMin, hru_extent.YMax]]
    ## Add points between corners
    hru_points = []
    for point_a, point_b in zip(hru_corners[:-1], hru_corners[1:]):
        steps = float(max(
            abs(point_b[0] - point_a[0]),
            abs(point_b[1] - point_a[1]))) / hru_cs
        for x,y in zip(np.linspace(point_a[0], point_b[0], steps + 1),
                       np.linspace(point_a[1], point_b[1], steps + 1)):
            hru_points.append(arcpy.Point(x,y))
    ## Project all points to output spatial reference and get projected extent
    transform = transform_func(hru_sr, target_sr)
    if transform:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(
                target_sr, transform).extent
    else:
        projected_extent = arcpy.Polygon(
            arcpy.Array(hru_points), hru_sr).projectAs(target_sr).extent
    logging.debug('  Projected Extent: {0}'.format(
        extent_string(projected_extent)))
    ## Adjust extent to match snap
    projected_extent = adjust_extent_to_snap(
        projected_extent, target_extent.lowerLeft, target_cs, 'EXPAND', False)
    logging.debug('  Snapped Extent: {0}'.format(
        extent_string(projected_extent)))
    ## Buffer extent 4 input cells
    ##projected_extent = buffer_extent_func(projected_extent, 4 * target_cs)
    projected_extent = buffer_extent_func(
        projected_extent, 4 * max(target_cs, hru_cs))
    logging.debug('  Buffered Extent: {0}'.format(
        extent_string(projected_extent)))
    return projected_extent

def project_raster_func(input_raster, output_raster, output_sr,
                        proj_method, input_cs, transform_str,
                        reg_point, input_sr, hru_param):
    ## Input raster can be a raster object or a raster path
    ##print isinstance(input_raster, Raster), isinstance(input_raster, str)
    try: input_extent = Raster(input_raster).extent
    except: input_extent = input_raster.extent
    ## DEADBEEF - Arc10.2 ProjectRaster does not honor extent
    ## Clip the input raster with the projected HRU extent first
    ## Project extent from "output" to "input" to get clipping extent
    proj_extent = project_hru_extent_func(
        hru_param.extent, hru_param.cs, output_sr,
        input_extent, input_cs, input_sr)
    ##clip_path = output_raster.replace('.img', '_clip.img')
    clip_path = os.path.join('in_memory', 'clip_raster')
    env.extent = proj_extent
    arcpy.Clip_management(
        input_raster, ' '.join(str(proj_extent).split()[:4]), clip_path)
    arcpy.ClearEnvironment('extent')
    ## Then project the clipped raster
    arcpy.ProjectRaster_management(
        clip_path, output_raster, output_sr, proj_method.upper(), input_cs,
        transform_str, reg_point, input_sr)
    ## Cleanup
    arcpy.Delete_management(clip_path)
   
def cell_area_func(hru_param_path, area_field):
    arcpy.CalculateField_management(
        hru_param_path, area_field, '!SHAPE.AREA@acres!', 'PYTHON')

## These two functions only set values that are in zone
## They don't reset values that are out of zone
def zone_by_area_func(
    zone_path, zone_field, zone_value, hru_param_path, hru_param,  
    hru_area_field='HRU_AREA', zone_area_field=None, area_pct=50):
    zone_value_field = 'ZONE_VALUE'
    int_area_field = 'INT_AREA'
    int_pct_field = 'INT_PCT'
    ## Need to set zone value into a field before intersect
    ## If zone_value is FID, add 1 so that only non-lake cells are 0
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    arcpy.AddField_management(zone_path, int_area_field, 'DOUBLE')
    arcpy.AddField_management(zone_path, int_pct_field, 'DOUBLE')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}! + 1'.format(zone_value), 'PYTHON')
    ## If zone value is an INT, save it into a field first
    elif type(zone_value) is int:
        ##zone_value = int(zone_value)
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
    ## Use zone_value field directly
    else:
        ##zone_value_field = zone_value
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}!'.format(zone_value), 'PYTHON')
    ## Calculate area of HRU cell if necessary
    ##if not arcpy.ListFields(zone_path, area_field):
    ##    arcpy.AddField_management(zone_path, area_field, 'DOUBLE')
    ##    cell_area_func(zone_path, area_field)
    ## Intersect the zone layer with the fishnet
    ##zone_int_path = os.path.join('in_memory', 'hru_lakes')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_param_path, zone_path), zone_int_path)
    ## Calculate using cell_area_func to force units to match
    cell_area_func(zone_int_path, int_area_field)
    ## Read in FID of selected cells
    hru_cell_dict = dict()
    fields = [
        hru_param.fid_field, hru_area_field,
        int_area_field, zone_value_field]
    with arcpy.da.SearchCursor(zone_int_path, fields) as s_cursor:
        for row in s_cursor:
            if (100 * float(row[2]) / float(row[1])) >= area_pct:
                hru_cell_dict[int(row[0])] = [float(row[3]), float(row[2])]
    ## Set value of selected HRU cells
    fields = [hru_param.fid_field, zone_field]
    if zone_area_field: fields.append(zone_area_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields) as u_cursor:
        for row in u_cursor:
            ## Remove items to speed up subsequent searches
            try:
                if len(fields) == 3:
                    row[1], row[2] = hru_cell_dict.pop(int(row[0]))
                elif len(fields) == 2:
                    row[1] = hru_cell_dict.pop(int(row[0]))
                u_cursor.updateRow(row)
            except KeyError: pass

def zone_by_centroid_func(
    zone_path, zone_field, zone_value, 
    hru_param_path, hru_point_path, hru_param):
    ## Need to set zone value into a field before intersect
    ## If zone_value is FID, add 1 so that only zone cells are 0
    zone_value_field = 'ZONE_VALUE'
    arcpy.AddField_management(zone_path, zone_value_field, 'LONG')
    if zone_value == arcpy.Describe(zone_path).OIDFieldName:
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}! + 1'.format(zone_value), 'PYTHON')
    ## If zone value is an INT, save it into a field first
    elif type(zone_value) is int:
        ##zone_value = int(zone_value)
        arcpy.CalculateField_management(
            zone_path, zone_value_field, zone_value, 'PYTHON')
    ## Use zone_value field directly
    else:
        ##zone_value_field = zone_value
        arcpy.CalculateField_management(
            zone_path, zone_value_field,
            '!{0}!'.format(zone_value), 'PYTHON')
    ## Intersect the zone layer with the fishnet
    ##zone_int_path = os.path.join('in_memory', 'hru_ppt_zones')
    zone_int_path = zone_path.replace('.shp', '_intersect.shp')
    arcpy.Intersect_analysis(
        (hru_point_path, zone_path), zone_int_path)
    zone_int_layer = 'zone_int_layer'
    arcpy.MakeFeatureLayer_management(zone_int_path, zone_int_layer)
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'CLEAR_SELECTION')
    arcpy.SelectLayerByAttribute_management(zone_int_layer, 'SWITCH_SELECTION')
    ## Read in FID of selected cells
    hru_point_dict = dict()
    fields = (hru_param.fid_field, zone_value_field)
    with arcpy.da.SearchCursor(zone_int_layer, fields) as s_cursor:
        for row in s_cursor:
            hru_point_dict[int(row[0])] = row[1]
    ## Set value of selected HRU cells
    fields = (hru_param.fid_field, zone_field)
    with arcpy.da.UpdateCursor(hru_param_path, fields) as u_cursor:
        for row in u_cursor:
            ## Remove items to speed up subsequent searches
            try:
                row[1] = hru_point_dict.pop(int(row[0]))
                u_cursor.updateRow(row)
            except KeyError:
                pass
    ## Cleanup
    arcpy.Delete_management(zone_int_layer)

def jensen_haise_func(
    hru_param_path, jh_coef_field, hru_elev_field,
    jh_tmin_field, jh_tmax_field):
    jh_cb = (
        'def ea(temp_c):\n'+
        '    return 6.1078 * math.exp((17.269 * temp_c) / (temp_c + 237.3))\n'+
        'def jensen_haise(elev, t_low, t_high):\n'+
        '    return 27.5 - 0.25 * (ea(t_high) - ea(t_low)) - (elev / 1000)\n')
    arcpy.CalculateField_management(
        hru_param_path, jh_coef_field,
        'jensen_haise(!{0}!, !{1}!, !{2}!)'.format(
            hru_elev_field, jh_tmin_field, jh_tmax_field),
        'PYTHON', jh_cb)    

################################################################################

#### Remap aspect
##logging.info('\nRemapping Aspect to HRU_ASPECT')
##arcpy.CalculateField_management(
##    polygon_path, hru_aspect_field,
##    'Reclass(!{0}!)'.format(dem_aspect_field),
##    'PYTHON', remap_code_block(aspect_remap_path))
####arcpy.DeleteField_management(polygon_path, dem_aspect_field)
def remap_code_block(remap_path):
    ## Read remap file into memory
    with open(remap_path) as remap_f: lines = remap_f.readlines()
    remap_f.close()
    remap_cb = ''
    for l in lines:
        ## Skip comment lines
        if '#' in l: continue
        ## Remove remap description
        l = l.strip().split('/*')[0]
        ## Split line on spaces and semi-colon
        l_split = [item.strip() for item in re.split('[ :]+',l)]
        ## Remap as a range if a min, max and value are all present
        if len(l_split) == 3: range_remap_flag = True
        ## Otherwise remap directly
        elif len(l_split) == 2: range_remap_flag = False
        ## Skip lines that don't match format
        else: continue
        ## Write remap code block
        if not range_remap_flag:
            if not remap_cb:
                remap_cb = '    if value == {0}: return {1}\n'.format(*l_split)
            else:
                remap_cb += '    elif value == {0}: return {1}\n'.format(*l_split)
        else:
            if not remap_cb:
                remap_cb = ('    if (value >= {0} and value <= {1}): '+
                             'return {2}\n').format(*l_split)
            else:
                remap_cb += ('    elif (value > {0} and value <= {1}): '+
                             'return {2}\n').format(*l_split)
    remap_cb = 'def Reclass(value):\n' + remap_cb
    return remap_cb

##def reclass_ascii_float_func(raster_path, remap_path):
##    ## Read remap file into memory
##    with open(remap_path) as remap_f: lines = remap_f.readlines()
##    remap_f.close()
##    first_line = True
##    raster_obj = Raster(raster_path)
##    for l in lines:
##        ## Skip comment lines
##        if '#' in l: continue
##        ## Remove remap description
##        l = l.split('/*')[0]
##        ## Split line on spaces and semi-colon
##        l_split = map(float, [item for item in re.split('[ :]+', l) if item])
##        ## Remap as a range if a min, max and value are all present
##        if len(l_split) == 3: range_remap_flag = True
##        ## Otherwise remap directly
##        elif len(l_split) == 2: range_remap_flag = False
##        ## Skip lines that don't match format
##        else: continue
##        ## Write remap code block
##        if not range_remap_flag:
##            raster_obj = Con(raster_obj == l_split[0], l_split[1], raster_obj)
##        elif first_line:
##            raster_obj = Con(
##                ((raster_obj >= l_split[0]) & (raster_obj <= l_split[1])),
##                l_split[2], raster_obj)
##        else:
##            raster_obj = Con(
##                ((raster_obj > l_split[0]) & (raster_obj <= l_split[1])),
##                l_split[2], raster_obj)
##        first_line = False
##    return raster_obj

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
