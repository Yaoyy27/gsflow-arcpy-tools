#--------------------------------
# Name:         crt_fill_parameters.py
# Purpose:      GSFLOW CRT fill parameters
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
import shutil
import subprocess
import sys
from time import clock

import arcpy
from arcpy import env
from arcpy.sa import *

##import numpy as np

from support_functions import *

################################################################################

def gsflow_crt_fill_parameters(workspace, config_path=None):
    """Calculate GSFLOW CRT Fill Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW CRT Fill Parameters')

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

        ## Parameters
        exit_seg = 0

        ## CRT Parameters
        calc_fill_work_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_fill_work_flag')
        crt_hruflg = 0
        crt_flowflg = 1
        crt_iprn = 1
        crt_dpit = 0.01
        crt_outitmax = 10000

        ## CRT Fill Parameters
        fill_ws_name = 'fill_work'
        fill_strmflg = 0
        fill_visflg = 0
        fill_ifill = 1

        ## CRT Executable
        crt_exe_path = inputs_cfg.get('INPUTS', 'crt_exe_path')
        ##crt_exe_name = 'CRT_1.1.1.exe'
        output_name = 'outputstat.txt'

        ## Check input paths
        if not arcpy.Exists(hru.polygon_path):
            logging.error(
                '\nERROR: Fishnet ({0}) does not exist\n'.format(
                    hru.polygon_path))
            raise SystemExit()
        flow_temp_ws = os.path.join(hru.param_ws, 'flow_rasters')
        if not os.path.isdir(flow_temp_ws): 
            logging.error(
                ('\nERROR: Flow_rasters folder does not exist'+
                 '\nERROR:   {0}'+
                 '\nERROR: Try re-running gsflow_dem_2_stream.py\n').format(
                     flow_temp_ws))
            raise SystemExit()
        ## Check that input fields exist and have data
        for f in [hru.irunbound_field, hru.iseg_field, hru.flow_dir_field]:
            if not arcpy.ListFields(hru.polygon_path, f): 
                logging.error(
                    ('\nERROR: Input field {0} is not present in fishnet'+
                     '\nERROR: Try re-running gsflow_dem_2_stream.py\n').format(
                         f))
                raise SystemExit()
            elif field_stat_func(hru.polygon_path, f, 'MAXIMUM') == 0:
                logging.error(
                    ('\nERROR: Input field {0} contains only 0'+
                     '\nERROR: Try re-running gsflow_dem_2_stream.py\n').format(
                         f))
                raise SystemExit()


        ## Build output folder if necessary
        fill_ws = os.path.join(hru.param_ws, fill_ws_name)
        if not os.path.isdir(fill_ws):
            os.mkdir(fill_ws)

        ## Copy CRT executable if necessary
        crt_exe_name = os.path.basename(crt_exe_path)
        if not os.path.isfile(os.path.join(fill_ws, crt_exe_name)):
            shutil.copy(crt_exe_path, fill_ws)
        if not os.path.isfile(os.path.join(fill_ws, crt_exe_name)):
            logging.error(
                '\nERROR: CRT executable ({0}) does not exist\n'.format(
                    os.path.join(fill_ws, crt_exe_name)))
            raise SystemExit()

        ## Fill files
        fill_hru_casc_path     = os.path.join(fill_ws, 'HRU_CASC.DAT')
        fill_outflow_hru_path  = os.path.join(fill_ws, 'OUTFLOW_HRU.DAT')
        fill_land_elev_path    = os.path.join(fill_ws, 'LAND_ELEV.DAT')
        fill_xy_path           = os.path.join(fill_ws, 'XY.DAT')

        ## Output names
        dem_adj_raster_name    = 'dem_adj'
        hru_type_raster_name   = 'hru_type'
        lakes_raster_name      = 'lakes'
        streams_raster_name    = 'streams'
        iseg_raster_name       = 'iseg'
        irunbound_raster_name  = 'irunbound'

        ## Output raster paths
        dem_adj_raster   = os.path.join(fill_ws, dem_adj_raster_name+'.img')
        hru_type_raster  = os.path.join(fill_ws, hru_type_raster_name+'.img')
        ## Output ascii paths
        a_fmt = '{0}_ascii.txt'
        dem_adj_ascii   = os.path.join(fill_ws, a_fmt.format(dem_adj_raster_name))
        hru_type_ascii  = os.path.join(fill_ws, a_fmt.format(hru_type_raster_name))


        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        ##env.pyramid = 'PYRAMIDS -1'
        env.pyramid = 'PYRAMIDS 0'
        env.workspace = fill_ws
        env.scratchWorkspace = hru.scratch_ws

        ## Add fields if necessary
        logging.info('\nAdding fields if necessary')
        add_field_func(hru.polygon_path, hru.krch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.irch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.jrch_field, 'LONG')
        add_field_func(hru.polygon_path, hru.iseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.reach_field, 'LONG')
        ##add_field_func(hru.polygon_path, hru.rchlen_field, 'LONG')
        add_field_func(hru.polygon_path, hru.maxreach_field, 'LONG')
        add_field_func(hru.polygon_path, hru.outseg_field, 'LONG')
        add_field_func(hru.polygon_path, hru.irunbound_field, 'LONG')
        add_field_func(hru.polygon_path, hru.crt_dem_field, 'DOUBLE')
        add_field_func(hru.polygon_path, hru.crt_fill_field, 'DOUBLE')

        ## Calculate KRCH, IRCH, JRCH for stream segments
        logging.info("\nKRCH, IRCH, & JRCH for streams")
        fields = [
            hru.type_field, hru.iseg_field, hru.row_field, hru.col_field,
            hru.krch_field, hru.irch_field, hru.jrch_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                if (int(row[0]) == 1 and int(row[1]) > 0):
                    row[4], row[5], row[6] = 1, int(row[2]), int(row[3])
                else:
                    row[4], row[5], row[6] = 0, 0, 0
                update_c.updateRow(row)

        ## Get list of segments and downstream cell for each stream/lake cell
        ## Downstream is calulated from flow direction
        ## Use IRUNBOUND instead of ISEG, since ISEG will be zeroed for lakes
        logging.info("Cell out-flow dictionary")
        cell_dict = dict()
        fields = [
            hru.type_field, hru.krch_field, hru.lake_id_field, hru.iseg_field,
            hru.irunbound_field, hru.dem_adj_field, hru.flow_dir_field,
            hru.col_field, hru.row_field, hru.id_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            ## Skip inactive cells
            if int(row[0]) == 0:
                continue
            ## Skip non-lake and non-stream cells
            if (int(row[1]) == 0 and int(row[2]) == 0):
                continue
            ## Read in parameters
            cell = (int(row[7]), int(row[8]))
            ## next_row_col(FLOW_DIR, CELL)
            ## HRU_ID, ISEG,  NEXT_CELL, DEM_ADJ, X, X, X
            cell_dict[cell] = [
                int(row[9]), int(row[4]), next_row_col(int(row[6]), cell),
                float(row[5]), 0, 0, 0]
            del cell
        ## Build list of unique segments
        iseg_list = sorted(list(set([v[1] for v in cell_dict.values()])))

        ## Calculate IREACH and OUTSEG
        logging.info("Calculate {0} and {1}".format(hru.reach_field, hru.outseg_field))
        outseg_dict = dict()
        for iseg in iseg_list:
            logging.debug("    Segment: {0}".format(iseg))
            ## Subset of cell_dict for current iseg
            iseg_dict = dict(
                [(k,v) for k,v in cell_dict.items() if v[1] == iseg])
            ## List of all cells in current iseg
            iseg_cells = iseg_dict.keys()
            ## List of out_cells for all cells in current iseg
            out_cells = [value[2] for value in iseg_dict.values()]
            ## Every iseg will (should?) have one out_cell
            out_cell = list(set(out_cells)-set(iseg_cells))[0]
            ## If not output cell, assume edge of domain
            try:
                outseg = cell_dict[out_cell][1]
            except KeyError:
                outseg = exit_seg
            ## Track sub-basin outseg
            outseg_dict[iseg] = outseg
            if iseg > 0: 
                ## Calculate reach number for each cell
                reach_dict = dict()
                start_cell = list(set(iseg_cells)-set(out_cells))[0]
                for i in xrange(len(out_cells)):
                    ##logging.debug("    Reach: {0}  Cell: {1}".format(i+1, start_cell))
                    reach_dict[start_cell] = i+1
                    start_cell = iseg_dict[start_cell][2]
                ## For each cell in iseg, save outseg, reach, & maxreach
                for iseg_cell in iseg_cells:
                    cell_dict[iseg_cell][4:] = [
                        outseg, reach_dict[iseg_cell], len(iseg_cells)]
                del reach_dict, start_cell
            else:
                ## For each lake segment cell, only save outseg
                ## All lake cells are routed directly to the outseg
                for iseg_cell in iseg_cells:
                    cell_dict[iseg_cell][4:] = [outseg, 0, 0]
            del iseg_dict, iseg_cells, iseg
            del out_cells, out_cell, outseg

        ## Saving ireach and outseg
        logging.info("Save {0} and {1}".format(hru.reach_field, hru.outseg_field))
        fields = [
            hru.type_field, hru.iseg_field, hru.col_field, hru.row_field,
            hru.outseg_field, hru.reach_field, hru.maxreach_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                ##if (int(row[0]) > 0 and int(row[1]) > 0):
                ###DEADBEEF - I'm not sure why only iseg > 0 in above line
                ## DEADBEEF - This should set outseg for streams and lakes
                if (int(row[0]) > 0 and int(row[1]) <> 0):
                    row[4:] = cell_dict[(int(row[2]), int(row[3]))][4:]
                else:
                    row[4:] = [0, 0, 0]
                update_c.updateRow(row)
      
        ## Set all lake iseg to 0
        logging.info("Lake {0}".format(hru.iseg_field))
        update_rows = arcpy.UpdateCursor(hru.polygon_path)
        for row in update_rows:
            if int(row.getValue(hru.type_field)) <> 2:
                continue
            iseg = int(row.getValue(hru.iseg_field))
            if iseg < 0:
                row.setValue(hru.iseg_field, 0)
            update_rows.updateRow(row)
            del row, iseg
        del update_rows

        ## Set environment parameters
        env.extent = hru.extent
        env.cellsize = hru.cs
        env.outputCoordinateSystem = hru.sr

        #### Build rasters
        ##logging.info("\nOutput model grid rasters")
        ##arcpy.PolygonToRaster_conversion(
        ##    hru.polygon_path, hru.type_field, hru_type_raster,
        ##    "CELL_CENTER", "", hru.cs)
        ##arcpy.PolygonToRaster_conversion(
        ##    hru.polygon_path, hru.dem_adj_field, dem_adj_raster,
        ##    "CELL_CENTER", "", hru.cs)
        ##
        #### Build rasters
        ##logging.info("Output model grid ascii")
        ##arcpy.RasterToASCII_conversion(hru_type_raster, hru_type_ascii)
        ##arcpy.RasterToASCII_conversion(dem_adj_raster, dem_adj_ascii)

        ## Input parameters files for Cascade Routing Tool (CRT)
        logging.info("\nOutput CRT fill files")
            
        ## Generate OUTFLOW_HRU.DAT for CRT
        logging.info("  {0}".format(
            os.path.basename(fill_outflow_hru_path)))
        outflow_hru_list = []
        fields = [
            hru.type_field, hru.iseg_field, hru.outseg_field, hru.reach_field,
            hru.maxreach_field, hru.col_field, hru.row_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            if int(row[0]) <> 1 or int(row[1]) == 0:
                continue
            if int(row[2]) == 0 and int(row[3]) == int(row[4]): 
                outflow_hru_list.append([int(row[6]), int(row[5])])
        if outflow_hru_list:
            with open(fill_outflow_hru_path, 'w+') as f:
                f.write('{0}    NUMOUTFLOWHRU\n'.format(
                    len(outflow_hru_list)))
                for i, outflow_hru in enumerate(outflow_hru_list):
                    f.write('{0} {1} {2}   OUTFLOW_ID ROW COL\n'.format(
                        i+1, outflow_hru[0], outflow_hru[1]))
            f.close()
        del outflow_hru_list

        ## Generate HRU_CASC.DAT for CRT from hru_polygon
        logging.info("  {0}".format(os.path.basename(fill_hru_casc_path)))
        hru_type_dict = defaultdict(dict)
        for row in sorted(arcpy.da.SearchCursor(
            hru.polygon_path,
            [hru.row_field, hru.col_field, hru.type_field, hru.dem_adj_field])):
            ## Calculate CRT fill for all non-lake and non-ocean (elev > 0) cells
            ##if row[3] > 0 and row[2] == 0:
            ##    hru_type_dict[int(row[0])][int(row[1])] = 1
            ##else: hru_type_dict[int(row[0])][int(row[1])] = row[2]
            ## Calculate CRT fill for all active cells
            hru_type_dict[int(row[0])][int(row[1])] = row[2]
        hru_casc_header = (
            '{0} {1} {2} {3} {4} {5} {6} {7}     '+
            'HRUFLG STRMFLG FLOWFLG VISFLG '+
            'IPRN IFILL DPIT OUTITMAX\n').format(
                crt_hruflg, fill_strmflg, crt_flowflg, fill_visflg,
                crt_iprn, fill_ifill, crt_dpit, crt_outitmax)
        with open(fill_hru_casc_path, 'w+') as f:
            f.write(hru_casc_header)
            for row, col_data in sorted(hru_type_dict.items()):
                f.write(' '.join([str(t) for c,t in sorted(col_data.items())])+'\n')
        f.close()
        del hru_casc_header, hru_type_dict
        #### Generate HRU_CASC.DATA for CRT from raster/ascii
        ##with open(hru_type_ascii, 'r') as f: ascii_data = f.readlines()
        ##f.close()
        ##hru_casc_header = (
        ##    '{0} {1} {2} {3} {4} {5} {6} {7}     '+
        ##    'HRUFLG STRMFLG FLOWFLG VISFLG '+
        ##    'IPRN IFILL DPIT OUTITMAX\n').format(
        ##        crt_hruflg, fill_strmflg, crt_flowflg, fill_visflg,
        ##        crt_iprn, fill_ifill, crt_dpit, crt_outitmax)
        ##with open(fill_hru_casc_path, 'w+') as f:
        ##    f.write(hru_casc_header)
        ##    for ascii_line in ascii_data[6:]: f.write(ascii_line)
        ##f.close()
        ##del hru_casc_header, ascii_data
        
        ## Generate LAND_ELEV.DAT for CRT from hru_polygon
        logging.info("  {0}".format(os.path.basename(fill_land_elev_path)))
        dem_adj_dict = defaultdict(dict)
        for row in sorted(arcpy.da.SearchCursor(
            hru.polygon_path, [hru.row_field, hru.col_field, hru.dem_adj_field])):
            dem_adj_dict[int(row[0])][int(row[1])] = row[2]
        with open(fill_land_elev_path, 'w+') as f:
            f.write('{0} {1}       NROW NCOL\n'.format(
                len(dem_adj_dict.keys()), len(dem_adj_dict[1].keys())))
            for row, col_data in sorted(dem_adj_dict.items()):
                f.write(' '.join(
                    ['{0:10.6f}'.format(t) for c,t in sorted(col_data.items())])+'\n')
        f.close()
        del dem_adj_dict
        #### Generate LAND_ELEV.DAT for CRT from raster/ascii
        ##logging.info("  {0}".format(os.path.basename(fill_land_elev_path)))
        ##with open(dem_adj_ascii, 'r') as f: ascii_data = f.readlines()
        ##f.close()
        ##with open(fill_land_elev_path, 'w+') as f:
        ##    f.write('{0} {1}       NROW NCOL\n'.format(
        ##        ascii_data[1].split()[1], ascii_data[0].split()[1]))
        ##    for ascii_line in ascii_data[6:]: f.write(ascii_line)
        ##f.close()
        ##del ascii_data

        ## Generate XY.DAT for CRT
        logging.info("  {0}".format(os.path.basename(fill_xy_path)))
        xy_list = [
            map(int, row)
            for row in sorted(arcpy.da.SearchCursor(
                hru.polygon_path, [hru.id_field, hru.x_field, hru.y_field]))]
        with open(fill_xy_path, 'w+') as f:
            for line in sorted(xy_list):
                f.write(' '.join(map(str, line))+'\n')
        f.close()
        del xy_list

        ## Run CRT
        logging.info('\nRunning CRT')
        os.chdir(fill_ws)
        subprocess.check_call(crt_exe_name)
        os.chdir(workspace)

        ## Read in outputstat.txt and get filled DEM
        logging.info("\nReading CRT {0}".format(output_name))
        output_path = os.path.join(fill_ws, output_name)
        with open(output_path, 'r') as f:
            output_data = [l.strip() for l in f.readlines()]
        f.close()

        ## Determine where filled data is in file
        crt_dem_i = output_data.index(
            'CRT FILLED LAND SURFACE MODEL USED TO GENERATE CASCADES')
        crt_fill_i = output_data.index(
            'DIFFERENCES BETWEEN FILLED AND UNFILLED LAND SURFACE MODELS')
        crt_type_i = output_data.index(
            'FINAL HRU CASCADE TYPE ARRAY USED TO COMPUTE CASCADES')
        logging.info('  Break indices: {0}, {1}, {2}'.format(
            crt_dem_i, crt_fill_i, crt_type_i))
        crt_dem_data = [r.split() for r in output_data[crt_dem_i+1:crt_fill_i-1]]
        crt_fill_data = [r.split() for r in output_data[crt_fill_i+1:crt_type_i-1]]
        logging.info('  ROWS/COLS: {0}/{1}'.format(
            len(crt_dem_data), len(crt_dem_data[0])))
        logging.info('  ROWS/COLS: {0}/{1}'.format(
            len(crt_fill_data), len(crt_fill_data[0])))

        ## Build dictionaries of the CRT data
        crt_dem_dict = defaultdict(dict)
        crt_fill_dict = defaultdict(dict)
        for i,r in enumerate(crt_dem_data):
            crt_dem_dict[i+1] = dict(
                [(j+1,c) for j,c in enumerate(crt_dem_data[i])])
        for i,r in enumerate(crt_fill_data):
            crt_fill_dict[i+1] = dict(
                [(j+1,c) for j,c in enumerate(crt_fill_data[i])])

        ## Write CRT values to hru_polygon
        logging.info("Writing CRT data to fishnet")
        fields = [hru.row_field, hru.col_field, hru.crt_dem_field, hru.crt_fill_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                row[2] = crt_dem_dict[int(row[0])][int(row[1])]
                row[3] = crt_fill_dict[int(row[0])][int(row[1])]
                update_c.updateRow(row)

    except:
        logging.exception('Unhandled Exception Error\n\n')
        raw_input('ENTER')

    finally:
        try: arcpy.CheckInExtension('Spatial')
        except: pass
        ##arcpy.ResetEnvironments()

################################################################################

def cell_distance(cell_a, cell_b, cs):
    ai, aj = cell_a
    bi, bj = cell_b
    return math.sqrt((ai - bi) ** 2 + (aj - bj) ** 2) * cs

##def calc_stream_width(flow_acc):
##    return -2E-6 * flow_acc ** 2 + 0.0092 * flow_acc + 1

################################################################################

if __name__ == '__main__':
    workspace = os.getcwd()
    log_ws = os.path.join(workspace, 'logs')
    if not os.path.isdir(log_ws):
        os.mkdir(log_ws)
    log_file_name = 'gsflow_streams_log.txt'

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
        ini_path = get_ini_file(workspace, ini_re, 'gsflow_dem_parameters')
    del ini_re

    ## Run Information
    logging.info('\n{0}'.format('#'*80))
    log_f = '{0:<20s} {1}'
    logging.info(log_f.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', workspace))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    logging.info(log_f.format('INI File:', os.path.basename(ini_path)))

    ## Calculate CRT Fill Parameters
    gsflow_crt_fill_parameters(workspace, ini_path)
