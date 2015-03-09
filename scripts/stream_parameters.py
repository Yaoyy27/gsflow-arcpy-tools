#--------------------------------
# Name:         stream_parameters.py
# Purpose:      GSFLOW stream parameters
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
from time import clock, sleep

import arcpy
from arcpy import env
from arcpy.sa import *

##import numpy as np

from support_functions import *

################################################################################

def gsflow_stream_parameters(workspace, config_path=None):
    """Calculate GSFLOW Stream Parameters

    keyword arguments:
    workspace -- the workspace (path) of the landsat scene folder
    config_path -- the config file (path)

    """

    try:
        logging.info('\nGSFLOW Stream Parameters')

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

        ## CRT Parameters
        calc_cascade_work_flag = inputs_cfg.getboolean(
            'INPUTS', 'calc_cascade_work_flag')
        crt_hruflg = 0
        crt_flowflg = 1
        crt_iprn = 1
        crt_dpit = 0.01
        crt_outitmax = 10000

        ## CRT Fill Parameters
        fill_strmflg = 0
        fill_visflg = 0
        fill_ifill = 1

        ## CRT Streams paramters
        crt_ws = os.path.join(hru.param_ws, 'cascade_work')
        crt_strmflg = 1
        crt_visflg = 1
        crt_ifill = 0

        ## CRT Executable
        crt_exe_path = inputs_cfg.get('INPUTS', 'crt_exe_path')
        ##crt_exe_name = 'CRT_1.1.1.exe'

        ## Override ascii and rasters flags to generate CRT inputs
        output_ascii_flag = True
        output_rasters_flag = True

        ## Parameters
        ##lake_seg_offset = fields_cfg.getint('INPUTS', 'lake_seg_offset')
        exit_seg = 0

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
        streams_path = os.path.join(flow_temp_ws, 'streams.shp')
        if not os.path.isfile(streams_path): 
            logging.error(
                ('\nERROR: Stream shapefiles does not exist'+
                 '\nERROR:   {0}'+
                 '\nERROR: Try re-running gsflow_dem_2_stream.py\n').format(
                     streams_path))
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
        stream_temp_ws = os.path.join(hru.param_ws, 'stream_rasters')
        if not os.path.isdir(stream_temp_ws):
            os.mkdir(stream_temp_ws)
        if not os.path.isdir(crt_ws):
            os.mkdir(crt_ws)

        ## Copy CRT executable if necessary
        crt_exe_name = os.path.basename(crt_exe_path)
        if not os.path.isfile(os.path.join(crt_ws, crt_exe_name)):
            shutil.copy(crt_exe_path, crt_ws)
        if not os.path.isfile(os.path.join(crt_ws, crt_exe_name)):
            logging.error(
                '\nERROR: CRT executable ({0}) does not exist\n'.format(
                    os.path.join(crt_ws, crt_exe_name)))
            raise SystemExit()

        ## Cascades files
        crt_hru_casc_path     = os.path.join(crt_ws, 'HRU_CASC.DAT')
        crt_outflow_hru_path  = os.path.join(crt_ws, 'OUTFLOW_HRU.DAT')
        crt_land_elev_path    = os.path.join(crt_ws, 'LAND_ELEV.DAT')
        crt_stream_cells_path = os.path.join(crt_ws, 'STREAM_CELLS.DAT')
        crt_xy_path           = os.path.join(crt_ws, 'XY.DAT')

        ## Output names
        dem_adj_raster_name    = 'dem_adj'
        hru_type_raster_name   = 'hru_type'
        lakes_raster_name      = 'lakes'
        streams_raster_name    = 'streams'
        iseg_raster_name       = 'iseg'
        irunbound_raster_name  = 'irunbound'
        subbasin_raster_name   = 'sub_basins'
        segbasin_raster_name   = 'seg_basins'

        ## Output raster paths
        dem_adj_raster   = os.path.join(stream_temp_ws, dem_adj_raster_name+'.img')
        hru_type_raster  = os.path.join(stream_temp_ws, hru_type_raster_name+'.img')
        iseg_raster      = os.path.join(stream_temp_ws, iseg_raster_name+'.img')
        irunbound_raster = os.path.join(stream_temp_ws, irunbound_raster_name+'.img')
        subbasin_raster  = os.path.join(stream_temp_ws, subbasin_raster_name+'.img')
        segbasin_raster  = os.path.join(stream_temp_ws, segbasin_raster_name+'.img')
        ## Output ascii paths
        a_fmt = '{0}_ascii.txt'
        dem_adj_ascii   = os.path.join(stream_temp_ws, a_fmt.format(dem_adj_raster_name))
        hru_type_ascii  = os.path.join(stream_temp_ws, a_fmt.format(hru_type_raster_name))
        iseg_ascii      = os.path.join(stream_temp_ws, a_fmt.format(iseg_raster_name))
        irunbound_ascii = os.path.join(stream_temp_ws, a_fmt.format(irunbound_raster_name))
        subbasin_ascii  = os.path.join(stream_temp_ws, a_fmt.format(subbasin_raster_name))
        segbasin_ascii  = os.path.join(stream_temp_ws, a_fmt.format(segbasin_raster_name))

        ## Layers
        hru_polygon_lyr = 'hru_polygon_lyr'

        ## Set ArcGIS environment variables
        arcpy.CheckOutExtension('Spatial')
        env.overwriteOutput = True
        ##env.pyramid = 'PYRAMIDS -1'
        env.pyramid = 'PYRAMIDS 0'
        env.workspace = stream_temp_ws
        env.scratchWorkspace = hru.scratch_ws

        ## Add fields if necessary
        logging.info('\nAdding fields if necessary')
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

        ## Check watershed and stream values from gsflow_dem_2_stream.py
        ## Lakes must be negative of LAKE_ID, not LAKE_ID + OFFSET
        ##logging.info("\nChecking watershed & stream values")
        ##fields = [hru.irunbound_field, hru.iseg_field]
        ##with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
        ##    for row in update_c:
        ##        irunbound, iseg = map(int, row)
        ##        if irunbound > lake_seg_offset:
        ##            row[0] = lake_seg_offset - irunbound
        ##        if iseg > lake_seg_offset:
        ##            row[1] = lake_seg_offset - iseg
        ##        del irunbound, iseg

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

        ## Get stream length for each cell
        logging.info("Stream length")
        arcpy.MakeFeatureLayer_management(hru.polygon_path, hru_polygon_lyr)
        arcpy.SelectLayerByAttribute_management(
            hru_polygon_lyr, "NEW_SELECTION",
            ' \"{0}\" = 1 And "{1}" <> 0'.format(hru.type_field, hru.iseg_field))
        length_path = os.path.join('in_memory', 'length')
        arcpy.Intersect_analysis(
            [hru_polygon_lyr, streams_path],
            length_path, "ALL", "", "LINE")
        arcpy.Delete_management(hru_polygon_lyr)
        length_field = 'LENGTH'
        arcpy.AddField_management(length_path, length_field, 'LONG')
        arcpy.CalculateField_management(
            length_path, length_field, '!shape.length@meters!', "PYTHON")
        length_dict = defaultdict(int)
        ## DEADBEEF - This probably needs a maximum limit
        for row in arcpy.da.SearchCursor(
            length_path, [hru.id_field, length_field]):
            length_dict[int(row[0])] += int(row[1])
        fields = [hru.type_field, hru.iseg_field, hru.rchlen_field, hru.id_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                if (int(row[0]) == 1 and int(row[1]) <> 0):
                    row[2] = length_dict[int(row[3])]
                else:
                    row[2] = 0
                update_c.updateRow(row)
        del length_dict, length_field, fields, hru_polygon_lyr

        ## Get list of segments and downstream cell for each stream/lake cell
        ## Downstream is calulated from flow direction
        ## Use IRUNBOUND instead of ISEG, since ISEG will be zeroed for lakes
        logging.info("Cell out-flow dictionary")
        cell_dict = dict()
        fields = [
            hru.type_field, hru.krch_field, hru.lake_id_field, hru.iseg_field,
            hru.irunbound_field, hru.subbasin_field, hru.dem_adj_field,
            hru.flow_dir_field, hru.col_field, hru.row_field, hru.id_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            ## Skip inactive cells
            if int(row[0]) == 0:
                continue
            ## Skip non-lake and non-stream cells
            elif (int(row[1]) == 0 and int(row[2]) == 0):
                continue
            ## Read in parameters
            cell = (int(row[8]), int(row[9]))
            ## next_row_col(FLOW_DIR, CELL)
            ## HRU_ID, ISEG,  NEXT_CELL, DEM_ADJ, X, X, X
            cell_dict[cell] = [
                int(row[10]), int(row[4]), next_row_col(int(row[7]), cell),
                float(row[6]), 0, 0, 0]
            del cell
        ## Build list of unique segments
        iseg_list = sorted(list(set([v[1] for v in cell_dict.values()])))

        ## Calculate IREACH and OUTSEG
        logging.info("Calculate IREACH and OUTSEG")
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
            out_cell = list(set(out_cells)-set(iseg_cells))
            ## If there is more than one out_cell
            ##   there is a problem with the stream network
            if len(out_cell) <> 1:
                logging.error(
                    ('\nERROR: ISEG {0} has more than one out put cell'+
                     '\n  Out cells: {1}'+
                     '\n  Check for streams exiting then re-entering a lake'+
                     '\n  Lake cell elevations may not be constant\n').format(
                         iseg, out_cell))
                raise SystemExit()
            ## If not output cell, assume edge of domain
            try:
                outseg = cell_dict[out_cell[0]][1]
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

        ## Calculate stream elevation
        logging.info("Stream elevation (DEM_ADJ - 1 for now)")
        fields = [hru.type_field, hru.iseg_field, hru.dem_adj_field, hru.strm_top_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                if int(row[0]) == 1 and int(row[1]) <> 0:
                    row[3] = float(row[2]) - 1
                else:
                    row[3] = 0
                update_c.updateRow(row)

        ## Saving ireach and outseg
        logging.info("Save IREACH and OUTSEG")
        fields = [
            hru.type_field, hru.iseg_field, hru.col_field, hru.row_field,
            hru.outseg_field, hru.reach_field, hru.maxreach_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                ##if (int(row[0]) > 0 and int(row[1]) > 0):
                ## DEADBEEF - I'm not sure why only iseg > 0 in above line
                ## DEADBEEF - This should set outseg for streams and lakes
                if (int(row[0]) > 0 and int(row[1]) <> 0):
                    row[4:] = cell_dict[(int(row[2]), int(row[3]))][4:]
                else:
                    row[4:] = [0, 0, 0]
                update_c.updateRow(row)

        ## Calculate IUPSEG for all segments flowing out of lakes
        logging.info("IUPSEG for streams flowing out of lakes")
        upseg_dict = dict(
            [(v, k) for k, v in outseg_dict.iteritems() if k < 0])
        fields = [hru.type_field, hru.iseg_field, hru.iupseg_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                if (int(row[0]) == 1 and int(row[1]) <> 0 and
                    int(row[1]) in upseg_dict.keys()):
                    row[2] = upseg_dict[int(row[1])]
                else:
                    row[2] = 0
                update_c.updateRow(row)

        ## Build dictionary of which segments flow into each segment
        ## Used to calculate seg-basins (sub watersheds) for major streams
        ## Also save list of all segments that pour to exit
        logging.info("Segment in/out-flow dictionary")
        inseg_dict = defaultdict(list)
        pourseg_dict = dict()
        pourseg_list = []
        for key, value in outseg_dict.iteritems():
            if key == exit_seg:
                continue
                ##inseg_dict[key].append(key)
            elif value == exit_seg:
                pourseg_list.append(key)
                inseg_dict[key].append(key)
            else:
                inseg_dict[value].append(key)
        
        ## Update pourseg for each segment, working up from initial pourseg
        ## Pourseg is the final exit segment for each upstream segment
        for pourseg in pourseg_list:
            testseg_list = inseg_dict[pourseg]
            while testseg_list:
                testseg = testseg_list.pop()
                if testseg == 316:
                    print 316
                pourseg_dict[testseg] = pourseg
                if pourseg == testseg:
                    continue
                testseg_list.extend(inseg_dict[testseg])
            del testseg_list

        ## Calculate SEG_BASIN for all active cells
        ## SEG_BASIN corresponds to the ISEG of the lowest segment
        logging.info("SEG_BASIN")
        fields = [hru.type_field, hru.irunbound_field, hru.segbasin_field]
        with arcpy.da.UpdateCursor(hru.polygon_path, fields) as update_c:
            for row in update_c:
                if int(row[0]) > 0 and int(row[1]) <> 0:
                    row[2] = pourseg_dict[int(row[1])]
                else:
                    row[2] = 0
                update_c.updateRow(row)
       
        ## Set all lake iseg to 0
        logging.info("Lake ISEG")
        update_rows = arcpy.UpdateCursor(hru.polygon_path)
        for row in update_rows:
            if int(row.getValue(hru.type_field)) <> 2: continue
            iseg = int(row.getValue(hru.iseg_field))
            if iseg < 0: row.setValue(hru.iseg_field, 0)
            update_rows.updateRow(row)
            del row, iseg
        del update_rows

        ## Set environment parameters
        env.extent = hru.extent
        env.cellsize = hru.cs
        env.outputCoordinateSystem = hru.sr

        ## Build rasters
        if output_rasters_flag:
            logging.info("\nOutput model grid rasters")
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.type_field, hru_type_raster,
                "CELL_CENTER", "", hru.cs)
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.dem_adj_field, dem_adj_raster,
                "CELL_CENTER", "", hru.cs)
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.iseg_field, iseg_raster,
                "CELL_CENTER", "", hru.cs)
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.irunbound_field, irunbound_raster,
                "CELL_CENTER", "", hru.cs)
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.segbasin_field, segbasin_raster,
                "CELL_CENTER", "", hru.cs)
            arcpy.PolygonToRaster_conversion(
                hru.polygon_path, hru.subbasin_field, subbasin_raster,
                "CELL_CENTER", "", hru.cs)

        ## Build rasters
        if output_ascii_flag:
            logging.info("Output model grid ascii")
            arcpy.RasterToASCII_conversion(hru_type_raster, hru_type_ascii)
            arcpy.RasterToASCII_conversion(dem_adj_raster, dem_adj_ascii)
            arcpy.RasterToASCII_conversion(iseg_raster, iseg_ascii)
            arcpy.RasterToASCII_conversion(irunbound_raster, irunbound_ascii)
            arcpy.RasterToASCII_conversion(segbasin_raster, segbasin_ascii)
            arcpy.RasterToASCII_conversion(subbasin_raster, subbasin_ascii)
            sleep(5)

        ## Input parameters files for Cascade Routing Tool (CRT)
        logging.info("\nOutput CRT files")

        ## Generate STREAM_CELLS.DAT file for CRT
        logging.info("  {0}".format(
            os.path.basename(crt_stream_cells_path)))
        stream_cells_list = []
        fields = [
            hru.type_field, hru.iseg_field, hru.reach_field,
            hru.col_field, hru.row_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            if int(row[0]) == 1 and int(row[1]) > 0: 
                stream_cells_list.append(
                    [int(row[4]), int(row[3]), int(row[1]), int(row[2]), 1])
        if stream_cells_list:
            with open(crt_stream_cells_path, 'w+') as f:
                f.write('{0}    NREACH\n'.format(len(stream_cells_list)))
                for stream_cells_l in sorted(stream_cells_list):
                    f.write(' '.join(map(str, stream_cells_l))+'\n')
            f.close
        del stream_cells_list
            
        ## Generate OUTFLOW_HRU.DAT for CRT
        logging.info("  {0}".format(
            os.path.basename(crt_outflow_hru_path)))
        outflow_hru_list = []
        fields = [
            hru.type_field, hru.iseg_field, hru.outseg_field, hru.reach_field,
            hru.maxreach_field, hru.col_field, hru.row_field]
        for row in arcpy.da.SearchCursor(hru.polygon_path, fields):
            if int(row[0]) <> 1 or int(row[1]) == 0: continue
            if int(row[2]) == 0 and int(row[3]) == int(row[4]): 
                outflow_hru_list.append([int(row[6]), int(row[5])])
        if outflow_hru_list:
            with open(crt_outflow_hru_path, 'w+') as f:
                f.write('{0}    NUMOUTFLOWHRU\n'.format(
                    len(outflow_hru_list)))
                for i, outflow_hru in enumerate(outflow_hru_list):
                    f.write('{0} {1} {2}   OUTFLOW_ID ROW COL\n'.format(
                        i+1, outflow_hru[0], outflow_hru[1]))
            f.close()
        del outflow_hru_list

        ## Generate HRU_CASC.DAT for CRT
        logging.info("  {0}".format(os.path.basename(crt_hru_casc_path)))
        with open(hru_type_ascii, 'r') as f: ascii_data = f.readlines()
        f.close()
        hru_casc_header = (
            '{0} {1} {2} {3} {4} {5} {6} {7}     '+
            'HRUFLG STRMFLG FLOWFLG VISFLG '+
            'IPRN IFILL DPIT OUTITMAX\n').format(
                crt_hruflg, crt_strmflg, crt_flowflg, crt_visflg,
                crt_iprn, crt_ifill, crt_dpit, crt_outitmax)
        with open(crt_hru_casc_path, 'w+') as f:
            f.write(hru_casc_header)
            for ascii_line in ascii_data[6:]: f.write(ascii_line)
        f.close()
        del hru_casc_header, ascii_data

        ## Generate LAND_ELEV.DAT for CRT
        logging.info("  {0}".format(os.path.basename(crt_land_elev_path)))
        with open(dem_adj_ascii, 'r') as f: ascii_data = f.readlines()
        f.close()
        with open(crt_land_elev_path, 'w+') as f:
            f.write('{0} {1}       NROW NCOL\n'.format(
                ascii_data[1].split()[1], ascii_data[0].split()[1]))
            for ascii_line in ascii_data[6:]: f.write(ascii_line)
        f.close()
        del ascii_data

        ## Generate XY.DAT for CRT
        logging.info("  {0}".format(os.path.basename(crt_xy_path)))
        xy_list = [
            map(int, row)
            for row in sorted(arcpy.da.SearchCursor(
                hru.polygon_path, [hru.id_field, hru.x_field, hru.y_field]))]
        with open(crt_xy_path, 'w+') as f:
            for line in sorted(xy_list):
                f.write(' '.join(map(str, line))+'\n')
        f.close()

        ## Run CRT
        logging.info('\nRunning CRT')
        os.chdir(crt_ws)
        subprocess.check_call(crt_exe_name)
        os.chdir(workspace)

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

    ## Calculate GSFLOW Stream Parameters
    gsflow_stream_parameters(workspace, ini_path)
