#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import time
import os
import math
import json
import csv
import copy
from StringIO import StringIO
from datetime import date, datetime, timedelta
from collections import defaultdict
import types
import sys
#print sys.path
import zmq
#print "pyzmq version: ", zmq.pyzmq_version(), " zmq version: ", zmq.zmq_version()
import re

import sqlite3
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from pyproj import Proj, transform

LOCAL_PRODUCER = True
LOCAL_YIELDSTAT = True

PATHS = {
    "berg": {
        "local_path_to_data_dir": "D:/",
        "cluster_path_to_data_dir": "/archiv-daten/md/data/",
        "local_path_to_output_dir": "out/",
        "cluster_path_to_output_dir": "out/"
    }
}

def run_producer(setup = None, custom_crop = None, server = {"server": None, "port": None, "nd-port": None}, shared_id = None):
    "main"

    context = zmq.Context()
    socket = context.socket(zmq.PUSH)

    config = {
        "user": "berg",
        "port": server["port"] if server["port"] else "6666",
        "server": server["server"] if server["server"] else "localhost",
        "shared_id": shared_id,
        "region": "quillow", #"ddr",
        "ref_mmk_type": "stt",
        "start_row": "0",
        "end_row": "-1",
        "start_year": "1991",
        "end_year": "2012",
        "crop_rotation": "1017pi,1013n", # "1017ci"
        "climate_scenario": "A1B",
        "return_corn_units": False,
        "use_dev_trend": False,
        "trend_base_year": "2005",
        "use_co2_increase": True,
        "get_dry_year_water_need": False,
        "debug_mode": False
    }
    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print "config:", config

    paths = PATHS[config["user"]]
    path_to_data_dir = paths["local_path_to_data_dir"] if LOCAL_PRODUCER else paths["cluster_path_to_data_dir"]
    path_to_yieldstat_climate_dir = paths["local_path_to_data_dir"] + "climate/" if LOCAL_YIELDSTAT else paths["cluster_path_to_data_dir"] + "climate/"
    
    socket.connect("tcp://" + config["server"] + ":" + str(config["port"]))

    def read_header(path_to_ascii_grid_file):
        "read metadata from esri ascii grid file"
        metadata = {}
        header_str = ""
        with open(path_to_ascii_grid_file) as _:
            for i in range(0, 6):
                line = _.readline()
                header_str += line
                sline = [x for x in line.split() if len(x) > 0]
                if len(sline) > 1:
                    metadata[sline[0].strip().lower()] = float(sline[1].strip())
        return metadata, header_str

    wgs84 = Proj(init="epsg:4326")
    #gk3 = Proj(init="epsg:3396")
    gk5 = Proj(init="epsg:31469")

    def create_ascii_grid_interpolator(arr, meta, ignore_nodata=True):
        "create interpolator from numpy array"

        rows, cols = arr.shape

        cellsize = int(meta["cellsize"])
        xll = int(meta["xllcorner"])
        yll = int(meta["yllcorner"])
        nodata_value = meta["nodata_value"]

        xll_center = xll + cellsize // 2
        yll_center = yll + cellsize // 2
        yul_center = yll_center + (rows - 1)*cellsize

        points = []
        values = []

        for row in range(rows):
            for col in range(cols):
                value = arr[row, col]
                if ignore_nodata and value == nodata_value:
                    continue
                r = xll_center + col * cellsize
                h = yul_center - row * cellsize
                points.append([r, h])
                values.append(value)

        return NearestNDInterpolator(np.array(points), np.array(values))

    def read_file_and_create_interpolator(path_to_grid, dtype=int, skiprows=6, confirm_creation=False):
        "read file and metadata and create interpolator"

        metadata, _ = read_header(path_to_grid)
        grid = np.loadtxt(path_to_grid, dtype=dtype, skiprows=skiprows)
        interpolate = create_ascii_grid_interpolator(grid, metadata)
        if confirm_creation: 
            print "created interpolator from:", path_to_grid
        return (interpolate, grid, metadata)
    
    gk5_ref_grid = None
    ref_metadata = None
    gk5_interpolators = {}
    gk5_grids = {}
    for mmk_type, dtype in [
        ("hft", int), ("nft", int), ("sft", int), ("steino", int), ("slope", int),
        ("stt", int), ("dgm", int), ("az", int), ("klz", int)
    ]:
        path_to_grid = path_to_data_dir + config["region"] + "/" + mmk_type + "_" + config["region"] + "_100_gk5.asc"
        inter, grid, metadata = read_file_and_create_interpolator(path_to_grid, confirm_creation=True)
        gk5_interpolators[mmk_type] = inter
        if mmk_type == config["ref_mmk_type"]:
            gk5_ref_grid = grid
            ref_metadata = metadata

    cdict = {}
    def create_climate_gk5_interpolator_from_json_file(path_to_latlon_to_rowcol_file, wgs84, gk5):
        "create interpolator from json list of lat/lon to row/col mappings"
        with open(path_to_latlon_to_rowcol_file) as _:
            points = []
            values = []

            for latlon, rowcol in json.load(_):
                row, col = rowcol
                clat, clon = latlon
                try:
                    cr_gk5, ch_gk5 = transform(wgs84, gk5, clon, clat)
                    cdict[(row, col)] = (round(clat, 4), round(clon, 4))
                    points.append([cr_gk5, ch_gk5])
                    values.append((row, col))
                    #print "row:", row, "col:", col, "clat:", clat, "clon:", clon, "h:", h, "r:", r, "val:", values[i]
                except:
                    continue

            return NearestNDInterpolator(np.array(points), np.array(values))

    climate_gk5_interpolator = create_climate_gk5_interpolator_from_json_file(path_to_data_dir + "climate/dwd/csvs/latlon_to_rowcol.json", wgs84, gk5)
    print "created climate gk5 interpolator:", path_to_data_dir + "climate/dwd/csvs/latlon_to_rowcol.json"

    sent_env_count = 1
    start_time = time.clock()

    env_template = {
        "type": "Yieldstat::Core::Env",
        "climateScenario": config["climate_scenario"],
        "startYear": int(config["start_year"]),
        "endYear": int(config["end_year"]),
        "trendBaseYear": int(config["trend_base_year"]),
        "useDevTrend": config["use_dev_trend"],
        "useCO2Increase": config["use_co2_increase"],
        "returnCornUnits": config["return_corn_units"],
        "getDryYearWaterNeed": config["get_dry_year_water_need"]
    }

    def parse_crop(crop_string):
        m = re.search(r"(\d{3,4})(p|n|c)(i?)", crop_string)
        try:
            return {
                "type": "Yieldstat::Core::YSCrop",
                "id": int(m.group(1)),
                "tillageType": ("conserving" if m.group(2) == "c" else ("noTillage" if m.group(2) == "n" else "plough")),
                "irrigate": True if m.group(3) == "i" else False    
            }
        except:
            print "The crop string:", crop_string, "doesn't resemble a valid crop! It will be ignored."
            return None

    # create crop rotation
    env_template["cropRotation"] = map(parse_crop, config["crop_rotation"].split(","))

    if config["shared_id"]:
        env_template["sharedId"] = config["shared_id"]

    rcols = int(ref_metadata["ncols"])
    rrows = int(ref_metadata["nrows"])
    rcellsize = int(ref_metadata["cellsize"])
    xllcorner = int(ref_metadata["xllcorner"])
    yllcorner = int(ref_metadata["yllcorner"])

    for rrow in xrange(0, rrows):
        print rrow,

        if rrow < int(config["start_row"]):
            continue
        elif int(config["end_row"]) > 0 and srow > int(config["end_row"]):
            break

        for rcol in xrange(0, rcols):

            stt = gk5_ref_grid[rrow, rcol]
            if stt == -9999:
                continue
            
            #get coordinate of clostest climate element of real soil-cell
            rh_gk5 = yllcorner + (rcellsize / 2) + (rrows - rrow - 1) * rcellsize
            rr_gk5 = xllcorner + (rcellsize / 2) + rcol * rcellsize

            crow, ccol = climate_gk5_interpolator(rr_gk5, rh_gk5)

            env_template["dgm"] = gk5_interpolators["dgm"](rr_gk5, rh_gk5)
            env_template["hft"] = gk5_interpolators["hft"](rr_gk5, rh_gk5)
            env_template["nft"] = gk5_interpolators["nft"](rr_gk5, rh_gk5)
            env_template["sft"] = gk5_interpolators["sft"](rr_gk5, rh_gk5)
            env_template["slope"] = gk5_interpolators["slope"](rr_gk5, rh_gk5)
            env_template["steino"] = gk5_interpolators["steino"](rr_gk5, rh_gk5)
            env_template["az"] = gk5_interpolators["az"](rr_gk5, rh_gk5)
            env_template["klz"] = gk5_interpolators["klz"](rr_gk5, rh_gk5)
            env_template["stt"] = gk5_interpolators["stt"](rr_gk5, rh_gk5)

            env_template["csvViaHeaderOptions"] = {
                "start-date": config["start_year"] + "-01-01",
		        "end-date": config["end_year"] + "-12-31",
        		"no-of-climate-file-header-lines": 2,
		        "csv-separator": ","#,
		        #"header-to-acd-names": {
			    #    "DE-date": "de-date",
			    #    "globrad": ["globrad", "/", 100]
		        #}
            }

            env_template["pathToClimateCSV"] = path_to_yieldstat_climate_dir + "dwd/csvs/germany/row-" + str(crow+1) + "/col-" + str(ccol+1) + ".csv"
            #print env_template["pathToClimateCSV"]

            env_template["customId"] = {
                "row": rrow, "col": rcol,
                "crow": crow, "ccol": ccol
            }

            socket.send_json(env_template)
            #print env_template
            print "sent env ", sent_env_count, " customId: ", env_template["customId"]
            #exit()
            sent_env_count += 1


    stop_time = time.clock()

    print "sending ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds"
    #print "ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
    print "exiting run_producer()"

if __name__ == "__main__":
    run_producer()