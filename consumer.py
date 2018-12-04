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

import sys
import gc
import csv
import types
import os
import json
from datetime import datetime
from collections import defaultdict, OrderedDict
import numpy as np

import zmq
#print "pyzmq version: ", zmq.pyzmq_version(), " zmq version: ", zmq.zmq_version()

LOCAL_CONSUMER = True

PATHS = {
    "berg": {
        "local_path_to_data_dir": "D:/",
        "cluster_path_to_data_dir": "/archiv-daten/md/data/",
        "local_path_to_output_dir": "out/",
        "cluster_path_to_output_dir": "out/"
    }
}

def run_consumer(path_to_output_dir = None, leave_after_finished_run = True, server = {"server": None, "port": None, "nd-port": None}, shared_id = None):
    "collect data from workers"

    config = {
        "user": "berg",
        "port": server["port"] if server["port"] else "7777",
        "server": server["server"] if server["server"] else "localhost", 
        "region": "quillow", #"ddr",
        "ref_mmk_type": "stt",
        "start_row": "0",
        "end_row": "0", #"-1",
        "shared_id": shared_id,
        "out": "out/" #None,
    }
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k,v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["user"]]
    path_to_data_dir = paths["local_path_to_data_dir"] if LOCAL_CONSUMER else paths["cluster_path_to_data_dir"]

    print "consumer config:", config

    context = zmq.Context()
    if config["shared_id"]:
        socket = context.socket(zmq.DEALER)
        socket.setsockopt(zmq.IDENTITY, config["shared_id"])
    else:
        socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])

    #socket.RCVTIMEO = 1000
    leave = False
    write_normal_output_files = False

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

    template_metadata, template_header = read_header(path_to_data_dir + config["region"] + "/" + config["ref_mmk_type"] + "_" + config["region"] + "_100_gk5.asc")
    print "read template metadata from:", path_to_data_dir + config["region"] + "/" + config["ref_mmk_type"] + "_" + config["region"] + "_100_gk5.asc"

    start_row = int(config["start_row"])
    end_row = int(config["end_row"])
    nrows = int(template_metadata["nrows"])
    ncols = int(template_metadata["ncols"])
    nodata_value = int(template_metadata["nodata_value"])
    
    new_grid = lambda: np.full((nrows, ncols), nodata_value)
    grids = defaultdict(lambda: defaultdict(new_grid))
    
    def process_message(msg):

        leave = False

        if msg["type"] == "finish":
            print "c: received finish message"
            leave = True

        elif not write_normal_output_files:

            process_message.received_env_count += 1
            
            custom_id = msg["customId"]
                        
            row = custom_id["row"]
            col = custom_id["col"]
            if not process_message.no_of_datacells:
                process_message.no_of_datacells = custom_id.get("ndatacells", None)
            
            leave = process_message.no_of_datacells == process_message.received_env_count

            print "env-count/no-datacells:", process_message.received_env_count, "/", process_message.no_of_datacells, ", leave:", leave

            if msg["runFailed"]:
                print "run with customId:", custom_id, "failed. Reason:", msg["reasonForRunFailed"]
                return leave

            for year, crop_result in msg["year2cropResult"].iteritems():
                if not crop_result["isNoData"]:
                    for res_id, value in crop_result["values"].iteritems():
                        grids[year][res_id][row, col] = value
            
            #if process_message.received_env_count % 10 == 0:
            #print process_message.received_env_count, "/", process_message.no_of_datacells,
            
        elif write_normal_output_files:
            process_message.received_env_count += 1
            
            #print "received work result ", process_message.received_env_count, " customId: ", str(msg.get("customId", "").values())
            print process_message.received_env_count,
            print msg

        return leave

    process_message.no_of_datacells = None
    process_message.received_env_count = 0

    while not leave:
        try:
            msg = socket.recv_json(encoding="latin-1")
            leave = process_message(msg)
        except Exception as e: 
            print(e)
            continue

    res_id_to_avgs = {}
    for year, res_id_to_grid in grids.iteritems():
        for res_id, grid in res_id_to_grid.iteritems():
            if res_id not in res_id_to_avgs:
                res_id_to_avgs[res_id] = (np.full((nrows, ncols), 0.0), 0)
            res_id_to_avgs[res_id] = (res_id_to_avgs[res_id][0] + grid, res_id_to_avgs[res_id][1] + 1)

            np.savetxt(config["out"] + res_id + "_" + str(year) + ".asc", grid, delimiter=" ", fmt="%.2f", header=template_header.strip(), comments="")

    for res_id, (avg_grid, count) in res_id_to_avgs.iteritems():
        avg = avg_grid / count

        np.savetxt(config["out"] + res_id + "_avg.asc", avg, delimiter=" ", fmt="%.2f", header=template_header.strip(), comments="")


    print "exiting run_consumer()"
    #debug_file.close()

if __name__ == "__main__":
    run_consumer()
#main()


