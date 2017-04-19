#! /usr/bin/env python3

# System
import mysql
import pickle
import sys
import json
from json import JSONEncoder
import argparse

# Data Structures
import collections
from collections import Counter
from collections import defaultdict
import os
from os import path

# CPLOP
sys.path.append(os.getcwd())
import fullsearch
import clusterEval
import dbscan
import config
from config import Config
import pyroprinting
from pyroprinting import Isolate
from pyroprinting import Region

import strainpurity
from strainpurity import MicrobialSourceTracking

DESCRIPTION='Performs MST on the isolates in the isolate configuration file.'
def get_arg_parser():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description=DESCRIPTION)
    parser.add_argument('isolate_config',
            metavar='<isolate configuration file>', 
            help='path (absolute or relative) to the configuration file containing isolates to perform MST on')
    return parser

def do_the_mst(min_neighbors, isolate):
    cfg = None
    cfg_json = None
    db_cfg = None

    # Open Configuration Files
    with open("clusterConfig.json", mode='r') as config_file:
        cfg = Config(config_file)
    with open("clusterConfig.json", mode='r') as config_file:
        cfg_json = json.load(config_file)
    assert(cfg)
    with open("mysqlConfig.json", mode='r') as file:
        db_cfg = json.load(file)
    # Query Database for Isolates
    assert(db_cfg)
    db_cnx = mysql.connector.connect(**db_cfg)

    # Do the MST
    for n in min_neighbors:
        cfg.minNeighbors = n
        mst = MicrobialSourceTracking(cfg, db_cnx, cfg_json, isolate)
    return

def main():
    # CONFIG AND DB LOADING
    parser = get_arg_parser()
    args = parser.parse_args()

    # Get Isolates
    isolate_config = args.isolate_config
    isolate_json = None
    with open(isolate_config, 'r') as file:
        isolate_json = json.load(file)
    isolates = isolate_json['isolates']

    # Run MST
    for i in isolates:
        do_the_mst([1, 2, 3, 4, 5, 6, 7, 8], i)
        #do_the_mst([3], i)
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)

