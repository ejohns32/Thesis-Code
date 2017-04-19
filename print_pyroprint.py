#! /usr/bin/env python3

import mysql
import pickle
import sys
import json
from json import JSONEncoder

import collections
from collections import Counter
from collections import defaultdict
import os
from os import path

sys.path.append(os.getcwd())
import fullsearch
import clusterEval
import dbscan
import config
from config import Config
import pyroprinting
from pyroprinting import Isolate
from pyroprinting import Region

def main():
    cfg = None
    cfg_json = None
    with open("clusterConfig.json", mode='r') as config_file:
        cfg = Config(config_file)
    with open("clusterConfig.json", mode='r') as config_file:
        cfg_json = json.load(config_file)
    db_cfg = None
    with open("mysqlConfig.json", mode='r') as file:
        db_cfg = json.load(file)
    assert(cfg)
    # Query Database for Isolates
    assert(db_cfg)
    db_cnx = mysql.connector.connect(**db_cfg)
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
