#! /usr/bin/python3

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
    # CONFIG AND DB LOADING
    do_the_thing([1, 2, 3, 4, 5, 6, 7])
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
