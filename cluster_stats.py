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

sys.path.append(os.getcwd())
import args
from args import get_stats_args

import config
from config import Config

import pyroprinting
from pyroprinting import Isolate
from pyroprinting import Region

import strain_purity
from strain_purity import ClusterPurity
from strain_purity import StrainPurities

class StrainPurityStats(object):
    def __init__(self, purities, min_size):
        min = None
        max = None
        num_clusters = 0
        sum_clust_size = 0
        sum_most_plural = 0
        most_plural = defaultdict(list)
        #most_plural = defaultdict(dict)
        above = {}
        above['40%'] = {}
        above['40%']['count'] = 0
        above['40%']['clusters'] = []
        above['40%']['species'] = set()
        above['60%'] = {}
        above['60%']['count'] = 0
        above['60%']['clusters'] = []
        above['60%']['species'] = set()
        above['80%'] = {}
        above['80%']['count'] = 0
        above['80%']['clusters'] = []
        above['80%']['species'] = set()

        purities.filter_low(min_size)
        for p in purities.purities:
            num_clusters += 1
            sum_clust_size += p.cluster_size
            sum_most_plural += p.most_plural['purity']
            label = p.most_plural['label']
            purity = p.most_plural['purity']
            count = p.most_plural['count']
            cluster = p.cluster_num
            most_plural[label].append({'cluster': cluster, 'purity ': purity, 'count  ': count})
            #most_plural[label][cluster] = (purity, count)
            if purity > .40:
                above['40%']['count'] += 1
                above['40%']['clusters'].append(cluster)
                above['40%']['species'].add(label)
            if purity > .60:
                above['60%']['count'] += 1
                above['60%']['clusters'].append(cluster)
                above['60%']['species'].add(label)
            if purity > .80:
                above['80%']['count'] += 1
                above['80%']['clusters'].append(cluster)
                above['80%']['species'].add(label)
            if min is None or p.cluster_size < min['size'] :
                min = {}
                min['size'] = p.cluster_size
                min['clusters'] = [p.cluster_num]
            elif p.cluster_size is min['size']:
                min['clusters'].append(p.cluster_num)
            if max is None or p.cluster_size > max['size'] :
                max = {}
                max['size'] = p.cluster_size
                max['clusters'] = [p.cluster_num]
            elif p.cluster_size is max['size']:
                max['clusters'].append(p.cluster_num)
        for species, cluster_list in most_plural.items():
            most_plural[species] = sorted(cluster_list, key=lambda k: k['purity '], reverse=True)
        above['40%']['species'] = sorted(above['40%']['species'])
        above['60%']['species'] = sorted(above['60%']['species'])
        above['80%']['species'] = sorted(above['80%']['species'])
        self.num = num_clusters
        self.avg_size = sum_clust_size / num_clusters
        self.avg_purity = sum_most_plural / num_clusters
        self.purest = most_plural
        self.above_percent = above
        self.min = min
        self.max = max

def main():
    # PARSE ARGS
    arg_parser = get_stats_args()
    args = arg_parser.parse_args()
    filename = args.filename
    k = args.min_size
    # PRINT ARGS
    print('Filename   : {}'.format(filename))
    print('Min Cluster: {}'.format(k))
    # OPEN STRAIN PURITIES
    purities = None
    with open(filename, 'r+b') as file:
        purities = pickle.load(file)
    stats = StrainPurityStats(purities, k)
    print(json.dumps(stats.__dict__, indent=2, sort_keys=True))
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
