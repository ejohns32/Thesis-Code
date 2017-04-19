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
from pprint import pformat

sys.path.append(os.getcwd())
import fullsearch
import clusterEval
import dbscan
import config
from config import Config
import pyroprinting
from pyroprinting import Isolate
from pyroprinting import Region

"""
    CLASSES
"""
"""
    ClusterPurity
    Take an identifier to label this cluster (a number, name, whatever) an
    iterable object representing a cluster, and a map that maps that object to
    the desired label to count.
"""
class ClusterPurity(object):
    def __init__(self, cluster_num, cluster, map):
        self.cluster_num = cluster_num
        #self.cluster = [c for c in cluster]
        self.cluster_size = len(cluster)
        self.counts = defaultdict(int)
        #self.members = {}
        self.purity = {}
        isol_set = set()
        for c in cluster:
            #self.members[c.name] = map[c]
            self.counts[map[c]] += 1
        self.most_plural = None
        for label, count in self.counts.items():
            self.purity[label] = count / self.cluster_size
            if self.most_plural is None or self.most_plural['count'] < count:
                self.most_plural = {}
                self.most_plural['label']  = label
                self.most_plural['purity'] = count / self.cluster_size
                self.most_plural['count']  = count
        self.sorted_purities = sorted(['%.3f %3d %s' % (self.purity[k],
            self.counts[k], k) for k in self.counts], reverse=True)
        self.cluster_purity = self.most_plural['count'] / self.cluster_size
    def __str__(self):
        return self.counts.__str__()

class MicrobialSourceTracking(object):
    def __init__(self, cfg, db_cnx, cfg_json, unknown_isolate_name):
        isolates = get_isolates(cfg)
        iso_map = get_isolate_species(db_cnx)
        print('Looking for {}...'.format(unknown_isolate_name))
        unknown_isolate = None
        for i in isolates:
            if i.name == unknown_isolate_name:
                print('Found {}.'.format(unknown_isolate_name))
                unknown_isolate = i
        if not unknown_isolate:
            print('FUCK can\'t find {}'.format(unknown_isolate_name))

        isolates = filter_isolates(isolates, iso_map)
        isolates.append(unknown_isolate)

        print('Caching DBScan Clusters...')
        cache_mst_clusters(cfg, isolates, unknown_isolate_name)

        print('Getting Clusters...')
        clusters = get_clusters(cfg)
        print('Done Getting Clusters')

        # Do Stats
        clustered_isolates = set()
        multiply_clustered_isolates = set()
        unclustered_isolates = set()
        self.species_counts = defaultdict(int)
        for cluster in clusters:
            if len(cluster) <= 1:
                print('NOISE:{}'.format(cluster))
            for i in cluster:
                if i in clustered_isolates:
                    #print('ALREADY SEEN:{}'.format(i))
                    multiply_clustered_isolates.add(i)
                clustered_isolates.add(i)
        for i in isolates:
            if i not in clustered_isolates:
                #print('DID NOT CLUSTER:{}'.format(i))
                unclustered_isolates.add(i)
            self.species_counts[iso_map[i]] += 1
        if unknown_isolate in clustered_isolates:
            print('{} was clustered!!'.format(unknown_isolate))
        if unknown_isolate in multiply_clustered_isolates:
            print('{} was clustered!!'.format(unknown_isolate))
        if unknown_isolate in unclustered_isolates:
            print('{} was NOT clustered!!'.format(unknown_isolate))
        self.num_clustered = len(clustered_isolates)
        self.num_unclustered = len(unclustered_isolates)
        self.num_multiply_clustered = len(multiply_clustered_isolates)
        self.num_isolates = len(isolates)
        print('Clustered:{}'.format(self.num_clustered))
        print('UnClustered:{}'.format(self.num_unclustered))
        print('SUM:{}'.format(self.num_clustered + self.num_unclustered))
        print('Isolates:{}'.format(self.num_isolates))
        print('Intersection:{}'.format([c for c in clustered_isolates.intersection(unclustered_isolates)]))
        isol_set = {i for i in isolates}
        counted = clustered_isolates.union(unclustered_isolates)
        print('WTF:{}'.format([c for c in counted if c not in isol_set]))
        assert(self.num_clustered + self.num_unclustered == self.num_isolates)
        self.cfg_json = cfg_json
        self.clusters = len(clusters)
        self.isolates = len(isolates)
        self.purities = [ClusterPurity(i, c, iso_map)
                         for i, c in zip(range(len(clusters)), clusters)]
        self.purities = sorted(self.purities, key=lambda k: k.cluster_size, reverse=True)
        self.clustering_purity = 0.0;
        self.clustering_size = 0
        for c in self.purities:
            if c.cluster_size <= 1:
                continue
            self.clustering_purity += c.cluster_size * c.cluster_purity
            self.clustering_size += c.cluster_size
        self.clustering_purity = self.clustering_purity / self.clustering_size
        return

def cache_mst_clusters(cfg, isolates, unknown_isolate_name):
    cacheFileName = dbscan.getDBscanClustersCacheFileName(cfg)
    cacheFileName = cacheFileName.split('.')
    cacheFileName[0] = cacheFileName[0] + unknown_isolate_name
    cacheFileName = '.'.join(cacheFileName)
    print('Checking for {}'.format(cacheFileName))
    if os.path.isfile(cacheFileName):
        return
    print('Clustering for {}'.format(unknown_isolate_name))
    clusters = dbscan.computeDBscanClusters(isolates, cfg)
    with open(cacheFileName, mode='w+b') as cacheFile:
        pickle.dump(clusters, cacheFile)

"""
    StrainPurities
    From an iterable object of isolates and clusters, it uses a map from
    isolates to species to count the number of each occurrence in each cluster.
"""
class StrainPurities(object):
    def __init__(self, cfg, db_cnx, cfg_json):
        isolates = get_isolates(cfg)
        iso_map = get_isolate_species(db_cnx)
        isolates = filter_isolates(isolates, iso_map)

        cache_dbscan_clusters(cfg, isolates)
        clusters = get_clusters(cfg)
        clustered_isolates = set()
        multiply_clustered_isolates = set()
        unclustered_isolates = set()
        self.species_counts = defaultdict(int)
        for cluster in clusters:
            if len(cluster) <= 1:
                print('NOISE:{}'.format(cluster))
            for i in cluster:
                if i in clustered_isolates:
                    print('ALREADY SEEN:{}'.format(i))
                    multiply_clustered_isolates.add(i)
                clustered_isolates.add(i)
        for i in isolates:
            if i not in clustered_isolates:
                print('DID NOT CLUSTER:{}'.format(i))
                unclustered_isolates.add(i)
            self.species_counts[iso_map[i]] += 1
        self.num_clustered = len(clustered_isolates)
        self.num_unclustered = len(unclustered_isolates)
        self.num_multiply_clustered = len(multiply_clustered_isolates)
        self.num_isolates = len(isolates)
        print('Clustered:{}'.format(self.num_clustered))
        print('UnClustered:{}'.format(self.num_unclustered))
        print('SUM:{}'.format(self.num_clustered + self.num_unclustered))
        print('Isolates:{}'.format(self.num_isolates))
        print('Intersection:{}'.format([c for c in clustered_isolates.intersection(unclustered_isolates)]))
        isol_set = {i for i in isolates}
        counted = clustered_isolates.union(unclustered_isolates)
        print('WTF:{}'.format([c for c in counted if c not in isol_set]))
        assert(self.num_clustered + self.num_unclustered == self.num_isolates)
        self.cfg_json = cfg_json
        self.clusters = len(clusters)
        self.isolates = len(isolates)
        self.purities = [ClusterPurity(i, c, iso_map)
                         for i, c in zip(range(len(clusters)), clusters)]
        self.purities = sorted(self.purities, key=lambda k: k.cluster_size, reverse=True)
        self.clustering_purity = 0.0;
        self.clustering_size = 0
        for c in self.purities:
            if c.cluster_size <= 1:
                continue
            self.clustering_purity += c.cluster_size * c.cluster_purity
            self.clustering_size += c.cluster_size
        self.clustering_purity = self.clustering_purity / self.clustering_size
    def __len__(self):
        return len(self.purities)
    def get_cluster_purity(self, i):
        return self.purities[i]
    def filter_low(self, min):
        filtered_purities = []
        filtered_isolates = 0
        for p in self.purities:
            if p.cluster_size >= min:
                filtered_purities.append(p)
                filtered_isolates += p.cluster_size
        self.purities = filtered_purities
        self.clusters = len(filtered_purities)
        self.isolates = filtered_isolates
        return self.purities
"""
   Encodes all of the data structures relevant to the StrainPurities object
   (StrainPurities, ClusterPurity, Isolates) into JSONable objects.
"""
class StrainPurityJSONEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return [x for x in obj]
        elif isinstance(obj, Isolate):
            return obj.name
        elif isinstance(obj, ClusterPurity):
            new_dict = {}
            for k in obj.__dict__:
                if k is not 'counts' and k is not 'purity' and k is not 'members':
                    new_dict[k] = obj.__dict__[k]
            return obj.__dict__
        elif isinstance(obj, StrainPurities):
            return obj.__dict__
        elif isinstance(obj, Config):
            return obj.__dict__
        elif isinstance(obj, Region):
            return obj.__repr__()
        return JSONEncoder.default(self, obj)

"""
   Reads in an Isolate's pickling in the way cachethings.py does and returns
   the isolates object initially pickled.
"""
def get_isolates(cfg):
    isolates = None
    cacheFileName = "isolatesAll.pickle"
    if os.path.isfile(cacheFileName):
        with open(cacheFileName, mode='rb') as cacheFile:
            isolates = pickle.load(cacheFile)
    else:
        isolates = pyroprinting.loadIsolates(cfg)
        with open(cacheFileName, mode='w+b') as cacheFile:
            pickle.dump(isolates, cacheFile)
    assert(isolates)
    return isolates

"""
    FUNCTIONS
"""
"""
   Reads in a dbscan pickling according to the config dictionary passed in
   Returns the resulting clusters object initially pickled.
"""
def get_clusters(cfg):
    clusters = None
    dbscan_filename = dbscan.getDBscanClustersCacheFileName(cfg)
    with open(dbscan_filename, 'r+b') as file:
        clusters = pickle.load(file)
    assert(clusters)
    return clusters

"""
    Caches the DBSCAN calculated on a collection of isolates according to the
    config dictionary passed in.
"""
def cache_dbscan_clusters(cfg, isolates):
    cacheFileName = dbscan.getDBscanClustersCacheFileName(cfg)
    if os.path.isfile(cacheFileName):
        return
    clusters = dbscan.computeDBscanClusters(isolates, cfg)
    with open(cacheFileName, mode='w+b') as cacheFile:
        pickle.dump(clusters, cacheFile)

"""
    Filters out the environmental isolates.
"""
environmental_filter = set([
    "SLO Creek Water",
    "unknown",
    "Unknown",
    "Unknown Avian",
    "unknown avian",
    "unknown Avian",
    "Unknown avian",
    "Pennington Creek Water",
    "Pacific Ocean Water",
    'Human and Cow Mix',
    'Human and Dog Mix',
    'Cow and Dog Mix',
    'Human, Cow, and Dog Mix'
    ])
filter = set(environmental_filter)

"""
    Takes an iterable collection of isolates, an isolate->whatever map, and a
    whatever filter
    Returns a list of isolates, filtering out the isolates who map to wahtever
    is in the filter
"""
def filter_isolates(isolates, iso_map, filter=filter):
    filtered_isolates = []
    for iso in isolates:
        if iso in iso_map:
            if iso_map[iso] not in filter:
                filtered_isolates.append(iso)
            #else:
                #print('Removing {} becasue {}'.format(iso, iso_map[iso]))
        else:
            print('Couldn\'t find {} in iso_map'.format(iso))
    return filtered_isolates

"""
   Queries the CPLOP database to get a table of the isoID's and commonName for
   its host species and builds a dictionary that maps isoID's to commonName's.
"""
isolate_species_query = \
    ' \
SELECT \
    iso.isoID, \
    iso.commonName \
    FROM \
        Isolates as iso \
;'
known_empty_isolates = {
    'Ck-001': 'Chicken',
    'Ck-002': 'Chicken',
    'Ck-003': 'Chicken',
    'Ck-004': 'Chicken',
    'Ck-005': 'Chicken',
    'Ct-102': 'Cat',
    'Ct-103': 'Cat',
    'Ct-104': 'Cat',
    'Ct-105': 'Cat',
    'CW-1868': 'Cow',
    'Ho-001': 'Horse',
    'Hu-077': 'Human',
    'Hu-078': 'Human'
}
remap = {
    'Cw': 'Cow',
    'Dg': 'Dog',
    'Hu': 'Human',
    'Hu and Cw': 'Human and Cow Mix',
    'Hu and Dg': 'Human and Dog Mix',
    'Cw and Dg': 'Cow and Dog Mix',
    'Hu, Cw and Dg': 'Human, Cow, and Dog Mix',
    'cougar': 'Cougar',
    'orangutan': 'Orangutan',
    'Human ': 'Human',
    'Neonatal Human': 'Human',
    'Human UTI': 'Human',
    'Barn Owl': 'Owl',
    'Long-eared Owl': 'Owl',
    'Wild Pig': 'Pig',
    'Pig/Swine': 'Pig',
    'Bond tail Pigeon': 'Pigeon'
}
"""
    Takes in a database connection
    Returns a mapping of isolates to species
"""
def get_isolate_species(cnx):
    print('Building IsoID -> CommonName map...')
    isolate_species = {}
    # Map known empty isolates
    for iso_id, common_name in known_empty_isolates.items():
        #print("{}: {} (EMPTY)".format(iso_id, common_name))
        isolate_species[Isolate(iso_id, None)] = common_name
    # Query DB for Isolate:CommonName
    cursor = cnx.cursor()
    cursor.execute(isolate_species_query)
    for (iso_id, common_name) in cursor:
        #print("{}->{}".format(iso_id, common_name))
        if common_name in remap:   # Replace with known mappings
            common_name = remap[common_name]
        #print("{}: {}".format(iso_id, common_name))
        isolate_species[Isolate(iso_id, None)] = common_name
    return isolate_species

isolate_host_query = \
    ' \
SELECT \
    iso.isoID, \
    iso.hostID \
    FROM \
        Isolates as iso \
;'
def get_isolate_host(cnx):
    print('Building IsoID -> CommonName map...')
    isolate_host = {}
    # Map known empty isolates
    #for iso_id, host_id in known_empty_isolates.items():
    #    #print("{}: {} (EMPTY)".format(iso_id, common_name))
    #    isolate_species[Isolate(iso_id, None)] = common_name
    # Query DB for Isolate:CommonName
    cursor = cnx.cursor()
    cursor.execute(isolate_host_query)
    for (iso_id, hostID) in cursor:
        #print("{}->{}".format(iso_id, common_name))
        #if hostID in remap:   # Replace with known mappings
        #    hostID = remap[hostID]
        #print("{}: {}".format(iso_id, common_name))
        isolate_host[Isolate(iso_id, None)] = hostID
    return isolate_host

def clusters_to_csv(min_neighbors, cfg, cnx, file=sys.stdout):
    iso_spec_map = get_isolate_species(cnx)
    iso_host_map = get_isolate_host(cnx)
    isolates = get_isolates(cfg)
    filtered_isolates = filter_isolates(isolates, iso_spec_map)

    clusters = get_clusters(cfg)
    clustered = set()
    unclustered = set()

    clusters = sorted(clusters, key = lambda x : len(x), reverse=True)
    print('{},{},{},{}'.format('cluster','hostID','isoID','species'), file=file)
    num_clust = 1
    for cluster in clusters:
        clust_name = 'clust_{}'.format(num_clust)
        cluster = sorted(list(cluster), key = lambda x : x.__str__())
        for i in cluster:
            clustered.add(i)
            host = iso_host_map[i] if i in iso_host_map else 'unknown?'
            species = iso_spec_map[i] if i in iso_spec_map else 'unknown??'
            host = host.replace(',','_')
            species = species.replace(',','_')
            print('{},{},{},{}'.format(clust_name, host, i, species), file=file)
        num_clust += 1
    for i in filtered_isolates:
        if i not in clustered:
            host = iso_host_map[i] if i in iso_host_map else 'unknown?'
            species = iso_spec_map[i] if i in iso_spec_map else 'unknown??'
            host = host.replace(',','_')
            species = species.replace(',','_')
            print('{},{},{},{}'.format('unclustered', host, i, species), file=file)
    filtered_isolates = set(filtered_isolates)
    for i in isolates:
        if i not in filtered_isolates:
            host = iso_host_map[i] if i in iso_host_map else 'unknown?'
            species = iso_spec_map[i] if i in iso_spec_map else 'unknown??'
            host = host.replace(',','_')
            species = species.replace(',','_')
            print('{},{},{},{}'.format('filtered', host, i, species), file=file)
    return


def main():
    print("Nothing implemented just yet.")
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
