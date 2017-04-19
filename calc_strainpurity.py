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

import strainpurity
from strainpurity import ClusterPurity
from strainpurity import StrainPurities
from strainpurity import StrainPurityJSONEncoder
from strainpurity import clusters_to_csv

import graphstrainpurity
from graphstrainpurity import table_purities
from graphstrainpurity import table_species
from graphstrainpurity import plot_hist_species
from graphstrainpurity import plot_cluster_scatter_unfilled
from graphstrainpurity import plot_cluster_scatter_filled_below
from graphstrainpurity import plot_cluster_scatter_filled_between
from graphstrainpurity import plot_cluster_scatter_filled_stack
from graphstrainpurity import plot_lower_upper_scatter_filled_stack
from graphstrainpurity import plot_lower_upper_unpure_scatter_filled_stack
from graphstrainpurity import plot_scatter
from graphstrainpurity import plot_two_scatter
from graphstrainpurity import plot_dist_clust
from graphstrainpurity import plot_dist_size
from graphstrainpurity import plot_dist_data

"""
    Does the thing.
"""
def do_the_thing(min_neighbors):
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
    #isolates = get_isolates(cfg)
    # Query Database for Isolates -> CommonName mappings
    #iso_map = get_isolate_species(db_cnx)
    # Filter isolates accordingly
    #isolates = filter_isolates(isolates, iso_map, filter)

    clustering_purities = []
    all_cluster_purities = []
    all_cluster_sizes = []
    all_cluster_unique_spec = []
    all_datapoint_purities = []

    all_num_clustered = []
    all_num_unclustered = []
    all_num_multiply = []
    all_num_isolates = []
    all_num_minor_isolates = []
    all_num_major_isolates = []
    all_num_major_pure_isolates = []
    all_num_major_unpure_isolates = []
    all_cluster_counts = []
    species_counts = None
    # Run Each Min Neighbor Clustering
    for n in min_neighbors:
        print('Clustering for min_neighbors of {}'.format(n))
        cfg.minNeighbors = n

        # Build Filename Strings
        json_filename = "purities{}_{}_{}.json".format(cfg.isolateSubsetSize, "_".join(
            str(region.clusterThreshold) for region in cfg.regions), cfg.minNeighbors)
        pickle_filename = "purities{}_{}_{}.pickle".format(cfg.isolateSubsetSize, "_".join(
            str(region.clusterThreshold) for region in cfg.regions), cfg.minNeighbors)

        # Run DBSCAN
        #cache_dbscan_clusters(cfg, isolates)
        # Get Clusters
        #clusters = get_clusters(cfg)
        # Calc Purities
        purities = StrainPurities(cfg, db_cnx, cfg_json)
        print('Writing JSON...')
        with open(json_filename, mode='w') as file:
            json.dump(purities, file, cls=StrainPurityJSONEncoder, sort_keys=True)
        print('Pickling...')
        with open(pickle_filename, mode='w+b') as file:
            pickle.dump(purities, file)
        clustering_purities.append(purities.clustering_purity);
        cluster_purities = []
        cluster_sizes = []
        cluster_unique_spec = []
        cluster_counts = []
        datapoint_purities = []
        num_minor_isolates = 0
        num_major_isolates = 0
        num_major_pure_isolates = 0
        num_major_unpure_isolates = 0
        for c in purities.purities:
            if c.cluster_size <= 1:
                print("Skipping: {}".format(c))
                continue
            cluster_purities.append(c.cluster_purity)
            cluster_sizes.append(c.cluster_size)
            cluster_unique_spec.append(len(c.counts))
            datapoint_purities.extend([c.cluster_purity for i in range(c.cluster_size)])
            cluster_counts.append(c.counts)
            for spec, count in c.counts.items():
                if spec is c.most_plural['label']:
                    num_major_isolates += count
                    if c.cluster_purity < 1.0:
                        num_major_unpure_isolates += count
                    else:
                        num_major_pure_isolates += count
                else:
                    num_minor_isolates += count
        all_cluster_purities.append(cluster_purities)
        all_cluster_sizes.append(cluster_sizes)
        all_cluster_unique_spec.append(cluster_unique_spec)
        all_datapoint_purities.append(datapoint_purities)

        all_num_clustered.append(purities.num_clustered)
        all_num_unclustered.append(purities.num_unclustered)
        all_num_multiply.append(purities.num_multiply_clustered)
        all_num_isolates.append(purities.num_isolates)
        all_num_minor_isolates.append(num_minor_isolates)
        all_num_major_isolates.append(num_major_isolates)
        all_num_major_pure_isolates.append(num_major_pure_isolates)
        all_num_major_unpure_isolates.append(num_major_unpure_isolates)
        species_counts = purities.species_counts
        all_cluster_counts.append(cluster_counts)
        with open('neigh_all_clusters_{}.csv'.format(n), mode='w') as file:
            clusters_to_csv(n, cfg, db_cnx, file)
    data_clustering_purities = zip(min_neighbors, clustering_purities)
    print("Making Species Table....")
    #for n,c in zip(min_neighbors, all_cluster_counts):
    #    table_purities(n,c)
    #return
    #table_species(species_counts)
    print('# SCATTER -- Minor vs Major Impure')
    plot_lower_upper_unpure_scatter_filled_stack(
            min_neighbors,
            all_num_minor_isolates,
            all_num_major_unpure_isolates,
            all_num_major_pure_isolates,
            all_num_isolates,
            'Number of Minor and Major Pure and Impure Isolates',
            'neigh_minor_major_impure_filled_stack.pdf'
            )
    return -1
    print('SCATTER -- Overall Clustering Purity AND Accuracy')
    overall_accuracy = list(map(lambda x: x / max(all_num_isolates), all_num_major_isolates))
    plot_two_scatter(min_neighbors, clustering_purities, overall_accuracy,
            "Overall Accuracy Against Overall Clustering Purity",
            'Minimum Neighbors', 'Overall Accuracy / Clustering Purity', 'neigh_clust_accuracy.pdf',
            y_low = 0, y_high = 1.05, y_ticks = .05,
            y1_label = 'Overall Clustering Purity', y2_label = 'Overall Accuracy')
    #return -1
    print('SCATTER -- Overall Clustering Purity')
    plot_scatter(min_neighbors, clustering_purities,
            "Overall Clustering Purity as Minimum Neighbors Increases",
            'Minimum Neighbors', 'Overall Clustering Purity', 'neigh_clust.pdf',y_low = 0, y_high = 1.05, y_ticks = .05)
    #return -1
    print('# HISTOGRAM -- Num Datapoints of Cluster Purity')
    for n,c in zip(min_neighbors, all_datapoint_purities):
        # print('Cluster Purities for each Datapoint for n = {}:{}'.format(n, c))
        plot_dist_data(c, 'Number of Datapoints in Cluster Purity for Minimum Neighbors of {}'.format(n), 'neigh_dist_data_{}.pdf'.format(n))
    print('# HISTOGRAM -- Num Clusters of Cluster Size')
    for n,c in zip(min_neighbors, all_cluster_sizes):
        # print('Cluster Sizes for n = {}:{}'.format(n, c))
        plot_dist_size(c, 'Cluster Sizes for Minimum Neighbors of {}'.format(n), 'neigh_size_{}.pdf'.format(n))
    print('Done.')
    print('# SCATTER -- Minor vs Major')
    plot_lower_upper_scatter_filled_stack(
            min_neighbors,
            all_num_minor_isolates,
            all_num_major_isolates,
            all_num_isolates,
            'Number of Major and Minor Isolates',
            'neigh_minor_major_filled_stack.pdf'
            )
    print('# SCATTER -- Clustered vs Unclustered')
    plot_cluster_scatter_unfilled(min_neighbors, all_num_clustered,
            all_num_unclustered, all_num_multiply, all_num_isolates,
            "Number of Clustered and Unclustered Isolates as Minimum Neighbors Increases",
            'neigh_clust_unclust_unfilled.pdf')
    plot_cluster_scatter_filled_stack(min_neighbors, all_num_clustered,
            all_num_unclustered, all_num_multiply, all_num_isolates,
            "Number of Clustered and Unclustered Isolates as Minimum Neighbors Increases",
            'neigh_clust_unclust_filled_stack.pdf')
    plot_cluster_scatter_filled_between(min_neighbors, all_num_clustered,
            all_num_unclustered, all_num_multiply, all_num_isolates,
            "Number of Clustered and Unclustered Isolates as Minimum Neighbors Increases",
            'neigh_clust_unclust_filled_between.pdf')
    #plot_cluster_scatter_filled_below(min_neighbors, all_num_clustered,
    #        all_num_unclustered, all_num_multiply, all_num_isolates,
    #        "Number of Clustered and unclustered Isolates as Minimum Neighbors Increases",
    #        'neigh_clust_unclust_filled_below.pdf')
    print('# HISTOGRAM -- Species in CPLOP')
    assert(species_counts)
    return
    plot_hist_species(species_counts)
    print('# SCATTER -- Sizes vs Purities')
    for n,s,p in zip(min_neighbors, all_cluster_sizes, all_cluster_purities):
        max_size = 375
        marker_size = [2*max_size*( (x) / max_size) for x in s]
        plot_scatter(s, p,
                "Clustering Purity as it Relates to Cluster Size for Neighbors of {}".format(n),
                'Cluster Size (Cluster Size Affects Dot Size)', 'Individual Cluster Purity', 'neigh_size_purity_{}.pdf'.format(n),
                x_low=0, x_high=375,
                y_low=0, y_high=1.05, y_ticks=.05,
                marker_size = marker_size
                )
    print('# SCATTER -- Unique vs Purities')
    for n,s,p,u in zip(min_neighbors, all_cluster_sizes, all_cluster_purities, all_cluster_unique_spec):
        max_spec = 16
        max_size = 375
        marker_size = [2*max_size*( (x) / max_spec) for x in u]
        plot_scatter(u, p,
                "Clustering Purity as it Relates to Unique Species Neighbors of {}".format(n),
                'Unique Species (Num Unique Affects Dot Size)', 'Individual Cluster Purity', 'neigh_unique_purity_{}.pdf'.format(n),
                x_low=0, x_high=max_spec,
                y_low=0, y_high=1.05, y_ticks=.05,
                marker_size=marker_size
                )
        plot_scatter(p, u,
                "Clustering Purity as it Relates to Unique Species Neighbors of {}".format(n),
                'Individual Cluster Purity', 'Unique Species (Num Unique Affects Dot Size)', 'neigh_purity_unique_{}.pdf'.format(n),
                y_low=0, y_high=max_spec,
                x_low=0, x_high=1.05, x_ticks=.05,
                marker_size=marker_size
                )
        plot_scatter(s, p,
                "Clustering Purity as it Relates to Cluster Size for Neighbors of {}".format(n),
                'Cluster Size (Num Unique Affects Dot Size)', 'Individual Cluster Purity', 'neigh_size_purity_unique_{}.pdf'.format(n),
                x_low=0, x_high=max_size,
                y_low=0, y_high=1.05, y_ticks=.05,
                marker_size = marker_size
                )
        marker_size = [2*max_size*( (x) / max_size) for x in s]
        plot_scatter(u, p,
                "Clustering Purity as it Relates to Unique Species Neighbors of {}".format(n),
                'Unique Species (Cluster Size Affects Dot Size)', 'Individual Cluster Purity', 'neigh_unique_purity_{}_size.pdf'.format(n),
                x_low=0, x_high=max_spec,
                y_low=0, y_high=1.05, y_ticks=.05,
                marker_size=marker_size
                )
        plot_scatter(p, u,
                "Unique Species as it Relates to Clustering Purity Neighbors of {}".format(n),
                'Individual Cluster Purity', 'Unique Species (Cluster Size Affects Dot Size)', 'neigh_purity_unique_{}_size.pdf'.format(n),
                y_low=0, y_high=max_spec,
                x_low=0, x_high=1.05, x_ticks=.05,
                marker_size=marker_size
                )
    print('# HISTOGRAM -- Num Clusters of Cluster Purity', )
    for n,c in zip(min_neighbors, all_cluster_purities):
        # print('Cluster Purities for n = {}:{}'.format(n, c))
        plot_dist_clust(c, 'Cluster Purities for Minimum Neighbors of {}'.format(n), 'neigh_dist_{}.pdf'.format(n))
    return None

def main():
    # CONFIG AND DB LOADING
    do_the_thing([1, 2, 3, 4, 5, 6, 7, 8])
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
