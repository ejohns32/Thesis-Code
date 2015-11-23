#!/bin/bash

for count in 1000 2000 3000 4000 5000 6000
do
    for index in spatial fullsearch
    do
    	for cluster in popDBSCAN dbscan 
    	do
    		/usr/bin/time -l python dbscan.py $count $index $cluster
    	done
    done
done