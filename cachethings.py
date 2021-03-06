import pickle
import os.path

import config
import pyroprinting
import fullsearch
import clusterEval
import dbscan




def cacheAllIsolates(cfg):
	cacheFileName = "isolatesAll.pickle"
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolatesFromDB(cfg)

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(isolates, cacheFile)

# def sliceIsolates(isolates, sliceSize):
# 	with open("isolatesAll.pickle", mode='w+b') as isolatesFile:
# 		pickle.dump(isolates, isolatesFile)

# 	for count in range(sliceSize, len(isolates), sliceSize):
# 		slicedIsolates = list(random.sample(isolates, count))
# 		print("{}".format(count))
# 		with open("isolates{}.pickle".format(min(count, len(slicedIsolates))), mode='w+b') as slicedFile:
# 			pickle.dump(slicedIsolates, slicedFile)

def cacheIsolateSubset(cfg):
	cacheFileName = pyroprinting.getIsolatesCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolatesFromFile("isolatesAll.pickle")

	isolateSubset = pyroprinting.getRandomSubset(isolates, cfg)
	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(isolateSubset, cacheFile)



def cacheNeighbors(cfg):
	cacheFileName = fullsearch.getNeighborsMapCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates(cfg)

	neighbors = fullsearch.computeNeighborsMap(isolates, cfg)

	avgMatches = sum(len(matches) for _, matches in neighbors.items()) / len(neighbors)
	print("avgMatches: {}".format(avgMatches))

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(neighbors, cacheFile)



def cacheReplicatePearsons(cfg):
	cacheFileName = clusterEval.getReplicatePearsonsCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return

	pearsons = clusterEval.loadReplicatePearsonsFromDB(cfg)

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(pearsons, cacheFile)



def cacheDBscanClusters(cfg):
	cacheFileName = dbscan.getDBscanClustersCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates(cfg)

	clusters = dbscan.computeDBscanClusters(isolates, cfg)

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(clusters, cacheFile)


if __name__ == '__main__':
	cfg = config.loadConfig()
	# isolates = loadIsolatesFromDB(cfg)
	# sliceIsolates(isolates, 1000)
	cacheAllIsolates(cfg)
	cacheIsolateSubset(cfg)
	cacheNeighbors(cfg)
	# cacheReplicatePearsons(cfg)
	cacheDBscanClusters(cfg)
