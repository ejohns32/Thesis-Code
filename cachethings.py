import random
import pickle
import os.path

import pyroprinting




def cacheAllIsolates(cfg):
	cacheFileName = "isolatesAll.pickle"
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolatesFromDB(cfg.regions)

	print("loading all isolates from DB...")
	with open(cacheFileName, mode='w+b') as isolatesFile:
		pickle.dump(isolates, isolatesFile)

# def sliceIsolates(isolates, sliceSize):
# 	with open("isolatesAll.pickle", mode='w+b') as isolatesFile:
# 		pickle.dump(isolates, isolatesFile)

# 	for count in range(sliceSize, len(isolates), sliceSize):
# 		slicedIsolates = list(random.sample(isolates, count))
# 		print("{}".format(count))
# 		with open("isolates{}.pickle".format(min(count, len(slicedIsolates))), mode='w+b') as slicedFile:
# 			pickle.dump(slicedIsolates, slicedFile)

def cacheIsolateSubset(cfg):
	cacheFileName = "isolates{}.pickle".format(cfg.isolateSubsetSize)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates()

	print("sampling subset of isolates of size {}...".format(cfg.isolateSubsetSize))
	isolateSubset = list(random.sample(isolates, cfg.isolateSubsetSize))
	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(isolateSubset, cacheFile)



def cacheNeighbors(cfg):
	cacheFileName = "neighbors{}T{}.pickle".format(cfg.isolateSubsetSize, cfg.threshold)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates(cfg.isolateSubsetSize)

	print("calculating isolate neighbors at threshold {} for a subset of size {}...".format(cfg.threshold, cfg.isolateSubsetSize))
	neighbors = fullsearch.computeNeighborsMap(isolates, cfg.radii)

	avgMatches = sum(len(matches) for _, matches in neighbors.items()) / len(neighbors)
	print("avgMatches: {}".format(avgMatches))

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(neighbors, cacheFile)



def cacheRegionsPearsonMap(cfg):
	cacheFileName = "pearsons{}.pickle".format(cfg.isolateSubsetSize)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates(cfg.isolateSubsetSize)

	print("calculating pairwise pearson correlations for a subset of size {}...".format(cfg.threshold, cfg.isolateSubsetSize))
	regionsPearsonMap = getRegionsPearsonMap(isolates, cfg.regions)

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(regionsPearsonMap, cacheFile)



def cacheDBscanClusters(cfg):
	cacheFileName = "dbscan{}.pickle".format(cfg.isolateSubsetSize)
	if os.path.isfile(cacheFileName):
		return

	isolates = pyroprinting.loadIsolates(cfg)
	correctNeighbors = fullsearch.getNeighborsMap(isolates, cfg)

	precomputedSearcher = fullsearch.PrecomputedSearcher(correctNeighbors)
	clusters = dbscan(precomputedSearcher, cfg.radii, cfg.minNeighbors)

	with open(cacheFileName, mode='w+b') as cacheFile:
		pickle.dump(clusters, cacheFile)


if __name__ == '__main__':
	cfg = config.loadConfig()
	# isolates = loadIsolatesFromDB(cfg)
	# sliceIsolates(isolates, 1000)
	cacheAllIsolates(cfg)
	cacheIsolateSubset(cfg)
	cacheNeighbors(cfg)
	cacheRegionsPearsonMap(cfg)
	cacheDBscanClusters(cfg)
