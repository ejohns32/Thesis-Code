import pickle
import cProfile
import os.path

import config
import pyroprinting


class FullSearchIndex:
	def __init__(self, isolates):
		self.isolates = set(isolates)
	def __iter__(self):
		return iter(self.isolates)

	def __len__(self):
		return len(self.isolates)

	def pop(self):
		return self.isolates.pop()


	def getNeighborsOf(self, queryIsolate, radii):
		result = set()
		for iso in self.isolates:
			if iso.isWithinRadiiOf(queryIsolate, radii):
				result.add(iso)

		return result - {queryIsolate}

	def popNeighborsOf(self, queryIsolate, radii):
		result = set()
		for iso in self.isolates:
			if iso.isWithinRadiiOf(queryIsolate, radii):
				result.add(iso)

		self.isolates -= result

		return result - {queryIsolate}


class PrecomputedIndex:
	def __init__(self, neighborMap):
		self.neighborMap = neighborMap

	def __iter__(self):
		return iter(self.neighborMap.keys())

	def __len__(self):
		return len(self.neighborMap)

	def getNeighborsOf(self, queryIsolate, radii):
		return self.neighborMap[queryIsolate]


# TODO: calculate average pairwise dist
def computeNeighborsMap(isolates, cfg):
	print("calculating isolate neighbors for a subset of size {}...".format(cfg.isolateSubsetSize))

	neighbors = {isolate: set() for isolate in isolates}
	for i, iso1 in enumerate(isolates):
		print("{}/{}".format(i+1, len(isolates)))
		for iso2 in isolates[i+1:]:
			if iso1.isWithinRadiiOf(iso2, cfg.radii):
				neighbors[iso1].add(iso2)
				neighbors[iso2].add(iso1)

	return neighbors

def getNeighborsMapCacheFileName(cfg):
	return "neighbors{}_{}.pickle".format(cfg.isolateSubsetSize, "_".join(str(region.clusterThreshold) for region in cfg.regions))
	# return "neighbors{}T{}.pickle".format(cfg.isolateSubsetSize, '_'.join(str(region.pSimThresholdAlpha) for region in  cfg.regions))

def loadNeighborsMapFromFile(cacheFileName):
	with open(cacheFileName, mode='r+b') as cacheFile:
		neighborMap = pickle.load(cacheFile)
	return neighborMap


def getNeighborsMap(isolates, cfg):
	cacheFileName = getNeighborsMapCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return loadNeighborsMapFromFile(cacheFileName)
	else:
		return computeNeighborsMap(isolates, cfg)

