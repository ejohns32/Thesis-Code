import pickle
import cProfile

import config
import pyroprinting


class FullSearcher:
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


class PrecomputedSearcher:
	def __init__(self, neighborMap):
		self.neighborMap = neighborMap

	def __iter__(self):
		return iter(self.neighborMap.keys())

	def __len__(self):
		return len(self.neighborMap)

	def getNeighborsOf(self, queryIsolate, radii):
		return self.neighborMap[queryIsolate]



def precomputeNeighbors(isolates, radii):
	neighbors = {isolate: set() for isolate in isolates}
	for i, iso1 in enumerate(isolates):
		print("{}/{}".format(i, len(isolates)))
		for iso2 in isolates[i+1:]:
			if iso1.isWithinRadiiOf(iso2, radii):
				neighbors[iso1].add(iso2)
				neighbors[iso2].add(iso1)
	return neighbors


def loadNeighborMapFromFile(isolateSubsetSize, threshold):
	with open("neighbors{}T{}.pickle".format(isolateSubsetSize, threshold), mode='r+b') as correctFile:
		neighborMap = pickle.load(correctFile)
	return neighborMap




if __name__ == '__main__':
	cfg = config.loadConfig()
	isolates = pyroprinting.loadIsolatesFromFile(cfg.isolateSubsetSize)

	neighbors = precomputeNeighbors(isolates, cfg.radii)
	# cProfile.run("neighbors = precomputeNeighbors(isolates, cfg.radii)")

	print("size: {}".format(len(neighbors)))
	avgMatches = sum(len(matches) for _, matches in neighbors.items()) / len(neighbors)
	print("avgMatches: {}".format(avgMatches))

	# if len(sys.argv) == 1:
	with open("neighbors{}T{}.pickle".format(cfg.isolateSubsetSize, cfg.threshold), mode='w+b') as correctFile:
		pickle.dump(neighbors, correctFile)
