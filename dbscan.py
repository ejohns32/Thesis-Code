import collections
import itertools
import os.path
import pickle

import config
import pyroprinting
import fullsearch
import spatial

def dbscan(points, radii, minNeighbors):
	clusters = []

	unclassifiedPoints = set(points)
	while len(unclassifiedPoints) > 0:
		point = unclassifiedPoints.pop()
		neighbors = points.getNeighborsOf(point, radii)
		if len(neighbors) >= minNeighbors:
			clusters.append(expandCluster(points, radii, minNeighbors, unclassifiedPoints, neighbors) | {point})

	return clusters


def expandCluster(points, radii, minNeighbors, unclassifiedPoints, reachable):
	cluster = set(reachable)
	unclassifiedPoints -= reachable

	while len(reachable) > 0:
		point = reachable.pop()
		neighbors = points.getNeighborsOf(point, radii)
		if len(neighbors) >= minNeighbors:
			unclassifiedNeighbors = set(neighbor for neighbor in neighbors if neighbor in unclassifiedPoints)
			unclassifiedPoints -= unclassifiedNeighbors # or -= neighbors
			reachable |= unclassifiedNeighbors
			cluster |= neighbors

	return cluster






def myDBscan(points, radii, minNeighbors):
	clusters = []

	unclassifiedPoints = points #TODO: copy spatial tree instead of assuming they are ok with us consuming it
	while len(unclassifiedPoints) > 0:
		point = unclassifiedPoints.pop()
		# print("popped {}".format(point))
		neighbors = unclassifiedPoints.popNeighborsOf(point, radii)
		clusters.extend(investigateNeighbors(unclassifiedPoints, radii, minNeighbors, point, neighbors))

	return clusters


def investigateNeighbors(unseenPoints, radii, minNeighbors, startPoint, startNeighbors):
	clusters = []
	seen = collections.defaultdict(set)

	for neighbor in startNeighbors:
		seen[neighbor].add(startPoint)

	#TODO: probably better for reachable and noise to be ordered. either bf or df?
	# perhaps to the point of combining them and not doing a single cluster at a time.
	# though we want to do more dense areas first right?
	# or do we want to find edges first to keep the queues short?
	if len(startNeighbors) >= minNeighbors:
		reachable = startNeighbors
		noise = set()
		currentCluster = {startPoint} | set(startNeighbors)
		clusters.append(currentCluster)
	else:
		noise = startNeighbors
		reachable = set()
		currentCluster = None

	while len(reachable) + len(noise) > 0:
		#exhaust reachable before investigating noise
		if len(reachable) > 0:
			point = reachable.pop()
		else:
			currentCluster = None
			point = noise.pop()

		unseenNeighbors = unseenPoints.popNeighborsOf(point, radii)

		# need to check all points in queue because they have been removed from unseen and havent ran their own queries yet
		# in the degererate case (where all noise and clusters are connected by radii) this could make whole algorithm n^2
		for queuedPoint in itertools.chain(reachable, noise):
			#TODO: try to shortcut this comparison
			# maybe if unseen = 0?
			# maybe if unseen + seen is > minNeighbors? (for both)
			# maybe can only skip points in reachable
			# maybe already have distances to neighbors so combine for bestcase distance and skip if bestcase doesnt pass threshold?
			# maybe these are stored in spatial trees?
			# maybe depends if point came from noise/reachable
			if point.isWithinRadiiOf(queuedPoint, radii):
				assert(queuedPoint not in seen[point])
				seen[point].add(queuedPoint)
				assert(point not in seen[queuedPoint])
				seen[queuedPoint].add(point)

		for neighbor in unseenNeighbors:
			assert(point not in seen[neighbor])
			seen[neighbor].add(point)

		if len(unseenNeighbors) + len(seen[point]) >= minNeighbors:
			if currentCluster is None:
				currentCluster = {point}
				clusters.append(currentCluster)

			for newlyReachable in seen[point] & noise:
				noise.remove(newlyReachable)
				reachable.add(newlyReachable)

			reachable |= unseenNeighbors
			currentCluster |= unseenNeighbors | seen[point]
		else:
			noise |= unseenNeighbors

	return clusters






def verifyClusters(pointsList, radii, minNeighbors, clusters):
	correct = True

	for i, cluster in enumerate(clusters):
		print("verifying cluster {}/{}...".format(i, len(clusters)))
		correct &= verifyClusterConnectivity(radii, minNeighbors, cluster)
		correct &= verifyClusterMaximality(pointsList, radii, minNeighbors, cluster)

	print("verifying clusters existence...")
	correct = verifyClustersExistence(pointsList, radii, minNeighbors, clusters)

	return correct


def verifyClustersExistence(pointsList, neighborMap, minNeighbors, clusters):
	correct = True
	clusterMap = collections.defaultdict(list)

	for cluster in clusters:
		for point in cluster:
			clusterMap[point].append(cluster)

	for i, point1 in enumerate(pointsList):
		print("{}/{} - {}".format(i+1, len(pointsList), point1))
		for j, point2 in enumerate(pointsList):
			if i != j:
				if isDensityReachable(neighborMap, minNeighbors, point1, point2):
					thisCorrect = len(clusterMap[point1]) == 1
					correct &= thisCorrect
					if not thisCorrect:
						print("\t\tbad {}, {}".format(point1, point2))
					continue

	return correct


def verifyClusterMaximality(pointsList, neighborMap, minNeighbors, cluster):
	correct = True

	for point1 in cluster:
		for point2 in pointsList:
			if point1 != point2:
				if isDensityReachable(neighborMap, minNeighbors, point1, point2):
					thisCorrect = point2 in cluster
					correct &= thisCorrect
					if not thisCorrect:
						print("\t\tbad {}, {}".format(point1, point2))

	return correct


def verifyClusterConnectivity(neighborMap, minNeighbors, cluster):
	correct = True
	clusterPointsList = list(cluster)

	for i, point1 in enumerate(clusterPointsList):
		for point2 in clusterPointsList[i+1:]:
			thisCorrect = areDensityConnected(neighborMap, minNeighbors, point1, point2)
			correct &= thisCorrect
			if not thisCorrect:
				print("\t\tbad {}, {}".format(point1, point2))

	return correct


def isDensityReachable(neighborMap, minNeighbors, fromPoint, point):
	assert(fromPoint != point)

	queue = collections.deque()
	if len(neighborMap[fromPoint]) >= minNeighbors:
		queue.append(fromPoint)

	seen = {fromPoint}

	#breadth first search all reachable points
	while len(queue) > 0:
		neighbors = neighborMap[queue.popleft()]

		if len(neighbors) >= minNeighbors:
			if point in neighbors:
				return True

			queue.extend({neighbor for neighbor in neighbors if neighbor not in seen})
			seen |= neighbors

	return False


def areDensityConnected(neighborMap, minNeighbors, point1, point2):
	assert(point1 != point2)

	if point2 in neighborMap[point1]:
		if len(neighborMap[point1]) >= minNeighbors or len(neighborMap[point2]) >= minNeighbors:
			return True
	#look for a point from which both point1 and point2 are reachable
	for point in neighborMap[point1]:
		if point is not point2 and isDensityReachable(neighborMap, minNeighbors, point, point1) and isDensityReachable(neighborMap, minNeighbors, point, point2):
			return True

	return False


def getMultipleClusterPoints(points, clusters):
	clusterMap = collections.defaultdict(list)
	multipleClusterPoints = []

	for cluster in clusters:
		for point in cluster:
			clusterMap[point].append(cluster)

	for point in points:
		if len(clusterMap[point]) > 1:
			multipleClusterPoints.append(point)
			print(point)

	return multipleClusterPoints

def printClusters(clusters):
	for i, cluster in enumerate(clusters):
		print("{}/{}: {}".foramt(i+1, len(clusters), cluster))
		print()



def computeDBscanClusters(isolates, cfg):
	print("clustering for a subset of size {}...".format(cfg.isolateSubsetSize))
	correctNeighbors = fullsearch.getNeighborsMap(isolates, cfg)
	precomputedSearcher = fullsearch.PrecomputedSearcher(correctNeighbors)
	return dbscan(precomputedSearcher, cfg.radii, cfg.minNeighbors)

def getDBscanClustersCacheFileName(cfg):
	return "dbscan{}.pickle".format(cfg.isolateSubsetSize)

def loadDBscanClustersFromFile(cacheFileName):
	with open(cacheFileName, mode='r+b') as cacheFile:
		dbscanClusters = pickle.load(cacheFile)
	return dbscanClusters

def getDBscanClusters(isolates, cfg):
	cacheFileName = getDBscanClustersCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return loadDBscanClustersFromFile(cacheFileName)
	else:
		return computeDBscanClusters(isolates, cfg)



if __name__ == '__main__':
	cfg = config.loadConfig()
	isolates = pyroprinting.loadIsolates(cfg)
	correctNeighbors = fullsearch.getNeighborsMap(isolates, cfg)

	tree = spatial.Tree(isolates, cfg)
	# fullSearcher = fullsearch.FullSearcher(isolates)
	precomputedSearcher = fullsearch.PrecomputedSearcher(correctNeighbors)

	# spatialGetClusters = dbscan(tree, cfg.radii, cfg.minNeighbors)
	spatialPopClusters = myDBscan(tree, cfg.radii, cfg.minNeighbors)
	# fullGetClusters = dbscan(fullSearcher, cfg.radii, cfg.minNeighbors)
	# fullPopClusters = myDBscan(fullSearcher, cfg.radii, cfg.minNeighbors)
	preFullClusters = dbscan(precomputedSearcher, cfg.radii, cfg.minNeighbors)

	# spatialGetClusters = {frozenset(cluster) for cluster in spatialGetClusters}
	spatialPopClusters = {frozenset(cluster) for cluster in spatialPopClusters}
	# fullGetClusters = {frozenset(cluster) for cluster in fullGetClusters}
	# fullPopClusters = {frozenset(cluster) for cluster in fullPopClusters}
	preFullClusters = {frozenset(cluster) for cluster in preFullClusters}

	# printClusters(preFullClusters)
	# assert len(getMultipleClusterPoints(isolates, preFullClusters)) == 0
	assert spatialPopClusters == preFullClusters
	# assert spatialGetClusters == spatialPopClusters
	# assert fullGetClusters == fullPopClusters
	# assert verifyClusters(isolates, correctNeighbors, cfg.minNeighbors, clusters)
