import collections
import itertools

def dbscan(points, radii, minNeighbors):
	clusters = []

	unclassifiedPoints = set(points)
	while len(unclassifiedPoints) > 0:
		point = unclassifiedPoints.pop()
		neighbors = points.getPointsInSphere(point, radii) - set(point)
		if len(neighbors) >= minNeighbors:
			clusters.append(expandCluster(points, radii, minNeighbors, unclassifiedPoints, neighbors) | set(point))

	return clusters

def expandCluster(points, radii, minNeighbors, unclassifiedPoints, reachable):
	cluster = set(reachable)
	unclassifiedPoints -= reachable

	while len(reachable) > 0:
		point = reachable.pop()
		neighbors = points.getPointsInSphere(point, radii) - set(point)
		if len(neighbors) >= minNeighbors:
			unclassifiedNeighbors = neighbor for neighbor in neighbors if neighbor in unclassifiedPoints
			unclassifiedPoints -= unclassifiedNeighbors # or -= neighbors
			reachable |= unclassifiedNeighbors
			cluster |= neighbors

	return cluster



def myDBscan(points, radii, minNeighbors):
	clusters = []

	unclassifiedPoints = points.copy() #spatial tree
	while len(unclassifiedPoints) > 0:
		point = unclassifiedPoints.pop()
		neighbors = unclassifiedPoints.popPointsInSphere(point, radii)
		clusters.extend(investigateNeighbors(unclassifiedPoints, radii, minNeighbors, point, neighbors))

	return clusters

def investigateNeighbors(unseenPoints, radii, minNeighbors, startPoint, startNeighbors):
	clusters = []
	seen = collections.defaultdict(set)

	for neighbor in startNeighbors:
		seen[neighbor] |= startPoint

	#TODO(5): probably better for reachable and noise to be ordered. either bf or df?
	#perhaps to the point of combining them and not doing a single cluster at a time. though i want to do more dense areas first right? or do I want to find edges first to keep the queues short?
	if len(startNeighbors) >= minNeighbors:
		reachable = startNeighbors
		noise = set()
		currentCluster = set(startPoint) | set(startNeighbors)
		clusters.append(currentCluster)
	else:
		noise = startNeighbors
		reachable = set()
		currentCluster = None

	while len(reachable) + len(noise) > 0:
		if len(reachable) > 0: #exhaust reachable before investigating noise
			point = reachable.pop()
		else: #if len(noise) > 0
			currentCluster = None
			point = noise.pop()

		unseenNeighbors = unseenPoints.popPointsInSphere(point, radii)
		#TODO(1): maybe don't check noise? not compatible with (3)
		for queuedPoint in itertools.chain(reachable, noise):
			#TODO(4): maybe have a way of shortcutting this comparison
			#maybe if unseen = 0?
			#maybe already have distances to neighbors so combine?
			#maybe these are stored in spatial trees?
			#maybe depends if point came from noise/reachable
			if point.isWithinRadiiOf(queuedPoint, radii):
				assert(queuedPoint not in seen[point])
				seen[point] |= queuedPoint
				assert(point not in seen[queuedPoint])
				seen[queuedPoint] |= point

		for neighbor in unseenNeighbors:
			assert(point not in seen[neighbor])
			seen[neighbor] |= point

		if len(unseenNeighbors) + len(seen[point]) >= minNeighbors:
			if currentCluster is None:
				#TODO(2): check if point is already in a cluster. not necessary with (3)
				currentCluster = set(point)
				clusters.append(currentCluster)

			#TODO(3): move seen[point] & noise from noise to reachable?
			reachable |= unseenNeighbors
			currentCluster |= unseenNeighbors | seen[point]
		else:
			noise |= unseenNeighbors

	return clusters



def verifyClusters(points, radii, minNeighbors, clusters):
	correct = verifyClustersMaximality(points, radii, minNeighbors, clusters)

	for cluster in clusters:
		corrent &= verifyClusterConnectivity(points, radii, minNeighbors, cluster)

	return correct

def verifyClustersMaximality(points, neighborMap, minNeighbors, clusters):
	correct = True
	pointsList = list(points)
	clusterMap = collections.defaultdict(set)

	for cluster in clusters:
		for point in cluster:
			clusterMap[point] |= cluster

	for i, point1 in enumerate(pointsList):
		for j, point2 in enumerate(pointsList):
			if i != j:
				for cluster in clusterMap[point1]:
					if isDensityReachable(neighborMap, minNeighbors, point1, point2):
						correct &= point2 in cluster

	return correct

def verifyClusterConnectivity(points, neighborMap, minNeighbors, cluster):
	correct = True
	clusterPointsList = list(cluster)

	for i, point1 in enumerate(clusterPointsList):
		for point2 in clusterPointsList[i:]:
			correct &= areDensityConnected(neighborMap, minNeighbors, point1, point2)

	return correct


def isDensityReachable(neighborMap, minNeighbors, fromPoint, point):
	assert(fromPoint != point)

	queue = collections.deque()
	if len(neighborMap[fromPoint]) >= minNeighbors:
		queue.append(fromPoint)

	#breadth first search all reachable points
	while len(queue) > 0:
		neighbors = neighborMap[queue.popleft()]

		if len(neighbors) >= minNeighbors:
			queue.extend(neighbors)

			if point in neighbors:
				return True

	return False

def areDensityConnected(neighborMap, minNeighbors, point1, point2):
	assert(point1 != point2)

	#look for a point from which both point1 and point2 are reachable
	for point in neighborMap[point1]:
		if isDensityReachable(neighborMap, minNeighbors, point, point1) and isDensityReachable(neighborMap, minNeighbors, point, point2)
			return True

	return False
