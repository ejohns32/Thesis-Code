import collections
import itertools

import protospatial #TODO: would be cool if not dependant on this

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
			if arePointsWithinRadii(point, queuedPoint, radii):
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