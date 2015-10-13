import pickle
import numpy
import itertools
import random
import math
import cProfile
from operator import itemgetter

class Box:
	def __init__(self, lowerCorner, upperCorner):
		self.lowerCorner = lowerCorner
		self.upperCorner = upperCorner
	def combine(box1, box2):
		return Box(numpy.minimum(box1.lowerCorner, box2.lowerCorner), numpy.maximum(box1.upperCorner, box2.upperCorner))
	def combineAll(boxes):
		combined = next(boxes)
		for box in boxes:
			combined = Box.combine(combined, box)
		return combined

	def bound(points, dispCount):
		dimMin = numpy.full(dispCount, numpy.finfo(numpy.float64).max)
		dimMax = numpy.full(dispCount, numpy.finfo(numpy.float64).min)
		for point in points:
			dimMin = numpy.minimum(dimMin, point)
			dimMax = numpy.maximum(dimMax, point)
		return Box(dimMin, dimMax)


def ballContainsBox(point, radius, box):
	return numpy.linalg.norm(numpy.maximum(numpy.absolute(point-box.lowerCorner), numpy.absolute(point-box.upperCorner))) < radius

def pointInBox(point, box):
	return numpy.all(numpy.greater_equal(point, box.lowerCorner)) and numpy.all(numpy.less_equal(point, box.upperCorner))
def pointDistFromBox(point, box, stats):
	# dist = min(numpy.linalg.norm(point-numpy.maximum(point, box.lowerCorner)), numpy.linalg.norm(point-numpy.minimum(point, box.upperCorner)))
	dist = numpy.linalg.norm(point-numpy.maximum(point, box.lowerCorner)), point-numpy.minimum(point, box.upperCorner)
	stats.levelStats[stats.level].totalDists += 1
	stats.levelStats[stats.level].avgFullDist += dist
	return dist
def ballIntersectsBox(point, radius, box, stats):
	return pointDistFromBox(point, box, stats) <= radius

def multiDimSplitVal(point, dims, coeffs):
	return sum(point[dims[i]] * coeffs[i] for i in range(len(dims)))

def multiDimSplitValDiff(val1, val2, coeffNorm):
	return max(0, (val1 - val2)/math.sqrt(coeffNorm))

class LeafNode:
	def __init__(self, pyroIDs):
		self.pyroIDs = pyroIDs

class InnerNode:
	def __init__(self, children):
		self.children = children



class MultiRegionSplitNode:
	def newLeaf(isoIDs, regions, zScores):
		return MultiRegionSplitLeafNode(isoIDs)
	def newInner(children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues):
		return MultiRegionSplitInnerNode(children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues)
	def isContainedBy(self, center, radius):
		return False
	def doesContainPoint(self, point):
		return False

class MultiRegionSplitInnerNode(MultiRegionSplitNode, InnerNode):
	def __init__(self, children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues):
		InnerNode.__init__(self, children)
		self.parent = None
		for childIndex, child in enumerate(self.children):
			child.parent = self
			child.childIndex = childIndex
		self.regions = regions
		self.regionSplitDims = regionSplitDims
		self.regionDimCoeffs = regionDimCoeffs
		self.regionCoeffNorms = [sum(regionDimCoeffs[i][j]**2 for j in range(len(regionDimCoeffs[i])))for i in range(len(regions))]
		self.regionSplitValues = [[-666.666] + splitValues + [666.666] for splitValues in regionSplitValues]
	def getRegionSplitDist(self, childIndex, region, center):
		dimVal = multiDimSplitVal(center, self.regionSplitDims[region], self.regionDimCoeffs[region])
		lowerDist = multiDimSplitValDiff(self.regionSplitValues[region][childIndex], dimVal, self.regionCoeffNorms[region])
		upperDist = multiDimSplitValDiff(dimVal, self.regionSplitValues[region][childIndex+1], self.regionCoeffNorms[region])
		# print("dist({}, {}) = {}".format(dimVal, self.splitValues[childIndex], lowerDist))
		# print("dist({}, {}) = {}".format(dimVal, self.splitValues[childIndex+1], upperDist))
		# at least one should always be 0, both if inside
		return max(lowerDist, upperDist)
	# def getChildSplitDist(self, childIndex, center):
	# 	return max(self.getRegionSplitDist(childIndex >> i % 2, region, center) for i, region in enumerate(self.regions))

	def getIntersectingChildren(self, center, radius, stats):
		anscDist = {region: 0.0 for region in self.regions}
		if self.parent:
			anscDist = self.parent.memoDists[self.childIndex]

		result = []
		self.memoDists = [{region: 0.0 for region in self.regions} for _ in range(len(self.children))]
		for childIndex, child in enumerate(self.children):
			splitDists = []
			regionDists = []
			for i, region in enumerate(self.regions):
				splitDist = self.getRegionSplitDist((childIndex >> i) % 2, i, center[region])
				dist = math.sqrt(anscDist[region]**2 + splitDist**2) if splitDist > 0 else anscDist[region]
				self.memoDists[childIndex][region] = dist
				splitDists.append(splitDist)
				regionDists.append(dist)

			stats.levelStats[stats.level].totalDists += 1
			stats.levelStats[stats.level].avgSplitDist += max(splitDists)
			stats.levelStats[stats.level].avgFullDist += max(regionDists)

			if max(regionDists) <= radius:
				result.append(child)

		return result
	def getDistSortedChildren(self, point, stats):
		anscDist = {region: 0.0 for region in self.regions}
		if self.parent:
			anscDist = self.parent.memoDists[self.childIndex]

		result = []
		self.memoDists = [{region: 0.0 for region in self.regions} for _ in range(len(self.children))]
		for childIndex, child in enumerate(self.children):
			splitDists = []
			regionDists = []
			for i, region in enumerate(self.regions):
				splitDist = self.getRegionSplitDist((childIndex >> i) % 2, i, point[region])
				dist = math.sqrt(anscDist[region]**2 + splitDist**2) if splitDist > 0 else anscDist[region]
				self.memoDists[childIndex][region] = dist
				splitDists.append(splitDist)
				regionDists.append(dist)

			stats.levelStats[stats.level].totalDists += 1
			stats.levelStats[stats.level].avgSplitDist += max(splitDists)
			stats.levelStats[stats.level].avgFullDist += max(regionDists)

			result.append((max(regionDists), child))

		return sorted(result, key=itemgetter(0))

class MultiRegionSplitLeafNode(MultiRegionSplitNode, LeafNode):
	def __init__(self, pyroIDs):
		LeafNode.__init__(self, pyroIDs)


class MultiRegionBoundedNode:
	def newLeaf(isoIDs, regions, zScores):
		return MultiRegionBoundedLeafNode(isoIDs, regions, zScores)
	def newInner(children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues):
		return MultiRegionBoundedInnerNode(children, regions)

	def __init__(self, regions, regionBoundingBoxes):
		self.regions = regions
		self.regionBoundingBoxes = regionBoundingBoxes
	def isContainedBy(self, center, radius):
		return all(ballContainsBox(center[region], radius, self.regionBoundingBoxes[i]) for i, region in enumerate(self.regions))
	def doesContainPoint(self, point):
		return all(pointInBox(point[region], self.regionBoundingBoxes[i]) for i, region in enumerate(self.regions))

class MultiRegionBoundedInnerNode(MultiRegionBoundedNode, InnerNode):
	def __init__(self, children, regions):
		MultiRegionBoundedNode.__init__(self, regions, [Box.combineAll(child.regionBoundingBoxes[i] for child in children) for i, region in enumerate(regions)])
		InnerNode.__init__(self, children)
	def getIntersectingChildren(self, center, radius, stats):
		return (child for child in self.children if all(ballIntersectsBox(center[region], radius, child.regionBoundingBoxes[i], stats) for i, region in enumerate(self.regions)))
	def getDistSortedChildren(self, point, stats):
		return sorted(((max(pointDistFromBox(point[region], child.regionBoundingBoxes[i], stats) for i, region in enumerate(self.regions)), child) for child in self.children), key=itemgetter(0))

class MultiRegionBoundedLeafNode(MultiRegionBoundedNode, LeafNode):
	def __init__(self, pyroIDs, regions, zScores):
		MultiRegionBoundedNode.__init__(self, [region for region, _ in regions], [Box.bound((zScores[pyroID][region] for pyroID in pyroIDs), dispCount) for region, dispCount in regions])
		LeafNode.__init__(self, pyroIDs)


class MultiRegionCombinedNode(MultiRegionBoundedNode, MultiRegionSplitNode):
	def newLeaf(isoIDs, regions, zScores):
		return MultiRegionCombinedLeafNode(isoIDs, regions, zScores)
	def newInner(children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues):
		return MultiRegionCombinedInnerNode(children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues)

	def isContainedBy(self, center, radius):
		return MultiRegionBoundedNode.isContainedBy(self, center, radius) or MultiRegionSplitNode.isContainedBy(self, center, radius)
	def doesContainPoint(self, point):
		return MultiRegionBoundedNode.doesContainPoint(self, point) # and MultiRegionSplitNode.doesContainPoint(self, point)

class MultiRegionCombinedInnerNode(MultiRegionCombinedNode, MultiRegionBoundedInnerNode, MultiRegionSplitInnerNode):
	def __init__(self, children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues):
		MultiRegionBoundedInnerNode.__init__(self, children, regions)
		MultiRegionSplitInnerNode.__init__(self, children, regions, regionSplitDims, regionDimCoeffs, regionSplitValues)
	def getIntersectingChildren(self, center, radius, stats):
		return set(MultiRegionBoundedInnerNode.getIntersectingChildren(self, center, radius, stats)) & set(MultiRegionSplitInnerNode.getIntersectingChildren(self, center, radius, stats))
	def getDistSortedChildren(self, point, stats):
		boundedDists = {child: dist for dist, child in MultiRegionBoundedInnerNode.getDistSortedChildren(self, point, stats)}
		splitDists = {child: dist for dist, child in MultiRegionSplitInnerNode.getDistSortedChildren(self, point, stats)}

		for child in self.children:
			stats.levelStats[stats.level].totalDists -= 1
			stats.levelStats[stats.level].avgFullDist -= min(boundedDists[child], splitDists[child])

		return sorted(((max(boundedDists[child], splitDists[child]), child) for child in self.children), key=itemgetter(0))

class MultiRegionCombinedLeafNode(MultiRegionCombinedNode, LeafNode):
	def __init__(self, pyroIDs, regions, zScores):
		MultiRegionBoundedLeafNode.__init__(self, pyroIDs, regions, zScores)
		MultiRegionSplitLeafNode.__init__(self, pyroIDs)



def splitMultiRegion(zScores, isoIDs, regions, nodeType, pointsPerLeaf):
	if len(isoIDs) <= pointsPerLeaf:
		return nodeType.newLeaf(set(isoIDs), regions, zScores)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(regions))]
		# print(truthTable)
		childrenIDs = {truth: [] for truth in truthTable}

		regionSplitDims = []
		regionDimCoeffs = []
		regionSplitValue = []
		for region, dispCount in regions:
			dimSum = sum((zScores[isoID][region] for isoID in isoIDs), numpy.zeros(dispCount))
			dimAvg = dimSum / len(isoIDs)
			dimDevSum = sum(((dimAvg-zScores[isoID][region]) ** 2 for isoID in isoIDs), numpy.zeros(dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isoIDs))
			
			splitDims = [max(range(dispCount), key=lambda dim: dimStdDev[dim])]
			dimCoeffs = [1.0]
			isoIDs.sort(key=lambda isoID: multiDimSplitVal(zScores[isoID][region], splitDims, dimCoeffs))
			splitValue = [multiDimSplitVal(zScores[isoIDs[len(isoIDs)//2]][region], splitDims, dimCoeffs)]

			regionSplitDims.append(splitDims)
			regionDimCoeffs.append(dimCoeffs)
			regionSplitValue.append(splitValue)

		for i, isoID in enumerate(isoIDs):
			truth = tuple(bool(multiDimSplitVal(zScores[isoID][region], regionSplitDims[j], regionDimCoeffs[j]) < regionSplitValue[j]) for j, (region, _) in enumerate(regions))
			childrenIDs[truth].append(isoID)

		children = [splitMultiRegion(zScores, childrenIDs[truth], regions, nodeType, pointsPerLeaf) for truth in truthTable]

		return nodeType.newInner(children, [region for region, _ in regions], regionSplitDims, regionDimCoeffs, regionSplitValue)

def splitMultiRegionCorrelatedDims(zScores, isoIDs, regions, nodeType, pointsPerLeaf, dimsPerSplit):
	if len(isoIDs) <= pointsPerLeaf:
		return nodeType.newLeaf(set(isoIDs), regions, zScores)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(regions))]
		# print(truthTable)
		childrenIDs = {truth: [] for truth in truthTable}

		regionSplitDims = []
		regionDimCoeffs = []
		regionSplitValue = []
		for region, dispCount in regions:
			dimSum = sum((zScores[isoID][region] for isoID in isoIDs), numpy.zeros(dispCount))
			dimAvg = dimSum / len(isoIDs)
			dimDevSum = sum(((dimAvg-zScores[isoID][region]) ** 2 for isoID in isoIDs), numpy.zeros(dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isoIDs))
			
			mainDim = max(range(93), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(zScores[isoID][region] * zScores[isoID][region][mainDim] for isoID in isoIDs) - dimSum * dimSum[mainDim] / len(isoIDs), ((len(isoIDs) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(range(93), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:dimsPerSplit]

			dimCoeffs = [1.0] * dimsPerSplit

			for i in range(1, dimsPerSplit):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]

			isoIDs.sort(key=lambda isoID: multiDimSplitVal(zScores[isoID][region], splitDims, dimCoeffs))
			splitValue = [multiDimSplitVal(zScores[isoIDs[len(isoIDs)//2]][region], splitDims, dimCoeffs)]

			regionSplitDims.append(splitDims)
			regionDimCoeffs.append(dimCoeffs)
			regionSplitValue.append(splitValue)

		for i, isoID in enumerate(isoIDs):
			truth = tuple(bool(multiDimSplitVal(zScores[isoID][region], regionSplitDims[j], regionDimCoeffs[j]) < regionSplitValue[j]) for j, (region, _) in enumerate(regions))
			childrenIDs[truth].append(isoID)

		children = [splitMultiRegionCorrelatedDims(zScores, childrenIDs[truth], regions, nodeType, pointsPerLeaf, dimsPerSplit) for truth in truthTable]

		return nodeType.newInner(children, [region for region, _ in regions], regionSplitDims, regionDimCoeffs, regionSplitValue)

# def printTree(tree, depth):
# 	if isinstance(tree, LeafNode):
# 		print(" "*depth+"<{}>".format(len(tree.pyroIDs)))
# 	elif isinstance(tree, SplitInnerNode):
# 		print(" "*depth+"{}:{}".format(tree.splitDims, tree.splitValues))
# 		for child in tree.children:
# 			printTree(child, depth+1)
# 	else:
# 		print(" "*depth+"BoundedNode")

def countNodes(tree):
	if isinstance(tree, LeafNode):
		return (0, 1)
	else:
		inner, leaves = (1, 0)
		for child in tree.children:
			childInner, childLeaves = countNodes(child)
			inner, leaves = inner + childInner, leaves + childLeaves
		return inner, leaves

def getMinMaxDepth(tree, depth):
	if isinstance(tree, LeafNode):
		return (depth, depth)
	else:
		minDepth, maxDepth = (666, -666)
		for child in tree.children:
			childMin, childMax = getMinMaxDepth(child, depth+1)
			minDepth, maxDepth = min(childMin, minDepth), max(childMax, maxDepth)
		return minDepth, maxDepth


class LevelStats:
	def __init__(self):
		self.calls = 0
		self.totalDists = 0
		self.avgSplitDist = 0
		self.avgFullDist = 0
		self.inner = 0
		self.branches = 0
	def average(self, count):
		self.calls /= count
		self.totalDists /= count
		self.avgSplitDist /= count
		self.avgFullDist /= count
		self.inner /= count
		self.branches /= count
	def calculate(self):
		if self.totalDists > 0:
			self.avgSplitDist /= self.totalDists
			self.avgFullDist /= self.totalDists
		if self.inner > 0:
			self.branches /= self.inner
	def printStats(self):
		print("calls at this level: {}".format(self.calls))
		print("avgSplitDist, avgFullDist, totalDists: {}, {}, {}".format(self.avgSplitDist, self.avgFullDist, self.totalDists))
		print("branchesPerNode, innerNodes: {}, {}".format(self.branches, self.inner))

class QueryStats:
	def __init__(self):
		self.level = -1
		self.levelStats = []
		self.correct = 0
		self.returnCount = 0
		# self.contained = 0
		# self.inBox = 0
		# self.intersected = 0
		self.innerChecked = 0
		self.leavesChecked = 0
		self.pointsChecked = 0
		self.pointsFetched = 0
		self.extraNodes = 0
		self.innerUsed = 0
		self.leavesUsed = 0
		self.alternateCount = 0
		self.branchPerInner = 0
		self.alternateDepth = 0
		self.altDeathDepth = 0
		self.altDeaths = 0
	def increaseLevel(self):
		self.level += 1
		if len(self.levelStats) <= self.level:
			self.levelStats.append(LevelStats())
		self.levelStats[self.level].calls += 1
	def decreaseLevel(self):
		self.level -= 1
	def average(self, count):
		self.correct /= count
		self.returnCount /= count
		# self.contained /= count
		# self.inBox /= count
		# self.intersected /= count
		self.innerChecked /= count
		self.leavesChecked /= count
		self.pointsChecked /= count
		self.pointsFetched /= count
		self.extraNodes /= count
		self.innerUsed /= count
		self.leavesUsed /= count
		self.alternateCount /= count
		self.branchPerInner /= count
		self.alternateDepth /= count
		self.altDeathDepth /= count
		self.altDeaths /= count
		for levelStat in self.levelStats:
			levelStat.average(count)
	def calculate(self):
		self.extraPoints = self.pointsChecked - (self.returnCount-self.pointsFetched)
		if self.innerChecked > 0:
			self.branchPerInner /= self.innerChecked
		if self.alternateCount > 0:
			self.alternateDepth /= self.alternateCount
		if self.altDeaths > 0:
			self.altDeathDepth /= self.altDeaths
		else:
			self.altDeathDepth = -1
		for levelStat in self.levelStats:
			levelStat.calculate()
	def printStats(self):
		print("\n\ttotals:")
		print("correct queries: {}".format(self.correct))
		print("returned = fetched + checked : {} = {} + {}".format(self.returnCount, self.pointsFetched, self.returnCount - self.pointsFetched))
		print("nodes checked = leaves + inner: {} = {} + {}".format(self.leavesChecked + self.innerChecked, self.leavesChecked, self.innerChecked))
		print("alternates checked: {}".format(self.alternateCount))
		print("extras = points + nodes: {} = {} + {}".format(self.extraPoints + self.extraNodes, self.extraPoints, self.extraNodes))
		print("points checked / used points = {}".format(self.pointsChecked/self.returnCount))
		print("points checked / used checked points = {}".format(self.pointsChecked/(self.returnCount-self.pointsFetched)))
		print("nodes checked / used nodes = {}".format((self.leavesChecked+self.innerChecked)/(self.leavesChecked+self.innerChecked-self.extraNodes)))
		print("leaves checked / used leaves = {}".format(self.leavesChecked/self.leavesUsed))
		print("inner checked / used inner = {}".format(self.innerChecked/self.innerUsed))
		print("avgBranchPerNode: {}".format(self.branchPerInner))
		print("avgAltDepth: {}".format(self.alternateDepth))
		print("altDeathsBeforeLeaves: {}".format(self.altDeaths))
		print("avgAltDeathDepth: {}".format(self.altDeathDepth))
		# print("nodes contained: {}".format(self.contained))
		# print("nodes containing point: {}".format(self.inBox))
		# print("nodes intersected: {}".format(self.intersected))

	def printLevelStats(self):
		for i, levelStat in enumerate(self.levelStats):
			print("\tlevel: {}".format(i))
			levelStat.printStats()

def fetchAll(tree, stats):
	if isinstance(tree, LeafNode):
		stats.pointsFetched += len(tree.pyroIDs)
		stats.leavesUsed += 1
		return tree.pyroIDs
	else:
		stats.innerUsed += 1
		result = set()
		for child in tree.children:
			result |= fetchAll(child, stats)
		return result


def regionDist(p1, p2):
	return numpy.linalg.norm(p1-p2)

def pointDist(regions, p1, p2):
	return max(regionDist(p1[region], p2[region]) for region in regions)

def rangeQuery(zScores, regions, tree, center, radius, stats, depth):
	stats.increaseLevel()

	alt = False
	if not tree.doesContainPoint(center):
		alt = True
		stats.alternateCount += 1
		stats.alternateDepth += depth

	result = set()
	if isinstance(tree, LeafNode):
		stats.leavesChecked += 1
		for pyroID in tree.pyroIDs:
			stats.pointsChecked += 1
			if pointDist(regions, center, zScores[pyroID]) <= radius:
				result.add(pyroID)

		if len(result) > 0:
			stats.leavesUsed += 1
	else:
		stats.innerChecked += 1
		stats.levelStats[stats.level].inner += 1

		i = 0
		for child in tree.getIntersectingChildren(center, radius, stats):
			if child.isContainedBy(center, radius):
				result |= fetchAll(child, stats)
			else:
				result |= rangeQuery(zScores, regions, child, center, radius, stats, depth+1)
			i += 1

		stats.branchPerInner += i
		stats.levelStats[stats.level].branches += i
		if alt and i == 0:
			stats.altDeaths += 1
			stats.altDeathDepth += depth

		if len(result) > 0:
			stats.innerUsed += 1

	if len(result) == 0:
		stats.extraNodes += 1

	stats.decreaseLevel()
	return result

# def nearestNeighbors(zScores, regions, tree, center, maxRadius, count, stats, depth):
# 	stats.increaseLevel()

# 	alt = False
# 	if not tree.doesContainPoint(center):
# 		alt = True
# 		stats.alternateCount += 1
# 		stats.alternateDepth += depth

# 	if isinstance(tree, LeafNode):
# 		stats.pointsChecked += len(tree.pyroIDs)
# 		stats.leavesChecked += 1

# 		result = sorted([(pointDist(regions, center, zScores[pyroID]), pyroID) for pyroID in tree.pyroIDs], key=itemgetter(0))
# 		result = [(dist, pyroID) for dist, pyroID in result if dist <= maxRadius]

# 		# result = sorted(tree.pyroIDs, key=lambda pyroID: pointDist(regions, center, zScores[pyroID]))[:count]
# 		# result = list(filter(lambda pyroID: pointDist(regions, center, zScores[pyroID]) <= maxRadius, result))

# 		if len(result) > 0:
# 			stats.leavesUsed += 1
# 	else:
# 		stats.innerChecked += 1
# 		stats.levelStats[stats.level].inner += 1

# 		result = []
# 		bestRadius = maxRadius
# 		distSortedChildren = tree.getDistSortedChildren(center, stats)
# 		i = 0
# 		for dist, child in distSortedChildren:
# 			if dist >= bestRadius:
# 				break

# 			i += 1
# 			# result = sorted(result + nearestNeighbors(zScores, regions, child, center, bestRadius, count, stats, depth+1), key=lambda pyroID: pointDist(regions, center, zScores[pyroID]))[:count]
# 			result += nearestNeighbors(zScores, regions, child, center, bestRadius, count, stats, depth+1)
# 			result = sorted(result, key=itemgetter(0))[:count]
# 			# if len(result) == count and pointDist(regions, center, zScores[result[-1]]) < bestRadius:
# 				# bestRadius = pointDist(regions, center, zScores[result[-1]])
# 			if len(result) == count and result[-1][0] < bestRadius:
# 				bestRadius = result[-1][0]

# 		# result = list(filter(lambda pyroID: pointDist(regions, center, zScores[pyroID]) <= maxRadius, result))
# 		result = [(dist, pyroID) for dist, pyroID in result if dist <= maxRadius]

# 		stats.branchPerInner += i
# 		stats.levelStats[stats.level].branches += i
# 		if alt and i == 0:
# 			stats.altDeaths += 1
# 			stats.altDeathDepth += depth

# 		if len(result) > 0:
# 			stats.innerUsed += 1

# 	if len(result) == 0:
# 		stats.extraNodes += 1

# 	stats.decreaseLevel()
# 	return result


def loadZScores():
	with open("zScores.pickle", mode='r+b') as zScoresFile:
		return pickle.load(zScoresFile)

def makeTrees(zScores, regions):
	regionDispCounts = [93, 95]
	pointsPerLeaf = 4
	# childrenPerInner = 3
	dimsPerSplit = 5

	trees = []
	# trees.append(splitMultiRegion(zScores, list(zScores.keys()), list(zip(regions, regionDispCounts)), MultiRegionSplitNode, pointsPerLeaf))
	# trees.append(splitMultiRegionCorrelatedDims(zScores, list(zScores.keys()), list(zip(regions, regionDispCounts)), MultiRegionSplitNode, pointsPerLeaf, 4))
	# trees.append(splitMultiRegionCorrelatedDims(zScores, list(zScores.keys()), list(zip(regions, regionDispCounts)), MultiRegionBoundedNode, pointsPerLeaf, 4))
	trees.append(splitMultiRegionCorrelatedDims(zScores, list(zScores.keys()), list(zip(regions, regionDispCounts)), MultiRegionCombinedNode, pointsPerLeaf, 4))

	return trees

def loadCorrect(radius):
	# with open("correctN100R{}.pickle".format(radius), mode='r+b') as correctFile:
	# 	correctNearest = pickle.load(correctFile)
	with open("correctR{}.pickle".format(radius), mode='r+b') as correctFile:
		correctRange = pickle.load(correctFile)
	return correctNearest, correctRange

def test(zScores, trees, correctNearest, correctRange, regions, radius, nearestCount):
	queryCount = 100
	queryIDs = random.sample(zScores.keys(), queryCount)
	# queryIDs = list(zScores.keys())
	# queryCount = len(queryIDs)

	for i, tree in enumerate(trees):
		print("\n----------")
		print("| tree {} |".format(i))
		print("----------\n")

		# printTree(tree, 0)
		print("inner, leaf nodes: {}".format(countNodes(tree)))
		print("min, max depth: {}".format(getMinMaxDepth(tree, 0)))

		# nearestStats = QueryStats()
		rangeStats = QueryStats()

		for i, pyroID in enumerate(queryIDs):
			print("{}/{}".format(i, len(queryIDs)))
			# resultN = set(nearestNeighbors(zScores, regions, tree, zScores[pyroID], radius, nearestCount, nearestStats, 0))
			# correctN = set(correctNearest[pyroID][:nearestCount])
			# if resultN == correctN:
			# 	nearestStats.correct += 1
			# else:
			# 	print("\n\tN {} --- {} / {} : {} / {}".format(pyroID, len(resultN - correctN), len(resultN), len(correctN - resultN), len(correctN)))
				# print("{} ----- {}".format({(pID, pointDist(regions, zScores[pyroID], zScores[pID])) for pID in resultN - correctN}, {(pID, pointDist(regions, zScores[pyroID], zScores[pID])) for pID in correctN - resultN}))
			# nearestStats.returnCount += len(resultN)

			resultR = rangeQuery(zScores, regions, tree, zScores[pyroID], radius, rangeStats, 0)
			# correctR = correctRange[pyroID]
			# if resultR == correctR:
			# 	rangeStats.correct += 1
			# else:
			# 	print("\tR {} --- {} / {} : {} / {}".format(pyroID, len(resultR - correctR), len(resultR), len(correctR - resultR), len(correctR)))
				# print("{} ----- {}".format({(pID, pointDist(regions, zScores[pyroID], zScores[pID])) for pID in resultR - correctR}, {(pID, pointDist(regions, zScores[pyroID], zScores[pID])) for pID in correctR -resultR}))
			rangeStats.returnCount += len(resultR)

		# nearestStats.average(queryCount)
		rangeStats.average(queryCount)

		# nearestStats.calculate()
		rangeStats.calculate()

		# print("\n\tNearest Stats")
		# nearestStats.printLevelStats()
		# nearestStats.printStats()
		print("\n\tRange Stats")
		rangeStats.printLevelStats()
		rangeStats.printStats()

if __name__ == '__main__':
	radius = 1.37
	nearestCount = 5
	regions = ['23-5', '16-23']

	zScores = loadZScores()
	trees = makeTrees(zScores, regions)
	correctNearest, correctRange = loadCorrect(radius)

	cProfile.run("test(zScores, trees, correctNearest, correctRange, regions, radius, nearestCount)")
