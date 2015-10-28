import pickle
import numpy
import itertools
import math

import primatives
import pyroprints


class SpatialFilter:
	def __init__(self, region, children):
		self.region = region
		self.parent = None
		self.children = children
		for child in children:
			child.parent = self

	def deleteFromParent(self):
		if self.parent is not None:
			self.parent.children.remove(self)

	# def matches(self, other):
	# 	return type(self) is type(other) and self.region is other.region
		

class BoundingBoxFilter(SpatialFilter):
	def __init__(self, region, children, boundingBox):
		SpatialFilter.__init__(self, region, children)
		self.boundingBox = boundingBox

	def fromIsolates(region, isolates):
		boundingBox = primatives.BoundingBox.bound((isolate.regionsPyroprint[region] for isolate in isolates), region.dispCount)
		return BoundingBoxFilter(region, set(), boundingBox)

	def fromChildrenFilters(region, chidlrenFilters):
		boundingBox = primatives.BoundingBox.combineAll((spatialFilter.boundingBox for spatialFilter in chidlrenFilters), region.dispCount)
		return BoundingBoxFilter(region, chidlrenFilters, boundingBox)

	def update(self, isolates):
		self.boundingBox = primatives.BoundingBox.bound((isolate.regionsPyroprint[region] for isolate in isolates), self.region.dispCount)

	# def aggregate(self, others):
	# 	self.boundingBox = primatives.BoundingBox.combineAll((spatialFilter.boundingBox for spatialFilter in others), self.egion.dispCount)
	def aggregate(self):
		self.boundingBox = primatives.BoundingBox.combineAll((child.boundingBox for child in self.children), self.region.dispCount)


	def containedByQuery(self, center, radii):
		return self.boundingBox.containedByBall(center, radii[self.region])

	def intersectsQuery(self, center, radii):
		return self.boundingBox.intersectsBall(center, radii[self.region])


class PlanePartitionFilter(SpatialFilter):
	def __init__(self, region, children, plane, splitValue, isLeftOfParentSplit):
		SpatialFilter.__init__(self, chidlren, region)
		self.plane = plane
		self.childQueryValue = None
		self.splitValue = splitValue
		self.isLeftOfParentSplit = isLeftOfParentSplit
		self.queryDist = None

	def fromSide(region, isLeftOfParentSplit):
		return PlanePartitionFilter(region, set(), None, None, isLeftOfParentSplit)

	def fromSplitPlane(region, leftChildren, rightChildren, plane, splitValue, isLeftOfParentSplit):
		return PlanePartitionFilter(region, children, plane, splitValue, isLeftOfParentSplit)

	def update(self, isolates):
		pass # this would require looking at all isolates (unless we know which ones are were just removed...) and changing children

	def aggregate(self, others):
		pass # this would require looking at all isolates (unless we know which ones are were just removed...)

	def containedByQuery(self, center, radii):
		return False # plane partition has infite volume

	def intersectsQuery(self, center, radii):
		# return self.planePartition.intersectsBall(center, radii[self.region])
		self.childQueryValue = None

		if self.parent is None:
			self.queryDist = 0
		else:
			if self.parent.childQueryValue is None:
				self.parent.childQueryValue = self.parent.plane.valueOf(center)
			anscDist = self.parent.queryDist

			distSign = -1 if self.isLeftOfParentSplit else 1
			splitDist = distSign * self.parent.plane.signedValueDist(self.parent.splitValue, self.parent.childQueryValue)
			if self.isLeftOfParentSplit:
				splitDist = -splitDist
			self.queryDist = math.sqrt(anscDist**2 + splitDist**2) if splitDist > 0 else anscDist

		return self.queryDist < radii[self.region]


class Node:
	def __init__(self, spatialFilters):
		self.spatialFilters = spatialFilters

	def isContainedBy(self, center, radius):
		return any(spatialFilter.containedByQuery(center, radii) for spatialFilter in self.spatialFilters)
	# def doesContainPoint(self, point):
	# 	return all(spatialFilter.containsPoint(center, radii) for spatialFilter in self.spatialFilters)
	def intersectsQuery(self, center, radii):
		return all(spatialFilter.intersectsQuery(center, radii) for spatialFilter in self.spatialFilters)

class LeafNode(Node):
	def __init__(self, spatialFilters, isolates):
		Node.__init__(self, spatialFilters)
		self.isolates = set(isolates)
		self.count = len(self.isolates)

	def fetchAll(self):
		self.count = 0
		return self.isolates

	def rangeQuery(self, center, radii):
		# if self.isContainedBy(center, radii):
		# 	return self.fetchAll()

		result = set()

		for isolate in self.isolates:
			if isolate.isWithinRadiiOf(center, radii):
				result.add(isolate)

		if len(result) > 0:
			self.isolates -= result
			self.count = len(self.isolates)
			if self.count > 0:
				self.updateSpatialFilters()
			# else it will be deleted

		return result

	def updateSpatialFilters(self):
		for spatialFilter in self.spatialFilters:
			spatialFilter.update(self.isolates)

	# def getMatchingSpatialFilters(self, spatialFilterToMatch):
	# 	for spatialFilter in self.spatialFilters:
	# 		if spatialFilterToMatch.matches(spatialFilter):
	# 			return (spatialFilter,)

	# 	assert False, "Leaf doesn't have all spatial node types"


class InnerNode(Node):
	def __init__(self, spatialFilters, children):
		Node.__init__(self, spatialFilters)
		self.children = children
		self.count = len(self.children)

	def fetchAll(self):
		self.count = 0
		return set(itertools.chain(child.fetchAll() for child in self.children))

	def rangeQuery(self, center, radii):
		result = set()

		for child in self.children if child.intersectsQuery(center, radii):
			if child.isContainedBy(center, radii): #TODO: maybe only do this for leaves
				result |= child.fetchAll()
			else:
				result |= child.rangeQuery(center, radii)

		if len(result) > 0:
			for child in (child for child in self.children if child.count == 0):
				for spatialFilter in child.spatialFilters:
					spatialFilter.deleteFromParent()
				self.children.remove(child)
			self.count = len(self.children)

			if self.count > 0:
				self.updateSpatialFilters()
			#else node will be deleted by parent

		return result

	def updateSpatialFilters(self):
		for spatialFilter in self.spatialFilters:
			# #TODO: single traversal instead of one for each spatialFilter
			# spatialFilter.aggregate(itertools.chain(child.getMatchingSpatialFilters(spatialFilter) for child in self.children)
			spatialFilter.aggregate()

	# def getMatchingSpatialFilters(self, spatialFilterToMatch):
	# 	for spatialFilter in self.spatialFilters:
	# 		if spatialFilterToMatch.matches(spatialFilter):
	# 			return (spatialFilter,)

	# 	# none of ours matched, check children
	# 	return itertools.chain(child.getMatchingSpatialFilters(spatialFilterToMatch) for child in self.children


	# def rangeQuery(self, center, radii):
	# 	#TODO: maybe only do this for leaves
	# 	if self.isContainedBy(center, radii):
	# 		return self.fetchAll()

	# 	result = set()

	# 	if self.intersectsQuery(center, radii):
	# 		# traversal.push()
	# 		for child in self.children:
	# 			result |= child.rangeQuery(center, radii)
	# 		#traversal.pop()

	# 		if len(result) > 0:
	# 			self.children.removeAll(child for child in self.children if child.count == 0)
	# 			self.count = len(self.children)
	# 			if self.count > 0:
	# 				# update spatial
	# 			#else node will be deleted by parent

	# 	return result


	# def rangeQuery(self, center, radii):
	# 	result = set()

	# 	for child in self.getIntersectingChildren(center, radii):
	# 		if child.isContainedBy(center, radii):
	# 			result |= child.fetchAll()
	# 		else:
	# 			result |= child.rangeQuery(center, radii)

	# 	if len(result) > 0:
	# 		if len(self.children) == 0:
	# 			#delete node
	# 		else:
	# 			#merge children
	# 			self.regionBoundingBoxes = Box.combineRegions(self.regions, list(child.regionBoundingBoxes for child in self.children))

	# 	return result

	# def getIntersectingChildren(self, center, radii):
	# 	intersectingChildren = self.children
	# 	for region in self.regions:
	# 		for spatialFilter in self.spatialFilters[region]
	# 			# intersectingChildren = child for child in intersectingChildren if spatialFilter.doesChildIntersect(child, center, radii)
	# 			intersectingChildren = spatialFilter.filterIntersectingChildren(intersectingChildren, center, radii)
	# 	return intersectingChildren


class TreeConfig:
	def __init__(self, regions, pointsPerLeaf, dimsPerSplit):
		self.regions = regions
		self.pointsPerLeaf = pointsPerLeaf
		self.dimsPerSplit = dimsPerSplit

class Tree:
	def __init__(self, config, isolates):
		self.config = config

		regionsAllDims = {region: range(region.dispCount) for region in self.config.regions}
		dummyRegionsIsLeftOfParentSplit = tuple(None for region in self.config.regions)
		dummyFilterCollector = {(filterType, region): [] for region in regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}
		self.root = splitMultiRegionCorrelatedDims(isolates, self.config, regionsAllDims, dummyRegionsIsLeftOfParentSplit, dummyFilterCollector)

	def rangeQuery(self, center, radii):
		return self.root.rangeQuery(center, radii)


#TODO: verify dims not reused
def splitMultiRegionCorrelatedDims(isolates, config, regionsUnusedDims, regionsIsLeftOfParentSplit, parentFilterCollector):
	if len(isolates) <= config.pointsPerLeaf:
		spatialFilters = []
		for region, isLeft in zip(config.regions, regionsIsLeftOfParentSplit):
			planeFilter = PlanePartitionFilter.fromSide(region, isLeft)
			spatialFilters.append(planeFilter)
			parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		for region in config.regions:
			bboxFilter = BoundingBoxFilter.fromIsolates(region, isolates)
			spatialFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		return LeafNode(spatialFilters, isolates)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(config.regions))]
		# print(truthTable)
		childrenIsolates = {truth: [] for truth in truthTable}

		regionsSplitPlane = []
		regionsSplitValue = []
		for region in config.regions:
			dimSum = sum((isolate[region] for isolate in isolates), numpy.zeros(region.dispCount))
			dimAvg = dimSum / len(isolates)
			dimDevSum = sum(((dimAvg-isolate[region]) ** 2 for isolate in isolates), numpy.zeros(region.dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isolates))
			
			mainDim = max(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim])
			# mainDim = max(range(dispCount), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(isolate[region] * isolate[region][mainDim] for isolate in isolates) - dimSum * dimSum[mainDim] / len(isolates), ((len(isolates) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:config.dimsPerSplit]
			# splitDims = sorted(range(dispCount), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:config.dimsPerSplit]

			dimCoeffs = [1.0] * len(splitDims)

			for i in range(1, len(splitDims)):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]


			splitPlane = primatives.MultiDimPlane(splitDims, dimCoeffs)
			isolates.sort(key=lambda isolate: splitPlane.valueOf(isolate[region]))
			splitValue = splitPlane.valueOf(isolates[len(isolates)//2][region])

			regionsSplitPlane.append(splitPlane)
			regionsSplitValue.append(splitValue)

		for isolate in isolates:
			truth = tuple(bool(splitPlane.valueOf(isolate[region]) < splitValue) for region, splitPlane, splitValue in zip(config.regions, regionsSplitPlane, regionsSplitValue))
			childrenIsolates[truth].append(isolate)

		children = []

		childFilterCollector = {(filterType, region): [] for region in regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}

		childrenRegionsUnusedDims = {region: [dim for dim in regionsUnusedDims[region] if dim not in regionsSplitDims[region]] for region in regions}
		for truth in truthTable:
			children.append(splitMultiRegionCorrelatedDims(childrenIsolates[truth], config, childrenRegionsUnusedDims, truth, childFilterCollector))

		spatialFilters = []
		for region, splitPlane, splitValue, isLeft in zip(config.regions, regionsSplitPlane, regionsSplitValue, regionsIsLeftOfParentSplit):
			planeFilter = PlanePartitionFilter.fromSplitPlane(region, childFilterCollector[PlanePartitionFilter, region], splitPlane, splitValue, isLeft)
			spatialFilters.append(planeFilter)
			parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		for region in config.regions:
			bboxFilter = BoundingBoxFilter.fromChildrenFilters(region, childFilterCollector[BoundingBoxFilter, region])
			spatialFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		return InnerNode(spatialFilters, children)























# class MultiRegionSplitInnerNode(MultiRegionSplitNode, InnerNode):
# 	def __init__(self, children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
# 		InnerNode.__init__(self, children)
# 		self.parent = None
# 		for childIndex, child in enumerate(self.children):
# 			child.parent = self
# 			child.childIndex = childIndex
# 		self.regions = regions
# 		self.regionsSplitDims = regionsSplitDims
# 		self.regionsDimCoeffs = regionsDimCoeffs
# 		self.regionCoeffNorms = [sum(regionsDimCoeffs[i][j]**2 for j in range(len(regionsDimCoeffs[i])))for i in range(len(regions))]
# 		self.regionsSplitValues = [[-666.666] + splitValues + [666.666] for splitValues in regionsSplitValues]

# 	#TODO: reduce from 4 multipDimSplitValDiffs to 1 and from 2 multiDimSplitVal to 1
# 	def getRegionSplitDist(self, childIndex, region, center):
# 		dimVal = multiDimSplitVal(center, self.regionsSplitDims[region], self.regionsDimCoeffs[region])
# 		lowerDist = multiDimSplitValDiff(self.regionsSplitValues[region][childIndex], dimVal, self.regionCoeffNorms[region])
# 		upperDist = multiDimSplitValDiff(dimVal, self.regionsSplitValues[region][childIndex+1], self.regionCoeffNorms[region])
# 		# at least one should always be 0, both if inside
# 		return max(lowerDist, upperDist)



# class MultiRegionBoundedNode:
# 	def __init__(self, regions, regionBoundingBoxes):
# 		self.regions = regions
# 		self.regionBoundingBoxes = regionBoundingBoxes
# 		self.regionSeenBoxes = regionBoundingBoxes
# 	def isContainedBy(self, center, radii):
# 		return all(bbox.containedByBall(center[region], radius) for region, radius, bbox in zip(self.regions, radii, self.regionBoundingBoxes))
# 	def doesContainPoint(self, point):
# 		return all(bbox.containsPoint(point[region]) for region, bbox in zip(self.regions, self.regionBoundingBoxes))

# class MultiRegionBoundedInnerNode(MultiRegionBoundedNode, InnerNode):
# 	def __init__(self, children, regions):
# 		MultiRegionBoundedNode.__init__(self, regions, Box.combineRegions(regions, list(child.regionBoundingBoxes for child in children)))
# 		InnerNode.__init__(self, children)

# class MultiRegionBoundedLeafNode(MultiRegionBoundedNode, LeafNode):
# 	def __init__(self, isoIDs, regions, regionsDispCount, zScores):
# 		MultiRegionBoundedNode.__init__(self, [region for region in regions], Box.boundRegions(zScores, isoIDs, regions, regionsDispCount))
# 		LeafNode.__init__(self, isoIDs)



# class MultiRegionCombinedInnerNode(MultiRegionCombinedNode, MultiRegionBoundedInnerNode, MultiRegionSplitInnerNode):
# 	def getIntersectingChildren(self, center, radii):
# 		anscDist = {region: 0.0 for region in self.regions}
# 		if self.parent:
# 			anscDist = self.parent.memoDists[self.childIndex]

# 		result = set()
# 		self.memoDists = [{region: 0.0 for region in self.regions} for _ in range(len(self.children))]
# 		for childIndex, child in enumerate(self.children):
# 			if all(bbox.intersectsBall(center[region], radius) for region, radius, bbox in zip(self.regions, radii, child.regionBoundingBoxes)):

# 				splitDists = []
# 				regionDists = []
# 				for i, region in enumerate(self.regions):
# 					splitDist = self.getRegionSplitDist((childIndex >> i) % 2, i, center[region])
# 					dist = math.sqrt(anscDist[region]**2 + splitDist**2) if splitDist > 0 else anscDist[region]
# 					self.memoDists[childIndex][region] = dist
# 					splitDists.append(splitDist)
# 					regionDists.append(dist)

# 				if all(regionDist <= radius for regionDist, radius in zip(regionDists, radii)):
# 					result.add(child)

# 		return result
