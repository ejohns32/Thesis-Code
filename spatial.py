import pickle
import numpy
import itertools
import math

import pyroprinting



class BoundingBox:
	def __init__(self, lowerCorner, upperCorner):
		self.lowerCorner = lowerCorner
		self.upperCorner = upperCorner
	def combine(box1, box2):
		return BoundingBox(numpy.minimum(box1.lowerCorner, box2.lowerCorner), numpy.maximum(box1.upperCorner, box2.upperCorner))
	def combineAll(boxes, dispCount):
		combined = BoundingBox(numpy.full(dispCount, numpy.finfo(numpy.float64).max), numpy.full(dispCount, numpy.finfo(numpy.float64).min))
		# combined = next(boxes)
		for box in boxes:
			combined = BoundingBox.combine(combined, box)
		return combined

	def bound(points, dispCount):
		dimMin = numpy.full(dispCount, numpy.finfo(numpy.float64).max)
		dimMax = numpy.full(dispCount, numpy.finfo(numpy.float64).min)
		for point in points:
			dimMin = numpy.minimum(dimMin, point)
			dimMax = numpy.maximum(dimMax, point)
		return BoundingBox(dimMin, dimMax)

	#TODO: verify this
	def containedByBall(self, point, radius):
		return numpy.linalg.norm(numpy.maximum(numpy.absolute(point-self.lowerCorner), numpy.absolute(point-self.upperCorner))) <= radius
	def containsPoint(self, point):
		return numpy.all(numpy.greater_equal(point, self.lowerCorner)) and numpy.all(numpy.less_equal(point, self.upperCorner))
	#TODO: verify this
	def distToPoint(self, point):
		return numpy.linalg.norm(numpy.maximum(numpy.maximum(point, self.lowerCorner)-point, point-numpy.minimum(point, self.upperCorner)))
	def intersectsBall(self, point, radius):
		return self.distToPoint(point) <= radius



class MultiDimPlane:
	def __init__(self, dims, coeffs):
		self.numDims = len(dims)
		self.dims = dims
		self.coeffs = coeffs
		self.recipSqrtNormCoeff = 1/math.sqrt(sum(coeffs[i]**2 for i in range(self.numDims)))

	def valueOf(self, point):
		#TODO: maybe coeffs of 0 for not dims would be faster?
		return sum(point[self.dims[i]] * self.coeffs[i] for i in range(self.numDims))
	# def valueDiff(self, val1, val2):
	# 	return max(0, (val1 - val2)/self.sqrtCoeffNorm)
	def signedValueDist(self, val1, val2):
		return (val1 - val2)*self.recipSqrtNormCoeff

class PlanePartition:
	def __init__(self, plane, lowerValue, upperValue):
		self.plane = plane
		self.lowerValue = lowerValue
		self.upperValue = upperValue

	# def distToPoint(self, point):
	# 	#TODO: don't recalc dimVal for every child
	# 	dimVal = self.plane.valueOf(center)
	# 	lowerDist = self.plane.valueDiff(self.lowerValue, dimVal)
	# 	upperDist = self.plane.valueDiff(dimVal, self.upperValue)
	# 	# at least one should always be 0, both if inside
	# 	return max(lowerDist, upperDist)






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
		boundingBox = BoundingBox.bound((isolate.regionsPyroprint[region] for isolate in isolates), region.dispCount)
		return BoundingBoxFilter(region, set(), boundingBox)

	def fromChildrenFilters(region, chidlrenFilters):
		boundingBox = BoundingBox.combineAll((child.boundingBox for child in chidlrenFilters), region.dispCount)
		return BoundingBoxFilter(region, chidlrenFilters, boundingBox)

	def update(self, isolates):
		self.boundingBox = BoundingBox.bound((isolate.regionsPyroprint[self.region] for isolate in isolates), self.region.dispCount)

	def aggregate(self):
		self.boundingBox = BoundingBox.combineAll((child.boundingBox for child in self.children), self.region.dispCount)


	def containedByQuery(self, queryIsolate, radii):
		return self.boundingBox.containedByBall(queryIsolate.regionsPyroprint[self.region], radii[self.region])

	def intersectsQuery(self, queryIsolate, radii):
		return self.boundingBox.intersectsBall(queryIsolate.regionsPyroprint[self.region], radii[self.region])



class PlanePartitionFilter(SpatialFilter):
	def __init__(self, region, children, plane, splitValue, isLeftOfParentSplit):
		SpatialFilter.__init__(self, region, children)
		self.plane = plane
		self.planeQueryValue = None
		self.splitValue = splitValue
		self.isLeftOfParentSplit = isLeftOfParentSplit
		self.queryDist = None

	def fromSide(region, isLeftOfParentSplit):
		return PlanePartitionFilter(region, set(), None, None, isLeftOfParentSplit)

	def fromSplitPlane(region, children, plane, splitValue, isLeftOfParentSplit):
		return PlanePartitionFilter(region, children, plane, splitValue, isLeftOfParentSplit)

	def update(self, isolates):
		pass # this would require looking at all isolates (unless we know which ones are were just removed...) and changing children

	def aggregate(self):
		pass # this would require looking at all isolates (unless we know which ones are were just removed...)

	def containedByQuery(self, queryIsolate, radii):
		return False # plane partition has infite volume

	def intersectsQuery(self, queryIsolate, radii):
		# return self.planePartition.intersectsBall(queryIsolate, radii[self.region])
		self.planeQueryValue = None

		if self.parent is None:
			self.queryDist = 0
		else:
			if self.parent.planeQueryValue is None:
				self.parent.planeQueryValue = self.parent.plane.valueOf(queryIsolate.regionsPyroprint[self.region])
			anscDist = self.parent.queryDist
			assert anscDist is not None

			splitDist = self.parent.plane.signedValueDist(self.parent.splitValue, self.parent.planeQueryValue)
			if self.isLeftOfParentSplit:
				splitDist = -splitDist
			self.queryDist = math.sqrt(anscDist**2 + splitDist**2) if splitDist > 0 else anscDist

		return self.queryDist <= radii[self.region]






class Node:
	def __init__(self, spatialFilters):
		self.spatialFilters = spatialFilters

	def isContainedBy(self, queryIsolate, radii):
		# TODO: needs at least one FROM EACH REGION
		# return any(spatialFilter.containedByQuery(queryIsolate, radii) for spatialFilter in self.spatialFilters)
		return False
	# def doesContainPoint(self, point):
	# 	return all(spatialFilter.containsPoint(queryIsolate, radii) for spatialFilter in self.spatialFilters)
	def intersectsQuery(self, queryIsolate, radii):
		return all(spatialFilter.intersectsQuery(queryIsolate, radii) for spatialFilter in self.spatialFilters)



class LeafNode(Node):
	def __init__(self, spatialFilters, isolates):
		Node.__init__(self, spatialFilters)
		self.isolates = set(isolates)
		self.count = len(self.isolates)

	def fetchAll(self):
		self.count = 0
		return self.isolates

	def rangeQuery(self, queryIsolate, radii):
		# if self.isContainedBy(queryIsolate, radii):
		# 	return self.fetchAll()

		result = set()

		for isolate in self.isolates:
			if isolate.isWithinRadiiOf(queryIsolate, radii):
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



class InnerNode(Node):
	def __init__(self, spatialFilters, children):
		Node.__init__(self, spatialFilters)
		self.children = children
		self.count = len(self.children)

	def fetchAll(self):
		self.count = 0
		return set(itertools.chain.from_iterable(child.fetchAll() for child in self.children))

	def rangeQuery(self, queryIsolate, radii):
		result = set()

		for child in self.children:
			if child.intersectsQuery(queryIsolate, radii):
				if child.isContainedBy(queryIsolate, radii): #TODO: maybe only do this for leaves
					result |= child.fetchAll()
				else:
					result |= child.rangeQuery(queryIsolate, radii)

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
			# spatialFilter.aggregate(itertools.chain.from_iterable(child.getMatchingSpatialFilters(spatialFilter) for child in self.children)
			spatialFilter.aggregate()


	# def rangeQuery(self, queryIsolate, radii):
	# 	#TODO: maybe only do this for leaves
	# 	if self.isContainedBy(queryIsolate, radii):
	# 		return self.fetchAll()

	# 	result = set()

	# 	if self.intersectsQuery(queryIsolate, radii):
	# 		# traversal.push()
	# 		for child in self.children:
	# 			result |= child.rangeQuery(queryIsolate, radii)
	# 		#traversal.pop()

	# 		if len(result) > 0:
	# 			self.children.removeAll(child for child in self.children if child.count == 0)
	# 			self.count = len(self.children)
	# 			if self.count > 0:
	# 				# update spatial
	# 			#else node will be deleted by parent

	# 	return result


	# def rangeQuery(self, queryIsolate, radii):
	# 	result = set()

	# 	for child in self.getIntersectingChildren(queryIsolate, radii):
	# 		if child.isContainedBy(queryIsolate, radii):
	# 			result |= child.fetchAll()
	# 		else:
	# 			result |= child.rangeQuery(queryIsolate, radii)

	# 	if len(result) > 0:
	# 		if len(self.children) == 0:
	# 			#delete node
	# 		else:
	# 			#merge children
	# 			self.regionBoundingBoxes = Box.combineRegions(self.regions, list(child.regionBoundingBoxes for child in self.children))

	# 	return result

	# def getIntersectingChildren(self, queryIsolate, radii):
	# 	intersectingChildren = self.children
	# 	for region in self.regions:
	# 		for spatialFilter in self.spatialFilters[region]
	# 			# intersectingChildren = child for child in intersectingChildren if spatialFilter.doesChildIntersect(child, queryIsolate, radii)
	# 			intersectingChildren = spatialFilter.filterIntersectingChildren(intersectingChildren, queryIsolate, radii)
	# 	return intersectingChildren



class TreeConfig:
	def __init__(self, regions, pointsPerLeaf, dimsPerSplit):
		self.regions = regions
		self.pointsPerLeaf = pointsPerLeaf
		self.dimsPerSplit = dimsPerSplit



class Tree:
	def __init__(self, treeConfig, isolates):
		self.treeConfig = treeConfig

		regionsAllDims = {region: range(region.dispCount) for region in self.treeConfig.regions}
		dummyRegionsIsLeftOfParentSplit = tuple(None for region in self.treeConfig.regions)
		dummyFilterCollector = {(filterType, region): [] for region in self.treeConfig.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}
		self.root = splitMultiRegionCorrelatedDims(isolates, self.treeConfig, regionsAllDims, dummyRegionsIsLeftOfParentSplit, dummyFilterCollector, 0)

	def rangeQuery(self, queryIsolate, radii):
		if self.root.intersectsQuery(queryIsolate, radii):
			if self.root.isContainedBy(queryIsolate, radii): #TODO: maybe only do this for leaves
				return self.root.fetchAll()
			else:
				return self.root.rangeQuery(queryIsolate, radii)
		else:
			return set()



#TODO: verify dims not reused
def splitMultiRegionCorrelatedDims(isolates, treeConfig, regionsUnusedDims, regionsIsLeftOfParentSplit, parentFilterCollector, depth):
	if len(isolates) <= treeConfig.pointsPerLeaf:
		spatialFilters = []
		for region in treeConfig.regions:
			bboxFilter = BoundingBoxFilter.fromIsolates(region, isolates)
			spatialFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		for region, isLeft in zip(treeConfig.regions, regionsIsLeftOfParentSplit):
			planeFilter = PlanePartitionFilter.fromSide(region, isLeft)
			spatialFilters.append(planeFilter)
			parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		return LeafNode(spatialFilters, isolates)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(treeConfig.regions))]
		# print(truthTable)
		childrenIsolates = {truth: [] for truth in truthTable}

		regionsSplitPlane = []
		regionsSplitValue = []
		for region in treeConfig.regions:
			dimSum = sum((isolate.regionsPyroprint[region] for isolate in isolates), numpy.zeros(region.dispCount))
			dimAvg = dimSum / len(isolates)
			dimDevSum = sum(((dimAvg-isolate.regionsPyroprint[region]) ** 2 for isolate in isolates), numpy.zeros(region.dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isolates))
			
			mainDim = max(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim])
			# mainDim = max(range(dispCount), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(isolate.regionsPyroprint[region] * isolate.regionsPyroprint[region][mainDim] for isolate in isolates) - dimSum * dimSum[mainDim] / len(isolates), ((len(isolates) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:treeConfig.dimsPerSplit]
			# splitDims = sorted(range(dispCount), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:treeConfig.dimsPerSplit]

			dimCoeffs = [1.0] * len(splitDims)

			for i in range(1, len(splitDims)):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]


			splitPlane = MultiDimPlane(splitDims, dimCoeffs)
			isolates.sort(key=lambda isolate: splitPlane.valueOf(isolate.regionsPyroprint[region]))
			splitValue = (splitPlane.valueOf(isolates[len(isolates)//2].regionsPyroprint[region]) + splitPlane.valueOf(isolates[len(isolates)//2 + 1].regionsPyroprint[region])) / 2

			regionsSplitPlane.append(splitPlane)
			regionsSplitValue.append(splitValue)

		for isolate in isolates:
			truth = tuple(bool(splitPlane.valueOf(isolate.regionsPyroprint[region]) < splitValue) for region, splitPlane, splitValue in zip(treeConfig.regions, regionsSplitPlane, regionsSplitValue))
			childrenIsolates[truth].append(isolate)

		children = []

		childFilterCollector = {(filterType, region): [] for region in treeConfig.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}

		childrenRegionsUnusedDims = {region: [dim for dim in regionsUnusedDims[region] if dim not in splitPlane.dims] for region, splitPlane in zip(treeConfig.regions, regionsSplitPlane)}
		for truth in truthTable:
			children.append(splitMultiRegionCorrelatedDims(childrenIsolates[truth], treeConfig, childrenRegionsUnusedDims, truth, childFilterCollector, depth+1))

		planeFilters = []
		for region, splitPlane, splitValue, isLeft in zip(treeConfig.regions, regionsSplitPlane, regionsSplitValue, regionsIsLeftOfParentSplit):
			planeFilter = PlanePartitionFilter.fromSplitPlane(region, childFilterCollector[PlanePartitionFilter, region], splitPlane, splitValue, isLeft)
			planeFilters.append(planeFilter)
			parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		bboxFilters = []
		for region in treeConfig.regions:
			bboxFilter = BoundingBoxFilter.fromChildrenFilters(region, childFilterCollector[BoundingBoxFilter, region])
			bboxFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		spatialFilters = planeFilters + bboxFilters if depth <= 3 else bboxFilters + planeFilters

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
