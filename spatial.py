import pickle
import numpy
import itertools
import math
import cProfile

import config
import pyroprinting
import fullsearch




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
		#TODO: this could be a zip()
		#TODO: maybe coeffs of 0 for not dims would be faster?
		#TODO: maybe using numpy array indexing would be faster?
		return sum(point[self.dims[i]] * self.coeffs[i] for i in range(self.numDims))
	# def valueDiff(self, val1, val2):
	# 	return max(0, (val1 - val2)/self.sqrtCoeffNorm)
	def signedValueDist(self, val1, val2):
		return (val1 - val2)*self.recipSqrtNormCoeff

	def __repr__(self):
		return "\n"+" + ".join("{:.2f}*{}".format(coeff, dim) for dim, coeff in zip(self.dims, self.coeffs))


# class BoundedPlane:
# 	def __init__(self, plane, lowerValue, upperValue):
# 		self.plane = plane
# 		self.lowerValue = lowerValue
# 		self.upperValue = upperValue

# 	def combine(bp1, bp2):
# 		assert bp1.plane
# 		return BoundedPlane

	# def distToPoint(self, point):
	# 	#TODO: don't recalc dimVal for every child
	# 	dimVal = self.plane.valueOf(center)
	# 	lowerDist = self.plane.valueDiff(self.lowerValue, dimVal)
	# 	upperDist = self.plane.valueDiff(dimVal, self.upperValue)
	# 	# at least one should always be 0, both if inside
	# 	return max(lowerDist, upperDist)

class BoundingPlanes:
	# def __init__(self, boundedPlanes):
		# self.boundedPlanes = boundedPlanes
	def __init__(self, planes, lowerCorner, upperCorner):
		self.planes = planes
		self.lowerCorner = lowerCorner #numpy.array
		self.upperCorner = upperCorner #numpy.array
		self.valDiffToDistScalers = numpy.array([plane.recipSqrtNormCoeff for plane in planes])
	def combine(bps1, bps2):
		assert bps1 == bps2
		# return BoundingPlanes(bps1.planes, [min(lower1, lower2) for lower1, lower2 in zip(bps1.lowerCorner, bps2.lowerCorner)], [max(upper1, upper2) for upper1, upper2 in zip(bps1.upperCorner, bps2.upperCorner)])
		return BoundingPlanes(bps1.planes, numpy.minimum(bps1.lowerCorner, bps2.lowerCorner), numpy.maximum(bps1.upperCorner, bps2.upperCorner))
		# return BoundingPlanes(BoundedPlane.combine(bp1, bp2) for bp1, bp2 in zip(bps1, bps2))
	def combineAll(bpses, planes):
		combinedLower = numpy.full(len(planes), numpy.finfo(numpy.float64).max)
		combinedUpper = numpy.full(len(planes), numpy.finfo(numpy.float64).min)
		# combined = next(bpses)
		for bps in bpses:
			assert bps.planes == planes
			combinedLower = numpy.minimum(combinedLower, bps.lowerCorner)
			combinedUpper = numpy.maximum(combinedUpper, bps.upperCorner)

		return BoundingPlanes(planes, combinedLower, combinedUpper)

	def bound(points, planes):
		planesMin = numpy.full(len(planes), numpy.finfo(numpy.float64).max)
		planesMax = numpy.full(len(planes), numpy.finfo(numpy.float64).min)
		for point in points:
			pointValues = numpy.array([plane.valueOf(point) for plane in planes])
			planesMin = numpy.minimum(planesMin, pointValues)
			planesMax = numpy.maximum(planesMax, pointValues)
		return BoundingPlanes(planes, planesMin, planesMax)

	#TODO: verify this
	def containedByBall(self, pointValues, radius):
		# TODO: use numpy for dist instead of signedValueDist()
		# pointLowerDists = numpy.array(plane.signedValueDist(pointVal, lower) for plane, lower, pointVal in zip(self.planes, self.lowerCorner, pointValues))
		# pointUpperDists = numpy.array(plane.signedValueDist(pointVal, upper) for plane, upper, pointVal in zip(self.planes, self.upperCorner, pointValues))
		pointLowerDists = (pointValues - self.lowerCorner) * self.valDiffToDistScalers
		pointUpperDists = (pointValues - self.upperCorner) * self.valDiffToDistScalers
		return numpy.linalg.norm(numpy.maximum(numpy.absolute(pointLowerDists), numpy.absolute(pointUpperDists))) <= radius
	def containsPoint(self, pointValues):
		return numpy.all(numpy.greater_equal(pointValues, self.lowerCorner)) and numpy.all(numpy.less_equal(pointValues, self.upperCorner))
	#TODO: verify this
	def distToPoint(self, pointValues):
		# TODO: use numpy for dist instead of signedValueDist()
		# pointLowerDists = numpy.array(plane.signedValueDist(lowerOutDist, pointVal) for plane, lowerOutDist, pointVal in zip(self.planes, numpy.maximum(pointValues, self.lowerCorner), pointValues))
		# pointUpperDists = numpy.array(plane.signedValueDist(pointVal, upperOutDist) for plane, upperOutDist, pointVal in zip(self.planes, numpy.minimum(pointValues, self.upperCorner), pointValues))
		pointLowerDists = (numpy.maximum(pointValues, self.lowerCorner) - pointValues) * self.valDiffToDistScalers
		pointUpperDists = (pointValues - numpy.minimum(pointValues, self.upperCorner)) * self.valDiffToDistScalers
		return numpy.linalg.norm(numpy.maximum(pointLowerDists, pointUpperDists))
	def intersectsBall(self, pointValues, radius):
		return self.distToPoint(pointValues) <= radius




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
		boundingBox = BoundingBox.bound((isolate.regionsPyroprintZscores[region] for isolate in isolates), region.dispCount)
		return BoundingBoxFilter(region, set(), boundingBox)

	def fromChildrenFilters(region, chidlrenFilters):
		boundingBox = BoundingBox.combineAll((child.boundingBox for child in chidlrenFilters), region.dispCount)
		return BoundingBoxFilter(region, chidlrenFilters, boundingBox)

	def update(self, isolates):
		self.boundingBox = BoundingBox.bound((isolate.regionsPyroprintZscores[self.region] for isolate in isolates), self.region.dispCount)

	def aggregate(self):
		self.boundingBox = BoundingBox.combineAll((child.boundingBox for child in self.children), self.region.dispCount)


	def containedByQuery(self, queryIsolate, radii):
		return self.boundingBox.containedByBall(queryIsolate.regionsPyroprintZscores[self.region], radii[self.region])

	def intersectsQuery(self, queryIsolate, radii):
		return self.boundingBox.intersectsBall(queryIsolate.regionsPyroprintZscores[self.region], radii[self.region])



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
				self.parent.planeQueryValue = self.parent.plane.valueOf(queryIsolate.regionsPyroprintZscores[self.region])
			anscDist = self.parent.queryDist
			assert anscDist is not None

			splitDist = self.parent.plane.signedValueDist(self.parent.splitValue, self.parent.planeQueryValue)
			if self.isLeftOfParentSplit:
				splitDist = -splitDist
			self.queryDist = math.sqrt(anscDist**2 + splitDist**2) if splitDist > 0 else anscDist

		return self.queryDist <= radii[self.region]

class BoundingPlanesFilter(SpatialFilter):
	def __init__(self, region, children, boundingPlanes):
		SpatialFilter.__init__(self, region, children)
		self.boundingPlanes = boundingPlanes

	def fromIsolatesAndPlanes(region, isolates, planes):
		boundingPlanes = BoundingPlanes.bound((isolate.regionsPyroprintZscores[region] for isolate in isolates), planes)
		return BoundingPlanesFilter(region, set(), boundingPlanes)

	def fromChildrenFiltersAndPlanes(region, chidlrenFilters, planes):
		boundingPlanes = BoundingPlanes.combineAll((child.boundingPlanes for child in chidlrenFilters), planes)
		return BoundingPlanesFilter(region, chidlrenFilters, boundingPlanes)

	def update(self, isolates):
		self.boundingPlanes = BoundingPlanes.bound((isolate.regionsPyroprintZscores[self.region] for isolate in isolates), self.boundingPlanes.planes)

	def aggregate(self):
		self.boundingPlanes = BoundingPlanes.combineAll((child.boundingPlanes for child in self.children), self.boundingPlanes.planes)


	# TODO: fix double grab for parent or double pointValues computation at top
	def containedByQuery(self, queryIsolate, radii):
		point = queryIsolate.regionsPyroprintZscores[self.region]
		if self.parent is None:
			self.pointValues = numpy.array([plane.valueOf(point) for plane in self.boundingPlanes.planes])
		else:
			self.pointValues = self.parent.pointValues
		return self.boundingPlanes.containedByBall(qself.pointValues, radii[self.region])

	def intersectsQuery(self, queryIsolate, radii):
		point = queryIsolate.regionsPyroprintZscores[self.region]
		if self.parent is None:
			self.pointValues = numpy.array([plane.valueOf(point) for plane in self.boundingPlanes.planes])
		else:
			self.pointValues = self.parent.pointValues
		return self.boundingPlanes.intersectsBall(self.pointValues, radii[self.region])




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

	def pop(self):
		rtn = self.isolates.pop()
		self.update()
		return rtn

	# def fetchAll(self):
	# 	self.count = 0
	# 	return self.isolates

	def rangeQuery(self, queryIsolate, radii, deleteResults):
		# if self.isContainedBy(queryIsolate, radii):
		# 	return self.fetchAll()

		result = set()

		for isolate in self.isolates:
			if isolate.isWithinRadiiOf(queryIsolate, radii):
				result.add(isolate)

		if deleteResults and len(result) > 0:
				self.isolates -= result
				self.update()

		return result

	def update(self):
		self.count = len(self.isolates)
		if self.count > 0:
			self.updateSpatialFilters()
		#else node will be deleted by parent

	def updateSpatialFilters(self):
		for spatialFilter in self.spatialFilters:
			spatialFilter.update(self.isolates)



class InnerNode(Node):
	def __init__(self, spatialFilters, children):
		Node.__init__(self, spatialFilters)
		self.children = children
		self.count = sum(child.count for child in self.children)

	def pop(self):
		rtn = self.children[0].pop()
		self.update()
		return rtn

	# def fetchAll(self):
	# 	self.count = 0
	# 	return set(itertools.chain.from_iterable(child.fetchAll() for child in self.children))

	def rangeQuery(self, queryIsolate, radii, deleteResults):
		result = set()

		for child in self.children:
			if child.intersectsQuery(queryIsolate, radii):
				# if child.isContainedBy(queryIsolate, radii): #TODO: maybe only do this for leaves
				# 	result |= child.fetchAll()
				# else:
				result |= child.rangeQuery(queryIsolate, radii, deleteResults)

		if deleteResults and len(result) > 0:
			self.update()

		return result

	def update(self):
		for child in [child for child in self.children if child.count == 0]:
			for spatialFilter in child.spatialFilters:
				spatialFilter.deleteFromParent()
			self.children.remove(child)
		self.count = sum(child.count for child in self.children)

		if self.count > 0:
			self.updateSpatialFilters()
		#else node will be deleted by parent

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






# class TreeConfig:
# 	def __init__(self, regions, pointsPerLeaf, dimsPerSplit):
# 		self.regions = regions
# 		self.pointsPerLeaf = pointsPerLeaf
# 		self.dimsPerSplit = dimsPerSplit


class SpatialIndex:
	def __init__(self, isolates, cfg):
		# self.treeConfig = treeConfig
		self.isolates = isolates
		self.root = splitConsistentMultiRegionCorrelatedDims(isolates, cfg)

	def __iter__(self):
		return iter(self.isolates)

	def __len__(self):
		return self.root.count

	def pop(self):
		# TODO: also remove from self.isolates
		return self.root.pop()

	def getNeighborsOf(self, queryIsolate, radii):
		#TODO have children check themselves in range query instead of having parents do it
		if self.root.intersectsQuery(queryIsolate, radii):
			return self.root.rangeQuery(queryIsolate, radii, False) - {queryIsolate}
		else:
			return set()


	def popNeighborsOf(self, queryIsolate, radii):
		#TODO have children check themselves in range query instead of having parents do it
		if self.root.intersectsQuery(queryIsolate, radii):
			# TODO: also remove from self.isolates
			return self.root.rangeQuery(queryIsolate, radii, True) - {queryIsolate}
		else:
			return set()


def splitMultiRegionCorrelatedDims(isolates, cfg):
	regionsAllDims = {region: range(region.dispCount) for region in cfg.regions}
	dummyRegionsIsLeftOfParentSplit = tuple(None for region in cfg.regions)
	dummyFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingBoxFilter,)}
	# dummyFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}
	return splitMultiRegionCorrelatedDimsRecurse(isolates, cfg, regionsAllDims, dummyRegionsIsLeftOfParentSplit, dummyFilterCollector, 0)

# sorry about this function :(
def splitMultiRegionCorrelatedDimsRecurse(isolates, cfg, regionsUnusedDims, regionsIsLeftOfParentSplit, parentFilterCollector, depth):
	if len(isolates) <= cfg.pointsPerLeaf:
		spatialFilters = []
		for region in cfg.regions:
			bboxFilter = BoundingBoxFilter.fromIsolates(region, isolates)
			spatialFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		# for region, isLeft in zip(cfg.regions, regionsIsLeftOfParentSplit):
		# 	planeFilter = PlanePartitionFilter.fromSide(region, isLeft)
		# 	spatialFilters.append(planeFilter)
		# 	parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		return LeafNode(spatialFilters, isolates)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(cfg.regions))]
		# print(truthTable)
		childrenIsolates = {truth: [] for truth in truthTable}

		regionsSplitPlane = []
		regionsSplitValue = []
		for region in cfg.regions:
			dimSum = sum((isolate.regionsPyroprintZscores[region] for isolate in isolates), numpy.zeros(region.dispCount))
			dimAvg = dimSum / len(isolates)
			dimDevSum = sum(((dimAvg-isolate.regionsPyroprintZscores[region]) ** 2 for isolate in isolates), numpy.zeros(region.dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isolates))
			
			mainDim = max(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim])
			# mainDim = max(range(dispCount), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(isolate.regionsPyroprintZscores[region] * isolate.regionsPyroprintZscores[region][mainDim] for isolate in isolates) - dimSum * dimSum[mainDim] / len(isolates), ((len(isolates) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(regionsUnusedDims[region], key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:cfg.dimsPerSplit]
			# splitDims = sorted(range(dispCount), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:cfg.dimsPerSplit]

			dimCoeffs = [1.0] * len(splitDims)

			for i in range(1, len(splitDims)):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]


			splitPlane = MultiDimPlane(splitDims, dimCoeffs)
			isolates.sort(key=lambda isolate: splitPlane.valueOf(isolate.regionsPyroprintZscores[region]))
			splitValue = (splitPlane.valueOf(isolates[len(isolates)//2].regionsPyroprintZscores[region]) + splitPlane.valueOf(isolates[len(isolates)//2 + 1].regionsPyroprintZscores[region])) / 2

			regionsSplitPlane.append(splitPlane)
			regionsSplitValue.append(splitValue)

		for isolate in isolates:
			truth = tuple(bool(splitPlane.valueOf(isolate.regionsPyroprintZscores[region]) < splitValue) for region, splitPlane, splitValue in zip(cfg.regions, regionsSplitPlane, regionsSplitValue))
			childrenIsolates[truth].append(isolate)

		children = []

		childFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingBoxFilter,)}
		# childFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}

		childrenRegionsUnusedDims = {region: [dim for dim in regionsUnusedDims[region] if dim not in splitPlane.dims] for region, splitPlane in zip(cfg.regions, regionsSplitPlane)}
		for truth in truthTable:
			children.append(splitMultiRegionCorrelatedDimsRecurse(childrenIsolates[truth], cfg, childrenRegionsUnusedDims, truth, childFilterCollector, depth+1))

		# planeFilters = []
		# for region, splitPlane, splitValue, isLeft in zip(cfg.regions, regionsSplitPlane, regionsSplitValue, regionsIsLeftOfParentSplit):
		# 	planeFilter = PlanePartitionFilter.fromSplitPlane(region, childFilterCollector[PlanePartitionFilter, region], splitPlane, splitValue, isLeft)
		# 	planeFilters.append(planeFilter)
		# 	parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		bboxFilters = []
		for region in cfg.regions:
			bboxFilter = BoundingBoxFilter.fromChildrenFilters(region, childFilterCollector[BoundingBoxFilter, region])
			bboxFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		# spatialFilters = planeFilters + bboxFilters if depth <= 2 else bboxFilters + planeFilters
		spatialFilters = bboxFilters

		return InnerNode(spatialFilters, children)




def splitConsistentMultiRegionCorrelatedDims(isolates, cfg):
	# dummyFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingPlanesFilter,)}
	dummyFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingBoxFilter,)}
	# dummyFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}
	dummyRegionsIsLeftOfParentSplit = tuple(None for region in cfg.regions)

	regionsSplitPlanes = []
	for region in cfg.regions:
		dimSum = sum((isolate.regionsPyroprintZscores[region] for isolate in isolates), numpy.zeros(region.dispCount))
		dimAvg = dimSum / len(isolates)
		dimDevSum = sum(((dimAvg-isolate.regionsPyroprintZscores[region]) ** 2 for isolate in isolates), numpy.zeros(region.dispCount))
		dimStdDev = numpy.sqrt(dimDevSum / len(isolates))
		

		unusedDims = range(region.dispCount)
		splitPlanes = []
		while len(unusedDims) > 0:
			mainDim = max(unusedDims, key=lambda dim: dimStdDev[dim])
			# mainDim = max(range(dispCount), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(isolate.regionsPyroprintZscores[region] * isolate.regionsPyroprintZscores[region][mainDim] for isolate in isolates) - dimSum * dimSum[mainDim] / len(isolates), ((len(isolates) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(unusedDims, key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:cfg.dimsPerSplit]
			# splitDims = sorted(range(dispCount), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:cfg.dimsPerSplit]

			dimCoeffs = [1.0] * len(splitDims)

			for i in range(1, len(splitDims)):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]


			splitPlanes.append(MultiDimPlane(splitDims, dimCoeffs))
			unusedDims = list(set(unusedDims) - set(splitDims))

		regionsSplitPlanes.append(splitPlanes)

	# print(regionsSplitPlanes)

	return splitConsistentMultiRegionCorrelatedDimsRecurse(isolates, cfg, regionsSplitPlanes, dummyRegionsIsLeftOfParentSplit, dummyFilterCollector, 0)

# sorry about this function :(
def splitConsistentMultiRegionCorrelatedDimsRecurse(isolates, cfg, regionsSplitPlanes, regionsIsLeftOfParentSplit, parentFilterCollector, depth):
	if len(isolates) <= cfg.pointsPerLeaf:
		spatialFilters = []

		# for region, splitPlanes in zip(cfg.regions, regionsSplitPlanes):
		# 	planeFilter = BoundingPlanesFilter.fromIsolatesAndPlanes(region, isolates, splitPlanes)
		# 	spatialFilters.append(planeFilter)
		# 	parentFilterCollector[BoundingPlanesFilter, region].append(planeFilter)

		for region in cfg.regions:
			bboxFilter = BoundingBoxFilter.fromIsolates(region, isolates)
			spatialFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		# for region, isLeft in zip(cfg.regions, regionsIsLeftOfParentSplit):
		# 	planeFilter = PlanePartitionFilter.fromSide(region, isLeft)
		# 	spatialFilters.append(planeFilter)
		# 	parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		return LeafNode(spatialFilters, isolates)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(cfg.regions))]
		# print(truthTable)
		childrenIsolates = {truth: [] for truth in truthTable}

		regionsSplitPlane = []
		regionsSplitValue = []
		for region, splitPlanes in zip(cfg.regions, regionsSplitPlanes):
			splitPlane = splitPlanes[depth % len(splitPlanes)]
			isolates.sort(key=lambda isolate: splitPlane.valueOf(isolate.regionsPyroprintZscores[region]))
			splitValue = (splitPlane.valueOf(isolates[len(isolates)//2].regionsPyroprintZscores[region]) + splitPlane.valueOf(isolates[len(isolates)//2 + 1].regionsPyroprintZscores[region])) / 2

			regionsSplitPlane.append(splitPlane)
			regionsSplitValue.append(splitValue)

		# TODO: combine with sort and take slices (or denote the range somehow ) instead of copying
		for isolate in isolates:
			truth = tuple(bool(splitPlane.valueOf(isolate.regionsPyroprintZscores[region]) < splitValue) for region, splitPlane, splitValue in zip(cfg.regions, regionsSplitPlane, regionsSplitValue))
			childrenIsolates[truth].append(isolate)

		children = []

		# childFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingPlanesFilter,)}
		childFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (BoundingBoxFilter,)}
		# childFilterCollector = {(filterType, region): [] for region in cfg.regions for filterType in (PlanePartitionFilter, BoundingBoxFilter)}

		for truth in truthTable:
			children.append(splitConsistentMultiRegionCorrelatedDimsRecurse(childrenIsolates[truth], cfg, regionsSplitPlanes, truth, childFilterCollector, depth+1))

		# spatialFilters = []
		# for region, splitPlanes in zip(cfg.regions, regionsSplitPlanes):
		# 	planeFilter = BoundingPlanesFilter.fromChildrenFiltersAndPlanes(region, childFilterCollector[BoundingPlanesFilter, region], splitPlanes)
		# 	spatialFilters.append(planeFilter)
		# 	parentFilterCollector[BoundingPlanesFilter, region].append(planeFilter)


		# planeFilters = []
		# for region, splitPlane, splitValue, isLeft in zip(cfg.regions, regionsSplitPlane, regionsSplitValue, regionsIsLeftOfParentSplit):
		# 	planeFilter = PlanePartitionFilter.fromSplitPlane(region, childFilterCollector[PlanePartitionFilter, region], splitPlane, splitValue, isLeft)
		# 	planeFilters.append(planeFilter)
		# 	parentFilterCollector[PlanePartitionFilter, region].append(planeFilter)

		bboxFilters = []
		for region in cfg.regions:
			bboxFilter = BoundingBoxFilter.fromChildrenFilters(region, childFilterCollector[BoundingBoxFilter, region])
			bboxFilters.append(bboxFilter)
			parentFilterCollector[BoundingBoxFilter, region].append(bboxFilter)

		# spatialFilters = planeFilters + bboxFilters if depth <= 2 else bboxFilters + planeFilters
		spatialFilters = bboxFilters

		return InnerNode(spatialFilters, children)



def testSpatial(isolates, index, correctNeighbors, cfg):
	queryIsolates = list(isolates)
	queryCount = len(queryIsolates)

	correctCount = 0
	nonZeroCorrectCount = 0
	nonZeroCount = 0
	extraCount = 0
	missingCount = 0
	seen = set()

	for i, isolate in enumerate(queryIsolates):
		resultR = index.getNeighborsOf(isolate, cfg.radii)
		correctR = correctNeighbors[isolate]
		# seen |= {isolate}
		# correctR -= seen
		# seen |= correctR

		if len(resultR) > 0:
			nonZeroCount += 1
		extraCount += len(resultR - correctR)
		missingCount += len(correctR - resultR)

		# # print("{}/{} - {}".format(i, len(queryIsolates), isolate))
		# print("{}/{} - {} - {}:{}".format(i+1, len(queryIsolates), isolate, len(resultR), len(correctR)))
		# # print("\t{} --- {} / {} : {} / {}".format(isolate, len(resultR - correctR), len(resultR), len(correctR - resultR), len(correctR)))
		if resultR == correctR:
			correctCount += 1
			if len(resultR) > 0:
				nonZeroCorrectCount += 1
		else:
			print("\t{} ----- {}".format([(extraIsolate, [isolate.regionDist(extraIsolate, region) for region in cfg.regions]) for extraIsolate in resultR - correctR], [(missingIsolate, [isolate.regionDist(missingIsolate, region) for region in cfg.regions]) for missingIsolate in correctR -resultR]))

	print("total correct: {}/{}".format(correctCount, queryCount))
	print("nonEmpty correct: {}/{}".format(nonZeroCorrectCount, nonZeroCount))
	print("total extra: {}, total missing: {}".format(extraCount, missingCount))

if __name__ == '__main__':
	cfg = config.loadConfig()
	isolates = pyroprinting.loadIsolates(cfg)
	index = SpatialIndex(isolates, cfg)
	correctNeighbors = fullsearch.getNeighborsMap(isolates, cfg)

	# testSpatial(isolates, index, correctNeighbors, cfg)
	cProfile.run("testSpatial(isolates, index, correctNeighbors, cfg)")

