import pickle
import numpy
import itertools
import math
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
	def combineRegions(regions, childrenRegionBoxes):
		return [Box.combineAll(regionBoxes[i] for regionBoxes in childrenRegionBoxes) for i in range(len(regions))]

	def bound(points, dispCount):
		dimMin = numpy.full(dispCount, numpy.finfo(numpy.float64).max)
		dimMax = numpy.full(dispCount, numpy.finfo(numpy.float64).min)
		for point in points:
			dimMin = numpy.minimum(dimMin, point)
			dimMax = numpy.maximum(dimMax, point)
		return Box(dimMin, dimMax)

	def boundRegions(zScores, isoIDs, regions, regionsDispCount):
		return [Box.bound((zScores[isoID][region] for isoID in isoIDs), dispCount) for region, dispCount in zip(regions, regionsDispCount)]

	#TODO: verify this
	def containedByBall(self, point, radius):
		return numpy.linalg.norm(numpy.maximum(numpy.absolute(point-self.lowerCorner), numpy.absolute(point-self.upperCorner))) < radius
	def containsPoint(self, point):
		return numpy.all(numpy.greater_equal(point, self.lowerCorner)) and numpy.all(numpy.less_equal(point, self.upperCorner))
	#TODO: verify this
	def distToPoint(self, point):
		return numpy.linalg.norm(numpy.maximum(numpy.maximum(point, self.lowerCorner)-point, point-numpy.minimum(point, self.upperCorner)))
	def intersectsBall(self, point, radius):
		return self.distToPoint(point) <= radius


def multiDimSplitVal(point, dims, coeffs):
	return sum(point[dims[i]] * coeffs[i] for i in range(len(dims)))
def multiDimSplitValDiff(val1, val2, coeffNorm):
	return max(0, (val1 - val2)/math.sqrt(coeffNorm))


def regionDist(p1, p2):
	return numpy.linalg.norm(p1-p2)
def pointInQuery(regions, center, radii, isolate):
	return all((regionDist(center[region], isolate[region]) <= radius for region, radius in zip(regions, radii)))



class UnseenTree:
	def constructWithAll(isolates):
		return UnseenTree()
	def popPointsInRange(self, center, radii):
		return set()

class SeenTree:
	def addBatch(self, isolates):
		pass
	def atLeastNInRange(self, center, radii, n):
		return False


class Node:
	def fetchAll(self):
		raise NotImplementedError()
	def rangeQuery(zScores, regions, regionsDispCount, tree, center, radii):
		raise NotImplementedError()


class LeafNode(Node):
	def __init__(self, isoIDs):
		Node.__init__(self)
		self.isoIDs = isoIDs

	def fetchAll(self):
		return self.isoIDs

	#TODO: only delete from unseen tree
	def rangeQuery(zScores, regions, regionsDispCount, tree, center, radii):
		result = set()

		for isoID in tree.isoIDs:
			if pointInQuery(regions, center, radii, zScores[isoID]):
				result.add(isoID)

		tree.isoIDs -= result
		#add result to seenTree? or is that outside this function?

		if len(result) > 0:
			if len(tree.isoIDs) == 0:
				#delete node
			else:
				tree.regionBoundingBoxes = Box.boundRegions(zScores, tree.isoIDs, regions, regionsDispCount)

		return result


class InnerNode(Node):
	def __init__(self, children):
		Node.__init__(self)
		self.children = children

	def fetchAll(self):
		# result = set()
		# for child in self.children:
		# 	result |= fetchAll(child)
		# return result
		return itertools.chain(child.fetchAll() for child in self.children)

	#TODO: only delete from unseen tree
	def rangeQuery(self, zScores, regions, regionsDispCount, center, radii):
		result = set()

		for child in self.getIntersectingChildren(center, radii):
			if child.isContainedBy(center, radii):
				result |= child.fetchAll()
				#delete/ merge
			else:
				result |= child.rangeQuery(zScores, regions, regionsDispCount, center, radii)

		if len(result) > 0:
			if len(self.children) == 0:
				#delete node
			else:
				#potentially merge children
				self.regionBoundingBoxes = Box.combineRegions(regions, list(child.regionBoundingBoxes for child in self.children))

		return result



class MultiRegionSplitNode:
	def newLeaf(isoIDs, regions, regionsDispCount, zScores):
		return MultiRegionSplitLeafNode(isoIDs)
	def newInner(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
		return MultiRegionSplitInnerNode(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues)
	def isContainedBy(self, center, radius):
		return False
	def doesContainPoint(self, point):
		return True

class MultiRegionSplitInnerNode(MultiRegionSplitNode, InnerNode):
	def __init__(self, children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
		InnerNode.__init__(self, children)
		self.parent = None
		for childIndex, child in enumerate(self.children):
			child.parent = self
			child.childIndex = childIndex
		self.regions = regions
		self.regionsSplitDims = regionsSplitDims
		self.regionsDimCoeffs = regionsDimCoeffs
		self.regionCoeffNorms = [sum(regionsDimCoeffs[i][j]**2 for j in range(len(regionsDimCoeffs[i])))for i in range(len(regions))]
		self.regionsSplitValues = [[-666.666] + splitValues + [666.666] for splitValues in regionsSplitValues]
	def getRegionSplitDist(self, childIndex, region, center):
		dimVal = multiDimSplitVal(center, self.regionsSplitDims[region], self.regionsDimCoeffs[region])
		lowerDist = multiDimSplitValDiff(self.regionsSplitValues[region][childIndex], dimVal, self.regionCoeffNorms[region])
		upperDist = multiDimSplitValDiff(dimVal, self.regionsSplitValues[region][childIndex+1], self.regionCoeffNorms[region])
		# print("dist({}, {}) = {}".format(dimVal, self.splitValues[childIndex], lowerDist))
		# print("dist({}, {}) = {}".format(dimVal, self.splitValues[childIndex+1], upperDist))
		# at least one should always be 0, both if inside
		return max(lowerDist, upperDist)

	def getIntersectingChildren(self, center, radii, stats):
		anscDist = {region: 0.0 for region in self.regions}
		if self.parent:
			anscDist = self.parent.memoDists[self.childIndex]

		result = []
		self.memoDists = [{region: 0.0 for region in self.regions} for _ in range(len(self.children))]
		for childIndex, child in enumerate(self.children):
			if child.hasBeenSeen(stats):
				continue
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

			if all(regionDist <= radius for regionDist, radius in zip(regionDists, radii)):
				result.append(child)

		return result

class MultiRegionSplitLeafNode(MultiRegionSplitNode, LeafNode):
	def __init__(self, isoIDs):
		LeafNode.__init__(self, isoIDs)


class MultiRegionBoundedNode:
	def newLeaf(isoIDs, regions, regionsDispCount, zScores):
		return MultiRegionBoundedLeafNode(isoIDs, regions, regionsDispCount, zScores)
	def newInner(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
		return MultiRegionBoundedInnerNode(children, regions)

	def __init__(self, regions, regionBoundingBoxes):
		self.regions = regions
		self.regionBoundingBoxes = regionBoundingBoxes
		self.regionSeenBoxes = regionBoundingBoxes
	def isContainedBy(self, center, radii):
		return all(bbox.containedByBall(center[region], radius) for region, radius, bbox in zip(self.regions, radii, self.regionBoundingBoxes))
	def doesContainPoint(self, point):
		return all(bbox.containsPoint(point[region]) for region, bbox in zip(self.regions, self.regionBoundingBoxes))

class MultiRegionBoundedInnerNode(MultiRegionBoundedNode, InnerNode):
	def __init__(self, children, regions):
		MultiRegionBoundedNode.__init__(self, regions, Box.combineRegions(regions, list(child.regionBoundingBoxes for child in children)))
		InnerNode.__init__(self, children)
	def getIntersectingChildren(self, center, radii):
		return (child for child in self.children if not child.hasBeenSeen(stats) and all(bbox.intersectsBall(center[region], radius) for region, radius, bbox in zip(self.regions, radii, child.regionBoundingBoxes)))

class MultiRegionBoundedLeafNode(MultiRegionBoundedNode, LeafNode):
	def __init__(self, isoIDs, regions, regionsDispCount, zScores):
		MultiRegionBoundedNode.__init__(self, [region for region in regions], Box.boundRegions(zScores, isoIDs, regions, regionsDispCount))
		LeafNode.__init__(self, isoIDs)


class MultiRegionCombinedNode(MultiRegionBoundedNode, MultiRegionSplitNode):
	def newLeaf(isoIDs, regions, regionsDispCount, zScores):
		return MultiRegionCombinedLeafNode(isoIDs, regions, regionsDispCount, zScores)
	def newInner(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
		return MultiRegionCombinedInnerNode(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues)

	def isContainedBy(self, center, radii):
		return MultiRegionBoundedNode.isContainedBy(self, center, radii) # or MultiRegionSplitNode.isContainedBy(self, center, radii)
	def doesContainPoint(self, point):
		return MultiRegionBoundedNode.doesContainPoint(self, point) # and MultiRegionSplitNode.doesContainPoint(self, point)

class MultiRegionCombinedInnerNode(MultiRegionCombinedNode, MultiRegionBoundedInnerNode, MultiRegionSplitInnerNode):
	def __init__(self, children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues):
		MultiRegionBoundedInnerNode.__init__(self, children, regions)
		MultiRegionSplitInnerNode.__init__(self, children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValues)
	def getIntersectingChildren(self, center, radii):
		anscDist = {region: 0.0 for region in self.regions}
		if self.parent:
			anscDist = self.parent.memoDists[self.childIndex]

		result = set()
		self.memoDists = [{region: 0.0 for region in self.regions} for _ in range(len(self.children))]
		for childIndex, child in enumerate(self.children):
			if all(bbox.intersectsBall(center[region], radius) for region, radius, bbox in zip(self.regions, radii, child.regionBoundingBoxes)):

				splitDists = []
				regionDists = []
				for i, region in enumerate(self.regions):
					splitDist = self.getRegionSplitDist((childIndex >> i) % 2, i, center[region])
					dist = math.sqrt(anscDist[region]**2 + splitDist**2) if splitDist > 0 else anscDist[region]
					self.memoDists[childIndex][region] = dist
					splitDists.append(splitDist)
					regionDists.append(dist)

				if all(regionDist <= radius for regionDist, radius in zip(regionDists, radii)):
					result.add(child)

		return result


class MultiRegionCombinedLeafNode(MultiRegionCombinedNode, MultiRegionBoundedLeafNode, MultiRegionSplitLeafNode):
	def __init__(self, isoIDs, regions, regionsDispCount, zScores):
		MultiRegionBoundedLeafNode.__init__(self, isoIDs, regions, regionsDispCount, zScores)
		MultiRegionSplitLeafNode.__init__(self, isoIDs)



#TODO: verify dims not reused
def splitMultiRegionCorrelatedDims(zScores, isoIDs, regions, regionsDispCount, nodeType, pointsPerLeaf, dimsPerSplit, regionsUnusedDims):
	if len(isoIDs) <= pointsPerLeaf:
		return nodeType.newLeaf(set(isoIDs), regions, regionsDispCount, zScores)
	else:
		truthTable = [tuple(reversed(truth)) for truth in itertools.product((True, False), repeat=len(regions))]
		# print(truthTable)
		childrenIDs = {truth: [] for truth in truthTable}

		regionsSplitDims = []
		regionsDimCoeffs = []
		regionsSplitValue = []
		for region, dispCount, unusedDims in zip(regions, regionsDispCount, regionsUnusedDims):
			dimSum = sum((zScores[isoID][region] for isoID in isoIDs), numpy.zeros(dispCount))
			dimAvg = dimSum / len(isoIDs)
			dimDevSum = sum(((dimAvg-zScores[isoID][region]) ** 2 for isoID in isoIDs), numpy.zeros(dispCount))
			dimStdDev = numpy.sqrt(dimDevSum / len(isoIDs))
			
			mainDim = max(unusedDims, key=lambda dim: dimStdDev[dim])
			# mainDim = max(range(dispCount), key=lambda dim: dimStdDev[dim])
			correlations = numpy.divide(sum(zScores[isoID][region] * zScores[isoID][region][mainDim] for isoID in isoIDs) - dimSum * dimSum[mainDim] / len(isoIDs), ((len(isoIDs) - 1) * dimStdDev * dimStdDev[mainDim]))

			splitDims = sorted(unusedDims, key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:dimsPerSplit]
			# splitDims = sorted(range(dispCount), key=lambda dim: dimStdDev[dim] * abs(correlations[dim]), reverse=True)[:dimsPerSplit]

			dimCoeffs = [1.0] * dimsPerSplit

			for i in range(1, dimsPerSplit):
				dimCoeffs[i] = correlations[splitDims[i]] * dimStdDev[splitDims[0]] / dimStdDev[splitDims[i]]

			isoIDs.sort(key=lambda isoID: multiDimSplitVal(zScores[isoID][region], splitDims, dimCoeffs))
			splitValue = [multiDimSplitVal(zScores[isoIDs[len(isoIDs)//2]][region], splitDims, dimCoeffs)]

			regionsSplitDims.append(splitDims)
			regionsDimCoeffs.append(dimCoeffs)
			regionsSplitValue.append(splitValue)

		for i, isoID in enumerate(isoIDs):
			truth = tuple(bool(multiDimSplitVal(zScores[isoID][region], splitDims, dimCoeffs) < splitValue) for region, splitDims, dimCoeffs, splitValue in zip(regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValue))
			childrenIDs[truth].append(isoID)

		children = [splitMultiRegionCorrelatedDims(zScores, childrenIDs[truth], regions, regionsDispCount, nodeType, pointsPerLeaf, dimsPerSplit, [[dim for dim in regionsUnusedDims[regionI] if dim not in regionsSplitDims[regionI]] for regionI in range(len(regions))]) for truth in truthTable]

		return nodeType.newInner(children, regions, regionsSplitDims, regionsDimCoeffs, regionsSplitValue)
