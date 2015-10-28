import numpy
import math

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
		self.recipSqrtNormCoeff = 1/math.sqrt(sum(coeffs[i]**2 for j in range(self.numDims)))

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