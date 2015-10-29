import numpy



class Region:
	def __init__(self, name, dispCount, index):
		self.name = name
		self.dispCount = dispCount
		self.index = index
		# self.hash = hash((name, dispCount))

	def __repr__(self):
		return self.name

	def __eq__(self, other):
		return self.index == other.index

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return self.index

# class Pyroprint:
# 	def __init__(self, region, zScores):
# 		self.region = region
# 		self.zScores = zScores

class Isolate:
	def __init__(self, name, regionsPyroprint):
		self.name = name
		self.regionsPyroprint = regionsPyroprint

	def regionsDist(self, other):
		return {region: numpy.linalg.norm(other.regionsPyroprint[region] - self.regionsPyroprint[region]) for region in set(self.regionsPyroprint.keys()) & set(other.regionsPyroprint.keys())}

	def isWithinRadiiOf(self, other, radii):
		# regionsDist = self.regionsDist(other)
		# return all(regionsDist[region] <= radii[region] for region in radii.keys())
		# print([(region.name, region.dispCount) for region in other.regionsPyroprint.keys()], [(region.name, region.dispCount) for region in self.regionsPyroprint.keys()], [(region.name, region.dispCount) for region in radii.keys()])
		return all((numpy.linalg.norm(other.regionsPyroprint[region] - self.regionsPyroprint[region]) <= radii[region] for region in radii.keys()))

	def __repr__(self):
		return self.name

	def __eq__(self, other):
		return self.name == other.name

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash(self.name)



# def combineRegionsBoundingBoxes(regions, regionsBoxes):
# 	return [BoundingBox.combineAll((regionBoxes[region] for regionBoxes in regionsBoxes), region.dispCount) for region in regions]

# def boundRegionsInBoundingBoxes(isolates, regions):
# 	return [BoundingBox.bound((isolate.regionsPyroprint[region] for isolate in isolates), region.dispCount) for region in regions]
