import numpy



class Region:
	def __init__(self, name, dispCount):
		self.name = name
		self.dispCount = dispCount

# class Pyroprint:
# 	def __init__(self, region, zScores):
# 		self.region = region
# 		self.zScores = zScores

class Isolate:
	def __init__(self, isoID, regionPyroprintMap):
		self.isoID = isoID
		self.regionsPyroprint = {region, numpy.array(pyroprint) for region, pyroprint in regionPyroprintMap}

	def isWithinRadiiOf(self, other, radii):
		return all((numpy.linalg.norm(other.regionsPyroprint[region] - self.regionsPyroprint[region]) <= radii[region] for region in radii.keys()))



# def combineRegionsBoundingBoxes(regions, regionsBoxes):
# 	return [BoundingBox.combineAll((regionBoxes[region] for regionBoxes in regionsBoxes), region.dispCount) for region in regions]

# def boundRegionsInBoundingBoxes(isolates, regions):
# 	return [BoundingBox.bound((isolate.regionsPyroprint[region] for isolate in isolates), region.dispCount) for region in regions]
