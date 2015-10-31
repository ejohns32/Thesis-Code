import mysql.connector
import pickle
import json
import numpy

import config




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




def loadIsolatesFromDB():
	with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
		mysqlConfig = json.load(mysqlConfigJson)

	cfg = config.loadConfig()
	regionNameLookup = {region.name: region for region in cfg.regions}

	# regions = [('23-5', 93), ('16-23', 95)]
	data = dict()

	cnx = mysql.connector.connect(**mysqlConfig)
	# cursor = cnx.cursor(buffered=True)
	cursor = cnx.cursor()
	query = ("SELECT p1.isoID, p1.appliedRegion, position, zHeight FROM zScores INNER JOIN pyroprints p1 USING(pyroID) LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID WHERE p2.isoID IS NULL ORDER BY isoID, appliedRegion, position")

	# 	Double checked with:
	# SELECT p1.isoID, p1.appliedRegion, p1.pyroID FROM pyroprints p1 LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID WHERE p2.isoID IS NULL ORDER BY isoID, appliedRegion
	# select pyroID, appliedRegion from pyroprints where isoID = "Sp-079";

	cursor.execute(query)

	for (isoID, region, position, zHeight) in cursor:
		if isoID not in data:
			data[isoID] = {region: []}
			print("{}[{}]".format(isoID, region))
		elif region not in data[isoID]:
			data[isoID][region] = []
			print("{}[{}]".format(isoID, region))
		assert position == len(data[isoID][region])
		data[isoID][region].append(zHeight)

	cursor.close()
	cnx.close()

	toDelete = []

	for isoID in data.keys():
		for region in cfg.regions:
			if region.name not in data[isoID]:
				toDelete.append(isoID)
				break
			
			assert len(data[isoID][region.name]) == region.dispCount

	for isoID in toDelete:
		del data[isoID]

	print("{}/{}".format(len(data), len(data) + len(toDelete)))

	return [Isolate(isoID, {regionNameLookup[regionName]: numpy.array(pyroprint) for regionName, pyroprint in regionsPyroprintMap.items()}) for isoID, regionsPyroprintMap in data.items()]

def loadIsolatesFromFile(subsetSize = "All"):
	with open("isolates{}.pickle".format(subsetSize), mode='r+b') as isolatesFile:
		return pickle.load(isolatesFile)

