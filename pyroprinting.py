import pickle
import json
import os.path
import math
import random

import mysql.connector
import numpy

import config




class Region:
	nextIndex = 0
	def __init__(self, name, dispCount):
		self.name = name
		self.dispCount = dispCount
		self.index = nextIndex
		nextIndex += 1
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

def pearsonFromDist(dist, dispCount):
	return 1 - dist**2 / (2*dispCount)

def distFromPearson(pearson, dispCount):
	return math.sqrt(2*dispCount * (1 - pearson))

class Isolate:
	def __init__(self, name, regionsPyroprintZscores):
		self.name = name
		self.regionsPyroprintZscores = regionsPyroprintZscores

	def regionsDist(self, other):
		return [(region, numpy.linalg.norm(other.regionsPyroprintZscores[region] - self.regionsPyroprintZscores[region])) for region in set(self.regionsPyroprintZscores.keys()) & set(other.regionsPyroprintZscores.keys())]

	def regionsPearson(self, other):
		return [(region, numpy.sum(self.regionsPyroprintZscores[region] * other.regionsPyroprintZscores[region]) / region.dispCount) for region in set(self.regionsPyroprintZscores.keys()) & set(other.regionsPyroprintZscores.keys())]

	def isWithinRadiiOf(self, other, radii):
		return all((numpy.linalg.norm(other.regionsPyroprintZscores[region] - self.regionsPyroprintZscores[region]) <= radius for region, radius in radii.items()))

	def __repr__(self):
		return self.name

	def __eq__(self, other):
		return self.name == other.name

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash(self.name)


def loadIsolatesFromDB(regions):
	with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
		mysqlConfig = json.load(mysqlConfigJson)

	regionNameLookup = {region.name: region for region in regions}

	cnx = mysql.connector.connect(**mysqlConfig)
	# cursor = cnx.cursor(buffered=True)
	cursor = cnx.cursor()
	query = ("SELECT p1.isoID, p1.appliedRegion, position, zHeight FROM zScores INNER JOIN pyroprints p1 USING(pyroID) LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID WHERE p2.isoID IS NULL ORDER BY isoID, appliedRegion, position")

	# 	Double checked with:
	# SELECT p1.isoID, p1.appliedRegion, p1.pyroID FROM pyroprints p1 LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID WHERE p2.isoID IS NULL ORDER BY isoID, appliedRegion
	# select pyroID, appliedRegion from pyroprints where isoID = "Sp-079";

	cursor.execute(query)

	data = dict()
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
		for region in regions:
			if region.name not in data[isoID]:
				toDelete.append(isoID)
				break
			
			assert len(data[isoID][region.name]) == region.dispCount

	for isoID in toDelete:
		del data[isoID]

	print("{}/{} isolates had pyroprints for all regions".format(len(data), len(data) + len(toDelete)))

	return [Isolate(isoID, {regionNameLookup[regionName]: numpy.array(pyroprint) for regionName, pyroprint in regionsPyroprintMap.items()}) for isoID, regionsPyroprintMap in data.items()]

def getRandomSubset(isolates, isolateSubsetSize):
	list(random.sample(isolates, isolateSubsetSize))

def loadIsolatesFromFile(cfg):
	cacheFileName = "isolates{}.pickle".format(cfg.isolateSubsetSize)
	with open(cacheFileName, mode='r+b') as cacheFile:
		return pickle.load(cacheFile)

def loadIsolates(cfg):
	cacheFileName = "isolates{}.pickle".format(cfg.isolateSubsetSize)
	if os.path.isfile(cacheFileName):
		return loadIsolatesFromFile(cfg)
	else:
		isolates = loadIsolatesFromDB(cfg.regions)
		if cfg.isolateSubsetSize = "All":
			return isolates
		else:
			return getRandomSubset(isolates, cfg.isolateSubsetSize)

