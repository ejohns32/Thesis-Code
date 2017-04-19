import pickle
import json
import os.path
import math
import random

import mysql.connector
import numpy

import config



# TODO: use region.index to store regions in arrays instead of dicts
class Region:
	nextIndex = 0
	def __init__(self, name, dispCount, clusterThreshold, pSimThresholdAlpha, pSimThresholdBeta, betaDistributionAlpha, betaDistributionBeta):
		self.name = name
		self.dispCount = dispCount
		self.clusterThreshold = clusterThreshold
		self.pSimThresholdAlpha = pSimThresholdAlpha
		self.pSimThresholdBeta = pSimThresholdBeta
		self.betaDistributionAlpha = betaDistributionAlpha
		self.betaDistributionBeta = betaDistributionBeta

		self.index = Region.nextIndex
		Region.nextIndex += 1

	def fromDecodedJSON(decodedJSON):
		return Region(
			decodedJSON["name"],
			decodedJSON["dispCount"],
			decodedJSON["clusterThreshold"],
			decodedJSON["pearsonSimilarityThresholdAlpha"],
			decodedJSON["pearsonSimilarityThresholdBeta"],
			decodedJSON["betaDistributionAlpha"],
			decodedJSON["betaDistributionBeta"])

#	def __repr__(self):
#		return self.name

	def __eq__(self, other):
		return self.name == other.name

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash(self.name)

def pearsonFromDist(dist, dispCount):
	return 1 - dist**2 / (2*dispCount)

def distFromPearson(pearson, dispCount):
	return math.sqrt(2*dispCount * (1 - pearson))

class Isolate:
	def __init__(self, name, regionsPyroprintZscores):
		self.name = name
		self.regionsPyroprintZscores = regionsPyroprintZscores

	def regionDist(self, other, region):
		return numpy.linalg.norm(other.regionsPyroprintZscores[region] - self.regionsPyroprintZscores[region])

	def regionPearson(self, other, region):
		return numpy.sum(self.regionsPyroprintZscores[region] * other.regionsPyroprintZscores[region]) / region.dispCount

	def isWithinRadiiOf(self, other, radii):
		return all((numpy.linalg.norm(other.regionsPyroprintZscores[region] - self.regionsPyroprintZscores[region]) <= radius for region, radius in radii.items()))

	def __repr__(self):
		return self.name

	def __eq__(self, other):
		return self.name == other.name

	def __ne__(self, other):
		return not self.__eq__(other)

	def __lt__(self, other):
		return self.name < other.name

	def __hash__(self):
		return hash(self.name)


def loadIsolatesFromDB(cfg):
	print("loading all isolates from DB...")
	with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
		mysqlConfig = json.load(mysqlConfigJson)

	regionNameLookup = {region.name: region for region in cfg.regions}

	cnx = mysql.connector.connect(**mysqlConfig)
	# cursor = cnx.cursor(buffered=True)
	cursor = cnx.cursor()
	query = ("SELECT p1.isoID, p1.appliedRegion, position, zHeight FROM zScores INNER JOIN pyroprints p1 USING(pyroID) LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p2.isErroneous IS FALSE AND p1.pyroID < p2.pyroID WHERE p1.isErroneous IS FALSE AND p2.isoID IS NULL ORDER BY isoID, appliedRegion, position")

	# 		Double checked with:
	# SELECT p1.isoID, p1.appliedRegion, p1.pyroID FROM pyroprints p1 LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID AND p1.isErroneous IS FALSE AND p2.isErroneous IS FALSE WHERE p2.isoID IS NULL ORDER BY isoID, appliedRegion;
	# 		And
	# SELECT p1.isoID, p1.appliedRegion, p1.pyroID, p3.pyroID FROM pyroprints p1 LEFT JOIN pyroprints p2 ON p1.isoID = p2.isoID AND p1.appliedRegion = p2.appliedRegion AND p1.pyroID < p2.pyroID AND p1.isErroneous IS FALSE AND p2.isErroneous IS FALSE join pyroprints p3 ON p3.isoID = p1.isoID and p3.appliedRegion = p1.appliedRegion left join pyroprints p4 on p4.isoID = p3.isoID AND p4.appliedRegion = p3.appliedRegion and p3.pyroID < p4.pyroID WHERE p2.isoID IS NULL and p4.isoID IS NULL and p1.pyroID != p3.pyroID ORDER BY isoID, appliedRegion;
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
		# if position != len(data[isoID][region]):
		# 	print("\t{}: {} != {}".format(region, position, len(data[isoID][region])))
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
		print("Deleting {}", isoID);
		del data[isoID]

	print("{}/{} isolates had pyroprints for all regions".format(len(data), len(data) + len(toDelete)))

	return [Isolate(isoID, {regionNameLookup[regionName]: numpy.array(pyroprint) for regionName, pyroprint in regionsPyroprintMap.items()}) for isoID, regionsPyroprintMap in data.items()]

def getRandomSubset(isolates, cfg):
	print("sampling subset of isolates of size {}...".format(cfg.isolateSubsetSize))
	return list(random.sample(isolates, cfg.isolateSubsetSize))

def getIsolatesCacheFileName(cfg):
	return "isolates{}.pickle".format(cfg.isolateSubsetSize)

def loadIsolatesFromFile(cacheFileName):
	with open(cacheFileName, mode='r+b') as cacheFile:
		return pickle.load(cacheFile)

def loadIsolates(cfg):
	cacheFileName = getIsolatesCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return loadIsolatesFromFile(cacheFileName)
	elif os.path.isfile("isolatesAll.pickle"):
		isolates = loadIsolatesFromFile("isolatesAll.pickle")
		return getRandomSubset(isolates, cfg)
	else:
		isolates = loadIsolatesFromDB(cfg)
		if cfg.isolateSubsetSize == "All":
			return isolates
		else:
			return getRandomSubset(isolates, cfg)

