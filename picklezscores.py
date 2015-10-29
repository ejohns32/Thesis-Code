import mysql.connector
import pickle
import json
import numpy

import pyroprinting
import config


if __name__ == '__main__':
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

	isolates = [pyroprinting.Isolate(isoID, {regionNameLookup[regionName]: numpy.array(pyroprint) for regionName, pyroprint in regionsPyroprintMap.items()}) for isoID, regionsPyroprintMap in data.items()]

	print("{}/{}".format(len(data), len(data) + len(toDelete)))

	with open("isolates.pickle", mode='w+b') as zScoresFile:
		pickle.dump(isolates, zScoresFile)
