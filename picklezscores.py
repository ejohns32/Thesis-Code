import numpy
import mysql.connector
import pickle
import json

with open("mysqlconfig.json", mode='r') as mysqlconfig:
	config = json.load(mysqlconfig)

regions = [('23-5', 93), ('16-23', 95)]
data = dict()

cnx = mysql.connector.connect(**config)
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
	elif region not in data[isoID]:
		data[isoID][region] = []
	assert position == len(data[isoID][region])
	data[isoID][region].append(zHeight)

cursor.close()
cnx.close()

toDelete = []

for isoID in data.keys():
	for region, dispCount in regions:
		if region not in data[isoID]:
			toDelete.append(isoID)
			break
		
		assert len(data[isoID][region]) == dispCount
		data[isoID][region] = numpy.array(data[isoID][region])

for isoID in toDelete:
	del data[isoID]

print("{}/{}".format(len(data), len(data) + len(toDelete)))

with open("zScores.pickle", mode='w+b') as zScoresFile:
	pickle.dump(data, zScoresFile)
