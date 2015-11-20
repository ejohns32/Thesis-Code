import csv
import pickle

import pyroprinting
import config

def loadIsolateList(filename):
	cfg = config.loadConfig()
	isolates = pyroprinting.loadIsolatesFromFile("isolatesAll.pickle")
	isolateIdMap = {iso.name.strip(): iso for iso in isolates}

	with open(filename) as listFile:
		isoIds = {isoId.strip().strip("'").strip() for isoId in listFile.readline().split(',')}

	# print(isoIds)

	missingIsos = [iso.name for iso in isolates if iso.name.strip() not in isoIds]
	extraIsos = [isoId for isoId in isoIds if isoId not in isolateIdMap]
	print(extraIsos)
	print("extraIsoCount: {}".format(len(extraIsos)))
	print(missingIsos)
	print("missingIsoCount: {}".format(len(missingIsos)))

	sharedIsos = [iso for iso in isolates if iso.name.strip() in isoIds]

	print("{}/{} shared".format(len(sharedIsos), len(isolates)))
	with open("isolatesShared.pickle", mode='w+b') as cacheFile:
		pickle.dump(sharedIsos, cacheFile)

def loadFromCSV(filename, outfile):
	cfg = config.loadConfig()
	isolates = pyroprinting.loadIsolates(cfg)
	isolateIdMap = {iso.name.strip(): iso for iso in isolates}
	clusters = []

	with open(filename) as csvFile:
		csvLines = "".join(line for line in csvFile if line.strip()).splitlines(True) # remove blank lines
		# print(csvLines)
		csvReader = csv.reader(csvLines, delimiter=',')
		pastHeader = False # because Aldrin's csv files have some header rows
		currentClusterId = None
		currentCluster = None

		for i, row in enumerate(csvReader):
			# print("{}/{}".format(i+1, len(csvLines)))
			if row[0] == "Cluster Id":
				pastHeader = True
			elif pastHeader:
				if row[0].startswith("Threshold:") or row[0] == "******":
					print("Multiple clusterings detected in file. Skipping the rest.")
					break

				isoId = row[1].strip()
				if isoId in isolateIdMap:
					if row[0] != currentClusterId:
						currentClusterId = row[0]
						currentCluster = set()
						clusters.append(currentCluster)
					currentCluster.add(isolateIdMap[isoId])
				else:
					print("extra isolate: {}".format(isoId))

	# print(clusters)
	print(len(clusters))

	with open(outfile, mode='w+b') as cacheFile:
		pickle.dump(clusters, cacheFile)

def getOHClustClusters():
	with open("ohclust99.pickle", mode='r+b') as cacheFile:
		ohclustClusters = pickle.load(cacheFile)
	return ohclustClusters

def getAgglomerativeClusters():
	with open("agglomerative99.pickle", mode='r+b') as cacheFile:
		agglomerativeClusters = pickle.load(cacheFile)
	return agglomerativeClusters


if __name__ == '__main__':
	loadIsolateList("database_iso_ids")
	loadFromCSV("aldrins/database_oh_clusters_josh_ont_99.csv", "ohclust99.pickle")
	loadFromCSV("aldrins/database_oh_clusters_josh_ont_995.csv", "ohclust995.pickle")
	loadFromCSV("aldrins/database_oh_clusters_agglomerative_99.csv", "agglomerative99.pickle")
