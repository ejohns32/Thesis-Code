import os.path
import pickle
import json

import mysql.connector
from scipy import stats

import config
import pyroprinting
import dbscan
import importClusters


def getDistributionCorrectnessFunc(expectedDistribution):
	def func(actual):
		return stats.ks_2samp(actual, expectedDistribution)
	return func

def getBetaCorrectnessFunc(expectedBeta):
	def func(actual):
		return stats.kstest(actual, expectedBeta.cdf)
	return func

def getIntraClusterPearsons(cluster, region):
	pearsons = []

	clusterList = list(cluster)
	for i, iso1 in enumerate(clusterList):
		for iso2 in clusterList[i+1:]:
			pearsons.append(iso1.regionPearson(iso2, region))
			# pearsons.append(pearsonMap[iso1, iso2] if iso1 < iso2 else pearsonMap[iso2, iso1])

	return pearsons

def getInterClusterPearsons(cluster1, cluster2, region):
	pearsons = []

	for iso1 in cluster1:
		for iso2 in cluster2:
			pearsons.append(iso1.regionPearson(iso2, region))


	return pearsons

def getCombinedIntraClusterPearsons(clusters, region):
	pearsons = []

	for cluster in clusters:
		pearsons.extend(getIntraClusterPearsons(cluster, region))

	return pearsons

def distributionSimilarity(clusters1, clusters2, region):
	ksStat, pValue = stats.ks_2samp(getCombinedIntraClusterPearsons(clusters1, region), getCombinedIntraClusterPearsons(clusters2, region))
	return pValue



def distributionCorrectness(clusters, region, correctnessFunc):
	ksStat, pValue = correctnessFunc(getCombinedIntraClusterPearsons(clusters, region))
	# ksStat, pValue = stats.kstest(getCombinedIntraClusterPearsons(clusters, region), betaDistribution.cdf)
	return pValue

# TODO: try looking at average chance in addition to complete chance
def individualClusterDistributionCorrectness(clusters, region, correctnessFunc):
	# TODO: double check that this is statistically valid
	# null hypothesis: all individual clusters fit
	# reject if any individuals reject
	# want chance that none reject. given chance that each doesn't reject
	# so all pValues multiplied together

	# print()
	totalCount = 0
	chanceCount = 0
	averageP = 0.0
	chanceNoneReject = 1.0
	for cluster in clusters:
		if len(cluster) > 1: # otherwise it wont have any pairwise pearsons
			ksStat, pValue = correctnessFunc(getIntraClusterPearsons(cluster, region))
			# ksStat, pValue = stats.kstest(getIntraClusterPearsons(cluster, region), betaDistribution.cdf)
			# print(pValue)
			chanceNoneReject *= pValue
			averageP += pValue
			totalCount += 1
			if pValue >= .1:
				chanceCount += 1

	return chanceCount/totalCount, averageP/totalCount, chanceNoneReject

# TODO: try looking at average chance in addition to complete chance
# seems to need large clusters to work
def pairedClustersDistributionCorrectness(clusters, region, correctnessFunc):
	# TODO: double check that this is statistically valid
	# null hypothesis: no pairs of clusters fit
	# reject if any pairs fit
	# want chance that all reject. given chance that each doesn't reject
	# so all inverse pValues multiplied together

	intraPearsons = [getIntraClusterPearsons(cluster, region) for cluster in clusters]

	# print()
	totalCount = 0
	chanceCount = 0
	averageP = 0.0
	chanceAllReject = 1.0
	for i, cluster1 in enumerate(clusters):
		# print("{}/{}".format(i+1, len(clusters)))
		for  j, cluster2 in enumerate(clusters[i+1:]):
			# each cluster will have at least 1 iso, so combined will have at least 2 iso, so at least one pairwise pearson
			pearsons = intraPearsons[i] + intraPearsons[i+1+j] + getInterClusterPearsons(cluster1, cluster2, region)
			# ksStat, pValue = stats.kstest(pearsons, betaDistribution.cdf)
			ksStat, pValue = correctnessFunc(pearsons)
			# ksStat, pValue = stats.kstest(getIntraClusterPearsons(cluster1 | cluster2, region), betaDistribution.cdf)
			# print(1-pValue)
			# print(pValue)
			# print(cluster1)
			# print(cluster2)
			chanceAllReject *= 1 - pValue
			averageP += 1 - pValue
			totalCount += 1
			if pValue <= .1:
				chanceCount += 1

	return chanceCount/totalCount, averageP/totalCount, chanceAllReject





def getClusterMap(clusters):
	clustersMap = dict()

	for cluster in clusters:
		for isolate in cluster:
			# assert isolate not in clustersMap
			if isolate in clustersMap:
				print("\n\tisolate {} in multiple clusters".format(isolate))
				print(clustersMap[isolate])
				print(cluster)
			clustersMap[isolate] = cluster

	return clustersMap

# TODO: look at both regions at once and only count same if both regions are definitely similar
def thresholdCorrectness(isolates, clusters, regions):
	sameDS = diffDD = totDS = totDD = totSquishy = 0
	clusterMap = getClusterMap(clusters)

	for i, iso1 in enumerate(isolates):
		# print("{}/{}".format(i+1, len(isolates)))
		for iso2 in isolates[i+1:]:
			# print("{},{}".format(iso1, iso2))
			# pearson = pearsonMap[iso1, iso2] if iso1 < iso2 else pearsonMap[iso2, iso1]
			regionsPearson = [iso1.regionPearson(iso2, region) for region in regions]
			# print(regionsPearson)

			if all(pearson >= region.pSimThresholdAlpha for region, pearson in zip(regions, regionsPearson)):
			# if pearson >= dsThreshold:
				# print("\tDS")
				totDS += 1
				if iso1 in clusterMap and iso2 in clusterMap and clusterMap[iso1] == clusterMap[iso2]:
					# print("Same Cluster")
					sameDS += 1
				# elif iso1 not in clusterMap or iso2 not in clusterMap:
				# 	print("not in any cluster")

			elif any(pearson <= region.pSimThresholdBeta for region, pearson in zip(regions, regionsPearson)):
			# elif pearson <= ddThreshold:
				# print("\tDD")
				totDD += 1
				if iso1 not in clusterMap or iso2 not in clusterMap or clusterMap[iso1] != clusterMap[iso2]:
					# print("Diff Cluster")
					diffDD += 1
			else:
				totSquishy += 1

	correctDS = sameDS/totDS
	correctDD = diffDD/totDD
	# tcScore = (sameDS+diffDD)/(totDS+totDD)
	squishyRatio = totSquishy/(totSquishy+totDD+totDS)

	# print("correctDS = {}".format(correctDS))
	# print("correctDD = {}".format(correctDD))
	# print("tcScore = {}".format(tcScore))
	# print("squishyRatio = {}".format(squishyRatio))
	# print("totDS: {}, totSquishy: {}, totDD: {}".format(totDS, totSquishy, totDD))

	return correctDS, correctDD, squishyRatio



def getPairCounts(isolates, clusters1, clusters2):
	bothSame = 0
	oneOfEach = 0
	bothDiff = 0

	clusterMap1 = getClusterMap(clusters1)
	clusterMap2 = getClusterMap(clusters2)

	for i, iso1 in enumerate(isolates):
		# print("{}/{}".format(i+1, len(isolates)))
		for iso2 in isolates[i+1:]:
			sameClust1 = clusterMap1[iso1] == clusterMap1[iso2]
			sameClust2 = clusterMap2[iso1] == clusterMap2[iso2]

			if sameClust1 == sameClust2:
				if sameClust1:
					bothSame += 1
				else:
					bothDiff += 1
			else:
				oneOfEach += 1

	return bothSame, oneOfEach, bothDiff

def similarityIndexes(isolates, clusters1, clusters2):
	bothSame, oneOfEach, bothDiff = getPairCounts(isolates, clusters1, clusters2)

	randIndex = (bothSame + bothDiff) / (bothSame + oneOfEach + bothDiff)
	jaccardIndex = (bothSame) / (bothSame + oneOfEach)
	anotherIndex = (bothDiff) / (bothDiff + oneOfEach)

	# print("randIndex: {}".format(randIndex))
	# print("jaccardIndex: {}".format(jaccardIndex))
	# print("anotherIndex: {}".format(anotherIndex))

	return randIndex, jaccardIndex, anotherIndex



def loadReplicatePearsonsFromDB(cfg):
	print("loading replicate pearsons from DB...")

	regionNameLookup = {region.name: region for region in cfg.regions}

	with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
		mysqlConfig = json.load(mysqlConfigJson)

	cnx = mysql.connector.connect(**mysqlConfig)
	# cursor = cnx.cursor(buffered=True)
	cursor = cnx.cursor()

	maxReplicatesPerIsolateRegion = 10 # to prevent bias
	query = (
		"SELECT r.appliedRegion as region, PearsonMatch(p1.pyroID, p2.pyroID, r.pearsonDispLength) as pearson FROM Pyroprints AS p1 INNER JOIN Pyroprints AS p2 USING(`isoID`, `appliedRegion`) INNER JOIN Regions AS r USING(`appliedRegion`) WHERE (SELECT COUNT(*) FROM Pyroprints AS pyro WHERE pyro.isoID = p1.isoID AND pyro.appliedRegion = p1.appliedRegion AND pyro.pyroID > p1.pyroID AND pyro.isErroneous IS FALSE) < {} AND p1.pyroID < p2.pyroID AND p1.isErroneous IS FALSE AND p2.isErroneous IS FALSE ORDER BY p1.isoID, p1.appliedRegion, p1.pyroID, p2.pyroID".format(maxReplicatesPerIsolateRegion)
	)
	# Isolates checked with:
	# select isoID, appliedRegion, count(*) from pyroprints where isErroneous is false group by isoID, appliedRegion having count(*) > 1 order by isoID, appliedRegion;
	
	regionPearsons = {region.name: [] for region in cfg.regions}
	cursor.execute(query)
	for (region, pearson) in cursor:
		regionPearsons[region].append(pearson)

	cursor.close()
	cnx.close()

	for region in cfg.regions:
		print("region {} has {} replicate pairs".format(region.name, len(regionPearsons[region.name])))

	return {regionNameLookup[region]: pearsons for region, pearsons in regionPearsons.items()}

def getReplicatePearsonsCacheFileName(cfg):
	return "replicatePearsons.pickle"

def loadReplicatePearsonsFromFile(cacheFileName):
	with open(cacheFileName, mode='r+b') as cacheFile:
		return pickle.load(cacheFile)

def loadReplicatePearsons(cfg):
	cacheFileName = getReplicatePearsonsCacheFileName(cfg)
	if os.path.isfile(cacheFileName):
		return loadReplicatePearsonsFromFile(cacheFileName)
	else:
		return loadReplicatePearsonsFromDB(cfg)


# def computeRegionPearsonMap(isolates, region):
# 	# regionsPearsonMap = {region: dict() for region in regions}
# 	pearsonMap = dict()
# 	for i, iso1 in enumerate(isolates):
# 		print("{}/{}".format(i, len(isolates)))
# 		for iso2 in isolates[i+1:]:
# 			# for region in regions:
# 			pearson = iso1.regionPearson(iso2, region)
# 			if iso1 < iso2:
# 				# regionsPearsonMap[region][iso1, iso2] = pearson
# 				pearsonMap[iso1, iso2] = pearson
# 			else:
# 				# regionsPearsonMap[region][iso2, iso1] = pearson
# 				pearsonMap[iso2, iso1] = pearson

# 	return pearsonMap

# def getPearsonMapCacheFileName(region, cfg):
# 	return "pearsons{}R{}.pickle".format(cfg.isolateSubsetSize, region)

# def loadPearsonMapFromFile(cacheFileName):
# 	with open(cacheFileName, mode='r+b') as cacheFile:
# 		pearsonMap = pickle.load(cacheFile)
# 	return pearsonMap

# def getPearsonMap(isolates, region, cfg):
# 	cacheFileName = getPearsonMapCacheFileName(region, cfg)
# 	if os.path.isfile(cacheFileName):
# 		return loadPearsonMapFromFile(cacheFileName)
# 	else:
# 		return computeRegionPearsonMap(isolates, region)


if __name__ == '__main__':
	cfg = config.loadConfig()
	assert cfg.isolateSubsetSize == "Shared"
	isolates = pyroprinting.loadIsolates(cfg)

	dbscanClusters = dbscan.getDBscanClusters(isolates, cfg)
	ohclustClusters = importClusters.getOHClustClusters()
	agglomerativeClusters = importClusters.getAgglomerativeClusters()
	replicatePearsons = loadReplicatePearsons(cfg)

	for name, clusters in [("DBSCAN", dbscanClusters), ("OHCLUST", ohclustClusters)]:
		print("\n\n{}\n".format(name))

		for region in cfg.regions:
			print("\n\t{}".format(region))
			print("\t\tdistributionSimilarity: {}".format(distributionSimilarity(clusters, agglomerativeClusters, region)))

		print("\n\tsimilarityIndexes: {}".format(similarityIndexes(isolates, clusters, agglomerativeClusters)))

	for name, clusters in [("DBSCAN", dbscanClusters), ("OHCLUST", ohclustClusters), ("AGGLOMERATIVE", agglomerativeClusters)]:
		print("\n\n{}\n".format(name))

		for region in cfg.regions:
			print("\n\t{}".format(region))

			# print(stats.beta.fit(replicatePearsons[region]))
			betaCorrectnessFunc = getBetaCorrectnessFunc(stats.beta(region.betaDistributionAlpha, region.betaDistributionBeta))
			distributionCorrectnessFunc = getDistributionCorrectnessFunc(replicatePearsons[region])
			# betaDistribution = stats.beta(region.betaDistributionAlpha, region.betaDistributionBeta)
			# pearsonMap = getPearsonMap(isolates, region, cfg)

			print("\t\tdistributionCorrectness: {}".format(distributionCorrectness(clusters, region, distributionCorrectnessFunc)))
			print("\t\tbetaDistributionCorrectness: {}".format(distributionCorrectness(clusters, region, betaCorrectnessFunc)))
			print("\t\tindividualClusterDistributionCorrectness: {}".format(individualClusterDistributionCorrectness(clusters, region, distributionCorrectnessFunc)))
			print("\t\tindividualClusterBetaDistributionCorrectness: {}".format(individualClusterDistributionCorrectness(clusters, region, betaCorrectnessFunc)))
			print("\t\tpairedClustersDistributionCorrectness: {}".format(pairedClustersDistributionCorrectness(clusters, region, distributionCorrectnessFunc)))
			print("\t\tpairedClustersBetaDistributionCorrectness: {}".format(pairedClustersDistributionCorrectness(clusters, region, betaCorrectnessFunc)))

			# del pearsonMap # the garbage collector should free this now

		print("\n\tthresholdCorrectness: {}".format(
			thresholdCorrectness(isolates, clusters, cfg.regions)))

