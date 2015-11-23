import os.path
import pickle
import json
import math

import mysql.connector
from scipy import stats

import config
import pyroprinting
import dbscan
import importClusters


def getDistributionCorrectnessFunc(expectedDistribution):
	def func(actual):
		ksStat, pValue = stats.ks_2samp(actual, expectedDistribution)
		return ksStat, ksStat/(math.sqrt((len(actual) + len(expectedDistribution)) / (len(actual) * (len(expectedDistribution))))), pValue
		# return ksStat, pValue, len(actual), len(expectedDistribution)
	return func

# def getBetaCorrectnessFunc(expectedBeta):
# 	def func(actual):
# 		return stats.kstest(actual, expectedBeta.cdf)
# 	return func

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
	distribution1 = getCombinedIntraClusterPearsons(clusters1, region)
	distribution2 = getCombinedIntraClusterPearsons(clusters2, region)
	ksStat, pValue = stats.ks_2samp(distribution1, distribution2)
	# return pValue
	return ksStat, ksStat/(math.sqrt((len(distribution1) + len(distribution2)) / (len(distribution1) * (len(distribution2))))), pValue
	# return ksStat, pValue, len(distribution1), len(distribution2)



def distributionCorrectness(clusters, region, correctnessFunc):
	ksStat, weightedKS, pValue = correctnessFunc(getCombinedIntraClusterPearsons(clusters, region))
	# ksStat, pValue = stats.kstest(getCombinedIntraClusterPearsons(clusters, region), betaDistribution.cdf)
	# return pValue
	return ksStat, weightedKS, pValue

#seems to need small clusters properly fail to reject
def individualClusterDistributionCorrectness(clusters, region, correctnessFunc):
	# TODO: double check that this is statistically valid
	# null hypothesis: all individual clusters fit
	# reject if any individuals reject
	# want chance that none reject. given chance that each doesn't reject
	# so all pValues multiplied together

	# print()
	KSs = []
	wKSs = []
	pVs = []
	# totalCount = 0
	# chanceCount = 0
	# averageKS = 0.0
	# averageWeightedKS = 0.0
	# averageP = 0.0
	# chanceNoneReject = 1.0
	for cluster in clusters:
		if len(cluster) >= 8: # otherwise it wont have enough pairwise pearsons to be significant
			ksStat, weightedKS, pValue = correctnessFunc(getIntraClusterPearsons(cluster, region))
			# ksStat, pValue = stats.kstest(getIntraClusterPearsons(cluster, region), betaDistribution.cdf)
			# print(pValue)
			# chanceNoneReject *= pValue
			# averageP += pValue
			# averageKS += ksStat
			# averageWeightedKS += weightedKS
			# totalCount += 1
			# if pValue >= .1:
			# 	chanceCount += 1
			KSs.append(ksStat)
			wKSs.append(weightedKS)
			pVs.append(pValue)

	# return chanceCount/totalCount, averageP/totalCount, chanceNoneReject
	# return averageKS/totalCount, averageWeightedKS/totalCount, averageP/totalCount
	return stats.describe(KSs), stats.describe(wKSs), stats.describe(pVs)

# seems to need large clusters to properly reject
def pairedClustersDistributionCorrectness(clusters, region, correctnessFunc):
	# TODO: double check that this is statistically valid
	# null hypothesis: no pairs of clusters fit
	# reject if any pairs fit
	# want chance that all reject. given chance that each doesn't reject
	# so all inverse pValues multiplied together

	intraPearsons = [getIntraClusterPearsons(cluster, region) for cluster in clusters]

	# print()
	KSs = []
	wKSs = []
	pVs = []
	# totalCount = 0
	# chanceCount = 0
	# averageP = 0.0
	# averageKS = 0.0
	# chanceAllReject = 1.0
	for i, cluster1 in enumerate(clusters):
		# print("{}/{}".format(i+1, len(clusters)))
		for  j, cluster2 in enumerate(clusters[i+1:]):
			if len(cluster1) + len(cluster2) >= 8:
				pearsons = intraPearsons[i] + intraPearsons[i+1+j] + getInterClusterPearsons(cluster1, cluster2, region)
				# ksStat, pValue = stats.kstest(pearsons, betaDistribution.cdf)
				ksStat, weightedKS, pValue = correctnessFunc(pearsons)
				# print(1-pValue)
				# print(pValue)
				# print(cluster1)
				# print(cluster2)
				# chanceAllReject *= 1 - pValue
				# averageP += 1 - pValue
				# averageKS += ksStat
				# totalCount += 1
				# if pValue <= .1:
				# 	chanceCount += 1
				KSs.append(ksStat)
				wKSs.append(weightedKS)
				pVs.append(pValue)

	# return chanceCount/totalCount, averageP/totalCount, chanceAllReject
	# return chanceCount/totalCount, averageKS/totalCount
	return stats.describe(KSs), stats.describe(wKSs), stats.describe(pVs)





def getClusterMap(clusters):
	clustersMap = dict()
	multipleCount = 0

	for i, cluster in enumerate(clusters):
		for isolate in cluster:
			# assert isolate not in clustersMap
			if isolate in clustersMap:
				clustersMap[isolate].add(i)
				multipleCount += 1
				# print("\n\tisolate {} in multiple clusters".format(isolate))
				# print(clustersMap[isolate])
				# print(cluster)
			clustersMap[isolate] = {i}

	print("\t{} isolates in multiple clusters".format(multipleCount))

	return clustersMap

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
				if iso1 in clusterMap and iso2 in clusterMap and len(clusterMap[iso1] & clusterMap[iso2]) > 0:
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
	tcScore = (sameDS+diffDD)/(totDS+totDD)
	squishyRatio = totSquishy/(totSquishy+totDD+totDS)

	# print("correctDS = {}".format(correctDS))
	# print("correctDD = {}".format(correctDD))
	# print("tcScore = {}".format(tcScore))
	# print("squishyRatio = {}".format(squishyRatio))
	# print("totDS: {}, totSquishy: {}, totDD: {}".format(totDS, totSquishy, totDD))

	return correctDS, correctDD, (correctDS+correctDD)/2, tcScore, squishyRatio


# TODO: maybe add option of ignoring noise / clusters of size 1 / isolates not in clusterMap
def getPairCounts(isolates, clusters1, clusters2):
	bothSame = 0
	oneOfEach = 0
	bothDiff = 0

	clusterMap1 = getClusterMap(clusters1)
	clusterMap2 = getClusterMap(clusters2)

	for i, iso1 in enumerate(isolates):
		# print("{}/{}".format(i+1, len(isolates)))
		for iso2 in isolates[i+1:]:
			sameClust1 = iso1 in clusterMap1 and iso2 in clusterMap1 and len(clusterMap1[iso1] & clusterMap1[iso2]) > 0
			sameClust2 = iso1 in clusterMap2 and iso2 in clusterMap2 and len(clusterMap2[iso1] & clusterMap2[iso2]) > 0
			
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

	return jaccardIndex, anotherIndex, (jaccardIndex+anotherIndex)/2, randIndex 



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


def clustersStats(isolates, clusters):
	clusterMap = getClusterMap(clusters)
	clusterSizes = [len(cluster) for cluster in clusters]
	# sizeRange = min(clusterSizes), max(clusterSizes)
	count, minmax, mean, var, _, _ = stats.describe(clusterSizes)
	sizeHistogram = (
		len([size for size in clusterSizes if size < 4]),
		len([size for size in clusterSizes if size >= 4 and size < 8]),
		len([size for size in clusterSizes if size >= 8 and size < 16]),
		len([size for size in clusterSizes if size >= 16 and size < 32]),
		len([size for size in clusterSizes if size >= 32 and size < 64]),
		len([size for size in clusterSizes if size >= 64 and size < 128]),
		len([size for size in clusterSizes if size >= 128])
	)

	noise = [isolate for isolate in isolates if isolate not in clusterMap]

	return count, minmax, (mean, math.sqrt(var)), len(noise), sizeHistogram

def filterSingletonClusters(clusters):
	filtered = [cluster for cluster in clusters if len(cluster) > 1]
	print("{}/{}".format(len(filtered), len(clusters)))
	return filtered

if __name__ == '__main__':
	cfg = config.loadConfig()
	assert cfg.isolateSubsetSize == "Shared"
	isolates = pyroprinting.loadIsolates(cfg)

	# dbscanClusters = dbscan.getDBscanClusters(isolates, cfg)
	dbscan1Clusters = dbscan.loadDBscanClustersFromFile("dbscanShared_0.995_0.995_1.pickle")
	dbscan3Clusters = dbscan.loadDBscanClustersFromFile("dbscanShared_0.995_0.995_3.pickle")
	ohclust99Clusters = filterSingletonClusters(importClusters.getOHClustClusters(99))
	# ohclust995Clusters = filterSingletonClusters(importClusters.getOHClustClusters(995))
	agglomerativeClusters = filterSingletonClusters(importClusters.getAgglomerativeClusters())
	replicatePearsons = loadReplicatePearsons(cfg)

	db1Pair = ("DBSCAN 1", dbscan1Clusters)
	db3Pair = ("DBSCAN 3", dbscan3Clusters)
	oh99Pair = ("OHClust 99", ohclust99Clusters)
	# oh995Pair = ("OHClust 995", ohclust995Clusters)
	aggPair = ("AGGLOMERATIVE", agglomerativeClusters)

	for name, clusters in [db1Pair, db3Pair, oh99Pair, aggPair]:
		print("\n\n{}".format(name))
		print("\tclusterStats: {}".format(clustersStats(isolates, clusters)))

	# for name, clusters in [db1Pair, db3Pair, oh99Pair]:
	# 	print("\n\n{}".format(name))
	# 	print(len(clusters))

	# 	for region in cfg.regions:
	# 		print("\t{}".format(region))
	# 		print("\t\tdistributionSimilarity: {}".format(distributionSimilarity(clusters, agglomerativeClusters, region)))

	# 	print("\n\tsimilarityIndexes: {}".format(similarityIndexes(isolates, clusters, agglomerativeClusters)))

	# for name, clusters in [db1Pair, db3Pair, oh99Pair, aggPair]:
	# 	print("\n\n{}".format(name))
	# 	print(len(clusters))

	# 	for region in cfg.regions:
	# 		print("\t{}".format(region))

	# 		# print(stats.beta.fit(replicatePearsons[region]))
	# 		# betaCorrectnessFunc = getBetaCorrectnessFunc(stats.beta(region.betaDistributionAlpha, region.betaDistributionBeta))
	# 		distributionCorrectnessFunc = getDistributionCorrectnessFunc(replicatePearsons[region])
	# 		# betaDistribution = stats.beta(region.betaDistributionAlpha, region.betaDistributionBeta)
	# 		# pearsonMap = getPearsonMap(isolates, region, cfg)

	# 		print("\t\tdistributionCorrectness: {}".format(distributionCorrectness(clusters, region, distributionCorrectnessFunc)))
	# 		# print("\t\tbetaDistributionCorrectness: {}".format(distributionCorrectness(clusters, region, betaCorrectnessFunc)))
	# 		print("\t\tindividualClusterDistributionCorrectness: {}".format(individualClusterDistributionCorrectness(clusters, region, distributionCorrectnessFunc)))
	# 		# print("\t\tindividualClusterBetaDistributionCorrectness: {}".format(individualClusterDistributionCorrectness(clusters, region, betaCorrectnessFunc)))
	# 		print("\t\tpairedClustersDistributionCorrectness: {}".format(pairedClustersDistributionCorrectness(clusters, region, distributionCorrectnessFunc)))
	# 		# print("\t\tpairedClustersBetaDistributionCorrectness: {}".format(pairedClustersDistributionCorrectness(clusters, region, betaCorrectnessFunc)))

	# 		# del pearsonMap # the garbage collector should free this now

	# 	print("\n\tthresholdCorrectness: {}".format(
	# 		thresholdCorrectness(isolates, clusters, cfg.regions)))

