from scipy import stats


def getClusterPearsons(pearsonMap, cluster):
	pearsons = []

	clusterList = list(cluster)
	for i, iso1 in enumerate(clusterList):
		for iso2 in clusterList[i+1:]:
			pearsons.append(pearsonMap[(iso1, iso2)])

	return pearsons

def getInterClusterPearsons(pearsonMap, clusters):
	pearsons = []

	for cluster in clusters:
		pearsons.extend(getClusterPearsons(pearsonMap, cluster))

	return pearsons

# call twice with an pearsonMap for each region
def distributionSimilarity(pearsonMap, clusters1, clusters2):
	ksStat, pValue = stats.ks_2samp(getInterClusterPearsons(pearsonMap, clusters1), getInterClusterPearsons(pearsonMap, clusters2))
	return pValue



#betaDistribution = stats.beta(distAlpha, distBeta)
def distributionCorrectness(pearsonMap, clusters, betaDistribution):
	ksStat, pValue = stats.kstest(getInterClusterPearsons(pearsonMap, clusters), betaDistribution.cdf)
	return pValue

def individualClusterDistributionCorrectness(pearsonMap, clusters, betaDistribution):
	# TODO: double check that this is statistically valid
	# null hypothesis: all individual clusters fit
	# reject if any individuals reject
	# want chance that none reject. given chance that each doesn't reject
	# so all pValues multiplied together

	# print()
	chanceNoneReject = 1.0
	for cluster in clusters:
		if len(cluster) > 1: # otherwise it wont have any pairwise pearsons
			ksStat, pValue = stats.kstest(getClusterPearsons(pearsonMap, cluster), betaDistribution.cdf)
			# print(pValue)
			chanceNoneReject *= pValue

	return chanceNoneReject

# seems to need large clusters to work
def pairedClustersDistributionCorrectness(pearsonMap, clusters, betaDistribution):
	# TODO: double check that this is statistically valid
	# null hypothesis: no pairs of clusters fit
	# reject if any pairs fit
	# want chance that all reject. given chance that each doesn't reject
	# so all inverse pValues multiplied together

	# print()
	chanceAllReject = 1.0
	for i, cluster1 in enumerate(clusters):
		for cluster2 in clusters[i+1:]:
			# each cluster will have at least 1 iso, so combined will have at least 2 iso, so at least one pairwise pearson
			ksStat, pValue = stats.kstest(getClusterPearsons(pearsonMap, cluster1 | cluster2), betaDistribution.cdf)
			# print(1-pValue)
			# print(pValue)
			# print(cluster1)
			# print(cluster2)
			chanceAllReject *= 1 - pValue

	return chanceAllReject





def getClusterMap(clusters):
	clustersMap = dict()

	for cluster in clusters:
		for isolate in cluster:
			assert isolate not in clustersMap
			clustersMap[isolate] = cluster

	return clustersMap

# call twice with an pearsonMap for each region
def thresholdCorrectness(isolates, pearsonMap, clusters, dsThreshold, ddThreshold):
	sameDS = diffDD = totDS = totDD = totSquishy = 0
	clusterMap = getClusterMap(clusters)

	for i, iso1 in enumerate(isolates):
		for iso2 in isolates[i+1:]:
			pearson = pearsonMap[(iso1, iso2)]

			if pearson >= dsThreshold:
				totDS += 1
				if clusterMap[iso1] == clusterMap[iso2]:
					sameDS += 1
			elif pearson <= ddThreshold:
				totDD += 1
				if clusterMap[iso1] != clusterMap[iso2]:
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

	return correctDS, correctDD, tcScore, squishyRatio



def getPairCounts(isolates, clusters1, clusters2):
	bothSame = 0
	oneOfEach = 0
	bothDiff = 0

	clusterMap1 = getClusterMap(clusters1)
	clusterMap2 = getClusterMap(clusters2)

	for i, iso1 in enumerate(isolates):
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





def getRegionsPearsonMap(isolates, regions):
	regionsPearsonMap = {region: dict() for region in regions}

	for i, iso1 in enumerate(isolates):
		for iso2 in isolates[i+1:]:
			for region, pearson in iso1.regionsPearson(iso2).items():
				regionsPearsonMap[region][(iso1, iso2)] = regionsPearsonMap[region][(iso2, iso1)] = pearson

	return regionsPearsonMap




if __name__ == '__main__':
	pass
