import unittest
import random

from scipy import stats
import numpy

import clusterEval

class MockIsolate:
	def __init__(self, mockID):
		self.mockID = mockID

	def regionPearson(self, other, region):
		if self.mockID < other.mockID:
			return MockIsolate.pearsonMap[region][(self.mockID, other.mockID)]
		else:
			return MockIsolate.pearsonMap[region][(other.mockID, self.mockID)]

	def __repr(self):
		return self.mockID

	def __eq__(self, other):
		return self.mockID == other.mockID

	def __ne__(self, other):
		return not self.__eq__(other)

	def __hash__(self):
		return self.mockID

class MockRegion:
	def __init__(self, pSimThresholdAlpha, pSimThresholdBeta):
		self.pSimThresholdAlpha = pSimThresholdAlpha
		self.pSimThresholdBeta = pSimThresholdBeta



class TestClusterDistributions(unittest.TestCase):
	def setUp(self):
		self.mockRegion = MockRegion(.995, .99)
		MockIsolate.pearsonMap = {self.mockRegion: {
			(0, 1): .990, (0, 2): .991, (0, 3): .992, (0, 4): .993,
			              (1, 2): .994, (1, 3): .995, (1, 4): .996,
			                            (2, 3): .997, (2, 4): .998,
			                                          (3, 4): .999}
		}

		self.clusters = [{MockIsolate(0), MockIsolate(4)}, {MockIsolate(1), MockIsolate(2), MockIsolate(3)}]

	def testGetIntraClusterPearsons(self):
		self.assertEqual(set(clusterEval.getIntraClusterPearsons(self.clusters[0], self.mockRegion)), {.993})
		self.assertEqual(set(clusterEval.getIntraClusterPearsons(self.clusters[1], self.mockRegion)), {.994, .995, .997})

	def testGetInterClusterPearsons(self):
		self.assertEqual(set(clusterEval.getInterClusterPearsons(self.clusters[0], self.clusters[1], self.mockRegion)), {.990, .991, .992, .996, .998, .999})

	def testGetCombinedIntraClusterPearsons(self):
		self.assertEqual(set(clusterEval.getCombinedIntraClusterPearsons(self.clusters, self.mockRegion)), {.993, .994, .995, .997})

def getMockCorrectnessFunc():
	def func():
		return
	return

class TestDistributionMatch(unittest.TestCase):
	def setUp(self):
		self.betaDistribution = stats.beta(668.1768, 1.532336)

		strainsStartIndex = 0
		strainsSwitchIndex = 40
		strainsEndIndex = 60
		# self.clusters1 = [{1, 2, 3, 4}, {5, 6}]
		# self.clusters2 = [{1, 3}, {2, 4}, {5}, {6}]
		# self.clusters3 = [{1, 3, 5}, {2, 4, 6}]
		self.clusters1 = [set(MockIsolate(i) for i in range(strainsStartIndex, strainsSwitchIndex)),
		                  set(MockIsolate(i) for i in range(strainsSwitchIndex, strainsEndIndex))]
		self.clusters2 = [set(MockIsolate(i) for i in range(strainsStartIndex, strainsSwitchIndex, 2)),
		                  set(MockIsolate(i) for i in range(strainsStartIndex + 1, strainsSwitchIndex, 2)),
		                 # [set(range(strainsStartIndex, strainsSwitchIndex//2)), set(range(strainsSwitchIndex//2, strainsSwitchIndex)),
		                  set(MockIsolate(i) for i in range(strainsSwitchIndex, strainsEndIndex, 2)),
		                  set(MockIsolate(i) for i in range(strainsSwitchIndex + 1, strainsEndIndex, 2))]
		self.clusters3 = [self.clusters2[0] | self.clusters2[2], self.clusters2[1] | self.clusters2[3]]

		# print(betaValues)
		# random.shuffle(betaValues)
		pearsonMap = dict()

		# 1-4 match each other and 5-6 match each other
		betaValues = list(self.betaDistribution.ppf(numpy.linspace(0, 1, num=(strainsSwitchIndex)**2 / 2 + 2)))[1:-1]
		for i in range(strainsStartIndex, strainsSwitchIndex):
			for j in range(i+1, strainsSwitchIndex):
				pearsonMap[i, j] = betaValues.pop()
			for j in range(strainsSwitchIndex, strainsEndIndex):
				pearsonMap[i, j] = random.uniform(0.9, 0.99)

		betaValues = list(self.betaDistribution.ppf(numpy.linspace(0, 1, num=(strainsEndIndex - strainsSwitchIndex)**2 / 2 + 2)))[1:-1]
		for i in range(strainsSwitchIndex, strainsEndIndex):
			for j in range(i+1, strainsEndIndex):
				pearsonMap[i, j] = betaValues.pop()

		self.mockRegion = MockRegion(.995, .99)
		MockIsolate.pearsonMap = {self.mockRegion: pearsonMap}

		# self.pearsonMap = {(1,2): .996, (1,3): .997, (1,4): .998, (1,5): .980, (1,6): .900,
		#                                 (2,3): .997, (2,4): .998, (2,5): .970, (2,6): .910,
		#                                              (3,4): .998, (3,5): .960, (3,6): .920,
		#                                                           (4,5): .950, (4,6): .930,
		#                                                                        (5,6): .999}
		# for (first, second), pearson in list(self.pearsonMap.items()):
		# 	self.pearsonMap[second, first] = pearson

	def testBetaDistribution(self):
		# confirm the values found in Diana's paper
		self.assertAlmostEqual(self.betaDistribution.cdf(0.9953), .10, places=2)
		self.assertAlmostEqual(self.betaDistribution.cdf(0.9941), .05, places=3) # 3 places, because the first 0 doesn't count
		self.assertAlmostEqual(self.betaDistribution.cdf(0.9915), .01, places=3) # 3 places, because the first 0 doesn't count

	def testDistributionSimilarity(self):
		self.assertAlmostEqual(clusterEval.distributionSimilarity(self.clusters1, self.clusters1, self.mockRegion), 1.0)
		# should not reject null hypothesis
		self.assertGreater(clusterEval.distributionSimilarity(self.clusters1, self.clusters2, self.mockRegion), .10)
		# should reject null hypothesis
		self.assertLess(clusterEval.distributionSimilarity(self.clusters1, self.clusters3, self.mockRegion), .10)

	def testDistributionCorrectness(self):
		# should not reject null hypothesis
		self.assertGreater(clusterEval.distributionCorrectness(self.clusters1, self.mockRegion, self.betaDistribution), .10)
		# should not reject null hypothesis
		self.assertGreater(clusterEval.distributionCorrectness(self.clusters2, self.mockRegion, self.betaDistribution), .10)
		# should reject null hypothesis
		self.assertLess(clusterEval.distributionCorrectness(self.clusters3, self.mockRegion, self.betaDistribution), .10)

	def testIndividualClusterDistributionCorrectness(self):
		# should not reject null hypothesis
		correctRatio, averageP, completeP = clusterEval.individualClusterDistributionCorrectness(self.clusters1, self.mockRegion, self.betaDistribution)
		self.assertGreater(correctRatio, .90)
		self.assertGreater(averageP, .50)
		self.assertGreater(completeP, .10)
		# should not reject null hypothesis
		correctRatio, averageP, completeP = clusterEval.individualClusterDistributionCorrectness(self.clusters2, self.mockRegion, self.betaDistribution)
		self.assertGreater(correctRatio, .90)
		self.assertGreater(averageP, .50)
		self.assertGreater(completeP, .10)
		# should reject null hypothesis
		correctRatio, averageP, completeP = clusterEval.individualClusterDistributionCorrectness(self.clusters3, self.mockRegion, self.betaDistribution)
		self.assertLess(correctRatio, .10)
		self.assertLess(averageP, .50)
		self.assertLess(completeP, .10)

	def testPairedClustersDistributionCorrectness(self):
		# should not reject null hypothesis
		correctRatio, averageP, completeP = clusterEval.pairedClustersDistributionCorrectness(self.clusters1, self.mockRegion, self.betaDistribution)
		self.assertGreater(correctRatio, .90)
		self.assertGreater(averageP, .50)
		self.assertGreater(completeP, .10)
		# should reject null hypothesis, but its not very strong
		correctRatio, averageP, completeP = clusterEval.pairedClustersDistributionCorrectness(self.clusters2, self.mockRegion, self.betaDistribution)
		self.assertLess(correctRatio, .75)
		self.assertLess(completeP, .25)
		# should not reject null hypothesis
		correctRatio, averageP, completeP = clusterEval.pairedClustersDistributionCorrectness(self.clusters3, self.mockRegion, self.betaDistribution)
		self.assertGreater(correctRatio, .90)
		self.assertGreater(averageP, .50)
		self.assertGreater(completeP, .10)



class TestThresholdCorrectness(unittest.TestCase):
	def testThresholdCorrectness(self):
		dsThreshold = .995
		ddThreshold = .99
		regions = [MockRegion(.995, .99), MockRegion(.995, .99)]

		MockIsolate.pearsonMap = dict()
		# totSquishy = 3, sameDS = 3, diffDS = 1, sameDD = 1, diffDD = 2
		MockIsolate.pearsonMap[regions[0]] = {
			(0, 1): .993, (0, 2): .999, (0, 3): .993, (0, 4): .999,
			              (1, 2): .999, (1, 3): .999, (1, 4): .999,
			                            (2, 3): .987, (2, 4): .987,
			                                          (3, 4): .987}

		# totSquishy = 3, sameDS = 3, diffDS = 1, sameDD = 0, diffDD = 3
		MockIsolate.pearsonMap[regions[1]] = {
			(0, 1): .999, (0, 2): .999, (0, 3): .993, (0, 4): .999,
			              (1, 2): .993, (1, 3): .999, (1, 4): .987,
			                            (2, 3): .993, (2, 4): .987,
			                                          (3, 4): .987}
		# for (first, second), pearson in list(pearsonMap.items()):
		# 	pearsonMap[second, first] = pearson

		# but combined:
		totSquishy = 3
		sameDS = 2
		diffDS = 1
		sameDD = 1
		diffDD = 3
		clusters = [{MockIsolate(0), MockIsolate(4)}, {MockIsolate(1), MockIsolate(2), MockIsolate(3)}]

		correctDS, correctDD, squishyRatio = clusterEval.thresholdCorrectness([MockIsolate(i) for i in range(5)], clusters, regions)
		self.assertAlmostEqual(correctDS, sameDS/(sameDS+diffDS))
		self.assertAlmostEqual(correctDD, diffDD/(sameDD+diffDD))
		# self.assertAlmostEqual(tcScore, (3+2)/(4+3))
		self.assertAlmostEqual(squishyRatio, totSquishy/(totSquishy+sameDS+diffDS+sameDD+diffDD))





class TestSimilarityIndexes(unittest.TestCase):
	def setUp(self):
		# self.clusters1 = [{0, 1, 2, 3, 4}, {5, 6, 7, 8, 9}]
		# self.clusters2 = [{0, 1, 2, 3, 5}, {4, 6, 7, 8, 9}]
		# self.clusters3 = [{0, 1, 2, 3, 4}, {5, 6, 7}, {8, 9}]
		# self.clusters3 = [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}]
		# self.clusters4 = [{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}]

		# with 5 points, there are 10 pairwise comparisons
		self.points = range(5)
		self.clusters1 = [{0, 1, 2}, {3, 4}]
		self.clusters2 = [{0, 1}, {2, 3, 4}]
		self.clusters3 = [{0, 1, 2, 3, 4}]
		self.clusters4 = [{0}, {1}, {2}, {3}, {4}]

	def testGetPairCounts(self):
		bothSame, oneOfEeach, bothDiff = clusterEval.getPairCounts(self.points, self.clusters1, self.clusters1)
		self.assertEqual(bothSame, 4)
		self.assertEqual(oneOfEeach, 0)
		self.assertEqual(bothDiff, 6)

		bothSame, oneOfEeach, bothDiff = clusterEval.getPairCounts(self.points, self.clusters1, self.clusters2)
		self.assertEqual(bothSame, 2)
		self.assertEqual(oneOfEeach, 4)
		self.assertEqual(bothDiff, 4)

		bothSame, oneOfEeach, bothDiff = clusterEval.getPairCounts(self.points, self.clusters1, self.clusters3)
		self.assertEqual(bothSame, 4)
		self.assertEqual(oneOfEeach, 6)
		self.assertEqual(bothDiff, 0)

		bothSame, oneOfEeach, bothDiff = clusterEval.getPairCounts(self.points, self.clusters1, self.clusters4)
		self.assertEqual(bothSame, 0)
		self.assertEqual(oneOfEeach, 4)
		self.assertEqual(bothDiff, 6)

	def testSimilarityIndexes(self):
		rand, jaccard, another = clusterEval.similarityIndexes(self.points, self.clusters1, self.clusters1)
		self.assertAlmostEqual(rand, 1)
		self.assertAlmostEqual(jaccard, 1)
		self.assertAlmostEqual(another, 1)

		rand, jaccard, another = clusterEval.similarityIndexes(self.points, self.clusters1, self.clusters2)
		self.assertAlmostEqual(rand, 6/10)
		self.assertAlmostEqual(jaccard, 2/6)
		self.assertAlmostEqual(another, 4/8)

		rand, jaccard, another = clusterEval.similarityIndexes(self.points, self.clusters1, self.clusters3)
		self.assertAlmostEqual(rand, 4/10)
		self.assertAlmostEqual(jaccard, 4/10)
		self.assertAlmostEqual(another, 0/6)

		rand, jaccard, another = clusterEval.similarityIndexes(self.points, self.clusters1, self.clusters4)
		self.assertAlmostEqual(rand, 6/10)
		self.assertAlmostEqual(jaccard, 0/4)
		self.assertAlmostEqual(another, 6/10)






if __name__ == '__main__':
    unittest.main()