import unittest
import json
import math

import mysql.connector
import numpy

import pyroprinting


class TestIsolate(unittest.TestCase):
	def setUp(self):
		self.region1 = pyroprinting.Region("region1", 6, 0, 0, 0, 0)
		self.region2 = pyroprinting.Region("region2", 4, 0, 0, 0, 0)

		# Note: not all arrays are valid zScore arrays. Pearson from zscores requires such a valid array.
		self.iso1 = pyroprinting.Isolate("iso1", {self.region1: numpy.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0]), self.region2: numpy.array([1.0, -1.0, 1.0, -1.0])})
		self.iso2 = pyroprinting.Isolate("iso2", {self.region1: numpy.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0]), self.region2: numpy.array([-1.0, 1.0, -1.0, 1.0])})

		# with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
		# 	mysqlConfig = json.load(mysqlConfigJson)

		# self.cnx = mysql.connector.connect(**mysqlConfig)

	# def tearDown(self):
	# 	self.cnx.close()

	def testDistAndPearsonConversions(self):
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 95), 95))
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 93), 93))
		self.assertNotAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 95), 93))
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(1.3784, 95), places=2)
		self.assertAlmostEqual(1.3784, pyroprinting.distFromPearson(.99, 95), places=4)

	# def assertRegionsFloatAlmostEqual(self, expRegionsFloat, actRegionsFloat):
	# 	for (expRegion, expFloat), (actRegion, actFloat) in zip(expRegionsFloat, actRegionsFloat):
	# 		self.assertEqual(expRegion, actRegion)
	# 		self.assertAlmostEqual(expFloat, actFloat)

	def testRegionsDistAndPearson(self):
		self.assertAlmostEqual((1 + -1 + 1 + 1 + -1 + 1)/6, self.iso1.regionPearson(self.iso2, self.region1))
		self.assertAlmostEqual((-1 + -1 + -1 + -1)/4, self.iso1.regionPearson(self.iso2, self.region2))

		self.assertAlmostEqual(math.sqrt(0**2 + 2**2 + 0**2 + 0**2 + 2**2 + 0**2), self.iso1.regionDist(self.iso2, self.region1))
		self.assertAlmostEqual(math.sqrt(2**2 + 2**2 + 2**2 + 2**2), self.iso1.regionDist(self.iso2, self.region2))

		for region in [self.region1, self.region2]:
			self.assertAlmostEqual(self.iso1.regionPearson(self.iso2, region), pyroprinting.pearsonFromDist(self.iso1.regionDist(self.iso2, region), region.dispCount)) 
			self.assertAlmostEqual(self.iso1.regionPearson(self.iso2, region), pyroprinting.pearsonFromDist(self.iso1.regionDist(self.iso2, region), region.dispCount)) 
			self.assertAlmostEqual(self.iso1.regionDist(self.iso2, region), pyroprinting.distFromPearson(self.iso1.regionPearson(self.iso2, region), region.dispCount)) 

			self.assertAlmostEqual(self.iso1.regionPearson(self.iso2, region), self.iso2.regionPearson(self.iso1, region))
			self.assertAlmostEqual(1, self.iso1.regionPearson(self.iso1, region))

			self.assertAlmostEqual(self.iso1.regionDist(self.iso2, region), self.iso2.regionDist(self.iso1, region))
			self.assertAlmostEqual(0, self.iso1.regionDist(self.iso1, region))



		# self.assertRegionsFloatAlmostEqual(self.iso1.regionsPearson(self.iso2), [(region, pyroprinting.pearsonFromDist(dist, region.dispCount)) for region, dist in self.iso1.regionsDist(self.iso2)])
		# self.assertAlmostEqual(self.iso1.regionsPearson(self.iso2), [(region, pyroprinting.pearsonFromDist(dist, region.dispCount)) for region, dist in self.iso1.regionsDist(self.iso2)])
		# self.assertRegionsFloatAlmostEqual(self.iso1.regionsDist(self.iso2), [(region, pyroprinting.distFromPearson(dist, region.dispCount)) for region, dist in self.iso1.regionsPearson(self.iso2)])

		# self.assertRegionsFloatAlmostEqual(self.iso1.regionsPearson(self.iso2), self.iso2.regionsPearson(self.iso1))
		# self.assertRegionsFloatAlmostEqual([(self.region1, (1 + -1 + 1 + 1 + -1 + 1)/6), (self.region2, (-1 + -1 + -1 + -1)/4)], self.iso1.regionsPearson(self.iso2))
		# self.assertRegionsFloatAlmostEqual([(self.region1, 1), (self.region2, 1)], self.iso1.regionsPearson(self.iso1))

		# self.assertRegionsFloatAlmostEqual(self.iso1.regionsDist(self.iso2), self.iso2.regionsDist(self.iso1))
		# self.assertRegionsFloatAlmostEqual([(self.region1, math.sqrt(0**2 + 2**2 + 0**2 + 0**2 + 2**2 + 0**2)), (self.region2, math.sqrt(2**2 + 2**2 + 2**2 + 2**2))], self.iso1.regionsDist(self.iso2))
		# self.assertRegionsFloatAlmostEqual([(self.region1, 0), (self.region2, 0)], self.iso1.regionsDist(self.iso1))



		# query = ("SELECT PearsonMatch(p1.pyroID, p2.pyroID, 95) FROM Pyroprints p1 JOIN Pyroprints p2 WHERE p1.appliedRegion = '16-23' AND p2.appliedRegion = '16-23' AND p1.pyroID < p2.pyroID AND p1.pyroID >= 167 AND p2.pyroID <= 170 ORDER BY p1.pyroID, p2.pyroID;")

		# cursor = self.cnx.cursor()
		# cursor.execute(query)

		# for (pearson) in cursor:
		# 	pass

		# cursor.close()

	def testIsWithinRadiiOf(self):
		# test if both regions are inside or outside radius
		self.assertFalse(self.iso1.isWithinRadiiOf(self.iso2, {self.region1: 2, self.region2: 3}))
		self.assertTrue(self.iso1.isWithinRadiiOf(self.iso1, {self.region1: 2, self.region2: 3}))
		self.assertTrue(self.iso2.isWithinRadiiOf(self.iso2, {self.region1: 2, self.region2: 3}))

		# test if one region is inside radius but not other
		self.assertFalse(self.iso1.isWithinRadiiOf(self.iso2, {self.region1: 20, self.region2: 3})) 
		self.assertFalse(self.iso1.isWithinRadiiOf(self.iso2, {self.region1: 2, self.region2: 30}))
		



if __name__ == '__main__':
    unittest.main()