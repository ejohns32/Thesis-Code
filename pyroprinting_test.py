import unittest
import json
import math

import mysql.connector
import numpy

import pyroprinting


class TestIsolate(unittest.TestCase):
	# def setUp(self):
	# 	with open("mysqlConfig.json", mode='r') as mysqlConfigJson:
	# 		mysqlConfig = json.load(mysqlConfigJson)

	# 	self.cnx = mysql.connector.connect(**mysqlConfig)

	# def tearDown(self):
	# 	self.cnx.close()

	def testDistAndPearsonConversions(self):
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 95), 95))
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 93), 93))
		self.assertNotAlmostEqual(.99, pyroprinting.pearsonFromDist(pyroprinting.distFromPearson(.99, 95), 93))
		self.assertAlmostEqual(.99, pyroprinting.pearsonFromDist(1.3784, 95), places=2)
		self.assertAlmostEqual(1.3784, pyroprinting.distFromPearson(.99, 95), places=4)

	def assertRegionsFloatAlmostEqual(self, expRegionsFloat, actRegionsFloat):
		for (expRegion, expFloat), (actRegion, actFloat) in zip(expRegionsFloat, actRegionsFloat):
			self.assertEqual(expRegion, actRegion)
			self.assertAlmostEqual(expFloat, actFloat)

	def testRegionsDistAndPearson(self):
		region1 = pyroprinting.Region("region1", 6)
		region2 = pyroprinting.Region("region2", 4)

		# Note: not all arrays are valid zScore arrays. Pearson from zscores requires such a valid array.
		iso1 = pyroprinting.Isolate("iso1", {region1: numpy.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0]), region2: numpy.array([1.0, -1.0, 1.0, -1.0])})
		iso2 = pyroprinting.Isolate("iso2", {region1: numpy.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0]), region2: numpy.array([-1.0, 1.0, -1.0, 1.0])})

		self.assertRegionsFloatAlmostEqual(iso1.regionsPearson(iso2), [(region, pyroprinting.pearsonFromDist(dist, region.dispCount)) for region, dist in iso1.regionsDist(iso2)])
		self.assertRegionsFloatAlmostEqual(iso1.regionsDist(iso2), [(region, pyroprinting.distFromPearson(dist, region.dispCount)) for region, dist in iso1.regionsPearson(iso2)])

		self.assertRegionsFloatAlmostEqual(iso1.regionsPearson(iso2), iso2.regionsPearson(iso1))
		self.assertRegionsFloatAlmostEqual([(region1, (1 + -1 + 1 + 1 + -1 + 1)/6), (region2, (-1 + -1 + -1 + -1)/4)], iso1.regionsPearson(iso2))
		self.assertRegionsFloatAlmostEqual([(region1, 1), (region2, 1)], iso1.regionsPearson(iso1))

		self.assertRegionsFloatAlmostEqual(iso1.regionsDist(iso2), iso2.regionsDist(iso1))
		self.assertRegionsFloatAlmostEqual([(region1, math.sqrt(0**2 + 2**2 + 0**2 + 0**2 + 2**2 + 0**2)), (region2, math.sqrt(2**2 + 2**2 + 2**2 + 2**2))], iso1.regionsDist(iso2))
		self.assertRegionsFloatAlmostEqual([(region1, 0), (region2, 0)], iso1.regionsDist(iso1))


		# query = ("SELECT PearsonMatch(p1.pyroID, p2.pyroID, 95) FROM Pyroprints p1 JOIN Pyroprints p2 WHERE p1.appliedRegion = '16-23' AND p2.appliedRegion = '16-23' AND p1.pyroID < p2.pyroID AND p1.pyroID >= 167 AND p2.pyroID <= 170 ORDER BY p1.pyroID, p2.pyroID;")

		# cursor = self.cnx.cursor()
		# cursor.execute(query)

		# for (pearson) in cursor:
		# 	pass

		# cursor.close()

	def testIsWithinRadiiOf(self):
		region1 = pyroprinting.Region("region1", 6)
		region2 = pyroprinting.Region("region2", 4)

		# Note: not all arrays are valid zScore arrays.
		iso1 = pyroprinting.Isolate("iso1", {region1: numpy.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0]), region2: numpy.array([1.0, -1.0, 1.0, -1.0])})
		iso2 = pyroprinting.Isolate("iso2", {region1: numpy.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0]), region2: numpy.array([-1.0, 1.0, -1.0, 1.0])})

		# test if both regions are inside or outside radius
		self.assertFalse(iso1.isWithinRadiiOf(iso2, {region1: 2, region2: 3}))
		self.assertTrue(iso1.isWithinRadiiOf(iso1, {region1: 2, region2: 3}))
		self.assertTrue(iso2.isWithinRadiiOf(iso2, {region1: 2, region2: 3}))

		# test if one region is inside radius but not other
		self.assertFalse(iso1.isWithinRadiiOf(iso2, {region1: 20, region2: 3})) 
		self.assertFalse(iso1.isWithinRadiiOf(iso2, {region1: 2, region2: 30}))
		



if __name__ == '__main__':
    unittest.main()