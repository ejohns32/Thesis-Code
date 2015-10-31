import json
import math
import pickle
import os

import pyroprinting

def loadConfig():
	with open("clusterConfig.json", mode='r') as clusterConfigJson:
		return Config(clusterConfigJson)

regionsPickleFileName = "regions.pickle"

class Config:
	def __init__(self, clusterConfigJson):
		decoded = json.load(clusterConfigJson)

		self.regions = [pyroprinting.Region(region["name"], region["dispCount"], i) for i, region in enumerate(decoded["regions"])]
		# if os.path.isfile(regionsPickleFileName):
		# 	with open(regionsPickleFileName, mode='r+b') as regionsPickleFile:
		# 		self.regions = pickle.load(regionsPickleFile)
		# 	assert(self.regions == configRegions)
		# else:
		# 	with open(regionsPickleFileName, mode='w+b') as regionsPickleFile:
		# 		pickle.dump(configRegions, regionsPickleFile)
		# 	self.regions = configRegions

		self.threshold = decoded["threshold"]
		self.radii = {region: math.sqrt(2*region.dispCount * (1 - self.threshold)) for region in self.regions}
		print(self.threshold, tuple((region.name, radius) for region, radius in self.radii.items()))
		self.minNeighbors = decoded["minNeighbors"]
		self.pointsPerLeaf = decoded["pointsPerLeaf"]
		self.dimsPerSplit = decoded["dimsPerSplit"]
		self.isolateSubsetSize = decoded["isolateSubsetSize"]
