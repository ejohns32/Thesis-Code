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

		self.regions = [pyroprinting.Region.fromDecodedJSON(decodedRegion) for decodedRegion in decoded["regions"]]

		# self.threshold = decoded["threshold"]
		# self.radii = {region: math.sqrt(2*region.dispCount * (1 - self.threshold)) for region in self.regions}
		# self.radii = {region: pyroprinting.distFromPearson(self.threshold, region.dispCount) for region in self.regions}
		# print(self.threshold, tuple((region.name, radius) for region, radius in self.radii.items()))

		self.radii = {region: pyroprinting.distFromPearson(region.pSimThresholdAlpha, region.dispCount) for region in self.regions}
		# print(tuple((region.name, region.pSimThresholdAlpha, radius) for region, radius in self.radii.items()))

		self.minNeighbors = decoded["minNeighbors"]
		self.pointsPerLeaf = decoded["pointsPerLeaf"]
		self.dimsPerSplit = decoded["dimsPerSplit"]
		self.isolateSubsetSize = decoded["isolateSubsetSize"]
