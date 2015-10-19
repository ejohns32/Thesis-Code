import pickle
import numpy
import itertools
import random
import math
import sys

fileName = sys.argv[1] if len(sys.argv) > 1 else "zScores.pickle"
with open(fileName, mode='r+b') as zScoresFile:
	zScores = pickle.load(zScoresFile)

isoIDs = list(zScores.keys())
regions = ['23-5', '16-23']
regionDispCounts = [93, 95]
thresholds = [0.995]

for threshold in thresholds:
	radii = [math.sqrt(2*dispCount*(1-threshold)) for dispCount in regionDispCounts]
	print(threshold, radii)

	correct = {isoID: set() for isoID in isoIDs}
	for i, p1 in enumerate(isoIDs):
		print("{}/{}".format(i, len(isoIDs)))
		for p2 in isoIDs[i:]:
			if all(numpy.linalg.norm(zScores[p1][region] - zScores[p2][region]) < radius for region, radius in zip(regions, radii)):
				correct[p1].add(p2)
				correct[p2].add(p1)

	print("size: {}".format(len(correct)))
	avgMatches = sum(len(matches) for _, matches in correct.items()) / len(correct)
	print("avgMatches: {}".format(avgMatches))

	if len(sys.argv) == 1:
		with open("correctR{}.pickle".format(threshold), mode='w+b') as correctFile:
			pickle.dump(correct, correctFile)
