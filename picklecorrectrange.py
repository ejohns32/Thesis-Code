import pickle
import numpy
import itertools
import random

with open("zScores.pickle", mode='r+b') as zScoresFile:
	zScores = pickle.load(zScoresFile)

isoIDs = list(zScores.keys())
regions = ['23-5', '16-23']
radii = [0.964, 1.37]
for radius in radii:
	correct = {isoID: set() for isoID in isoIDs}
	for i, p1 in enumerate(isoIDs):
		print("{}/{}".format(i, len(isoIDs)))
		for p2 in isoIDs[i:]:
			if all(numpy.linalg.norm(zScores[p1][region] - zScores[p2][region]) < radius for region in regions):
				correct[p1].add(p2)
				correct[p2].add(p1)

	print("size: {}".format(len(correct)))
	avgMatches = sum(len(matches) for _, matches in correct.items()) / len(correct)
	print("avgMatches: {}".format(avgMatches))

	with open("correctR{}.pickle".format(radius), mode='w+b') as correctFile:
		pickle.dump(correct, correctFile)
