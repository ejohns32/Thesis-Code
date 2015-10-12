import pickle
import numpy
import itertools
import random

count = 100

with open("zScores.pickle", mode='r+b') as zScoresFile:
	zScores = pickle.load(zScoresFile)
with open("correctN{}.pickle".format(count), mode='r+b') as correctFile:
	correctNearest = pickle.load(correctFile)

regions = ['23-5', '16-23']
radii = [0.964, 1.37]

for radius in radii:
	correct = {isoID: [iID for iID in correctNearest[isoID] if max(numpy.linalg.norm(zScores[isoID][region] - zScores[iID][region]) for region in regions) < radius] for isoID in correctNearest.keys()}

	print("size: {}".format(len(correct)))
	avgMatches = sum(len(matches) for _, matches in correct.items()) / len(correct)
	print("avgMatches: {}".format(avgMatches))

	with open("correctN{}R{}.pickle".format(count, radius), mode='w+b') as correctFile:
		pickle.dump(correct, correctFile)
