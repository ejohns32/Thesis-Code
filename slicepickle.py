import pickle
import random

with open("zScores.pickle", mode='r+b') as zScoresFile:
	zScores = pickle.load(zScoresFile)

for count in range(1000, len(zScores), 1000):
	slicedZScores = dict(random.sample(zScores.items(), count))
	print("{}".format(count))
	with open("zScores{}.pickle".format(min(count, len(slicedZScores))), mode='w+b') as slicedFile:
		pickle.dump(slicedZScores, slicedFile)
