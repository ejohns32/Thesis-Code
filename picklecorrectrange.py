import pickle
import sys
import cProfile

import config

def run(isolates, cfg):
	correct = {isolate: set() for isolate in isolates}
	for i, iso1 in enumerate(isolates):
		print("{}/{}".format(i, len(isolates)))
		for iso2 in isolates[i:]:
			if iso1.isWithinRadiiOf(iso2, cfg.radii):
				correct[iso1].add(iso2)
				correct[iso2].add(iso1)
	return correct

if __name__ == '__main__':
	fileName = sys.argv[1] if len(sys.argv) > 1 else "isolates.pickle"
	with open(fileName, mode='r+b') as zScoresFile:
		isolates = pickle.load(zScoresFile)

	cfg = config.loadConfig()

	# cProfile.run("correct = run(isolates, cfg)")
	correct = run(isolates, cfg)

	print("size: {}".format(len(correct)))
	avgMatches = sum(len(matches) for _, matches in correct.items()) / len(correct)
	print("avgMatches: {}".format(avgMatches))

	if len(sys.argv) == 1:
		with open("correctT{}.pickle".format(cfg.threshold), mode='w+b') as correctFile:
			pickle.dump(correct, correctFile)
