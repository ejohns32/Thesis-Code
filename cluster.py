import sys
import pickle
import cProfile

import config
import spatial

def loadIsolates():
	fileName = sys.argv[1] if len(sys.argv) > 1 else "isolates.pickle"
	with open(fileName, mode='r+b') as zScoresFile:
		return pickle.load(zScoresFile)

def makeTree(isolates, cfg):
	# printTree(tree, 0)
	# print("inner, leaf nodes: {}".format(countNodes(tree)))
	# print("min, max depth: {}".format(getMinMaxDepth(tree, 0)))
	return spatial.Tree(spatial.TreeConfig(cfg.regions, cfg.pointsPerLeaf, cfg.dimsPerSplit), isolates)

def loadCorrect(cfg):
	with open("correctT{}.pickle".format(cfg.threshold), mode='r+b') as correctFile:
		correctRange = pickle.load(correctFile)
	return correctRange

def testSpatial(isolates, tree, correctRange, cfg):
	# queryCount = 1000
	# queryIsolates = random.sample(isolates, queryCount)
	queryIsolates = set(isolates)
	queryCount = len(queryIsolates)

	correctCount = 0
	nonZeroCorrectCount = 0
	nonZeroCount = 0
	extraCount = 0
	missingCount = 0

	seenIsolates = set()
	unseenCore = set()

	i = 0
	while len(unseenCore) + len(queryIsolates) > 0:
		i += 1
		if len(unseenCore) > 0:
			print("*" ,end="")
			isolate = unseenCore.pop()
		else:
			isolate = queryIsolates.pop()

		# print("{}/{} - {}".format(i, len(queryIsolates), isolate.name))
		resultR = tree.rangeQuery(isolate, cfg.radii)
		# correctR = correctRange[isolate] - seenIsolates

		# if len(resultR) > 0:
		# 	nonZeroCount += 1
		# extraCount += len(resultR - correctR)
		# missingCount += len(correctR - resultR)

		# print("{}/{} - {} - {}:{}".format(i, len(queryIsolates), isolate.name, len(resultR), len(correctR)))
		# # print("\t{} --- {} / {} : {} / {}".format(isolate, len(resultR - correctR), len(resultR), len(correctR - resultR), len(correctR)))
		# if resultR == correctR:
		# 	correctCount += 1
		# 	if len(resultR) > 0:
		# 		nonZeroCorrectCount += 1
		# else:
		# 	print("\t{} ----- {}".format([(extraIsolate, isolate.regionsDist(extraIsolate)) for extraIsolate in resultR - correctR], [(missingIsolate, isolate.regionsDist(missingIsolate)) for missingIsolate in correctR -resultR]))

		# seenIsolates |= correctR
		# unseenCore |= correctR - {isolate}
		# queryIsolates -= correctR

		seenIsolates |= resultR
		unseenCore |= resultR - {isolate}
		queryIsolates -= resultR

	print("{}/{}".format(nonZeroCorrectCount, nonZeroCount))
	print("{}/{}".format(correctCount, queryCount))
	print("{},{}".format(extraCount, missingCount))


if __name__ == '__main__':
	cfg = config.loadConfig()
	isolates = loadIsolates()
	tree = makeTree(isolates, cfg)
	correctRange = loadCorrect(cfg)

	# testSpatial(isolates, tree, correctRange, cfg)
	cProfile.run("testSpatial(isolates, tree, correctRange, cfg)")