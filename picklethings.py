import random
import pickle

import pyroprinting

if __name__ == '__main__':
	isolates = pyroprinting.loadIsolatesFromDB()

	with open("isolatesAll.pickle", mode='w+b') as isolatesFile:
		pickle.dump(isolates, isolatesFile)

	for count in range(1000, len(isolates), 1000):
		slicedIsolates = list(random.sample(isolates, count))
		print("{}".format(count))
		with open("isolates{}.pickle".format(min(count, len(slicedIsolates))), mode='w+b') as slicedFile:
			pickle.dump(slicedIsolates, slicedFile)
