# Thesis-Code



	MySQL configuration

You need to a supply a file called mysqlconfig.json with the information described in
https://dev.mysql.com/doc/connector-python/en/connector-python-example-connecting.html

Example mysqlconfig.json:
{
  "user": "myusername",
  "password": "mypassword",
  "host": "127.0.0.1",
  "database": "CPLOP",
}





	Running
To just get clusters:
TODO

Most files can be run standalone for verification / timing purposes.
Some files have unit tests while others just check equality with alternate implementations.

use clusterEval.py to evaluate the quality of clusters
The output isn't verbose so you'll need to look at the code to see what the outputted numbers are





	Caching

THIS FEATURE IS ONLY INTENDED FOR FASTER DEBUGGING, TESTING and TUNING

You can run cachethings.py to cache calculations things to disk for performance or consistency between multiple runs. The CACHE WILL BE INVALIDATED if the database is updated or you change the config. cachethings.py only caches things that aren't already cached which means YOU HAVE TO MANUALLY DELETE INVALIDATED CACHED FILES and run cachethings.py again.





	Cluster Configuration

Most parameters that you might want to change can be found in clusterConfig.json

Here is a commented example of clusterConfig.json:
{
	"regions": [				# if you change the order or the contents you'll need to recache everything
		{
			"name": "23-5",								# needs to match the values in CPLOP
			"dispCount": 93,							# becauses zScores were precomputed with it

			"pearsonSimilarityThresholdAlpha": 0.995,	# eps for dbscan. changes
			"pearsonSimilarityThresholdBeta": 0.99,		# invalidate cached clusters

			"betaDistributionAlpha": 386.65,			# these two are used
			"betaDistributionBeta": 1.36074				# for cluster evaluation
		}, {
			"name": "16-23",
			"dispCount": 95,

			"pearsonSimilarityThresholdAlpha": 0.995,
			"pearsonSimilarityThresholdBeta": 0.99,

			"betaDistributionAlpha": 668.1768,
			"betaDistributionBeta": 1.532336
		}
	],

	"minNeighbors": 3,			# dbscan param: values of 1 or 2 make dbscan
								# degenerate to single link agglomerative
								# (without the hierarchy)
								# invalidates cached clusters

	"pointsPerLeaf": 8,			# spatial index param

	"dimsPerSplit": 10,			# spatial index param

	"isolateSubsetSize": 3000	# size of randomly sampled subset to use.
								# set to "All" to use all isolates
								# run cachethings.py for consistency between timing runs.
								# changing this does not invalidate other cached subsets
}
