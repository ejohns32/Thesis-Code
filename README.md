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





	Caching

THIS FEATURE IS ONLY INTENDED FOR FASTER DEBUGGING, TESTING and TUNING

You can run cachethings.py to cache calculations things to disk for performance or consistency between multiple runs. If the database is updated or you change the config, it'll invalidate the cache. For now that means you have to manually delete the cached files and run cachethings.py again.





	Cluster Configuration

Most parameters that you might want to change can be found in clusterConfig.json

Here is a commented example of clusterConfig.json:
{
	"regions": [				# if you change the order or the contents you'll need to recache everything
		{
			"name": "23-5",
			"dispCount": 93
		}, {
			"name": "16-23",
			"dispCount": 95
		}
	],
	"threshold": 0.995,			# dbscan param: defined as pearson correlation
								# but internally is converted to zScore distance

	"minNeighbors": 5,			# dbscan param: values of 1 or 2 make dbscan
								# degenerate to single link agglomerative
								# (without the hierarchy)

	"pointsPerLeaf": 8,			# spatial index param

	"dimsPerSplit": 10,			# spatial index param

	"isolateSubsetSize": 3000	# size of randomly sampled subset to use.
								# set to "All" to use all isolates
								# run cachethings.py for consistency between timing runs.
								# changing this does not invalidate other cached subsets
}
