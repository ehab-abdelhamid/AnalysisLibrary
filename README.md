Introduction
------------
Distributed data analysis library based on MPI
This library contains a number of clustering and clasification libraries;
Clustering systems:
1- dbscan
2- OPTICS
3- KMeans
4- PINK
5- HDBSCAN
6- Canopy
Classification systems:
1- SVM
2- Random Forest

Compilation
-----------
Each systems has its own compilation instructions. Check the details in their corresponding folders.

Running
-------
Use the cluster.sh script to run clustering techniques. script parameters are:
1- input file
2- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)
3- output file
4- clustering technique:
	1- dbscan, 2- optics, 3- kmeans, 4- pink, 5-HDBSCAN, 6- Canopy
5- clustering parameters
Example to run pink:
./cluster.sh ../datasets/Credit/credit-card.bin 4 out-pink.txt 4

Use the classify.sh script to run classification techniques. script parameters are:
1- input file
2- operation type; train or classify: 1-SVM train, 2-SVM classify, 3- RF train and classify
3- used clasification model
	1- SVM
	2- Random Forest
4- model (path for SVM, training file for RF)
5- output file
6- application specific parameters
7- MPI specific parameters
8- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)
9- classification parameter index (from)
10- classification parameter index (to)
Example for Random Forest train then test:
./classify.sh ../datasets/CLData/test/dev-first10_.arff 3 0 ../datasets/CLData/train/train-first10_.arff ./RF-out.out "" "-n 4" 3 10 10
