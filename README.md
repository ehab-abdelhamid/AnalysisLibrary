Introduction
------------
Distributed data analysis library based on MPI<br />
This library contains a number of clustering and clasification libraries;<br />
Clustering systems:<br />
1- dbscan<br />
2- OPTICS<br />
3- KMeans<br />
4- PINK<br />
5- HDBSCAN<br />
6- Canopy<br />
Classification systems:<br />
1- SVM<br />
2- Random Forest<br />

Compilation
-----------
Each systems has its own compilation instructions. Check the details in their corresponding folders.<br />

Running
-------
Use the cluster.sh script to run clustering techniques. script parameters are:<br />
1- input file<br />
2- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)<br />
3- output file<br />
4- clustering technique:<br />
	1- dbscan, 2- optics, 3- kmeans, 4- pink, 5-HDBSCAN, 6- Canopy<br />
5- clustering parameters<br />
Example to run pink:<br />
./cluster.sh ../datasets/Credit/credit-card.bin 4 out-pink.txt 4<br />

Use the classify.sh script to run classification techniques. script parameters are:<br />
1- input file<br />
2- operation type; train or classify: 1-SVM train, 2-SVM classify, 3- RF train and classify<br />
3- used clasification model<br />
	1- SVM<br />
	2- Random Forest<br />
4- model (path for SVM, training file for RF)<br />
5- output file<br />
6- application specific parameters<br />
7- MPI specific parameters<br />
8- file type: 0-csv, 1-csv with header, 2-  sparse csv, 3- ARFF file, 4- binary (4-bytes), 5- binary (8-bytes)<br />
9- classification parameter index (from)<br />
10- classification parameter index (to)<br />
Example for Random Forest train then test:<br />
./classify.sh ../datasets/CLData/test/dev-first10_.arff 3 0 ../datasets/CLData/train/train-first10_.arff ./RF-out.out "" "-n 4" 3 10 10<br />
