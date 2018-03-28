/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   Files: omp_main.cpp clusters.cpp  clusters.h utils.h utils.cpp          */
/*                      optics.cpp optics.h kdtree2.cpp kdtree2.hpp          */
/*			mutable_priority_queue.h			     */
/*                                                                           */
/*   Description: an openmp implementation of OPTICS clustering algorithm    */
/*                              using Graph Algorithmic Techniques           */
/*                                                                           */
/*   Author:  Md. Mostofa Ali Patwary                                        */
/*            EECS Department, Northwestern University                       */
/*            Email: mostofa.patwary@gmail.com                               */
/*                                                                           */
/*   Copyright, 2013, Northwestern University                                */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/*   Please cite the following publication if you use this package           */
/*                                                                           */
/*   Md. Mostofa Ali Patwary, Diana Palsetia, Ankit Agrawal, Wei-keng Liao,  */
/*   Fredrik Manne, and Alok Choudhary, "Scalable Parallel OPTICS Data 	     */
/*   Clustering Using Graph Algorithmic Techniques", Proceedings of the      */
/*   International Conference on High Performance Computing, Networking,     */
/*   Storage and Analysis (Supercomputing, SC'13), pp.49:1-49:12, 2013.      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "optics.h"
#include "utils.h"
#include "kdtree2.hpp"


static void usage(char *argv0) 
{
    const char *params =
	"Usage: %s [switches] [-i filename -b] [-d dimensions -p points] -m min_points_cluster -e epsilon -o output -t threads\n"
	"	USE -i and -b option to read data from file\n"
	"	ELSE use -d and -p option to generate points randomly\n"
    	"	-i filename					: file containing data to be clustered\n"
    	"	-b						: input file is in binary format (default no)\n"
	"	-m min_points_cluster				: min points to form a cluster, e.g. 2\n"
	"	-e epsilon					: radius or threshold on neighbourhoods retrieved, e.g. 0.8\n"
	"	-o output					: save the MST into the output file (optional)\n"
	"	-t threads					: number of threads to be employed\n\n";

    fprintf(stderr, params, argv0);
    exit(-1);
}

int main(int argc, char** argv)
{
	double	seconds;
	int 	opt;

	int 	dims, threads, minPts;
	double 	eps;
	char* 	outfilename;
	int     isBinaryFile;
	char*   infilename;

	// some default values
	dims 		= -1;
	eps		= -1;
	isBinaryFile 	= 0;
	outfilename 	= NULL;
	infilename	= NULL;
	threads 	= 1;

	while ((opt=getopt(argc,argv,"i:t:m:e:o:b"))!= EOF)
    	{
		switch (opt)
        	{
        		case 'i':
 	        	   	infilename = optarg;
        	       		break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'b':
            			isBinaryFile = 1;
            			break;
    		        case 'm':
            			minPts = atoi(optarg);
            			break;
            		case 'e':
            			eps  = atof(optarg);
            			break;
			case 'o':
				outfilename = optarg;
				break;
            		case '?':
           			usage(argv[0]);
     		      		break;
           		default:
           			usage(argv[0]);
           			break;
    		}
	}

	if(infilename == NULL || minPts < 0 || eps < 0 || threads <= 0)
	{
		usage(argv[0]);
		exit(-1);
	}
	
	omp_set_num_threads(threads);

	NWUClustering::ClusteringAlgo dbs;
	dbs.set_optics_params(eps, minPts);

	cout << "Input Parameters:" << " minPts " << minPts << " eps " << eps << endl;
	cout << "Number of threads " << threads << endl;

	double start = omp_get_wtime();
	cout << "Reading points from file: " << infilename << endl;

	// reading points from the input file
	if(dbs.read_file(infilename, isBinaryFile) == -1) // reading data files
	{
		cout << "Invalid file: " << infilename << endl;
		exit(-1);
	}
	cout << "Reading points took " << omp_get_wtime() - start << " seconds."<< endl;

	// build kdtree for the points
	start = omp_get_wtime();
	dbs.build_kdtree();
	cout << "Building kdtree took " << omp_get_wtime() - start << " seconds."<< endl;

	//MST based OPTICS algorithm	
	start = omp_get_wtime();
	compute_mst_parallel(dbs);
	cout << "Total time taken by parallel MST based OPTICS: " << omp_get_wtime() - start << " seconds."<< endl;
	cout << "Edges in global MST " << dbs.m_mst.size() << endl;

	if(outfilename != NULL)
	{
		start = omp_get_wtime();
		ofstream outputfile;
		outputfile.open(outfilename);
		saveMST(dbs, outputfile);
		outputfile.close();
		cout << "Saving MST to file took: " << omp_get_wtime() - start << " seconds."<< endl;
	}
	return 0;
}
