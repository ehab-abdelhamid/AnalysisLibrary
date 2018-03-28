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


#include "utils.h"
#include "mst2clusters.h"

using namespace NWUClustering;

static void usage(char *argv0) 
{
    const char *params =
	"Usage: %s [switches] -i filename -o filename\n"
    	"	-i filename					: input file containing the MSTs\n"
	"	-e eps_prime					: input parameter eps_prime to compute clusters\n"	
	"	-t threads					: number of threads to perform the computation\n"
	"	-o output					: output file to store the clusters\n\n";

    fprintf(stderr, params, argv0);
    exit(-1);
}

int main(int argc, char** argv)
{
	double	seconds;
	int 	opt;

	float 	eps_prime;
	int threads;	
	char* 	outfilename;
	char*   infilename;

	// some default values
	outfilename 	= NULL;
	infilename	= NULL;
	eps_prime	= -1;
	threads 	= -1;	

	while ((opt=getopt(argc,argv,"i:o:e:t:"))!= EOF)
    	{
		switch (opt)
        	{
        		case 'i':
 	        	   	infilename = optarg;
        	       		break;
			case 'o':
				outfilename = optarg;
				break;
			case 't':
				threads = atoi(optarg);		
				break;
			case 'e':
				eps_prime = atof(optarg);		
				break;
            		case '?':
           			usage(argv[0]);
     		      		break;
           		default:
           			usage(argv[0]);
           			break;
    		}
	}

	if(infilename == NULL || outfilename == NULL || eps_prime < 0 || threads <= 0)
	{
		usage(argv[0]);
		exit(-1);
	}
	
	omp_set_num_threads(threads);	

	MSTtoClusteringAlgo mca;	
	
	double start = omp_get_wtime();
	mca.read_file(infilename);
	cout << "Reading the input file took: " << omp_get_wtime() - start << " seconds."<< endl;

	start = omp_get_wtime();
	compute_clusters(mca, eps_prime);
	cout << "Computing the clusters took: " << omp_get_wtime() - start << " seconds."<< endl;
	
	start = omp_get_wtime();
	mca.saveClusters(outfilename);
	cout << "Saving MST to file took: " << omp_get_wtime() - start << " seconds."<< endl;

	return 0;
}
