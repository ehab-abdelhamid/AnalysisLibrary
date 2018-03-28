/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   Files: omp_main.cpp clusters.cpp  clusters.h utils.h utils.cpp          */
/*             	optics.cpp optics.h kdtree2.cpp kdtree2.hpp          	     */
/*             	geometric_partitioning.cpp geometric_partitioning.h 	     */
/*             	mutable_priority_queue.h                             	     */
/*                                                                           */
/*   Description: an mpi implementation of OPTICS clustering algorithm    */
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
/*   Fredrik Manne, and Alok Choudhary, "Scalable Parallel OPTICS Data       */
/*   Clustering Using Graph Algorithmic Techniques", Proceedings of the      */
/*   International Conference on High Performance Computing, Networking,     */
/*   Storage and Analysis (Supercomputing, SC'13), pp.49:1-49:12, 2013.      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



#include "optics.h"
#include "utils.h"
#include "kdtree2.hpp"
#include "geometric_partitioning.h"

static void usage(char *argv0) 
{
    const char *params =
	"Usage: %s [switches] [-i filename -b] -m min_points_cluster -e epsilon -o output\n"
    	"	-i filename		: file containing data to be clustered\n"
    	"	-b			: input file is in binary format (default no)\n"
	"	-m min_points_cluster	: min points to form a cluster, e.g. 2\n"
	"	-e epsilon		: radius or threshold on neighbourhoods retrieved, e.g. 0.8\n"
	"	-o output		: save the MST (optional)\n\n";
	
    fprintf(stderr, params, argv0);
    exit(-1);
}


int main(int argc, char** argv)
{
	double	seconds;
	int 	opt;
	int 	minPts, procs;
	double 	eps;
	char* 	outfilename;
	int     isBinaryFile;
	char*   infilename;
	int rank, nproc;

    	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	// some default values
	minPts 		= -1;
	eps		= -1;
	isBinaryFile 	= 0;
	outfilename 	= NULL;
	infilename	= NULL;

    	while ((opt=getopt(argc,argv,"i:m:e:o:b"))!= EOF)
    	{
		switch (opt)
        	{
        		case 'i':
            			infilename = optarg;
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

	if((infilename == NULL) || (minPts < 0 || eps < 0))
	{
		if(rank == proc_of_interest)
			usage(argv[0]);
		MPI_Finalize();
		exit(-1);
	}

	if(rank == proc_of_interest) cout << "\n\nNumber of processes " << nproc << endl;

	// check if nproc is NOT multiple of TWO
	unsigned int proc_count = nproc;
 	while (((proc_count % 2) == 0) && proc_count > 1)
   		proc_count /= 2;
 	
	if(proc_count != 1)
	{
		if(rank == proc_of_interest) cout << "\n\nNumber of processes (" << nproc << ") is NOT a multiple of TWO" << endl;
		MPI_Finalize();
		return 0;
	}

	// set the input parameters
	NWUClustering::ClusteringAlgo dbs;	
	dbs.set_dbscan_params(eps, minPts);

    	if(rank == proc_of_interest) cout << "Input paramters: Epsilon: " << eps << " MinPts: " << minPts << endl;
	
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();

	// read the input file, each process reads its own part ONLY
	if(rank == proc_of_interest) cout << "Reading points from file: " << infilename << endl;
	if(dbs.read_file(infilename, isBinaryFile) == -1)
	{
		MPI_Finalize();
		exit(-1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == proc_of_interest) cout << "Reading points took: " << MPI_Wtime() - start << " seconds." << endl;

	start = MPI_Wtime();
	
	// partion the points geometrically and exchange among the processes accordingly
	start_partitioning(dbs);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == proc_of_interest) cout << "Partitioning the points geometrically among the processes took: " << MPI_Wtime() - start << " seconds." << endl;

	start = MPI_Wtime();
	get_extra_points(dbs);
	MPI_Barrier(MPI_COMM_WORLD);	

	if(rank == proc_of_interest) cout << "Gathering extra points from neighbor processes took: " << MPI_Wtime() - start << " seconds." << endl;

	start = MPI_Wtime();
	// compute the lower and upper IDs of the points each process owns
	dbs.compute_upper_lower_limit_gids();
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank == proc_of_interest) cout << "Computing upper and lower limit of global IDs took: " << MPI_Wtime() - start << " seconds."<< endl;

	start = MPI_Wtime();
	
	// each process builds kdtree on its own points
	dbs.build_kdtree();
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == proc_of_interest) cout << "Building kdtree took: " << MPI_Wtime() - start << " seconds." << endl;

	// MST based parallel OPTICS algorithm	
	compute_mst_parallel(dbs);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == proc_of_interest) cout << "MST based parallel OPTICS algorithm took: " << MPI_Wtime() - start << " seconds."<< endl;

	if(outfilename != NULL)
        {
                start = MPI_Wtime();
		if(rank == proc_of_interest)
		{
                	ofstream outputfile;
                	outputfile.open(outfilename);
                	saveMST(dbs, outputfile);
                	outputfile.close();
		}
		double end = MPI_Wtime();
		if(rank == proc_of_interest) cout << "Saving MST to file took: " << end - start << " seconds."<< endl;
        }

	MPI_Finalize();
	return 0;
}
