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

#ifndef _MST_TO_CLUSTERS_
#define _MST_TO_CLUSTERS_

#include "utils.h"

namespace NWUClustering
{
	struct MST_Node
	{
		int m_reachable_to;
		int m_reachable_from;
		float m_f_between_distance;
	};

	class MSTtoClusteringAlgo
	{
	public:
		MSTtoClusteringAlgo(){ }
		virtual ~MSTtoClusteringAlgo();

		int read_file(char* infilename);
		void saveClusters(char* outfilename);

	public:
		vector<float>		m_vec_f_core_distance;
		vector <MST_Node> 	m_mst;

		vector<int>     	m_vec_clusterID;
		int             	m_clusters_found;
                int             	m_points_in_clusters;
                int             	m_noise_count;

		
		int 			m_i_num_points;
	};

	int compute_clusters(MSTtoClusteringAlgo& mca, float eps_prime); 
	//void saveMST(ClusteringAlgo& dbs, ostream& o);
};

#endif
