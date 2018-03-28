/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   Files: omp_main.cpp clusters.cpp  clusters.h utils.h utils.cpp          */
/*              optics.cpp optics.h kdtree2.cpp kdtree2.hpp                  */
/*              geometric_partitioning.cpp geometric_partitioning.h          */
/*              mutable_priority_queue.h                                     */
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

#ifndef _OPTICS_
#define _OPTICS_

#include "utils.h"
#include "clusters.h"

namespace NWUClustering
{
	struct cluster_profile
	{
		int             m_clusters_found;
                int             m_points_in_clusters;
                int             m_noise_count;
	};

	class ClusteringAlgo : public Clusters
	{
	public:
                ClusteringAlgo(){ }
                virtual ~ClusteringAlgo();

		void 	set_dbscan_params(double eps, int minPts); 

	public:
		double  m_epsSquare;
                int     m_minPts;
                int     m_messages_per_round;
                int     m_compression;
                int     m_collective_comm;	

		int 	m_connecting_edges_count;

		vector <float>   	m_vec_f_core_distance;
		vector <bool>    	m_visited;
		vector <MST_Node>       m_mst;

		map <int, float>	m_global_point_ID_to_core_distance;

                vector <float>          	m_vec_f_core_distance_root;		
		vector <vector <MST_Edge> >	m_vec_vec_adjacency_list_root;
		int 				m_global_owned_local_point_count_root;
		vector <bool>           	m_visited_root;
		vector <MST_Node>       	m_mst_root;		
		vector<int>     		m_vec_clusterID_root;
		cluster_profile 		m_cp_root;

		priority_queue<MST_Node>*	m_pq_mst_node;
		vector <int>			m_vec_global_point_ID_gathered;
		vector <float>                  m_vec_f_core_distance_gathered;
		int				m_all_gid_gathered_upper;
		int                             m_all_gid_gathered_lower;

		vector <int>			m_vec_gid_parent_gid_pair;
	};

	void saveMST(ClusteringAlgo& dbs, ostream& o);

	void compute_mst_parallel(ClusteringAlgo& dbs); // compute the MST from the point set in parallel
	void updateNeighbors_mst_based_parallel(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in, int tid, map<int,float>& reach_so_far_local_in);
	void compute_mst_global_from_mst_local(ClusteringAlgo& dbs);
	void gather_local_msts_and_compute_adj_root(ClusteringAlgo& dbs);
	void updateNeighbors_adjacency_list_based_global_mst_from_local_mst(vector <bool>& visited_root, int pid, vector <MST_Edge>* vec_edges_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in);	

	// kruskals merging
	void pairwise_local_mst_merging_kruskal(ClusteringAlgo& dbs);
	void pairwise_local_mst_merging_kruskal_smart(ClusteringAlgo& dbs);
	void sequential_local_mst_merging_kruskal(ClusteringAlgo& dbs);

	void local_mst_merging_boruvka_filter_kruskal(ClusteringAlgo& dbs);
};

#endif
