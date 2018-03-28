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

#ifndef _OPTICS_
#define _OPTICS_

#include "utils.h"
#include "clusters.h"

namespace NWUClustering
{
	class ClusteringAlgo : public Clusters
	{
	public:
		ClusteringAlgo():m_extra_dis_counter(0){ }
		virtual ~ClusteringAlgo();

		// functions for dbscan algorithm
		void set_optics_params(double eps, int minPts); //, int parcent_of_data);
	
	public:
		
		// parameters to run dbscan algorithm
		double 	m_epsSquare;
		int 	m_minPts;

		vector<float>	m_vec_f_core_distance;
		vector<float> 	m_vec_f_reachability_distance;
		vector<int>	m_vec_i_closest_point;
 		vector<int>	m_vec_i_ordered_points;
		vector<bool> 	m_visited;	

		//sequential
		vector <MST_Node> 		m_mst;
		vector <vector <MST_Edge> > 	m_vec_vec_adjacency_list;
		vector <int>			m_vec_parents;

		//parallel
		vector <vector <MST_Node> >	m_mst_per_thread;
		vector <vector <bool> >     	m_visited_per_thread;	

		vector <int> 	m_connecting_edges_count;
		vector <int>    m_initial_prID;

		int 	m_extra_dis_counter;
	
		priority_queue <MST_Node> m_mst_pq;
		vector < priority_queue <MST_Node> > m_mst_pq_per_thread;
	};
	
	// classcial
	void computeCoreDistance(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in);

	//parallel
	void compute_mst_parallel(ClusteringAlgo& dbs); // compute the MST from the point set in parallel
	void updateNeighbors_mst_based_parallel(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in, int tid, map<int,float>& reach_so_far_local_in);
	void updateNeighbors_adjacency_list_based_global_mst_from_local_mst(ClusteringAlgo& dbs, int pid, vector <MST_Edge>* vec_edges_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in);
	void compute_adjacency_list_from_mst_local_parallel(ClusteringAlgo& dbs); // internal function, get adjacency list from mst
	void compute_mst_global_from_mst_local(ClusteringAlgo& dbs);
	void compute_mst_global_from_mst_local_kruskal(ClusteringAlgo& dbs);
	void saveMST(ClusteringAlgo& dbs, ostream& o);
};

#endif
