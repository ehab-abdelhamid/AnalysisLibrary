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


#ifndef _CLUSTER_
#define _CLUSTER_

#include "utils.h"
#include "kdtree2.hpp"

namespace NWUClustering
{
	struct Points
	{
		array2dfloat 	m_points;
		int 		m_i_dims;
		int 		m_i_num_points;
		int 		m_i_num_local_points;
		interval* 	m_box;
		vector<int> 	m_vec_global_point_ID;
		vector <int> 	m_vec_processor_ID;
		vector <int> 	m_vec_local_ind;

		int 		m_owned_local_point_gID_upper;
		int		m_owned_local_point_gID_lower;
		int             m_not_owned_local_point_gID_upper;
                int             m_not_owned_local_point_gID_lower;
		int             m_all_local_point_gID_upper;
                int             m_all_local_point_gID_lower;

		int 		m_total_point_count;
	};

	struct Points_Outer
    	{
        	array2dfloat m_points;
		vector <int> m_prIDs;
		vector <int> m_ind;

        	int m_i_dims;
        	int m_i_num_points;
		interval* m_box;
		vector<int> m_vec_global_point_ID;
    	};
	
	struct MST_Node
	{
		int m_reachable_to;
		int m_reachable_from;
		//float m_reachability;
		//float m_f_best_distance;
		float m_f_between_distance;
	};

	inline bool operator<(const MST_Node& n1, const MST_Node& n2) {
 		return (n1.m_f_between_distance > n2.m_f_between_distance); // 
	//return (e1.dis < e2.dis) || ((e1.dis == e2.dis) && (e1.idx < e2.idx)); // RECHECK
	}

	
	struct MST_Edge
	{
		int 	m_other_end;
		float 	m_f_edge_weight;
	};

	struct Point_Index_Pair
	{
		int m_reachable_to;
		int m_reachable_from;
	};

	inline bool operator<(const Point_Index_Pair& e1, const Point_Index_Pair& e2) {
  		return (e1.m_reachable_to < e2.m_reachable_to);
	}

	class Clusters
	{
	public:
		Clusters():m_pts(NULL),m_kdtree(NULL),m_pts_outer(NULL),m_kdtree_outer(NULL), m_parcent_of_data(-1){ }
		virtual ~Clusters();
		
		bool 	allocate_outer(int dims);
		//bool 	addPoints(int source, int buf_size, int dims, vector<float>& raw_data);
		//bool 	updatePoints(vector< vector<int> >& raw_ind);

		bool    addPoints_NEW(int source, int buf_size, int dims, vector<float>& raw_data);
		bool    updatePoints_NEW(vector< vector<int> >& raw_ind, vector< vector<int> >& raw_gid);

		int     read_file(char* infilename, int isBinaryFile);

		void    writePoints(ostream& o);
		int     build_kdtree();
		int     build_kdtree_outer();	

		void assign_prID_ind_local_points();	
		void compute_upper_lower_limit_gids();
	public:
		Points* 	m_pts;
		kdtree2* 	m_kdtree;
		int     	m_parcent_of_data;

		Points_Outer*	m_pts_outer;
		kdtree2*    	m_kdtree_outer;

		vector <int> 	m_pid_to_cid; // point id to cluster id
		vector <vector <int> > m_clusters;
	};
};

#endif

