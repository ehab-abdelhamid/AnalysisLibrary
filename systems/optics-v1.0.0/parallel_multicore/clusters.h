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

#ifndef _CLUSTER_
#define _CLUSTER_

#include "utils.h"
#include "kdtree2.hpp"

namespace NWUClustering
{
    	struct record_anl
        {
        	float d[3];
        	int   cls;
        };

	struct Points
	{
		array2dfloat m_points;
		int	m_i_dims;
		int 	m_i_num_points;
		//float 	m_f_core_distance;
		//float	m_f_reachability_distance;
	};

	struct MST_Node
	{
		int m_reachable_to;
		int m_reachable_from;
		//float m_reachability;
		float m_f_best_distance;
	};

	inline bool operator<(const MST_Node& n1, const MST_Node& n2) {
 		return (n1.m_f_best_distance > n2.m_f_best_distance); // 
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
		Clusters():m_pts(NULL),m_kdtree(NULL){ }
		virtual ~Clusters();

		int     read_file(char* infilename, int isBinaryFile);
		//void 	randomInit (int dims = 5, int num_points = 10);

		void    writePoints(ostream& o);
		int     build_kdtree();
		
	public:
		Points* 	m_pts;
		kdtree2* 	m_kdtree;
		vector <int> 	m_pid_to_cid; // point id to cluster id
		vector <vector <int> > m_clusters;
		int     m_parcent_of_data;
	};
};

#endif

