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

#include "mst2clusters.h"

namespace NWUClustering
{
	MSTtoClusteringAlgo::~MSTtoClusteringAlgo()
	{
		m_vec_f_core_distance.clear();
		m_mst.clear();
		m_vec_clusterID.clear();
	}

	int MSTtoClusteringAlgo::read_file(char* infilename)
	{
		int num_points, edge_count, i;

		string line, line2, buf;
                ifstream file(infilename);
                stringstream ss;
	
		if (file.is_open())
		{
			// get the first line and get the point count and number of edges in MST
			getline(file, line);
			line2 = line;
			ss.clear();
			ss << line2;

			ss >> buf;
			num_points = atoi(buf.c_str());

			ss >> buf;
			edge_count = atoi(buf.c_str());
				
			cout << "Number of points: " << num_points << " edges: " << edge_count << endl;
		
			m_i_num_points = num_points;	

			m_vec_f_core_distance.resize(num_points, FLT_MAX);
			MST_Node node;
			m_mst.reserve(edge_count);

			int u, v;
			float w, core_dis;	
			for (i = 0; i < edge_count; i++)				
			{
				
				getline(file, line);	
				if(line.length() == 0)
					continue;

				ss.clear();
				ss << line;

				ss >> buf;
				node.m_reachable_to = atoi(buf.c_str());
				
				ss >> buf;	
				node.m_reachable_from = atoi(buf.c_str());

				ss >> buf;	
				node.m_f_between_distance = atof(buf.c_str());
				
				ss >> buf;	
				core_dis = atof(buf.c_str());


				m_mst.push_back(node);
				m_vec_f_core_distance[node.m_reachable_from] = core_dis;


				//if(i <= 3)
				//	cout << node.m_reachable_to << " " << node.m_reachable_from << " " << node.m_f_between_distance << " " << m_vec_f_core_distance[node.m_reachable_from] << endl;
				//if(i >= edge_count - 3)
				//	cout << node.m_reachable_to << " " << node.m_reachable_from << " " << node.m_f_between_distance << " " << m_vec_f_core_distance[node.m_reachable_from] << endl;
			}
		}
		else
		{
			cout << "Error: no such file: " << infilename << endl;
			return -1;
		}
		
		return 0;
	}


	void MSTtoClusteringAlgo::saveClusters(char* outfilename)
	{
		ofstream clusterfile;
		clusterfile.open(outfilename);	
		int i;
		//ostream& o;
		for(i = 0; i < m_i_num_points; i++)
		{
			clusterfile << i << " " << m_vec_clusterID[i] << endl;	
		}

		clusterfile.close();
	}

	int compute_clusters(MSTtoClusteringAlgo& mca, float eps_prime)
	{
		cout << "Computing clusters from MSTs in parallel using the union find algorithms" << endl;
		
		vector <int> vec_parents;

		vec_parents.reserve(mca.m_i_num_points);	
		vec_parents.resize(mca.m_i_num_points, 0);
		
		omp_lock_t *nlocks;
		nlocks = (omp_lock_t *) malloc(mca.m_i_num_points * sizeof(omp_lock_t));

		int i, u, v, root, root_u, root_v, x, z;
                float w;

		#pragma omp parallel for private(i) shared(mca)
                for(i = 0; i < mca.m_i_num_points; i++)
                {
                        vec_parents[i] = i;
                        omp_init_lock(&nlocks[i]);
                }

		// computing the trees for each cluster
		#pragma omp parallel for private(i, u, v, w, root_u, root_v, z) shared(mca, eps_prime, nlocks, vec_parents)
                for(i = 0; i < mca.m_mst.size(); i++)
                {
                        u = mca.m_mst[i].m_reachable_to;
                        v = mca.m_mst[i].m_reachable_from;
                        w = mca.m_mst[i].m_f_between_distance;

                        if(w <= eps_prime && mca.m_vec_f_core_distance[v] != FLT_MAX)
                        {
                                // union u and v
                                root_u = u;
                                root_v = v;
                                while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                        if(vec_parents[root_u] < vec_parents[root_v])
                                        {
                                                if(root_u == vec_parents[root_u])
                                                {
                                                        omp_set_lock(&nlocks[root_u]);
                                                        bool p_set = false;
                                                        if(root_u == vec_parents[root_u])
                                                        {
                                                                vec_parents[root_u] = vec_parents[root_v];
                                                                p_set = true;
                                                        }
                                                        omp_unset_lock(&nlocks[root_u]);
                                                        if(p_set == true)
                                                                break;
                                                }

                                                z = vec_parents[root_u];
                                                vec_parents[root_u] = vec_parents[root_v];
                                                root_u = z;
                                        }
                                        else
                                        {
                                                if(root_v == vec_parents[root_v])
                                                {
                                                        omp_set_lock(&nlocks[root_v]);
                                                        bool p_set = false;
                                                        if(root_v == vec_parents[root_v])
                                                        {
                                                                vec_parents[root_v] = vec_parents[root_u];
                                                                p_set = true;
                                                        }
                                                        omp_unset_lock(&nlocks[root_v]);
                                                        if(p_set == true)
                                                                break;
                                                }

                                                z = vec_parents[root_v];
                                                vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;
                                        }
                                }
                        }
                }
	
		// assign cluster IDs to each point
		// shorten the path so that each point points to its parents

		#pragma omp parallel for private(i, root, x, z) shared(mca)
                for(i = 0; i < mca.m_i_num_points; i++)
                {
                        if(vec_parents[i] == i)
                                continue;

                        root = i;

                        while(root != vec_parents[root])
                                root = vec_parents[root];

                        // compress the path
                        if(vec_parents[i] != root)
                        {
                                x = i;
                                while(root != vec_parents[x])
                                {
                                        z = vec_parents[x];
                                        vec_parents[x] = root;
                                        x = z;
                                }
                        }
                }

		// now the hight of the tree is only 1
                // find the noise and clusters
                int root_count = 0, cluster_count = 0, noise_count = 0, points_in_clusters = 0;
                vector <int> child_count;
                child_count.resize(mca.m_i_num_points, 0);

                #pragma omp parallel for private(i) shared(mca, nlocks, vec_parents)
                for(i = 0; i < mca.m_i_num_points; i++)
                {
                        if(vec_parents[i] == i)
                        {
                                #pragma omp atomic
                                root_count++;
                        }
                        //else
                        omp_set_lock(&nlocks[vec_parents[i]]);
                        child_count[vec_parents[i]]++; // height is one
                        omp_unset_lock(&nlocks[vec_parents[i]]);
                }


		mca.m_vec_clusterID.clear();
		mca.m_vec_clusterID.resize(mca.m_i_num_points, INT_MIN);

		// assign ID for the roots including the noise points
                for(i = 0; i < mca.m_i_num_points; i++)
                {
                        if(child_count[i] == 1)
                        {
                                noise_count++;
                                mca.m_vec_clusterID[i] = NOISE_ID; // 0 is for noise
                        }
                        else if(child_count[i] > 1)
                        {
                                cluster_count++;

                                points_in_clusters += child_count[i];
                                mca.m_vec_clusterID[i] = cluster_count;
                        }
                }

                // assign IDs for the rest of the points part of the clusters
                #pragma omp parallel for private(i) shared(mca)
                for(i = 0; i < mca.m_i_num_points; i++)
                {
                        if(mca.m_vec_clusterID[i] == INT_MIN)
                        //if(child_count[i] == 0)
                                mca.m_vec_clusterID[i] = mca.m_vec_clusterID[vec_parents[i]];
                }

                mca.m_clusters_found = cluster_count;
                mca.m_points_in_clusters = points_in_clusters;
                mca.m_noise_count = noise_count;

		cout << "Numbers of clusters: " << cluster_count << " points in clusters " << points_in_clusters << " noise " << noise_count << endl;
		
		child_count.clear();
		free(nlocks);
		vec_parents.clear();

		return 0;
	}
};

