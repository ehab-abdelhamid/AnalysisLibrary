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

namespace NWUClustering
{
	void ClusteringAlgo::set_optics_params(double eps, int minPts)//,int parcent_of_data)
	{
		m_epsSquare =  eps * eps;
		m_minPts =  minPts;
		m_parcent_of_data = 100; //parcent_of_data;
	}

	ClusteringAlgo::~ClusteringAlgo()
	{
		m_vec_f_core_distance.clear();
                m_vec_f_reachability_distance.clear();    
                m_vec_i_closest_point.clear();	
		m_vec_i_ordered_points.clear();		
		m_visited.clear();		

		m_mst.clear();	
		m_vec_vec_adjacency_list.clear();

		m_mst_per_thread.clear();
		m_visited_per_thread.clear();
		m_connecting_edges_count.clear();
		m_initial_prID.clear();
		m_mst_pq_per_thread.clear();		
	}
/////////////////////////////////////////////////////////////////////////
	void saveMST(ClusteringAlgo& dbs, ostream& o)
	{
		int i, u, v;
		float w, core_dis;
		o << dbs.m_pts->m_i_num_points << " " << dbs.m_mst.size() << endl; 
		for(i = 0; i < dbs.m_mst.size(); i++)
		{
			u = dbs.m_mst[i].m_reachable_to;
			v = dbs.m_mst[i].m_reachable_from;
			w = dbs.m_mst[i].m_f_best_distance;
			core_dis = dbs.m_vec_f_core_distance[v];

			o << u << " " << v << " " << sqrt(w) << " " << sqrt(core_dis) << endl;
		}
	}

/////////////////////////////////////////////////////////////////////////////////////
	
	void computeCoreDistance(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in)
	{
		if(ne_in.size() >= dbs.m_minPts)
		{
			int i, count_neighbors = ne_in.size();
			priority_queue<kdtree2_result> pq;
			for(i = 0; i < ne_in.size(); i++)
			{
				pq.push(ne_in[i]);
			}
			
			for(i = count_neighbors - 1; i >= 0; i--)
			{
				ne_in[i] = pq.top();
				if(i == dbs.m_minPts - 1)
					dbs.m_vec_f_core_distance[pid] = pq.top().dis;
				pq.pop();
			}

			//dbs.m_vec_f_core_distance[pid] = pq.top().dis;
			
			pq = priority_queue <kdtree2_result>();
		}
	}
//////////////////////////////////////////////////////////////////////////////
	void updateNeighbors_mst_based_parallel(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in, int tid, map<int,float>& reach_so_far_local_in)
	{
		map<int,int>::iterator other_ends_it_in;		
		map<int,float>::iterator reach_so_far_local_it_in;
		int i;
	
                float newRdist, oldRdist, new_dist_x_y, old_dist_x_y;
		for(i = 0; i < ne_in.size(); i++)
		{
			if (dbs.m_visited_per_thread[tid][ne_in[i].idx] == false)
			{
				newRdist = max(dbs.m_vec_f_core_distance[pid], ne_in[i].dis);
					
				if(orderSeeds_in.find_and_return_value(ne_in[i].idx, oldRdist))
				{
					if(newRdist < oldRdist)
					{
						orderSeeds_in.update(ne_in[i].idx, newRdist);
						other_ends_it_in = other_ends_in.find(ne_in[i].idx);

						other_ends_it_in->second = pid;
					}
					else if(newRdist == oldRdist)
					{
						other_ends_it_in = other_ends_in.find(ne_in[i].idx);	
						old_dist_x_y = dbs.m_kdtree->getDistance(other_ends_it_in->second, ne_in[i].idx);

						new_dist_x_y = dbs.m_kdtree->getDistance(pid, ne_in[i].idx);

						if(new_dist_x_y < old_dist_x_y || (new_dist_x_y == old_dist_x_y && pid < other_ends_it_in->second))
						{
							orderSeeds_in.update(ne_in[i].idx, newRdist);
							other_ends_it_in->second = pid;
						}
					}
				}
				else
				{
					orderSeeds_in.update(ne_in[i].idx, newRdist);
					other_ends_in.insert(pair<int,int>(ne_in[i].idx, pid));
				}

				
			}
			else
			{
			}
		}
	

	}

	void updateNeighbors_adjacency_list_based(ClusteringAlgo& dbs, int pid, vector <MST_Edge>* vec_edges_in, mutable_priority_queue<float,int>& orderSeeds_in)
	{

		// RECHECK: IF THE MST CONTAINS ANY CYCLE, THIS MIGHT NOT GIVE CRRECT ERESU
		int i;
		float newRdist, dis_x_y;
		MST_Edge* p_ed;

		for(i = 0; i < vec_edges_in->size(); i++)
		{
			p_ed = &(*vec_edges_in)[i];
			if(!dbs.m_visited[p_ed->m_other_end])
			{
				if(p_ed->m_f_edge_weight == FLT_MAX)
				{
					dis_x_y = dbs.m_kdtree->getDistance(pid, p_ed->m_other_end);
					if(dbs.m_vec_f_core_distance[pid] > dis_x_y)
						newRdist = dbs.m_vec_f_core_distance[pid];
					else
						newRdist = dis_x_y;
				}
				else
				{	
					newRdist = p_ed->m_f_edge_weight;
				}

				if(dbs.m_vec_f_reachability_distance[p_ed->m_other_end] == FLT_MAX)
				{
					dbs.m_vec_f_reachability_distance[p_ed->m_other_end] = newRdist;
					orderSeeds_in.insert(p_ed->m_other_end, newRdist);
				}
				else if(newRdist < dbs.m_vec_f_reachability_distance[p_ed->m_other_end]) // THIS SHOULDN"T BE THE CASE
				{
					dbs.m_vec_f_reachability_distance[p_ed->m_other_end] = newRdist;
					orderSeeds_in.update(p_ed->m_other_end, newRdist);
				}
			}
		}		
	}
////////////////////////////////////////////////////////////////////////////////////
	void compute_mst_parallel(ClusteringAlgo& dbs)
	{
		double start = omp_get_wtime();

		int sch, maxthreads = omp_get_max_threads();

		if(dbs.m_pts->m_i_num_points % maxthreads == 0)	
			sch = dbs.m_pts->m_i_num_points/maxthreads;	
		else
			sch = dbs.m_pts->m_i_num_points/maxthreads + 1;
		
		dbs.m_initial_prID.resize(dbs.m_pts->m_i_num_points, -1);

		dbs.m_vec_f_core_distance.resize(dbs.m_pts->m_i_num_points, FLT_MAX);
		dbs.m_visited.resize(dbs.m_pts->m_i_num_points, false);
		
		int i, j, pid, npid, cid, tid;
		vector<int>* ind = dbs.m_kdtree->getIndex();
		dbs.m_mst_per_thread.resize(maxthreads);
		
		dbs.m_connecting_edges_count.clear();
		dbs.m_connecting_edges_count.resize(maxthreads, 0);

		dbs.m_visited_per_thread.resize(maxthreads);
		
		#ifdef _BORUVKA_FILTER_AND_KRUSKAL
		dbs.m_mst_pq_per_thread.resize(maxthreads);
		#endif

		#pragma omp parallel private(i, pid, npid, tid) shared(sch, ind, dbs)
		{
			int lower, upper;
			tid = omp_get_thread_num();
			lower = sch * tid;
			upper = sch * (tid + 1);

			if(upper > dbs.m_pts->m_i_num_points)
				upper = dbs.m_pts->m_i_num_points;
	
			kdtree2_result_vector ne;
                	kdtree2_result_vector ne2;
                	ne.reserve(dbs.m_pts->m_i_num_points);
                	ne2.reserve(dbs.m_pts->m_i_num_points);
			
			mutable_priority_queue<float,int> orderSeeds;
			mutable_priority_queue<float,int> ::iterator it;
			map<int,int> other_ends;
			map<int,int>::iterator other_ends_it;
			
			map<int,float> reach_so_far_local;	

			MST_Node mst_node;

			dbs.m_visited_per_thread[tid].resize(dbs.m_pts->m_i_num_points, false);

			for(i = lower; i < upper; i++)
			{
				pid = (*ind)[i];
				dbs.m_initial_prID[pid] = tid;
			}
			
			#pragma omp barrier
			
			for(i = lower; i < upper; i++)
			{
				pid = (*ind)[i];
				if(!dbs.m_visited[pid])
				{
					dbs.m_visited[pid] = true;
                                	dbs.m_visited_per_thread[tid][pid] = true;
					ne.clear();
					dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne);
					computeCoreDistance(dbs, pid, ne);
					reach_so_far_local.insert(pair<int,float>(pid, FLT_MAX));

					if(dbs.m_vec_f_core_distance[pid] != FLT_MAX)
					{
						updateNeighbors_mst_based_parallel(dbs, pid, ne, orderSeeds, other_ends, tid, reach_so_far_local);
						while(orderSeeds.size() > 0)
						{
							it = orderSeeds.begin();
							npid = it->second;							
							reach_so_far_local.insert(pair<int,float>(npid, it->first));

							mst_node.m_reachable_to = npid;
							mst_node.m_f_best_distance = it->first;
                                        		other_ends_it = other_ends.find(npid);
                                        		if(other_ends_it != other_ends.end())
							{
                                                		mst_node.m_reachable_from = other_ends_it->second;
								other_ends.erase(other_ends_it);
							}
						
							#ifdef _KRUSKAL
								#ifdef _BORUVKA_FILTER_AND_KRUSKAL
								if(dbs.m_initial_prID[npid] == tid)
									dbs.m_mst_pq_per_thread[tid].push(mst_node);
									//dbs.m_mst_per_thread[tid].push_back(mst_node);
								else
								{
									#pragma omp critical   
                                                        		dbs.m_mst_pq.push(mst_node);
								}
								#else
								#pragma omp critical				
								dbs.m_mst_pq.push(mst_node);
								#endif
							#else
							dbs.m_mst_per_thread[tid].push_back(mst_node);
							#endif

							orderSeeds.erase(it->second);
						
							dbs.m_visited_per_thread[tid][npid] = true;

							if(dbs.m_initial_prID[npid] == tid)
							{
								dbs.m_visited[npid] = true;
								ne2.clear();
								dbs.m_kdtree->r_nearest_around_point(npid, 0, dbs.m_epsSquare, ne2);
								computeCoreDistance(dbs, npid, ne2);
								if(dbs.m_vec_f_core_distance[npid] != FLT_MAX)
									updateNeighbors_mst_based_parallel(dbs, npid, ne2, orderSeeds, other_ends, tid, reach_so_far_local);
							}
							else
							{
								dbs.m_connecting_edges_count[tid]++;
							}
						}
					}
				}
			}

			orderSeeds.clear();
			other_ends.clear();
			ne.clear();
			ne2.clear();
		}	

                cout << "Computing local MSTs in parallel took: " << omp_get_wtime() - start<< " seconds." << endl;

		#ifdef _DEBUG
		int edges_count = 0, total_connecting_edges_count = 0;

		
                for(i = 0; i < maxthreads; i++)
		{
			#ifdef _BORUVKA_FILTER_AND_KRUSKAL
			edges_count += dbs.m_mst_pq_per_thread[i].size();
			#else
			edges_count += dbs.m_mst_per_thread[i].size();
                        #endif
			total_connecting_edges_count += dbs.m_connecting_edges_count[i];
		}
		#ifdef _KRUSKAL
		#ifdef _BORUVKA_FILTER_AND_KRUSKAL
		edges_count += dbs.m_mst_pq.size();
		#else
		edges_count = dbs.m_mst_pq.size();
		#endif
                cout << "Total edges using Kruskal " << edges_count << " connecting " << total_connecting_edges_count << " rest " << edges_count - total_connecting_edges_count << endl;
		#else
		cout << "Total edges " << edges_count << " connecting " << total_connecting_edges_count << " rest " << edges_count - total_connecting_edges_count << endl;
		#endif
		#endif

		#ifdef _KRUSKAL
		compute_mst_global_from_mst_local_kruskal(dbs);
		#else
		compute_mst_global_from_mst_local(dbs);
		#endif		
		
		#ifdef _DEBUG
		cout << "Total Time taken to compute the final MST parallel " << omp_get_wtime() - start << endl;
		#endif

		ind = NULL;
	}

	void compute_mst_global_from_mst_local(ClusteringAlgo& dbs)
	{
		double start = omp_get_wtime();
                dbs.m_visited.clear();
		dbs.m_visited.resize(dbs.m_pts->m_i_num_points, false);
		dbs.m_visited_per_thread.clear(); // RECHECK IF NEEDED		

		// convert the mst to an adjacency list
		compute_adjacency_list_from_mst_local_parallel(dbs);
		#ifdef _DEBUG_DETAILS
		print_adjacency_list_build_from_mst(dbs);
		#endif		

		int i, j, pid, npid;
		int notvisited = 0;
		mutable_priority_queue<float,int> orderSeeds;
                mutable_priority_queue<float,int> ::iterator it;
		vector <MST_Edge>* vec_edges;
		vector <MST_Edge>* vec_edges_next;
		map<int,int> other_ends;
                map<int,int>::iterator other_ends_it;

		vector<int>* ind = dbs.m_kdtree->getIndex();
		MST_Node mst_node;

		for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
		{
			pid = (*ind)[i];
			if (!dbs.m_visited[pid])
			{
                                #ifdef _DEBUG
                                notvisited++;
                                #endif                               
                                dbs.m_visited[pid] = true;			
				vec_edges = &dbs.m_vec_vec_adjacency_list[pid];
				if(vec_edges->size() > 0)
				{
				
					updateNeighbors_adjacency_list_based_global_mst_from_local_mst(dbs, pid, vec_edges, orderSeeds, other_ends);
					while(orderSeeds.size() > 0)
					{
						it = orderSeeds.begin();
						npid = it->second;

						mst_node.m_reachable_to = npid;
                                                mst_node.m_f_best_distance = it->first;
                                                
						orderSeeds.erase(it->second);

					        other_ends_it = other_ends.find(npid);
                                                if(other_ends_it != other_ends.end())
                                                {
                                                	mst_node.m_reachable_from = other_ends_it->second;
                                                	other_ends.erase(other_ends_it);
                                                }

						dbs.m_mst.push_back(mst_node);
						
						dbs.m_visited[npid] = true;
						vec_edges_next = &dbs.m_vec_vec_adjacency_list[npid];
						if(vec_edges_next->size() > 0)
						{
							updateNeighbors_adjacency_list_based_global_mst_from_local_mst(dbs, npid, vec_edges_next, orderSeeds, other_ends);
						}		
					}
				}
			}
		}

		orderSeeds.clear();
		other_ends.clear();
		vec_edges_next = NULL;
		vec_edges = NULL;
		ind = NULL;

                cout << "Computing global MST from local MSTs took: " << omp_get_wtime() - start << " seconds."<< endl;
		//cout << "Edges in global MST " << dbs.m_mst.size() << endl;
	}

	void compute_mst_global_from_mst_local_kruskal(ClusteringAlgo& dbs)
	{
		double start = omp_get_wtime();

		#ifdef _BORUVKA_FILTER_AND_KRUSKAL
		int maxthreads = omp_get_max_threads();
		#endif

		vector <int> parents;
		parents.resize(dbs.m_pts->m_i_num_points);
		vector <bool> touched;
		touched.resize(dbs.m_pts->m_i_num_points, false);

		int i;
		#ifdef _KRUSKAL_PARALLEL_MERGING
		omp_lock_t *nlocks;
		nlocks = (omp_lock_t *) malloc(dbs.m_pts->m_i_num_points * sizeof(omp_lock_t));
		#endif

		#ifdef _KRUSKAL_PARALLEL_MERGING
		#pragma omp parallel for private(i) shared(dbs, parents, nlocks)
		#else
		#pragma omp parallel for private(i) shared(dbs, parents)
		#endif
		for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
		{
			parents[i] = i;
			#ifdef _KRUSKAL_PARALLEL_MERGING
			omp_init_lock(&nlocks[i]);
			#endif
		}
		
		int tid;
		MST_Node mst_node;
		int u, v, w, root_u, root_v, z, elements_pq;
		
		#ifdef _BORUVKA_FILTER_AND_KRUSKAL
		int filtered_edges = 0;
		#pragma omp parallel for private(z, tid, i, mst_node, u, v, w, root_u, root_v) shared(touched, filtered_edges, dbs, parents, maxthreads)
		for(tid = 0; tid < maxthreads; tid++)
		{
			while(dbs.m_mst_pq_per_thread[tid].size() > 0)
			{
				mst_node = dbs.m_mst_pq_per_thread[tid].top();
				dbs.m_mst_pq_per_thread[tid].pop();

				u = mst_node.m_reachable_to;
                        	v = mst_node.m_reachable_from;
                        	w = mst_node.m_f_best_distance;
				
				if(touched[u] == true || touched[v] == true)
                        	{
					#pragma omp critical
                                	dbs.m_mst_pq.push(mst_node);
                                	continue;
                        	}
				
				touched[u] = true;
				touched[v] = true;

				root_u = u;
                        	root_v = v;
		
				bool success = false;
				while(parents[root_u] != parents[root_v])
				{
					if(parents[root_u] < parents[root_v])
					{
						if(root_u == parents[root_u])
						{
                                                        parents[root_u] = parents[root_v];
							#pragma omp critical
                                                        dbs.m_mst.push_back(mst_node);
							#pragma omp atomic
							filtered_edges++;
							success = true;
                                                }
						
						z = parents[root_u];
						parents[root_u] = parents[root_v];
						root_u = z;
					}
					else
					{
						if(root_v == parents[root_v])
						{
        	                                	parents[root_v] = parents[root_u];
                	                                #pragma omp critical
							dbs.m_mst.push_back(mst_node);
							#pragma omp atomic
                                                        filtered_edges++;
							success = true;
						}
							
						z = parents[root_v];
						parents[root_v] = parents[root_u];
						root_v = z;	
					}
				}

				if(success == false)
				{
					#pragma omp critical
					dbs.m_mst_pq.push(mst_node);
				}
			}
		}
		cout << "Filterted edges " << filtered_edges << endl;
		#endif

		elements_pq = dbs.m_mst_pq.size();
	
		#ifdef _KRUSKAL_PARALLEL_MERGING
		#pragma omp parallel for private(i, mst_node, u, v, w, z, root_u, root_v) shared(dbs, parents, nlocks)
		#endif
                for(i = 0; i < elements_pq; i++)
		{
			#ifdef _KRUSKAL_PARALLEL_MERGING
			#pragma omp critical
			#endif
			{
				mst_node = dbs.m_mst_pq.top();
				dbs.m_mst_pq.pop();
			}

			u = mst_node.m_reachable_to;
                        v = mst_node.m_reachable_from;
                        w = mst_node.m_f_best_distance;

			root_u = u;
                        root_v = v;
			
			while(parents[root_u] != parents[root_v])
			{
				if(parents[root_u] < parents[root_v])
				{
					if(root_u == parents[root_u])
					{
						#ifdef _KRUSKAL_PARALLEL_MERGING
						omp_set_lock(&nlocks[root_u]);
						bool p_set = false;
						if(root_u == parents[root_u])
						{		
							parents[root_u] = parents[root_v];
							p_set = true;
						}
						
						omp_unset_lock(&nlocks[root_u]);
						if(p_set == true)
						{
							#pragma omp critical
							dbs.m_mst.push_back(mst_node);
							break;
						}
						#else
						if(root_u == parents[root_u])
                                                {               
                                                        parents[root_u] = parents[root_v];
                                                        dbs.m_mst.push_back(mst_node);
                                                }
						#endif
					}
						
					z = parents[root_u];
					parents[root_u] = parents[root_v];
					root_u = z;
				}
				else
				{
					if(root_v == parents[root_v])
					{
						#ifdef _KRUSKAL_PARALLEL_MERGING
						omp_set_lock(&nlocks[root_v]);
						bool p_set = false;
						if(root_v == parents[root_v])
						{		
							parents[root_v] = parents[root_u];
							p_set = true;
						}
						omp_unset_lock(&nlocks[root_v]);
						if(p_set == true)
						{
							#pragma omp critical
                                                        dbs.m_mst.push_back(mst_node);
							break;
						}
						#else
						if(root_v == parents[root_v])
                                                {
                                                        parents[root_v] = parents[root_u];
                                                        dbs.m_mst.push_back(mst_node);
                                                }
						#endif
					}
							
					z = parents[root_v];
					parents[root_v] = parents[root_u];
					root_v = z;	
				}
			}						
		}

		#ifdef _KRUSKAL_PARALLEL_MERGING		
		free(nlocks);
		#endif

		parents.clear();
                cout << "Computing global MST from local MSTs in parallel took: " << omp_get_wtime() - start << " seconds."<< endl;
		//cout << "Edges in global MST " << dbs.m_mst.size() << endl;
	}

	void compute_adjacency_list_from_mst_local_parallel(ClusteringAlgo& dbs)
	{
		#ifdef _DEBUG
		double start = omp_get_wtime();
		#endif

		// the mst might contain cycle
		dbs.m_vec_vec_adjacency_list.clear();
		// alocate for the first dimention
		dbs.m_vec_vec_adjacency_list.resize(dbs.m_pts->m_i_num_points);
		
		MST_Edge me, me_other;
		int i, tid, j;
		omp_lock_t *nlocks;
  		nlocks = (omp_lock_t *) malloc(dbs.m_pts->m_i_num_points * sizeof(omp_lock_t));
		
		#pragma omp parallel for private(i) shared(dbs, nlocks)
		for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
			omp_init_lock(&nlocks[i]); // Initialize locks

		#pragma omp parallel private(i, j, tid, me, me_other) shared(dbs, nlocks)
		{
			tid = omp_get_thread_num();
			for(j = 0; j < dbs.m_mst_per_thread[tid].size(); j++)
			{
				me.m_other_end = dbs.m_mst_per_thread[tid][j].m_reachable_to;
				me.m_f_edge_weight = dbs.m_mst_per_thread[tid][j].m_f_best_distance;
				omp_set_lock(&nlocks[dbs.m_mst_per_thread[tid][j].m_reachable_from]); // locking
				dbs.m_vec_vec_adjacency_list[dbs.m_mst_per_thread[tid][j].m_reachable_from].push_back(me);
				omp_unset_lock(&nlocks[dbs.m_mst_per_thread[tid][j].m_reachable_from]); // un-locking

				me_other.m_other_end = dbs.m_mst_per_thread[tid][j].m_reachable_from;
				me_other.m_f_edge_weight = FLT_MAX;
				omp_set_lock(&nlocks[dbs.m_mst_per_thread[tid][j].m_reachable_to]); // locking
				dbs.m_vec_vec_adjacency_list[dbs.m_mst_per_thread[tid][j].m_reachable_to].push_back(me_other);
				omp_unset_lock(&nlocks[dbs.m_mst_per_thread[tid][j].m_reachable_to]); // un-locking
			}
		}
		
		free(nlocks);
             
		// filter out same edge
		map<int, float> filter_identical;
		map<int, float>::iterator filter_identical_it;

		#pragma omp parallel for private(i, j, filter_identical, filter_identical_it, me) shared(dbs)
		for(i = 0; i < dbs.m_vec_vec_adjacency_list.size(); i++)
		{
			for(j = 0; j < dbs.m_vec_vec_adjacency_list[i].size(); j++)
			{
				filter_identical_it = filter_identical.find(dbs.m_vec_vec_adjacency_list[i][j].m_other_end);
				if(filter_identical_it != filter_identical.end())
				{
					if(filter_identical_it->second > dbs.m_vec_vec_adjacency_list[i][j].m_f_edge_weight)
                                        	filter_identical_it->second = dbs.m_vec_vec_adjacency_list[i][j].m_f_edge_weight;
				}
				else
					filter_identical.insert(pair<int,float>(dbs.m_vec_vec_adjacency_list[i][j].m_other_end, dbs.m_vec_vec_adjacency_list[i][j].m_f_edge_weight));
			}

			dbs.m_vec_vec_adjacency_list[i].clear();

			// sort based on the closest distance
			priority_queue<kdtree2_result> pq;
			kdtree2_result kr;

			for(filter_identical_it = filter_identical.begin(); filter_identical_it != filter_identical.end(); filter_identical_it++)
			{
				kr.dis = dbs.m_kdtree->getDistance(i, filter_identical_it->first);
				kr.idx = filter_identical_it->first;	
				pq.push(kr);
			}

			dbs.m_vec_vec_adjacency_list[i].resize(pq.size(), me);
			int count_neighbors_local = pq.size();
			for(j = count_neighbors_local - 1; j >= 0; j--)
			{
				filter_identical_it = filter_identical.find(pq.top().idx);
				me.m_other_end = filter_identical_it->first;
                                me.m_f_edge_weight = filter_identical_it->second;
                                dbs.m_vec_vec_adjacency_list[i][j] = me;
				
				pq.pop();
			}	
			
			filter_identical.clear();
		}	
   
		#ifdef _DEBUG
                cout << "Time taken to compute adjacency_list from local MSTs parallel " << omp_get_wtime() - start << endl;
                #endif		
	}

	void updateNeighbors_adjacency_list_based_global_mst_from_local_mst(ClusteringAlgo& dbs, int pid, vector <MST_Edge>* vec_edges_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in)
	{
		map<int,int>::iterator other_ends_it_in;		
		int i;
		float newRdist, dis_x_y, oldRdist, new_dist_x_y, old_dist_x_y;
                MST_Edge* p_ed;

		for(i = 0; i < vec_edges_in->size(); i++)
		{
			p_ed = &(*vec_edges_in)[i];
                        if(!dbs.m_visited[p_ed->m_other_end])
			{
				if(p_ed->m_f_edge_weight == FLT_MAX)
				{
					dis_x_y = dbs.m_kdtree->getDistance(pid, p_ed->m_other_end);

					if(dbs.m_vec_f_core_distance[pid] > dis_x_y)
						newRdist = dbs.m_vec_f_core_distance[pid];
					else
						newRdist = dis_x_y;
				}
				else
				{
					newRdist = p_ed->m_f_edge_weight;
				}

				if(orderSeeds_in.find_and_return_value(p_ed->m_other_end, oldRdist))
				{
					if(newRdist < oldRdist)
					{
						orderSeeds_in.update(p_ed->m_other_end, newRdist);
						other_ends_it_in = other_ends_in.find(p_ed->m_other_end);

						other_ends_it_in->second = pid;
                                                #ifdef _DEBUG
                                                if(other_ends_it_in->first != p_ed->m_other_end)
                                                        cout << "SOMETHING is Wrong at updateNeighbors_mst_based" << endl;
                                                #endif  
					}
					else if(newRdist == oldRdist)
					{
						other_ends_it_in = other_ends_in.find(p_ed->m_other_end);	
						old_dist_x_y = dbs.m_kdtree->getDistance(other_ends_it_in->second, p_ed->m_other_end);

						new_dist_x_y = dbs.m_kdtree->getDistance(pid, p_ed->m_other_end);

						if(new_dist_x_y < old_dist_x_y || (new_dist_x_y == old_dist_x_y && pid < other_ends_it_in->second))
						{
							orderSeeds_in.update(p_ed->m_other_end, newRdist);
							other_ends_it_in->second = pid;
							#ifdef _DEBUG
                                                	if(other_ends_it_in->first != p_ed->m_other_end)
                                                        	cout << "SOMETHING is Wrong at updateNeighbors_mst_based" << endl;
                                                	#endif
						}
					}
				}
				else
				{
					orderSeeds_in.update(p_ed->m_other_end, newRdist);
					other_ends_in.insert(pair<int,int>(p_ed->m_other_end, pid));
				}
				
			}
		}

		p_ed = NULL;		
	}
////////////////////////////////////////////////////////////////////////////////////	
};

