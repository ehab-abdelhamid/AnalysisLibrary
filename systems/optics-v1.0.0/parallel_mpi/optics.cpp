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


#include "optics.h"

namespace NWUClustering
{
    	void ClusteringAlgo::set_dbscan_params(double eps, int minPts) //, int messages_per_round, int compression, int collective_comm, int parcent_of_data)
	{
		m_epsSquare =  eps * eps;
		m_minPts =  minPts;
		m_messages_per_round = -1; //messages_per_round;
		m_compression = 0; //compression;
		m_collective_comm = 1; //collective_comm;
		m_parcent_of_data = 100; //parcent_of_data;
	}

	ClusteringAlgo::~ClusteringAlgo()
	{
		m_vec_f_core_distance.clear();
                m_visited.clear();
                m_mst.clear();
		
                m_vec_f_core_distance_root.clear();
                m_vec_vec_adjacency_list_root.clear();

		if(m_pq_mst_node != NULL)
			delete m_pq_mst_node;
		
		m_vec_global_point_ID_gathered.clear();
                m_vec_f_core_distance_gathered.clear();

		m_vec_gid_parent_gid_pair.clear();
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void saveMST(ClusteringAlgo& dbs, ostream& o)
        {
                int i, u, v;
		float w, core_dis; 
		o << dbs.m_pts->m_total_point_count << " " << dbs.m_mst_root.size() << endl; // number of points in the dataset
                for(i = 0; i < dbs.m_mst_root.size(); i++)
		{
			u = dbs.m_mst_root[i].m_reachable_to;
			v = dbs.m_mst_root[i].m_reachable_from;
			w = dbs.m_mst_root[i].m_f_between_distance;
			core_dis = dbs.m_vec_f_core_distance_root[v];

			o << u << " " << v << " " << sqrt(w) << " " << sqrt(core_dis) << endl; 
		}
        }

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
			
			pq = priority_queue <kdtree2_result>();
		}	
	}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void compute_mst_parallel(ClusteringAlgo& dbs)
	{
		//#ifdef _DEBUG
		double start = MPI_Wtime();
		//#endif
		int owned_local_points, total_local_points;

		owned_local_points = dbs.m_pts->m_i_num_local_points;
		total_local_points = dbs.m_pts->m_i_num_points;

		int rank, nproc;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		
		// we don't care about the core distance of the non local points
		dbs.m_vec_f_core_distance.clear();
		dbs.m_vec_f_core_distance.resize(owned_local_points, FLT_MAX);

		dbs.m_visited.clear();
		dbs.m_visited.resize(total_local_points, false);	
		
		int i, j, pid, npid, tid;
	
		kdtree2_result_vector ne;
		kdtree2_result_vector ne2;
		ne.reserve(total_local_points);
		ne2.reserve(total_local_points);

		vector<int>* ind = dbs.m_kdtree->getIndex();
		
		dbs.m_connecting_edges_count = 0;
		
		mutable_priority_queue <float, int> orderSeeds;
                mutable_priority_queue <float, int>::iterator it;
                map <int, int> other_ends;
                map <int, int>::iterator other_ends_it;
		map <int, float> reach_so_far_local;

		MST_Node mst_node;

		dbs.m_pq_mst_node = new priority_queue<MST_Node>; 

		for(i = 0; i < total_local_points; i++)
		{
			// process only the owned points
			pid = (*ind)[i];
			if(!dbs.m_visited[pid])
			{
				if(rank == dbs.m_pts->m_vec_processor_ID[pid])
				{
					//pid must be a local owned point
					dbs.m_visited[pid] = true;
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

                                                        mst_node.m_reachable_to = dbs.m_pts->m_vec_global_point_ID[npid];
                                                        //mst_node.m_f_best_distance = it->first;
                                                        other_ends_it = other_ends.find(npid);
                                                        if(other_ends_it != other_ends.end())
                                                        {
                                                                mst_node.m_reachable_from = dbs.m_pts->m_vec_global_point_ID[other_ends_it->second];
								mst_node.m_f_between_distance = dbs.m_kdtree->getDistance(other_ends_it->second, npid);
                                                                other_ends.erase(other_ends_it);
                                                        }
							
							//#ifdef _KRUSKAL
							#if defined(_KRUSKAL) || defined(_KRUSKAL_SEQUENTIAL) || defined(_KRUSKAL_SMART) || defined(_BORUVKA_FILTER_AND_KRUSKAL)
                        				dbs.m_pq_mst_node->push(mst_node);
							#else
							dbs.m_mst.push_back(mst_node);
							#endif

							orderSeeds.erase(it->second);
							dbs.m_visited[npid] = true;
                                			if(rank == dbs.m_pts->m_vec_processor_ID[npid])
							{
								ne2.clear();
								dbs.m_kdtree->r_nearest_around_point(npid, 0, dbs.m_epsSquare, ne2);
								computeCoreDistance(dbs, npid, ne2);
								if(dbs.m_vec_f_core_distance[npid] != FLT_MAX)
								{
									updateNeighbors_mst_based_parallel(dbs, npid, ne2, orderSeeds, other_ends, tid, reach_so_far_local);
								}
							}
							else
								dbs.m_connecting_edges_count++;
						}
					}
					
				}
			}
		}

		ne.clear();
                ne2.clear();
                orderSeeds.clear();
                other_ends.clear();
                reach_so_far_local.clear();	
		ind = NULL;
		
		MPI_Barrier(MPI_COMM_WORLD);
		
                if(rank == proc_of_interest) cout << "Computing the local MSTs in parallel took: " << MPI_Wtime() - start << " seconds."<< endl;

		#ifdef _BORUVKA_FILTER_AND_KRUSKAL
			#ifdef _CUT_OFF_ACTIVE
			#ifdef _USE_ROUND_TIME
 			//if(rank == proc_of_interest) cout << " bor_fil_kruskal_cut_time ";
			#endif
			#ifdef _USE_ROUND_COUNT
                        //if(rank == proc_of_interest) cout << " bor_fil_kruskal_cut_round ";
                        #endif
			#endif
		local_mst_merging_boruvka_filter_kruskal(dbs);
		#elif defined(_KRUSKAL)
		//if(rank == proc_of_interest) cout << " pairwise_kruskal ";
		pairwise_local_mst_merging_kruskal(dbs);
		#elif defined(_KRUSKAL_SMART)
			#ifdef _CUT_OFF_ACTIVE
			#ifdef _USE_ROUND_TIME
 			//if(rank == proc_of_interest) cout << " pairwise_kruskal_cut_time ";
			#endif
			#ifdef _USE_ROUND_COUNT
                        //if(rank == proc_of_interest) cout << " pairwise_kruskal_cut_round ";
                        #endif
			#endif	
		pairwise_local_mst_merging_kruskal_smart(dbs);
		#elif defined(_KRUSKAL_SEQUENTIAL)
		sequential_local_mst_merging_kruskal(dbs);
		#else
		gather_local_msts_and_compute_adj_root(dbs);
                compute_mst_global_from_mst_local(dbs);
		#endif
		
		#ifdef _DEBUG
		if(rank == proc_of_interest) cout << "Total time taken to compute the global MST in parallel " << MPI_Wtime() - start << endl;
		#endif
	}

	void local_mst_merging_boruvka_filter_kruskal(ClusteringAlgo& dbs)
	{
		double start = MPI_Wtime();
                double inter = start, inter2 = start;
                int rank, nproc, i, j, round_pos, round_count = 0, proc_count;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		map <int, int> gid_local_id;

		for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
		{
			gid_local_id.insert(pair<int,int>(dbs.m_pts->m_vec_global_point_ID[i], i));	
		}
		
		vector <int> vec_parents;
		vec_parents.reserve(dbs.m_pts->m_i_num_local_points);
		map <int, int> gid_parent_gid;

		for(i = 0; i < dbs.m_pts->m_i_num_local_points; i++)
		{
			vec_parents.push_back(i);
			gid_parent_gid.insert(pair<int,int>(dbs.m_pts->m_vec_global_point_ID[i], dbs.m_pts->m_vec_global_point_ID[i]));
		}
		
		int local_mst_edge_count_this = dbs.m_pq_mst_node->size();
                int global_mst_edge_count;

                MPI_Allreduce(&local_mst_edge_count_this, &global_mst_edge_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in local MSTs " << global_mst_edge_count << endl;
		#endif	
			
		priority_queue<MST_Node>* pq_mst_node;
                pq_mst_node = new priority_queue<MST_Node>;	
		MST_Node mst_node;

		int pid, npid;
		float distance;
		map <int, int>::iterator it_gid_local_id;
		map <int, int>::iterator it_gid_parent_gid1;
		map <int, int>::iterator it_gid_parent_gid2;	
	
		vector<bool> vec_touched;
		vec_touched.resize(dbs.m_pts->m_i_num_points, false);

		int filtered_edges = 0, root_u, root_v, root;

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Preprocessing on local MSTs (mappings) " << MPI_Wtime() -inter2 << endl;
		inter2 = MPI_Wtime();
		#else
		if(rank == proc_of_interest) cout << " prep_lo_msts_t " << MPI_Wtime() -inter2;
		inter2 = MPI_Wtime();
		#endif	
	
		while(dbs.m_pq_mst_node->size() > 0)
		{
			mst_node = dbs.m_pq_mst_node->top();
			dbs.m_pq_mst_node->pop();
			
			it_gid_local_id = gid_local_id.find(mst_node.m_reachable_to);
			if(it_gid_local_id == gid_local_id.end())
				cout << "SOMTHING IS WRONG Boruvka" << endl;
			npid = it_gid_local_id->second;
			
			it_gid_local_id = gid_local_id.find(mst_node.m_reachable_from);                        
                        if(it_gid_local_id == gid_local_id.end())
                                cout << "SOMTHING IS WRONG Boruvka" << endl;
                        pid = it_gid_local_id->second;
			
			if(vec_touched[pid] == true || vec_touched[npid] == true)
			{
				pq_mst_node->push(mst_node);
				continue;
			}
	
			if(rank != dbs.m_pts->m_vec_processor_ID[npid])
			{
				pq_mst_node->push(mst_node);
			}
			else
			{
				root_u = pid;
				while(root_u != vec_parents[root_u])
					root_u = vec_parents[root_u];

				root_v = npid;
                                while(root_v != vec_parents[root_v])
                                        root_v = vec_parents[root_v];

				if(root_u != root_v)
				{
					//this edge should directly go to the final solution
					//now perform union
					it_gid_parent_gid1 = gid_parent_gid.find(dbs.m_pts->m_vec_global_point_ID[root_u]);
					if(it_gid_parent_gid1 == gid_parent_gid.end() || it_gid_parent_gid1->second != dbs.m_pts->m_vec_global_point_ID[root_u])
						cout << "Something is wrong Boruvka" << endl;

					it_gid_parent_gid2 = gid_parent_gid.find(dbs.m_pts->m_vec_global_point_ID[root_v]);
                                        if(it_gid_parent_gid2 == gid_parent_gid.end() || it_gid_parent_gid2->second != dbs.m_pts->m_vec_global_point_ID[root_v])
                                                cout << "Something is wrong Boruvka" << endl;
					
					if(it_gid_parent_gid1->second < it_gid_parent_gid2->second)
					{
						it_gid_parent_gid1->second = it_gid_parent_gid2->second;
						vec_parents[root_u] = root_v;
						root = root_v;
					}
					else
					{
						it_gid_parent_gid2->second = it_gid_parent_gid1->second;
                                                vec_parents[root_v] = root_u;
                                                root = root_u;
					}

					vec_touched[pid] = true;
					vec_touched[npid] = true;

					root_u = pid;
                                	while(vec_parents[root_u] != root)
                                        	vec_parents[root_u] = root;

					root_v = npid;
                                	while(vec_parents[root_v] != root)
                                        	vec_parents[root_v] = root;
				
					filtered_edges++;
					dbs.m_mst_root.push_back(mst_node);
				}
				else
					pq_mst_node->push(mst_node);
			}	
		}

		delete dbs.m_pq_mst_node;
		dbs.m_pq_mst_node = pq_mst_node;
		pq_mst_node = NULL;

		int total_filtered_edges, rest_edges = dbs.m_pq_mst_node->size(), total_rest_edges, filtered_edges_per_thread[nproc]; 
		MPI_Allreduce(&filtered_edges, &total_filtered_edges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Gather(&filtered_edges, 1, MPI_INT, &filtered_edges_per_thread[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);	
		MPI_Allreduce(&rest_edges, &total_rest_edges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in local MSTs filtered " << total_filtered_edges << endl;
		if(rank == proc_of_interest) cout << "Total edges in local MSTs NOT filtered " << total_rest_edges << endl;
		if(rank == proc_of_interest) cout << "Taken time to filter edge " << MPI_Wtime() - inter2 << endl;
		if(rank == proc_of_interest) cout << "Total time by filtering edges " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
		#else
		if(rank == proc_of_interest) cout << " ed_lo_fil " << total_filtered_edges;
                if(rank == proc_of_interest) cout << " ed_lo_not_fil " << total_rest_edges;
                if(rank == proc_of_interest) cout << " fil_t " << MPI_Wtime() - inter2;
                if(rank == proc_of_interest) cout << " to_fil_t " << MPI_Wtime() - inter;
                inter = MPI_Wtime();
		#endif

		// now have to gather the edges and perform union operations...

		//vector <int> vec_gid_send;
		//vector <int> vec_parents_gid_send;
		//vec_gid_send.reserve(dbs.m_pts->m_i_num_local_points);
		//vec_parents_gid_send.reserve(dbs.m_pts->m_i_num_local_points);
		dbs.m_vec_gid_parent_gid_pair.reserve(gid_parent_gid.size());
		//i = 0;
		for(it_gid_parent_gid1 = gid_parent_gid.begin(); it_gid_parent_gid1 != gid_parent_gid.end(); it_gid_parent_gid1++)
		{
			if(it_gid_parent_gid1->first == it_gid_parent_gid1->second)
				continue;
			dbs.m_vec_gid_parent_gid_pair.push_back(it_gid_parent_gid1->first);
			dbs.m_vec_gid_parent_gid_pair.push_back(it_gid_parent_gid1->second);
			//i++;
		}

		//#ifdef _DEBUG
		//cout << "gid and parents gid pairs count " << i << endl;
		//#endif
		//if(i != dbs.m_pts->m_i_num_local_points)
		//	cout << "SOmething is wrong Boruvka" << endl;
		vec_parents.clear();
		gid_parent_gid.clear();
		gid_local_id.clear();
		vec_touched.clear();

		int global_owned_local_point_count;
                MPI_Allreduce(&dbs.m_pts->m_i_num_local_points, &global_owned_local_point_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                dbs.m_global_owned_local_point_count_root = global_owned_local_point_count;

		dbs.m_all_gid_gathered_upper = dbs.m_pts->m_all_local_point_gID_upper;
		dbs.m_all_gid_gathered_lower = dbs.m_pts->m_all_local_point_gID_lower;

		// copy the mst to root, but only for the first time;
                proc_count = nproc;

                while((proc_count % 2 == 0) && (proc_count > 1))
                {
                        round_count++;
                        proc_count = proc_count / 2;
                }

                #ifdef _DEBUG
                if(proc_count != 1)
                        cout << "Number of procs is not a multiple of 2" << endl;
                #endif

                vec_parents.reserve(global_owned_local_point_count);
                for(i = 0; i < global_owned_local_point_count; i++)
                        vec_parents.push_back(i);
				
		MPI_Datatype mst_node_type;
                MPI_Type_contiguous(sizeof(MST_Node), MPI_BYTE, &mst_node_type);
                MPI_Type_commit(&mst_node_type);

		int rank_from, rank_to, recv_size[5], send_size[5]; //mst_node_recv_size, mst_node_send_size;
		int tag = 200, z;
                MPI_Request req_send[3 * nproc], req_recv[3 * nproc];
                MPI_Status stat_send[3 * nproc], stat_recv[3 * nproc], stat_recv_single;
		int send_count, recv_count;
		//MST_Node mst_node;
		int rtag, rsource, rpos;
		double comp_time, comp_time_total = 0, comp_time2;
		
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "boruvka preprocessing time before starting pairwise merging " << MPI_Wtime() - inter << endl;
		inter = MPI_Wtime();
		inter2 = inter;
		#else
		if(rank == proc_of_interest) cout << " bor_pre_t " << MPI_Wtime() - inter;
                inter = MPI_Wtime();
                inter2 = inter;
		#endif

		//cout << "rank " << rank  << " MST size " << dbs.m_pq_mst_node->size() << endl;
		for(round_pos = 0; round_pos < round_count; round_pos++)
		{
			#ifdef _DEBUG
			int local_mst_edge_count_round = dbs.m_pq_mst_node->size(), local_mst_edge_count_round_total;			
			MPI_Allreduce(&local_mst_edge_count_round, &local_mst_edge_count_round_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	
			if(rank == proc_of_interest) cout << "Round " << round_pos << " local_mst_edge_count_round " << local_mst_edge_count_round_total << endl;	
			#endif

			//tag  = 200; // a random start
			if(rank % (int)(pow (2, round_pos)) == 0)
			{
				// rank is either a sender or a receiver
				if(rank % (int)(pow (2, round_pos + 1)) == 0)
				{
					// rank is receiver
					rank_from = rank + (int)(pow (2, round_pos));					
					//dbs.m_vec_gid_parent_gid_pair
					
					recv_size[1] = 0; // dbs.m_vec_gid_parent_gid_pair count
					recv_size[2] = 0; // corresponing global ID count
					recv_size[3] = 0; // global ID upper limit
					recv_size[4] = 0; // global ID lower limit
		
					MPI_Recv(&recv_size[0], 5, MPI_INT, rank_from, tag, MPI_COMM_WORLD, &stat_recv[0]);
					//cout << "rank " << rank << " round_pos " << round_pos << " tag " << tag << " rank_to " << rank_from << endl;

					recv_count = 0;
					int cur_vec_gid_parent_gid_pair_size = dbs.m_vec_gid_parent_gid_pair.size();
					if(recv_size[1] > 0)
					{
						dbs.m_vec_gid_parent_gid_pair.resize(cur_vec_gid_parent_gid_pair_size + recv_size[1]);
						MPI_Irecv(&dbs.m_vec_gid_parent_gid_pair[cur_vec_gid_parent_gid_pair_size], recv_size[1], MPI_INT, rank_from, tag + 1, MPI_COMM_WORLD, &req_recv[recv_count++]);
					}
					
					vector <MST_Node> vec_mst_node_recv;
					if(recv_size[2] > 0)
					{
                                        	vec_mst_node_recv.resize(recv_size[2]);
						MPI_Irecv(&vec_mst_node_recv[0], recv_size[2], mst_node_type, rank_from, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]); 
					}
						
					if(dbs.m_all_gid_gathered_upper < recv_size[3])
						dbs.m_all_gid_gathered_upper = recv_size[3];

					if(dbs.m_all_gid_gathered_lower > recv_size[4])
						dbs.m_all_gid_gathered_lower = recv_size[4];

					// reset that range from parents vector
					for(i = dbs.m_all_gid_gathered_lower; i < dbs.m_all_gid_gathered_upper; i++)
						vec_parents[i] = i;

					if(recv_count > 0)
                                                MPI_Waitall(recv_count, &req_recv[0], &stat_recv[0]);

					//#ifdef _DEBUG
					comp_time = MPI_Wtime();
					//#endif
					
					for(i = 0; i < dbs.m_vec_gid_parent_gid_pair.size(); i += 2)
						vec_parents[dbs.m_vec_gid_parent_gid_pair[i]] = dbs.m_vec_gid_parent_gid_pair[i + 1];

				
					for(j = 0; j < vec_mst_node_recv.size(); j++)
                                        	dbs.m_pq_mst_node->push(vec_mst_node_recv[j]);

					vec_mst_node_recv.clear();
					//cout << "rank= " << rank << " round_pos " << round_pos << " MST tree size " << dbs.m_pq_mst_node->size() << " vector " << vec_mst_node_recv.size() <<endl;
					// run Kruskals algorithm on the edges in the priority queue
					priority_queue<MST_Node>* pq_mst_node;
					pq_mst_node = new priority_queue<MST_Node>;

					if(round_pos == round_count - 1)
						dbs.m_mst_root.reserve(global_owned_local_point_count);

					while(dbs.m_pq_mst_node->size() > 0)
					{
						mst_node = dbs.m_pq_mst_node->top();
						dbs.m_pq_mst_node->pop();
								

						// now perform union_find operation
                        			root_u = mst_node.m_reachable_to;
                        			root_v = mst_node.m_reachable_from;
								
						while(vec_parents[root_u] != vec_parents[root_v])
                                		{
                                        		if(vec_parents[root_u] < vec_parents[root_v])
                                        		{
                                                		if(root_u == vec_parents[root_u])
                                                		{
                                                        		vec_parents[root_u] = vec_parents[root_v];
									if(round_pos == round_count - 1)
										dbs.m_mst_root.push_back(mst_node);
                                                        		else
										pq_mst_node->push(mst_node);
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
                                		                	vec_parents[root_v] = vec_parents[root_u];
									if(round_pos == round_count - 1)
										dbs.m_mst_root.push_back(mst_node);	
									else
										pq_mst_node->push(mst_node);
                                                        		break;
                                                		}

				                                z = vec_parents[root_v];
                                				vec_parents[root_v] = vec_parents[root_u];
                                                		root_v = z;                  
                                        		}
                                		}
					}

					delete dbs.m_pq_mst_node;

					dbs.m_pq_mst_node = pq_mst_node;
					pq_mst_node = NULL;
					//#ifdef _DEBUG
					comp_time2 = MPI_Wtime() - comp_time;
                                        comp_time_total += comp_time2;
                                        //#endif
					//cout << "rank==" << rank << " round_pos " << round_pos << " MST tree size " << dbs.m_pq_mst_node->size() << " vector " << vec_mst_node_recv.size() <<endl;
				}
				else
				{
					vec_parents.clear(); // this will reduce the required memory
					// rank is sender
					rank_to = rank - (int)(pow (2, round_pos));
					int cur_vec_gid_parent_gid_pair_size = dbs.m_vec_gid_parent_gid_pair.size();

					send_size[1] = cur_vec_gid_parent_gid_pair_size; // mst_node count
					send_size[2] = dbs.m_pq_mst_node->size(); // mst_node count 
					send_size[3] = dbs.m_all_gid_gathered_upper;
					send_size[4] = dbs.m_all_gid_gathered_lower;

					//cout << "rank " << rank << " round_pos " << round_pos << " tag " << tag << " rank_from " << rank_to << endl;	
					MPI_Send(&send_size[0], 5, MPI_INT, rank_to, tag, MPI_COMM_WORLD);

					send_count = 0;

					if(cur_vec_gid_parent_gid_pair_size > 0)
					{
						MPI_Isend(&dbs.m_vec_gid_parent_gid_pair[0], cur_vec_gid_parent_gid_pair_size, MPI_INT, rank_to, tag + 1, MPI_COMM_WORLD, &req_send[send_count++]);	
					}
					
					// send the MSTs in PQ
					vector <MST_Node> vec_mst_node_send;
					vec_mst_node_send.reserve(dbs.m_pq_mst_node->size());
					while(dbs.m_pq_mst_node->size() > 0)
					{
						vec_mst_node_send.push_back(dbs.m_pq_mst_node->top());
						dbs.m_pq_mst_node->pop();
					}
						
					if(send_size[2] > 0)
					{
						MPI_Isend(&vec_mst_node_send[0], send_size[2], mst_node_type, rank_to, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);		
					}
					
					if(send_count > 0)
                        			MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
					
					vec_mst_node_send.clear();
					dbs.m_vec_gid_parent_gid_pair.clear();
				}
			}
	
			tag += 3;
			//MPI_Barrier(MPI_COMM_WORLD);
			#ifdef _DEBUG
                	if(rank == proc_of_interest) cout << "Round " << round_pos << " total_time " << MPI_Wtime() - inter2 << " computation time " << comp_time2 << endl;
                	inter2 = MPI_Wtime();
                	#endif

			#ifdef _CUT_OFF_ACTIVE
				int stop_value = 0, stop_sum = 0;
				#ifdef _USE_ROUND_TIME
				if(comp_time2 > _CUT_OFF_TIME)
					stop_value = 1;
				#endif
			
				#ifdef _USE_ROUND_COUNT
				if(round_pos == _CUT_OFF_ROUND)
					stop_value = 1;	
				#endif
				MPI_Allreduce(&stop_value, &stop_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				if(stop_sum > 0)
					break;
			#endif
		}
					
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "boruvka Merging done: Taken time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif	

		vector <MST_Node> vec_mst_node_send_rest;
		vector <vector <MST_Node> > vec_vec_mst_node_recv;	
		//int rest_mst_edges_per_proc[nproc], rest_mst_edges_this = dbs.m_pq_mst_node->size();
		int pq_and_gid_parent_gid_pair_size[2], pq_and_gid_parent_gid_pair_size_per_proc[2 * nproc];
		int pos = 0, start_pos[nproc];

		//cout << "Rank " << rank << " entering cutoff " << endl;
		// smart code should go here
		#ifdef _CUT_OFF_ACTIVE
		if(round_pos < round_count)
		{
			pq_and_gid_parent_gid_pair_size[0] = dbs.m_pq_mst_node->size(); 
			pq_and_gid_parent_gid_pair_size[1] = dbs.m_vec_gid_parent_gid_pair.size();
			// send everyone to processor proc_of_interest then
                	MPI_Gather(&pq_and_gid_parent_gid_pair_size[0], 2, MPI_INT, &pq_and_gid_parent_gid_pair_size_per_proc[0], 2, MPI_INT, proc_of_interest, MPI_COMM_WORLD);

			//cout << "Rank " << rank << " starting communication" << endl;

			send_count = 0;
			recv_count = 0;
			if(rank != proc_of_interest)
			{
				vec_parents.clear();

				if(pq_and_gid_parent_gid_pair_size[0] > 0)
				{
					vec_mst_node_send_rest.clear();
                                	vec_mst_node_send_rest.reserve(pq_and_gid_parent_gid_pair_size[0]);
                                	while(dbs.m_pq_mst_node->size() > 0)
                                	{
                                        	vec_mst_node_send_rest.push_back(dbs.m_pq_mst_node->top());
                                        	dbs.m_pq_mst_node->pop();
                                	}

                                	MPI_Isend(&vec_mst_node_send_rest[0], pq_and_gid_parent_gid_pair_size[0], mst_node_type, proc_of_interest, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);
				}

				if(pq_and_gid_parent_gid_pair_size[1] > 0)
					MPI_Isend(&dbs.m_vec_gid_parent_gid_pair[0], pq_and_gid_parent_gid_pair_size[1], MPI_INT, proc_of_interest, tag + 4, MPI_COMM_WORLD, &req_send[send_count++]);

				if(send_count > 0)
                                                MPI_Waitall(send_count, &req_send[0], &stat_send[0]);

				//vec_mst_node_send.clear();
				vec_mst_node_send_rest.clear();
                                dbs.m_vec_gid_parent_gid_pair.clear();
			}
			else
			{
				vec_vec_mst_node_recv.clear();
				vec_vec_mst_node_recv.resize(nproc);

				if(rank != proc_of_interest || rank != 0)
                                        cout << "The following shouldNOT work" << endl;

				for(i = 0; i < nproc; i++)
                                {
                                        if(i == rank)
                                                continue;

					if(pq_and_gid_parent_gid_pair_size_per_proc[2 * i] > 0)
					{
						vec_vec_mst_node_recv[i].resize(pq_and_gid_parent_gid_pair_size_per_proc[2 * i]);
						MPI_Irecv(&vec_vec_mst_node_recv[i][0], pq_and_gid_parent_gid_pair_size_per_proc[2 * i], mst_node_type, i, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]);
					}	
				}	
				
				pos = dbs.m_vec_gid_parent_gid_pair.size();

				for(i = 0; i < nproc; i++)
				{
					if(i == rank)
                                                continue;

					start_pos[i] = pos;
					pos += pq_and_gid_parent_gid_pair_size_per_proc[2 * i + 1];
				}

				dbs.m_vec_gid_parent_gid_pair.resize(pos);
				
				for(i = 0; i < nproc; i++)
                                {
                                        if(i == rank)
                                                continue;

					if(pq_and_gid_parent_gid_pair_size_per_proc[2 * i + 1] > 0)
						MPI_Irecv(&dbs.m_vec_gid_parent_gid_pair[start_pos[i]], pq_and_gid_parent_gid_pair_size_per_proc[2 * i + 1], MPI_INT, i, tag + 4, MPI_COMM_WORLD, &req_recv[recv_count++]);	
				}

				// reset that range from parents vector
				for(i = 0; i < global_owned_local_point_count; i++)
					vec_parents[i] = i;

				for(i = 0; i < recv_count; i++)
				{
					MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv_single);
		
					rtag = stat_recv_single.MPI_TAG;
					rsource = stat_recv_single.MPI_SOURCE;
		
					if(rtag == tag + 3)
					{
						for(j = 0; j < vec_vec_mst_node_recv[rsource].size(); j++)
                                                	dbs.m_pq_mst_node->push(vec_vec_mst_node_recv[rsource][j]);

						vec_vec_mst_node_recv[rsource].clear();
					}
					else if(rtag == tag + 4)
					{
						// do nothing
					}
				}

				for(i = 0; i < dbs.m_vec_gid_parent_gid_pair.size(); i += 2)
					vec_parents[dbs.m_vec_gid_parent_gid_pair[i]] = dbs.m_vec_gid_parent_gid_pair[i + 1];
				
			}

			vec_vec_mst_node_recv.clear();
			vec_mst_node_send_rest.clear();

			//cout << "Rank " << rank << " communication done" << endl;
			if(rank == proc_of_interest)
			{
				//#ifdef _DEBUG
				comp_time = MPI_Wtime();
				//#endif
			
				dbs.m_mst_root.reserve(global_owned_local_point_count);

				while(dbs.m_pq_mst_node->size() > 0)
				{
					mst_node = dbs.m_pq_mst_node->top();
					dbs.m_pq_mst_node->pop();
								

					// now perform union_find operation
	                        	root_u = mst_node.m_reachable_to;
        	                	root_v = mst_node.m_reachable_from;
									
					while(vec_parents[root_u] != vec_parents[root_v])
                                	{
                                		if(vec_parents[root_u] < vec_parents[root_v])
	                                        {
        	                                	if(root_u == vec_parents[root_u])
                	                               	{
                        	                        	vec_parents[root_u] = vec_parents[root_v];
								dbs.m_mst_root.push_back(mst_node);
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
                                				vec_parents[root_v] = vec_parents[root_u];
								dbs.m_mst_root.push_back(mst_node);	
                                                        	break;
                                                	}

				                	z = vec_parents[root_v];
                                			vec_parents[root_v] = vec_parents[root_u];
                                                	root_v = z;                  
                                        	}
                                	}
				}

				//#ifdef _DEBUG
				comp_time2 = MPI_Wtime() - comp_time;
                                comp_time_total += comp_time2;
                               	//#endif
			}			
		
			#ifdef _DEBUG
                	if(rank == proc_of_interest) cout << "Cut_off_acitve: Round_count " << round_count << " total time " << MPI_Wtime() - inter << " computation " << comp_time2 << endl;
                	inter = MPI_Wtime();
                	#endif
		}
		#endif

		#ifndef _DEBUG
		if(rank == proc_of_interest) cout << " mer_t " << MPI_Wtime() - inter << " to_comp_t " << comp_time_total;
                inter = MPI_Wtime();
		#endif

		vec_parents.clear();
		gid_parent_gid.clear();			

		// now the filtered edges here edges
		//int pos = 0, start_pos[nproc], 
		int current_mst_root_size = dbs.m_mst_root.size();
		send_count = 0;
		recv_count = 0;
		if(proc_of_interest == rank)
		{
			if(proc_of_interest != 0)
				cout << "The following part need to be changed !!!" << endl;

			pos = current_mst_root_size;
			for(i = 1; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += filtered_edges_per_thread[i];
			}
			// reserve memory now at proc_of_interest

			dbs.m_mst_root.resize(pos);

			for(i = 1; i < nproc; i++)
			{
				if(filtered_edges_per_thread[i] > 0)
				{
					MPI_Irecv(&dbs.m_mst_root[start_pos[i]], filtered_edges_per_thread[i], mst_node_type, i, tag + 6, MPI_COMM_WORLD, &req_recv[recv_count++]);

				}
			}
			if(recv_count > 0)
                                MPI_Waitall(recv_count, &req_recv[0], &stat_recv[0]);
		}
		else
		{
			if(current_mst_root_size > 0)
                        {
                        	MPI_Isend(&dbs.m_mst_root[0], current_mst_root_size, mst_node_type, proc_of_interest, tag + 6, MPI_COMM_WORLD, &req_send[send_count++]);
                        }

                        if(send_count > 0)
                        	MPI_Waitall(send_count, &req_send[0], &stat_send[0]);

			dbs.m_mst_root.clear();
		}

		//if(rank == proc_of_interest) cout << "Gathering filtered edges DONE" << endl;

		//#ifdef _GATHER_CORE_DIS_AT_A_TIME
		int owned_local_points_count_per_proc[nproc];
                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		
		//int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_gathered.resize(global_owned_local_point_count, -1);
			dbs.m_vec_global_point_ID_gathered.resize(global_owned_local_point_count, -1);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &dbs.m_vec_f_core_distance_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &dbs.m_vec_global_point_ID_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		//#endif
		//if(rank == proc_of_interest) cout << "Gathering core distnaces DONE" << endl;

		if(rank == proc_of_interest)
		{
			//dbs.m_vec_f_core_distance_gathered_assigned.clear();
			//dbs.m_vec_f_core_distance_gathered_assigned.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root.clear();
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);

			for(i = 0; i < global_owned_local_point_count; i++)
				dbs.m_vec_f_core_distance_root[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
				//dbs.m_vec_f_core_distance_gathered_assigned[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
		}

		dbs.m_vec_global_point_ID_gathered.clear();
		dbs.m_vec_f_core_distance_gathered.clear();
		MPI_Type_free(&mst_node_type);
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "boruvka pairwise_local_mst_merging SMART using Kruskal - gathering filtered edges, assiging gid and core_distance time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
		#else
		if(rank == proc_of_interest) cout << " gat_fed_gid_cd_t " << MPI_Wtime() - inter;
                inter = MPI_Wtime();
                #endif

		/*	
		if(round_pos < round_count)
		{
			// send everyone to processor proc_of_interest then
                	MPI_Gather(&rest_mst_edges_this, 1, MPI_INT, &rest_mst_edges_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);

			send_count = 0;
			if(rest_mst_edges_this > 0 && rank != proc_of_interest)
			{
				// send the MSTs in PQ
				vec_mst_node_send_rest.clear();
				vec_mst_node_send_rest.reserve(rest_mst_edges_this);
				while(dbs.m_pq_mst_node->size() > 0)
				{
					vec_mst_node_send_rest.push_back(dbs.m_pq_mst_node->top());
					dbs.m_pq_mst_node->pop();
				}
						
				MPI_Isend(&vec_mst_node_send_rest[0], rest_mst_edges_this, mst_node_type, proc_of_interest, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);			
			}

			recv_count = 0;
			if(rank == proc_of_interest)
			{
				vec_vec_mst_node_recv.resize(nproc);
				for(i = 0; i < nproc; i++)
				{
					if(i == rank)
						continue;
					
					if(rest_mst_edges_per_proc[i] > 0)
					{
                                        	vec_vec_mst_node_recv[i].resize(rest_mst_edges_per_proc[i]);
						MPI_Irecv(&vec_vec_mst_node_recv[i][0], rest_mst_edges_per_proc[i], mst_node_type, i, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]);
					}					
				}
				
				#ifdef _DEBUG
				comp_time = MPI_Wtime();
				#endif
					
				// reset that range from parents vector
				for(i = 0; i < global_owned_local_point_count; i++)
					vec_parents[i] = i;

				#ifdef _DEBUG
				comp_time2 = MPI_Wtime() - comp_time;
                                comp_time_total += comp_time2;
                                #endif

				for(i = 0; i < recv_count; i++)
				{
					MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
					rtag = stat_recv.MPI_TAG;
					rsource = stat_recv.MPI_SOURCE;
		
					if(rtag == tag + 3)
					{
						for(j = 0; j < vec_vec_mst_node_recv[rsource].size(); j++)
                                                	dbs.m_pq_mst_node->push(vec_vec_mst_node_recv[rsource][j]);

						vec_vec_mst_node_recv[rsource].clear();
					}
				}
			}
			else
			{
				if(send_count > 0)
                        		MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
					
				vec_mst_node_send_rest.clear();	
			}

			vec_vec_mst_node_recv.clear();
			vec_mst_node_send_rest.clear();

			if(rank == proc_of_interest)
			{
				#ifdef _DEBUG
				comp_time = MPI_Wtime();
				#endif
			
				dbs.m_mst_root.reserve(global_owned_local_point_count);

				while(dbs.m_pq_mst_node->size() > 0)
				{
					mst_node = dbs.m_pq_mst_node->top();
					dbs.m_pq_mst_node->pop();
								

					// now perform union_find operation
	                        	root_u = mst_node.m_reachable_to;
        	                	root_v = mst_node.m_reachable_from;
									
					while(vec_parents[root_u] != vec_parents[root_v])
                                	{
                                		if(vec_parents[root_u] < vec_parents[root_v])
	                                        {
        	                                	if(root_u == vec_parents[root_u])
                	                               	{
                        	                        	vec_parents[root_u] = vec_parents[root_v];
								dbs.m_mst_root.push_back(mst_node);
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
                                				vec_parents[root_v] = vec_parents[root_u];
								dbs.m_mst_root.push_back(mst_node);	
                                                        	break;
                                                	}

				                	z = vec_parents[root_v];
                                			vec_parents[root_v] = vec_parents[root_u];
                                                	root_v = z;                  
                                        	}
                                	}
				}

				#ifdef _DEBUG
				comp_time2 = MPI_Wtime() - comp_time;
                                comp_time_total += comp_time2;
                               	#endif
			}			
		}

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Round " << round_count << " total_time " << MPI_Wtime() - inter2 << " computation time " << comp_time2 << endl;
                inter2 = MPI_Wtime();
                #endif
	
		vec_parents.clear();
		MPI_Type_free(&mst_node_type);

		#ifdef _DEBUG
		if(rank == proc_of_interest) cout << "Computation time while pairwise_local_mst_merging SMART " << comp_time_total << endl;
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging SMART using Kruskal - gathered and extracted global MST time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif

		#ifdef _GATHER_CORE_DIS_AT_A_TIME
		int owned_local_points_count_per_proc[nproc];
                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		
		int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_gathered.resize(global_owned_local_point_count, -1);
			dbs.m_vec_global_point_ID_gathered.resize(global_owned_local_point_count, -1);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &dbs.m_vec_f_core_distance_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &dbs.m_vec_global_point_ID_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		#endif

		if(rank == proc_of_interest)
		{
			//dbs.m_vec_f_core_distance_gathered_assigned.clear();
			//dbs.m_vec_f_core_distance_gathered_assigned.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root.clear();
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);

			for(i = 0; i < global_owned_local_point_count; i++)
				dbs.m_vec_f_core_distance_root[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
				//dbs.m_vec_f_core_distance_gathered_assigned[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
		}

		dbs.m_vec_global_point_ID_gathered.clear();
		dbs.m_vec_f_core_distance_gathered.clear();

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging SMART using Kruskal - assiging gid and core_distance time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif
	
		if(nproc == 1)
		{
			dbs.m_vec_f_core_distance_root.clear();
                        dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root = dbs.m_vec_f_core_distance;
			
			dbs.m_mst_root.clear();
			dbs.m_mst_root.reserve(dbs.m_pq_mst_node->size());
			
			while(dbs.m_pq_mst_node->size() > 0)
			{
				mst_node = dbs.m_pq_mst_node->top();
				dbs.m_pq_mst_node->pop();
								

				// now perform union_find operation
                        	root_u = mst_node.m_reachable_to;
                        	root_v = mst_node.m_reachable_from;
								
				while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                        if(vec_parents[root_u] < vec_parents[root_v])
                                	{
                                                if(root_u == vec_parents[root_u])
                                                {
                                                        vec_parents[root_u] = vec_parents[root_v];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);
                                                        //else
							//pq_mst_node->push(mst_node);
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
                                			vec_parents[root_v] = vec_parents[root_u];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);	
							//else
							//pq_mst_node->push(mst_node);
                                                        break;
                                                }

				                z = vec_parents[root_v];
                               			vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;                  
                                        }
                                }
			}
		}
	

		vec_gid_send.clear();
		vec_parents_gid_send.clear();
		*/

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in global MST " << dbs.m_mst_root.size() << endl;
		#else
		if(rank == proc_of_interest) cout << " ed_g_mst " << dbs.m_mst_root.size();
                #endif
	
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Time taken by pairwise_local_mst_merging using Kruskal  " << MPI_Wtime() - start << endl;
                #else
                if(rank == proc_of_interest) cout << " tot_mer_t " << MPI_Wtime() - start;
                #endif
	}

	void sequential_local_mst_merging_kruskal(ClusteringAlgo& dbs)
	{
		double start = MPI_Wtime();
		double inter = start;
                int rank, nproc, i, j, round_pos, round_count = 0, proc_count;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		int local_mst_edge_count_this = dbs.m_pq_mst_node->size();
                int global_mst_edge_count;
		int total_connecting_edges_count = 0;

                MPI_Allreduce(&local_mst_edge_count_this, &global_mst_edge_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&dbs.m_connecting_edges_count, &total_connecting_edges_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                #ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in local MSTs " << global_mst_edge_count << endl;
		if(rank == proc_of_interest) cout << "Total connecting edges " << total_connecting_edges_count << endl;
                #endif

		int global_owned_local_point_count;
		MPI_Allreduce(&dbs.m_pts->m_i_num_local_points, &global_owned_local_point_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		dbs.m_global_owned_local_point_count_root = global_owned_local_point_count;
		
		dbs.m_all_gid_gathered_upper = dbs.m_pts->m_all_local_point_gID_upper;
		dbs.m_all_gid_gathered_lower = dbs.m_pts->m_all_local_point_gID_lower;

		
		int local_mst_edge_count_per_proc[nproc];
		MPI_Gather(&local_mst_edge_count_this, 1, MPI_INT, &local_mst_edge_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		// copy the mst to root, but only for the first time;
                proc_count = nproc;

                while((proc_count % 2 == 0) && (proc_count > 1))
                {
                        round_count++;
                        proc_count = proc_count / 2;
                }

                #ifdef _DEBUG
                if(proc_count != 1)
                        cout << "Number of procs is not a multiple of 2" << endl;
                #endif


		MPI_Datatype mst_node_type;
                MPI_Type_contiguous(sizeof(MST_Node), MPI_BYTE, &mst_node_type);
                MPI_Type_commit(&mst_node_type);

		int rank_from, rank_to, recv_size[5], send_size[5]; //mst_node_recv_size, mst_node_send_size;
		int tag = 200, root_u, root_v, z;
                MPI_Request req_send[nproc], req_recv[nproc];
                MPI_Status stat_send[nproc], stat_recv;
		int send_count, recv_count;
		MST_Node mst_node;
		int rtag, rsource, rpos;
		double inter2, comp_time, comp_time_total = 0, comp_time2;

		// send the MSTs in PQ
		vector <MST_Node> vec_mst_node_send;
		vec_mst_node_send.reserve(dbs.m_pq_mst_node->size());
		while(dbs.m_pq_mst_node->size() > 0)
		{
			vec_mst_node_send.push_back(dbs.m_pq_mst_node->top());
			dbs.m_pq_mst_node->pop();
		}
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "sequential_local_mst_merging using Kruskal - preprocessing time " << MPI_Wtime() - inter << endl;
		inter = MPI_Wtime();
		inter2 = inter;
		#endif

		send_count = 0;
						
		if(local_mst_edge_count_this > 0)
		{
			MPI_Isend(&vec_mst_node_send[0], local_mst_edge_count_this, mst_node_type, proc_of_interest, tag, MPI_COMM_WORLD, &req_send[send_count++]);
		}

		vector <vector <MST_Node> > vec_vec_mst_node_recv;

		if(rank == proc_of_interest)
		{

			vec_vec_mst_node_recv.resize(nproc);
			
			recv_count = 0;
			cout << "starting receiving----" << endl;
			for(i = 0; i < nproc; i++)
			{

				if(local_mst_edge_count_per_proc[i] > 0)
				{
					cout << "allocating receiving---- from " << i << " size " << local_mst_edge_count_per_proc[i] << endl;
					vec_vec_mst_node_recv[i].resize(local_mst_edge_count_per_proc[i]);
					MPI_Irecv(&vec_vec_mst_node_recv[i][0], local_mst_edge_count_per_proc[i], mst_node_type, i, tag, MPI_COMM_WORLD, &req_recv[recv_count++]);
					cout << "starting receiving---- from " << i << " size " << local_mst_edge_count_per_proc[i] << endl;
				}
			}
			cout << "allocation and Irecv done----" << endl;	
		}
				
		if(send_count > 0)
                	MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
			
		
		vec_mst_node_send.clear(); 
	
		if(rank == proc_of_interest)
		{
			if(dbs.m_pq_mst_node->size() > 0)
				cout << "SOMETHING is wrong----" << endl;
			for(i = 0; i < recv_count; i++)
			{
				MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
				rtag = stat_recv.MPI_TAG;
                               	rsource = stat_recv.MPI_SOURCE;
				
				cout << "received from " << rsource << endl;

				if(rtag == tag)
				{					
					for(j = 0; j < vec_vec_mst_node_recv[rsource].size(); j++)
						dbs.m_pq_mst_node->push(vec_vec_mst_node_recv[rsource][j]);

					vec_vec_mst_node_recv[rsource].clear();

					cout << "inserted data of " << rsource << endl;
				}
			}	
	
		}

		vec_vec_mst_node_recv.clear();
	
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "sequential_local_mst_merging using Kruskal - gathered and inserted into the priority queue " << MPI_Wtime() - inter << endl;
		inter = MPI_Wtime();
		inter2 = inter;
		#endif

		dbs.m_mst_root.clear();
		vector <int> vec_parents;
		
		if(rank == proc_of_interest)
		{
                	vec_parents.reserve(global_owned_local_point_count);
                	for(i = 0; i < global_owned_local_point_count; i++)
                        	vec_parents.push_back(i);
		}
	
		if(rank == proc_of_interest)
		{
			dbs.m_mst_root.reserve(global_owned_local_point_count);
			while(dbs.m_pq_mst_node->size() > 0)
			{
				mst_node = dbs.m_pq_mst_node->top();
				dbs.m_pq_mst_node->pop();
								
				// now perform union_find operation
                        	root_u = mst_node.m_reachable_to;
                        	root_v = mst_node.m_reachable_from;
								
				while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                	if(vec_parents[root_u] < vec_parents[root_v])
                                        {
                                        	if(root_u == vec_parents[root_u])
                                               	{
                                                	vec_parents[root_u] = vec_parents[root_v];
							dbs.m_mst_root.push_back(mst_node);
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
                                			vec_parents[root_v] = vec_parents[root_u];
							dbs.m_mst_root.push_back(mst_node);	
                                                        break;
                                                }

				                z = vec_parents[root_v];
                                		vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;                  
                                        }
                                }
			}
		}
	
		vec_parents.clear();
		MPI_Type_free(&mst_node_type);

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "sequential_local_mst_merging using Kruskal - computing global MST from local MSTstime " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif

		#ifdef _GATHER_CORE_DIS_AT_A_TIME
		int owned_local_points_count_per_proc[nproc];
                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		
		int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_gathered.resize(global_owned_local_point_count, -1);
			dbs.m_vec_global_point_ID_gathered.resize(global_owned_local_point_count, -1);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &dbs.m_vec_f_core_distance_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &dbs.m_vec_global_point_ID_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		#endif

		if(rank == proc_of_interest)
		{
			//dbs.m_vec_f_core_distance_gathered_assigned.clear();
			//dbs.m_vec_f_core_distance_gathered_assigned.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root.clear();
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);

			for(i = 0; i < global_owned_local_point_count; i++)
				dbs.m_vec_f_core_distance_root[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
				//dbs.m_vec_f_core_distance_gathered_assigned[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
		}

		dbs.m_vec_global_point_ID_gathered.clear();
		dbs.m_vec_f_core_distance_gathered.clear();

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "sequential_local_mst_merging using Kruskal - assiging gid and core_distance time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif

		/*	
		if(nproc == 1)
		{
			dbs.m_vec_f_core_distance_root.clear();
                        dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root = dbs.m_vec_f_core_distance;
			
			dbs.m_mst_root.clear();
			dbs.m_mst_root.reserve(dbs.m_pq_mst_node->size());
		
			for(i = 0; i < global_owned_local_point_count; i++)
				vec_parents[i] = i;
	
			while(dbs.m_pq_mst_node->size() > 0)
			{
				mst_node = dbs.m_pq_mst_node->top();
				dbs.m_pq_mst_node->pop();
								

				// now perform union_find operation
                        	root_u = mst_node.m_reachable_to;
                        	root_v = mst_node.m_reachable_from;
								
				while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                        if(vec_parents[root_u] < vec_parents[root_v])
                                	{
                                                if(root_u == vec_parents[root_u])
                                                {
                                                        vec_parents[root_u] = vec_parents[root_v];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);
                                                        //else
							//pq_mst_node->push(mst_node);
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
                                			vec_parents[root_v] = vec_parents[root_u];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);	
							//else
							//pq_mst_node->push(mst_node);
                                                        break;
                                                }

				                z = vec_parents[root_v];
                               			vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;                  
                                        }
                                }
			}
		}
		*/

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in global MST " << dbs.m_mst_root.size() << endl;
                #endif
	
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Time taken by sequential_local_mst_merging using Kruskal  " << MPI_Wtime() - start << endl;
                #else
                if(rank == proc_of_interest) cout << " sequential_local_mst_merging_kruskal " << MPI_Wtime() - start;
                #endif
	}

	void pairwise_local_mst_merging_kruskal(ClusteringAlgo& dbs)
	{
		double start = MPI_Wtime();
		double inter = start;
                int rank, nproc, i, j, round_pos, round_count = 0, proc_count;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);


		int local_mst_edge_count_this = dbs.m_pq_mst_node->size();
                int global_mst_edge_count;
		int total_connecting_edges_count = 0;

                MPI_Allreduce(&local_mst_edge_count_this, &global_mst_edge_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&dbs.m_connecting_edges_count, &total_connecting_edges_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                #ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in local MSTs " << global_mst_edge_count << endl;
		if(rank == proc_of_interest) cout << "Total connecting edges " << total_connecting_edges_count << endl;
                #endif

		int global_owned_local_point_count;
		MPI_Allreduce(&dbs.m_pts->m_i_num_local_points, &global_owned_local_point_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		dbs.m_global_owned_local_point_count_root = global_owned_local_point_count;
		
		dbs.m_all_gid_gathered_upper = dbs.m_pts->m_all_local_point_gID_upper;
		dbs.m_all_gid_gathered_lower = dbs.m_pts->m_all_local_point_gID_lower;

		// copy the mst to root, but only for the first time;
                proc_count = nproc;

                while((proc_count % 2 == 0) && (proc_count > 1))
                {
                        round_count++;
                        proc_count = proc_count / 2;
                }

                #ifdef _DEBUG
                if(proc_count != 1)
                        cout << "Number of procs is not a multiple of 2" << endl;
                #endif

		#ifndef _GATHER_CORE_DIS_AT_A_TIME
		dbs.m_vec_global_point_ID_gathered.clear();
                dbs.m_vec_f_core_distance_gathered.clear();

		dbs.m_vec_f_core_distance_gathered.reserve(dbs.m_vec_f_core_distance.size());
		dbs.m_vec_global_point_ID_gathered.reserve(dbs.m_vec_f_core_distance.size());

		for(i = 0; i < dbs.m_vec_f_core_distance.size(); i++)
		{
			dbs.m_vec_f_core_distance_gathered.push_back(dbs.m_vec_f_core_distance[i]);
			dbs.m_vec_global_point_ID_gathered.push_back(dbs.m_pts->m_vec_global_point_ID[i]);
		}
		#endif

		vector <int> vec_parents;
                vec_parents.reserve(global_owned_local_point_count);
                for(i = 0; i < global_owned_local_point_count; i++)
                        vec_parents.push_back(i);
		
		
		MPI_Datatype mst_node_type;
                MPI_Type_contiguous(sizeof(MST_Node), MPI_BYTE, &mst_node_type);
                MPI_Type_commit(&mst_node_type);

		int rank_from, rank_to, recv_size[5], send_size[5]; //mst_node_recv_size, mst_node_send_size;
		int tag = 200, root_u, root_v, z;
                MPI_Request req_send[4 * nproc], req_recv[4 * nproc];
                MPI_Status stat_send[4 * nproc], stat_recv;
		int send_count, recv_count;
		MST_Node mst_node;
		int rtag, rsource, rpos;
		double inter2, comp_time, comp_time_total = 0, comp_time2;

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging using Kruskal - preprocessing time " << MPI_Wtime() - inter << endl;
		inter = MPI_Wtime();
		inter2 = inter;
		#endif

		//cout << "rank " << rank  << " MST size " << dbs.m_pq_mst_node->size() << endl;
		for(round_pos = 0; round_pos < round_count; round_pos++)
		{
			//tag  = 200; // a random start
			if(rank % (int)(pow (2, round_pos)) == 0)
			{
				// rank is either a sender or a receiver
				if(rank % (int)(pow (2, round_pos + 1)) == 0)
				{
					// rank is receiver
					rank_from = rank + (int)(pow (2, round_pos));					

					#ifndef _GATHER_CORE_DIS_AT_A_TIME
					recv_size[0] = 0; // mst_node count 
					recv_size[1] = 0; // core_distance count
					#endif
					recv_size[2] = 0; // corresponing global ID count
					recv_size[3] = 0; // global ID upper limit
					recv_size[4] = 0; // global ID lower limit
		
					MPI_Recv(&recv_size[0], 5, MPI_INT, rank_from, tag, MPI_COMM_WORLD, &stat_recv);
					//cout << "rank " << rank << " round_pos " << round_pos << " tag " << tag << " rank_to " << rank_from << endl;

					recv_count = 0;
					
					#ifndef _GATHER_CORE_DIS_AT_A_TIME	
					int cur_size = dbs.m_vec_f_core_distance_gathered.size();
					if(recv_size[0] > 0)
					{
						dbs.m_vec_f_core_distance_gathered.resize(cur_size + recv_size[0]);
						MPI_Irecv(&dbs.m_vec_f_core_distance_gathered[cur_size], recv_size[0], MPI_FLOAT, rank_from, tag + 1, MPI_COMM_WORLD, &req_recv[recv_count++]);
					}
										
					if(recv_size[1] > 0)
                                        {
						dbs.m_vec_global_point_ID_gathered.resize(cur_size + recv_size[1]);
                                                MPI_Irecv(&dbs.m_vec_global_point_ID_gathered[cur_size], recv_size[1], MPI_INT, rank_from, tag + 2, MPI_COMM_WORLD, &req_recv[recv_count++]);
                                        }
					#endif

					vector <MST_Node> vec_mst_node_recv;
					if(recv_size[2] > 0)
					{
                                        	vec_mst_node_recv.resize(recv_size[2]);
						MPI_Irecv(&vec_mst_node_recv[0], recv_size[2], mst_node_type, rank_from, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]); 
					}
						
					if(dbs.m_all_gid_gathered_upper < recv_size[3])
						dbs.m_all_gid_gathered_upper = recv_size[3];

					if(dbs.m_all_gid_gathered_lower > recv_size[4])
						dbs.m_all_gid_gathered_lower = recv_size[4];

					// reset that range from parents vector
					for(i = dbs.m_all_gid_gathered_lower; i < dbs.m_all_gid_gathered_upper; i++)
						vec_parents[i] = i;

					for(i = 0; i < recv_count; i++)
					{
						MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
						rtag = stat_recv.MPI_TAG;
						rsource = stat_recv.MPI_SOURCE;
		
						if(rtag == tag + 1)
						{
							// nothing to do with the core ditances for the time being		
						}
						else if(rtag == tag + 2)
						{
							// nothing to do for the time being
						}
						else if(rtag == tag + 3)
						{
							comp_time = MPI_Wtime();
							for(j = 0; j < vec_mst_node_recv.size(); j++)
                                                		dbs.m_pq_mst_node->push(vec_mst_node_recv[j]);

							// run Kruskals algorithm on the edges in the priority queue
							priority_queue<MST_Node>* pq_mst_node;
							pq_mst_node = new priority_queue<MST_Node>;

							if(round_pos == round_count - 1)
								dbs.m_mst_root.reserve(global_owned_local_point_count);

							while(dbs.m_pq_mst_node->size() > 0)
							{
								mst_node = dbs.m_pq_mst_node->top();
								dbs.m_pq_mst_node->pop();
								

								// now perform union_find operation
                        					root_u = mst_node.m_reachable_to;
                        					root_v = mst_node.m_reachable_from;
								
								while(vec_parents[root_u] != vec_parents[root_v])
                                				{
                                        				if(vec_parents[root_u] < vec_parents[root_v])
                                        				{
                                                				if(root_u == vec_parents[root_u])
                                                				{
                                                        				vec_parents[root_u] = vec_parents[root_v];
											if(round_pos == round_count - 1)
												dbs.m_mst_root.push_back(mst_node);
                                                        				else
												pq_mst_node->push(mst_node);
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
                                				                        vec_parents[root_v] = vec_parents[root_u];
											if(round_pos == round_count - 1)
												dbs.m_mst_root.push_back(mst_node);	
											else
												pq_mst_node->push(mst_node);
                                                        				break;
                                                				}

				                                                z = vec_parents[root_v];
                                				                vec_parents[root_v] = vec_parents[root_u];
                                                				root_v = z;                  
                                        				}
                                				}
							}

							delete dbs.m_pq_mst_node;

							dbs.m_pq_mst_node = pq_mst_node;
							pq_mst_node = NULL;
							comp_time2 = MPI_Wtime() - comp_time;
                                                        comp_time_total += comp_time2;
						}
					}
				
					vec_mst_node_recv.clear();
				}
				else
				{
					vec_parents.clear(); // this will reduce the required memory
					// rank is sender
					rank_to = rank - (int)(pow (2, round_pos));
					
					#ifndef _GATHER_CORE_DIS_AT_A_TIME
                                        send_size[0] = dbs.m_vec_f_core_distance_gathered.size(); // core_distance count
                                        send_size[1] = dbs.m_vec_global_point_ID_gathered.size(); // corresponing global ID count
					#endif

					send_size[2] = dbs.m_pq_mst_node->size(); // mst_node count 
					send_size[3] = dbs.m_all_gid_gathered_upper;
					send_size[4] = dbs.m_all_gid_gathered_lower;

					MPI_Send(&send_size[0], 5, MPI_INT, rank_to, tag, MPI_COMM_WORLD);

					send_count = 0;
					
					#ifndef _GATHER_CORE_DIS_AT_A_TIME
					if(send_size[0] > 0)
                                        {
                                                MPI_Isend(&dbs.m_vec_f_core_distance_gathered[0], send_size[0], MPI_FLOAT, rank_to, tag + 1, MPI_COMM_WORLD, &req_send[send_count++]);
                                        }

					if(send_size[1] > 0)
                                        {
                                                MPI_Isend(&dbs.m_vec_global_point_ID_gathered[0], send_size[1], MPI_INT, rank_to, tag + 2, MPI_COMM_WORLD, &req_send[send_count++]);
                                        }
					#endif
	
					// send the MSTs in PQ
					vector <MST_Node> vec_mst_node_send;
					vec_mst_node_send.reserve(dbs.m_pq_mst_node->size());
					while(dbs.m_pq_mst_node->size() > 0)
					{
						vec_mst_node_send.push_back(dbs.m_pq_mst_node->top());
						dbs.m_pq_mst_node->pop();
					}
						
					if(send_size[2] > 0)
					{
						MPI_Isend(&vec_mst_node_send[0], send_size[2], mst_node_type, rank_to, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);		
					}
					
					if(send_count > 0)
                        			MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
					
					vec_mst_node_send.clear(); 
				}
			}

			tag += 3;
			//MPI_Barrier(MPI_COMM_WORLD);
			#ifdef _DEBUG
                	if(rank == proc_of_interest) cout << "Round " << round_pos << " total_time " << MPI_Wtime() - inter2 << " computation time " << comp_time2 << endl;
                	inter2 = MPI_Wtime();
                	#endif
		}
	
		vec_parents.clear();
		MPI_Type_free(&mst_node_type);

		#ifdef _DEBUG
		if(rank == proc_of_interest) cout << "Computation time while pairwise_local_mst_merging " << comp_time_total << endl;
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging using Kruskal - gathered and extracted global MST time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif

		#ifdef _GATHER_CORE_DIS_AT_A_TIME
		int owned_local_points_count_per_proc[nproc];
                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		
		int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_gathered.resize(global_owned_local_point_count, -1);
			dbs.m_vec_global_point_ID_gathered.resize(global_owned_local_point_count, -1);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &dbs.m_vec_f_core_distance_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &dbs.m_vec_global_point_ID_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		#endif

		if(rank == proc_of_interest)
		{
			//dbs.m_vec_f_core_distance_gathered_assigned.clear();
			//dbs.m_vec_f_core_distance_gathered_assigned.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root.clear();
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);

			for(i = 0; i < global_owned_local_point_count; i++)
				dbs.m_vec_f_core_distance_root[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
				//dbs.m_vec_f_core_distance_gathered_assigned[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
		}

		dbs.m_vec_global_point_ID_gathered.clear();
		dbs.m_vec_f_core_distance_gathered.clear();

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging using Kruskal - assiging gid and core_distance time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
                #endif
	
		if(nproc == 1)
		{
			dbs.m_vec_f_core_distance_root.clear();
                        dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root = dbs.m_vec_f_core_distance;
			
			dbs.m_mst_root.clear();
			dbs.m_mst_root.reserve(dbs.m_pq_mst_node->size());
			
			while(dbs.m_pq_mst_node->size() > 0)
			{
				mst_node = dbs.m_pq_mst_node->top();
				dbs.m_pq_mst_node->pop();
								

				// now perform union_find operation
                        	root_u = mst_node.m_reachable_to;
                        	root_v = mst_node.m_reachable_from;
								
				while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                        if(vec_parents[root_u] < vec_parents[root_v])
                                	{
                                                if(root_u == vec_parents[root_u])
                                                {
                                                        vec_parents[root_u] = vec_parents[root_v];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);
                                                        //else
							//pq_mst_node->push(mst_node);
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
                                			vec_parents[root_v] = vec_parents[root_u];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);	
							//else
							//pq_mst_node->push(mst_node);
                                                        break;
                                                }

				                z = vec_parents[root_v];
                               			vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;                  
                                        }
                                }
			}
		}

		if(rank == proc_of_interest) cout << "Time taken by pairwise merging stage using Kruskal in parallel took: " << MPI_Wtime() - start << " seconds." << endl;
		if(rank == proc_of_interest) cout << "Total edges in global MST " << dbs.m_mst_root.size() << endl;
	}

	void pairwise_local_mst_merging_kruskal_smart(ClusteringAlgo& dbs)
	{
		double start = MPI_Wtime();
		double inter = start;
                int rank, nproc, i, j, round_pos, round_count = 0, proc_count;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		int local_mst_edge_count_this = dbs.m_pq_mst_node->size();
                int global_mst_edge_count;
		int total_connecting_edges_count = 0;

                MPI_Allreduce(&local_mst_edge_count_this, &global_mst_edge_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&dbs.m_connecting_edges_count, &total_connecting_edges_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                #ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Total edges in local MSTs " << global_mst_edge_count << endl;
		if(rank == proc_of_interest) cout << "Total connecting edges " << total_connecting_edges_count << endl;
		#else
		if(rank == proc_of_interest) cout << " ed_lo_msts " << global_mst_edge_count;
                if(rank == proc_of_interest) cout << " ed_co " << total_connecting_edges_count;
                #endif

		int global_owned_local_point_count;
		MPI_Allreduce(&dbs.m_pts->m_i_num_local_points, &global_owned_local_point_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		dbs.m_global_owned_local_point_count_root = global_owned_local_point_count;
		
		dbs.m_all_gid_gathered_upper = dbs.m_pts->m_all_local_point_gID_upper;
		dbs.m_all_gid_gathered_lower = dbs.m_pts->m_all_local_point_gID_lower;

		// copy the mst to root, but only for the first time;
                proc_count = nproc;

                while((proc_count % 2 == 0) && (proc_count > 1))
                {
                        round_count++;
                        proc_count = proc_count / 2;
                }

                #ifdef _DEBUG
                if(proc_count != 1)
                        cout << "Number of procs is not a multiple of 2" << endl;
                #endif


		vector <int> vec_parents;
                vec_parents.reserve(global_owned_local_point_count);
                for(i = 0; i < global_owned_local_point_count; i++)
                        vec_parents.push_back(i);
		
		
		MPI_Datatype mst_node_type;
                MPI_Type_contiguous(sizeof(MST_Node), MPI_BYTE, &mst_node_type);
                MPI_Type_commit(&mst_node_type);

		int rank_from, rank_to, recv_size[5], send_size[5]; //mst_node_recv_size, mst_node_send_size;
		int tag = 200, root_u, root_v, z;
                MPI_Request req_send[nproc], req_recv[nproc];
                MPI_Status stat_send[nproc], stat_recv;
		int send_count, recv_count;
		MST_Node mst_node;
		int rtag, rsource, rpos;
		double inter2, comp_time, comp_time_total = 0, comp_time2;
		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging SMART using Kruskal - preprocessing time " << MPI_Wtime() - inter << endl;
		inter = MPI_Wtime();
		inter2 = inter;
		#else
                if(rank == proc_of_interest) cout << " mer_pre " << MPI_Wtime() - inter;
                inter = MPI_Wtime();
                inter2 = inter;
		#endif

		//cout << "rank " << rank  << " MST size " << dbs.m_pq_mst_node->size() << endl;
		for(round_pos = 0; round_pos < round_count; round_pos++)
		{
			#ifdef _DEBUG
			int local_mst_edge_count_round = dbs.m_pq_mst_node->size(), local_mst_edge_count_round_total;			
			MPI_Allreduce(&local_mst_edge_count_round, &local_mst_edge_count_round_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);	
			if(rank == proc_of_interest) cout << "Round " << round_pos << " local_mst_edge_count_round " << local_mst_edge_count_round_total << endl;	
			#endif

			//tag  = 200; // a random start
			if(rank % (int)(pow (2, round_pos)) == 0)
			{
				// rank is either a sender or a receiver
				if(rank % (int)(pow (2, round_pos + 1)) == 0)
				{
					// rank is receiver
					rank_from = rank + (int)(pow (2, round_pos));					

					recv_size[2] = 0; // corresponing global ID count
					recv_size[3] = 0; // global ID upper limit
					recv_size[4] = 0; // global ID lower limit
		
					MPI_Recv(&recv_size[0], 5, MPI_INT, rank_from, tag, MPI_COMM_WORLD, &stat_recv);
					//cout << "rank " << rank << " round_pos " << round_pos << " tag " << tag << " rank_to " << rank_from << endl;

					recv_count = 0;
					
					vector <MST_Node> vec_mst_node_recv;
					if(recv_size[2] > 0)
					{
                                        	vec_mst_node_recv.resize(recv_size[2]);
						MPI_Irecv(&vec_mst_node_recv[0], recv_size[2], mst_node_type, rank_from, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]); 
					}
						
					if(dbs.m_all_gid_gathered_upper < recv_size[3])
						dbs.m_all_gid_gathered_upper = recv_size[3];

					if(dbs.m_all_gid_gathered_lower > recv_size[4])
						dbs.m_all_gid_gathered_lower = recv_size[4];

					// reset that range from parents vector
					for(i = dbs.m_all_gid_gathered_lower; i < dbs.m_all_gid_gathered_upper; i++)
						vec_parents[i] = i;

					for(i = 0; i < recv_count; i++)
					{
						MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
						rtag = stat_recv.MPI_TAG;
						rsource = stat_recv.MPI_SOURCE;
		
						if(rtag == tag + 1)
						{
							// nothing to do with the core ditances for the time being		
						}
						else if(rtag == tag + 2)
						{
							// nothing to do for the time being
						}
						else if(rtag == tag + 3)
						{
							//#ifdef _DEBUG
							comp_time = MPI_Wtime();
							//#endif
							//cout << "rank  " << rank << " round_pos " << round_pos << " MST tree size " << dbs.m_pq_mst_node->size() << " vector " << vec_mst_node_recv.size() <<endl;
							for(j = 0; j < vec_mst_node_recv.size(); j++)
                                                		dbs.m_pq_mst_node->push(vec_mst_node_recv[j]);

							//cout << "rank= " << rank << " round_pos " << round_pos << " MST tree size " << dbs.m_pq_mst_node->size() << " vector " << vec_mst_node_recv.size() <<endl;
							// run Kruskals algorithm on the edges in the priority queue
							priority_queue<MST_Node>* pq_mst_node;
							pq_mst_node = new priority_queue<MST_Node>;

							if(round_pos == round_count - 1)
								dbs.m_mst_root.reserve(global_owned_local_point_count);

							while(dbs.m_pq_mst_node->size() > 0)
							{
								mst_node = dbs.m_pq_mst_node->top();
								dbs.m_pq_mst_node->pop();
								

								// now perform union_find operation
                        					root_u = mst_node.m_reachable_to;
                        					root_v = mst_node.m_reachable_from;
								
								while(vec_parents[root_u] != vec_parents[root_v])
                                				{
                                        				if(vec_parents[root_u] < vec_parents[root_v])
                                        				{
                                                				if(root_u == vec_parents[root_u])
                                                				{
                                                        				vec_parents[root_u] = vec_parents[root_v];
											if(round_pos == round_count - 1)
												dbs.m_mst_root.push_back(mst_node);
                                                        				else
												pq_mst_node->push(mst_node);
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
                                				                        vec_parents[root_v] = vec_parents[root_u];
											if(round_pos == round_count - 1)
												dbs.m_mst_root.push_back(mst_node);	
											else
												pq_mst_node->push(mst_node);
                                                        				break;
                                                				}

				                                                z = vec_parents[root_v];
                                				                vec_parents[root_v] = vec_parents[root_u];
                                                				root_v = z;                  
                                        				}
                                				}
							}

							delete dbs.m_pq_mst_node;

							dbs.m_pq_mst_node = pq_mst_node;
							pq_mst_node = NULL;
							//#ifdef _DEBUG
							comp_time2 = MPI_Wtime() - comp_time;
                                                        comp_time_total += comp_time2;
                                                        //#endif
							//cout << "rank==" << rank << " round_pos " << round_pos << " MST tree size " << dbs.m_pq_mst_node->size() << " vector " << vec_mst_node_recv.size() <<endl;
						}
					}
				
					vec_mst_node_recv.clear();
				}
				else
				{
					vec_parents.clear(); // this will reduce the required memory
					// rank is sender
					rank_to = rank - (int)(pow (2, round_pos));
					
					send_size[2] = dbs.m_pq_mst_node->size(); // mst_node count 
					send_size[3] = dbs.m_all_gid_gathered_upper;
					send_size[4] = dbs.m_all_gid_gathered_lower;

					//cout << "rank " << rank << " round_pos " << round_pos << " tag " << tag << " rank_from " << rank_to << endl;	
					MPI_Send(&send_size[0], 5, MPI_INT, rank_to, tag, MPI_COMM_WORLD);

					send_count = 0;
					
					// send the MSTs in PQ
					vector <MST_Node> vec_mst_node_send;
					vec_mst_node_send.reserve(dbs.m_pq_mst_node->size());
					while(dbs.m_pq_mst_node->size() > 0)
					{
						vec_mst_node_send.push_back(dbs.m_pq_mst_node->top());
						dbs.m_pq_mst_node->pop();
					}
						
					if(send_size[2] > 0)
					{
						MPI_Isend(&vec_mst_node_send[0], send_size[2], mst_node_type, rank_to, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);		
					}
					
					if(send_count > 0)
                        			MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
					
					vec_mst_node_send.clear(); 
				}
			}

			tag += 3;
			//MPI_Barrier(MPI_COMM_WORLD);
			#ifdef _DEBUG
                	if(rank == proc_of_interest) cout << "Round " << round_pos << " total_time " << MPI_Wtime() - inter2 << " computation time " << comp_time2 << endl;
                	inter2 = MPI_Wtime();
                	#endif
			/*#ifdef _CUT_OFF_ACTIVE
			int stop_value = 0, stop_sum = 0;
			if(comp_time2 > _CUT_OFF_TIME)
				stop_value = 1;
			MPI_Allreduce(&stop_value, &stop_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if(stop_sum > 0)
				break;
			#endif
			*/
			#ifdef _CUT_OFF_ACTIVE
				int stop_value = 0, stop_sum = 0;
				#ifdef _USE_ROUND_TIME
				if(comp_time2 > _CUT_OFF_TIME)
					stop_value = 1;
				#endif
			
				#ifdef _USE_ROUND_COUNT
				if(round_pos == _CUT_OFF_ROUND)
					stop_value = 1;	
				#endif
				MPI_Allreduce(&stop_value, &stop_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				if(stop_sum > 0)
					break;
			#endif

		}
	
		vector <MST_Node> vec_mst_node_send_rest;
		vector <vector <MST_Node> > vec_vec_mst_node_recv;	
		int rest_mst_edges_per_proc[nproc], rest_mst_edges_this = dbs.m_pq_mst_node->size();	

		if(round_pos < round_count)
		{
			// send everyone to processor proc_of_interest then
                	MPI_Gather(&rest_mst_edges_this, 1, MPI_INT, &rest_mst_edges_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);

			send_count = 0;
			if(rest_mst_edges_this > 0 && rank != proc_of_interest)
			{
				vec_parents.clear();
				// send the MSTs in PQ
				vec_mst_node_send_rest.clear();
				vec_mst_node_send_rest.reserve(rest_mst_edges_this);
				while(dbs.m_pq_mst_node->size() > 0)
				{
					vec_mst_node_send_rest.push_back(dbs.m_pq_mst_node->top());
					dbs.m_pq_mst_node->pop();
				}
						
				MPI_Isend(&vec_mst_node_send_rest[0], rest_mst_edges_this, mst_node_type, proc_of_interest, tag + 3, MPI_COMM_WORLD, &req_send[send_count++]);			
			}

			recv_count = 0;
			if(rank == proc_of_interest)
			{
				vec_vec_mst_node_recv.resize(nproc);
				for(i = 0; i < nproc; i++)
				{
					if(i == rank)
						continue;
					
					if(rest_mst_edges_per_proc[i] > 0)
					{
                                        	vec_vec_mst_node_recv[i].resize(rest_mst_edges_per_proc[i]);
						MPI_Irecv(&vec_vec_mst_node_recv[i][0], rest_mst_edges_per_proc[i], mst_node_type, i, tag + 3, MPI_COMM_WORLD, &req_recv[recv_count++]);
					}					
				}
				
				//#ifdef _DEBUG
				comp_time = MPI_Wtime();
				//#endif
					
				// reset that range from parents vector
				for(i = 0; i < global_owned_local_point_count; i++)
					vec_parents[i] = i;

				//#ifdef _DEBUG
				comp_time2 = MPI_Wtime() - comp_time;
                                comp_time_total += comp_time2;
                                //#endif

				for(i = 0; i < recv_count; i++)
				{
					MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
					rtag = stat_recv.MPI_TAG;
					rsource = stat_recv.MPI_SOURCE;
		
					if(rtag == tag + 3)
					{
						for(j = 0; j < vec_vec_mst_node_recv[rsource].size(); j++)
                                                	dbs.m_pq_mst_node->push(vec_vec_mst_node_recv[rsource][j]);

						vec_vec_mst_node_recv[rsource].clear();
					}
				}
			}
			else
			{
				if(send_count > 0)
                        		MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
					
				vec_mst_node_send_rest.clear();	
			}

			vec_vec_mst_node_recv.clear();
			vec_mst_node_send_rest.clear();

			if(rank == proc_of_interest)
			{
				//#ifdef _DEBUG
				comp_time = MPI_Wtime();
				//#endif
			
				dbs.m_mst_root.reserve(global_owned_local_point_count);

				while(dbs.m_pq_mst_node->size() > 0)
				{
					mst_node = dbs.m_pq_mst_node->top();
					dbs.m_pq_mst_node->pop();
								

					// now perform union_find operation
	                        	root_u = mst_node.m_reachable_to;
        	                	root_v = mst_node.m_reachable_from;
									
					while(vec_parents[root_u] != vec_parents[root_v])
                                	{
                                		if(vec_parents[root_u] < vec_parents[root_v])
	                                        {
        	                                	if(root_u == vec_parents[root_u])
                	                               	{
                        	                        	vec_parents[root_u] = vec_parents[root_v];
								dbs.m_mst_root.push_back(mst_node);
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
                                				vec_parents[root_v] = vec_parents[root_u];
								dbs.m_mst_root.push_back(mst_node);	
                                                        	break;
                                                	}

				                	z = vec_parents[root_v];
                                			vec_parents[root_v] = vec_parents[root_u];
                                                	root_v = z;                  
                                        	}
                                	}
				}

				//#ifdef _DEBUG
				comp_time2 = MPI_Wtime() - comp_time;
                                comp_time_total += comp_time2;
                               	//#endif
			}			
		}

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Round " << round_count << " total_time " << MPI_Wtime() - inter2 << " computation time " << comp_time2 << endl;
                inter2 = MPI_Wtime();
                #endif
	
		vec_parents.clear();
		MPI_Type_free(&mst_node_type);

		#ifdef _DEBUG
		if(rank == proc_of_interest) cout << "Computation time while pairwise_local_mst_merging SMART " << comp_time_total << endl;
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging SMART using Kruskal - gathered and extracted global MST time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
		#else
                if(rank == proc_of_interest) cout << " mer_t " << MPI_Wtime() - inter;
                if(rank == proc_of_interest) cout << " comp_t " << comp_time_total;
                inter = MPI_Wtime();		
                #endif

		#ifdef _GATHER_CORE_DIS_AT_A_TIME
		int owned_local_points_count_per_proc[nproc];
                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		
		int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_gathered.resize(global_owned_local_point_count, -1);
			dbs.m_vec_global_point_ID_gathered.resize(global_owned_local_point_count, -1);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &dbs.m_vec_f_core_distance_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &dbs.m_vec_global_point_ID_gathered[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		#endif

		if(rank == proc_of_interest)
		{
			//dbs.m_vec_f_core_distance_gathered_assigned.clear();
			//dbs.m_vec_f_core_distance_gathered_assigned.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root.clear();
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);

			for(i = 0; i < global_owned_local_point_count; i++)
				dbs.m_vec_f_core_distance_root[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
				//dbs.m_vec_f_core_distance_gathered_assigned[dbs.m_vec_global_point_ID_gathered[i]] = dbs.m_vec_f_core_distance_gathered[i];
		}

		dbs.m_vec_global_point_ID_gathered.clear();
		dbs.m_vec_f_core_distance_gathered.clear();

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "pairwise_local_mst_merging SMART using Kruskal - assiging gid and core_distance time " << MPI_Wtime() - inter << endl;
                inter = MPI_Wtime();
		#else
                if(rank == proc_of_interest) cout << " ga_gid_cd " << MPI_Wtime() - inter;
                inter = MPI_Wtime();
                #endif
	
		if(nproc == 1)
		{
			dbs.m_vec_f_core_distance_root.clear();
                        dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			dbs.m_vec_f_core_distance_root = dbs.m_vec_f_core_distance;
			
			dbs.m_mst_root.clear();
			dbs.m_mst_root.reserve(dbs.m_pq_mst_node->size());
			
			while(dbs.m_pq_mst_node->size() > 0)
			{
				mst_node = dbs.m_pq_mst_node->top();
				dbs.m_pq_mst_node->pop();
								

				// now perform union_find operation
                        	root_u = mst_node.m_reachable_to;
                        	root_v = mst_node.m_reachable_from;
								
				while(vec_parents[root_u] != vec_parents[root_v])
                                {
                                        if(vec_parents[root_u] < vec_parents[root_v])
                                	{
                                                if(root_u == vec_parents[root_u])
                                                {
                                                        vec_parents[root_u] = vec_parents[root_v];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);
                                                        //else
							//pq_mst_node->push(mst_node);
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
                                			vec_parents[root_v] = vec_parents[root_u];
							//if(round_pos == round_count - 1)
							dbs.m_mst_root.push_back(mst_node);	
							//else
							//pq_mst_node->push(mst_node);
                                                        break;
                                                }

				                z = vec_parents[root_v];
                               			vec_parents[root_v] = vec_parents[root_u];
                                                root_v = z;                  
                                        }
                                }
			}
		}

                if(rank == proc_of_interest) cout << "Time taken by the merging stage using Kruskal (smart) " << MPI_Wtime() - start << endl;
                if(rank == proc_of_interest) cout << "Total edges in global MST " << dbs.m_mst_root.size() << endl;
	}

	void gather_local_msts_and_compute_adj_root(ClusteringAlgo& dbs)
	{
		double start = MPI_Wtime();

		int rank, nproc, i, j;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		int owned_local_points_count_per_proc[nproc];
                int global_owned_local_point_count;

                MPI_Gather(&dbs.m_pts->m_i_num_local_points, 1, MPI_INT, &owned_local_points_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
                MPI_Allreduce(&dbs.m_pts->m_i_num_local_points, &global_owned_local_point_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
		dbs.m_global_owned_local_point_count_root = global_owned_local_point_count;

		int local_mst_edge_count_this = dbs.m_mst.size();
		int local_mst_edge_count_per_proc[nproc];
		int global_mst_edge_count;

		MPI_Gather(&local_mst_edge_count_this, 1, MPI_INT, &local_mst_edge_count_per_proc[0], 1, MPI_INT, proc_of_interest, MPI_COMM_WORLD);
		MPI_Allreduce(&local_mst_edge_count_this, &global_mst_edge_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		#ifdef _DEBUG
		if(rank == proc_of_interest) cout << "Total edges in local MSTs " << global_mst_edge_count << endl;
		#endif

		//dbs.m_vec_local_msts_root.clear();
		//dbs.m_vec_local_msts_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv		

		int start_pos[nproc], pos = 0;
		dbs.m_vec_f_core_distance_root.clear();
		dbs.m_vec_f_core_distance_root.resize(1); // assign at least one element as we are accessing the 0th element at MPI_Gatherv 

		vector <float> 	vec_f_core_distance_root;
		vec_f_core_distance_root.resize(1);
		vector <int>	vec_global_point_ID_root;
		vec_global_point_ID_root.resize(1);

		if(proc_of_interest == rank)
		{
			dbs.m_vec_f_core_distance_root.resize(global_owned_local_point_count, -1);
			vec_f_core_distance_root.resize(global_owned_local_point_count);
			vec_global_point_ID_root.resize(global_owned_local_point_count);
			pos = 0;
			for(i = 0; i < nproc; i++)
			{
				start_pos[i] = pos;
				pos += owned_local_points_count_per_proc[i];
			}	
		}

		//MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_vec_f_core_distance.size(), MPI_FLOAT, &dbs.m_vec_f_core_distance_root[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_vec_f_core_distance[0], dbs.m_pts->m_i_num_local_points, MPI_FLOAT, &vec_f_core_distance_root[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_FLOAT, proc_of_interest, MPI_COMM_WORLD);

		MPI_Gatherv(&dbs.m_pts->m_vec_global_point_ID[0], dbs.m_pts->m_i_num_local_points, MPI_INT, &vec_global_point_ID_root[0], &owned_local_points_count_per_proc[0], &start_pos[0], MPI_INT, proc_of_interest, MPI_COMM_WORLD);

		vector < vector <MST_Node> > vec_vec_local_msts_root;
		if(proc_of_interest == rank)
		{
			vec_vec_local_msts_root.resize(nproc);
			for(i = 0; i < nproc; i++)
			{
				vec_vec_local_msts_root[i].resize(local_mst_edge_count_per_proc[i]);
			}
		}	
		
		MPI_Datatype mst_node_type;
		MPI_Type_contiguous(sizeof(MST_Node), MPI_BYTE, &mst_node_type); 
		MPI_Type_commit(&mst_node_type);
		
		int tag = 200, send_count, recv_count;
		MPI_Request req_send[nproc], req_recv[nproc];
		MPI_Status stat_send[nproc], stat_recv;

		recv_count = 0;
		
		if(proc_of_interest == rank)
		{
			for(i = 0; i < nproc; i++)
			{
				if(local_mst_edge_count_per_proc[i] > 0)
				{
					MPI_Irecv(&vec_vec_local_msts_root[i][0], local_mst_edge_count_per_proc[i], mst_node_type, i, tag, MPI_COMM_WORLD, &req_recv[recv_count++]);
				}
			}
		}

		send_count = 0;

            	if(local_mst_edge_count_this > 0)
            	{
                	MPI_Isend(&dbs.m_mst[0], local_mst_edge_count_this, mst_node_type, proc_of_interest, tag, MPI_COMM_WORLD, &req_send[send_count++]);
		}       	

		int rtag, rsource, rpos;

		// allocate mememry
		MST_Edge me, me_other;

		if(proc_of_interest == rank)
		{

			// assign the core distances now
			for(i = 0; i < global_owned_local_point_count; i++)
			{
				dbs.m_vec_f_core_distance_root[vec_global_point_ID_root[i]] = vec_f_core_distance_root[i];
			}

			vec_global_point_ID_root.clear();
			vec_f_core_distance_root.clear();

			// create the adjacency list
			dbs.m_vec_vec_adjacency_list_root.clear();
			dbs.m_vec_vec_adjacency_list_root.resize(dbs.m_global_owned_local_point_count_root);
			
			for(i = 0; i < recv_count; i++)
			{
				MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
				rtag = stat_recv.MPI_TAG;
				rsource = stat_recv.MPI_SOURCE;
		
				if(rtag == tag)
				{
					// memory efficient
					//for(i = 0; i < dbs.m_vec_local_msts_root.size(); i++)
					while(vec_vec_local_msts_root[rsource].size() > 0)
					{
						j = vec_vec_local_msts_root[rsource].size() - 1;
			
						me.m_other_end = vec_vec_local_msts_root[rsource][j].m_reachable_to;
						me.m_f_edge_weight = max(dbs.m_vec_f_core_distance_root[vec_vec_local_msts_root[rsource][j].m_reachable_from], vec_vec_local_msts_root[rsource][j].m_f_between_distance);
                        			dbs.m_vec_vec_adjacency_list_root[vec_vec_local_msts_root[rsource][j].m_reachable_from].push_back(me);
		
						me_other.m_other_end = vec_vec_local_msts_root[rsource][j].m_reachable_from;
                        			me.m_f_edge_weight = max(dbs.m_vec_f_core_distance_root[vec_vec_local_msts_root[rsource][j].m_reachable_to], vec_vec_local_msts_root[rsource][j].m_f_between_distance);
						dbs.m_vec_vec_adjacency_list_root[vec_vec_local_msts_root[rsource][j].m_reachable_to].push_back(me_other);

						// trim the last element as it has already been processed
						vec_vec_local_msts_root[rsource].resize(j);
					}
				}
			}
		}

		vec_global_point_ID_root.clear();
                vec_f_core_distance_root.clear();


		if(send_count > 0)
                        MPI_Waitall(send_count, &req_send[0], &stat_send[0]);
				
		//MPI_Gatherv(&dbs.m_mst[0], dbs.m_mst.size(), mst_node_type, &dbs.m_vec_local_msts_root[0], &local_mst_edge_count_per_proc[0], &start_pos[0], mst_node_type, proc_of_interest, MPI_COMM_WORLD);	
		// processor zero now has everything
		
		vec_vec_local_msts_root.clear();

		if(proc_of_interest == rank)
		{
			// filter out same edge
			map<int, float> filter_identical;
                	map<int, float>::iterator filter_identical_it;
		
			for(i = 0; i < dbs.m_vec_vec_adjacency_list_root.size(); i++)
			{
				// find uniqueness
				filter_identical.clear();
				for(j = 0; j < dbs.m_vec_vec_adjacency_list_root[i].size(); j++)
				{
					filter_identical_it = filter_identical.find(dbs.m_vec_vec_adjacency_list_root[i][j].m_other_end);
					if(filter_identical_it != filter_identical.end())
					{
						if(filter_identical_it->second > dbs.m_vec_vec_adjacency_list_root[i][j].m_f_edge_weight)
						{
							filter_identical_it->second = dbs.m_vec_vec_adjacency_list_root[i][j].m_f_edge_weight;		
						}
					}
					else
						filter_identical.insert(pair<int,float>(dbs.m_vec_vec_adjacency_list_root[i][j].m_other_end, dbs.m_vec_vec_adjacency_list_root[i][j].m_f_edge_weight));
				}

				dbs.m_vec_vec_adjacency_list_root[i].clear();

				// MUST THINK
				// HAVE TO SORT LATER
				// sort the edges now based on closeness
				//priority_queue<kdtree2_result> pq;
				//kdtree2_result kr;
				
				//for(filter_identical_it = filter_identical.begin(); filter_identical_it != filter_identical.end(); filter_identical_it++)
				//{
				//	kr.dis = 
				//}
				
				for(filter_identical_it = filter_identical.begin(); filter_identical_it != filter_identical.end(); filter_identical_it++)
				{
					me.m_other_end = filter_identical_it->first;
                                	me.m_f_edge_weight = filter_identical_it->second;
                                	dbs.m_vec_vec_adjacency_list_root[i].push_back(me);
				}

				filter_identical.clear();
			}
		}		
		
		//vector <float> vec_f_core_distance;
		// now we have all the core distance as well, hence proc_of_interest can process the data now sequentially
		MPI_Type_free(&mst_node_type);

		#ifdef _DEBUG
                if(rank == proc_of_interest) cout << "Time taken to gather_local_msts_and_compute_adj_root to root " << MPI_Wtime() - start << endl;
                #else
                if(rank == proc_of_interest) cout << " gather_local_msts_and_compute_adj_root " << MPI_Wtime() - start;
                #endif

	}

	void compute_mst_global_from_mst_local(ClusteringAlgo& dbs)
        {
		double start = MPI_Wtime();
		int rank, nproc, i, j;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		if(rank != proc_of_interest)
			return;
	
		// ONLY root is supposed to do the following
		dbs.m_visited_root.clear();
		dbs.m_visited_root.resize(dbs.m_global_owned_local_point_count_root);
		
		int pid, npid;
                mutable_priority_queue<float,int> orderSeeds;
                mutable_priority_queue<float,int> ::iterator it;
                vector <MST_Edge>* vec_edges;
                vector <MST_Edge>* vec_edges_next;
                map<int,int> other_ends;
                map<int,int>::iterator other_ends_it;

                //vector<int>* ind = dbs.m_kdtree->getIndex();
                MST_Node mst_node;

		for(i = 0; i < dbs.m_global_owned_local_point_count_root; i++)
		{
			pid = i; // could be (*ind)[i] but notice distributed kdt and gloabl ID
			if(!dbs.m_visited_root[pid])
			{
				dbs.m_visited_root[pid] = true;
		
				vec_edges = &dbs.m_vec_vec_adjacency_list_root[pid];
				if(vec_edges->size() > 0)
				{
					updateNeighbors_adjacency_list_based_global_mst_from_local_mst(dbs.m_visited_root, pid, vec_edges, orderSeeds, other_ends);
					while(orderSeeds.size() > 0)
					{
						it = orderSeeds.begin();
						npid = it->second;

						mst_node.m_reachable_to = npid;
						//mst_node.m_f_best_distance = it->first; ????????? BEST DISTANCE SHOULD BE REACHABLE DISTANCE
						mst_node.m_f_between_distance = it->first;
						orderSeeds.erase(it->second);

                                        	other_ends_it = other_ends.find(npid);
						#ifdef _DEBUG
                                                if(other_ends_it == other_ends.end())
                                                	cout << "SOMETHING is Wrong: not getting the other end" << endl;
                                                #endif
                                                if(other_ends_it != other_ends.end())
                                                {
                                                	mst_node.m_reachable_from = other_ends_it->second;
                                                        other_ends.erase(other_ends_it);
                                                }

                                                dbs.m_mst_root.push_back(mst_node);

                                                dbs.m_visited_root[npid] = true;

						vec_edges_next = &dbs.m_vec_vec_adjacency_list_root[npid];
                                                if(vec_edges_next->size() > 0)
							updateNeighbors_adjacency_list_based_global_mst_from_local_mst(dbs.m_visited_root, npid, vec_edges_next, orderSeeds, other_ends);
					}
				}
			}
		}

		orderSeeds.clear();
                vec_edges = NULL;
                vec_edges_next = NULL;
                other_ends.clear();

                if(rank == proc_of_interest) cout << "Time taken by the merging stage at root  " << MPI_Wtime() - start << endl;
                if(rank == proc_of_interest) cout << "Total edges in global MST " << dbs.m_mst_root.size() << endl;
	}

	void updateNeighbors_adjacency_list_based_global_mst_from_local_mst(vector <bool>& visited_root, int pid, vector <MST_Edge>* vec_edges_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in)
        {
		map<int,int>::iterator other_ends_it_in;
		
		int i;
                float newRdist, dis_x_y, oldRdist, new_dist_x_y, old_dist_x_y;
                MST_Edge* p_ed;

		for(i = 0; i < vec_edges_in->size(); i++)
                {
			
                        p_ed = &(*vec_edges_in)[i];
                        //if(!dbs.m_visited_root[p_ed->m_other_end])
                        if(!visited_root[p_ed->m_other_end])
                        {
				newRdist = p_ed->m_f_edge_weight;
				if(orderSeeds_in.find_and_return_value(p_ed->m_other_end, oldRdist))
				{
					if(newRdist < oldRdist)
					{
						orderSeeds_in.update(p_ed->m_other_end, newRdist);
                                                other_ends_it_in = other_ends_in.find(p_ed->m_other_end);
						other_ends_it_in->second = pid;
					}
					else if(newRdist == oldRdist)
                                        {
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


	void updateNeighbors_mst_based_parallel(ClusteringAlgo& dbs, int pid, kdtree2_result_vector& ne_in, mutable_priority_queue<float,int>& orderSeeds_in, map<int,int>& other_ends_in, int tid, map<int,float>& reach_so_far_local_in)
	{
		map<int,int>::iterator other_ends_it_in;
		map<int,float>::iterator reach_so_far_local_it_in;
		int i;

		float newRdist, oldRdist, new_dist_x_y, old_dist_x_y;
		for(i = 0; i < ne_in.size(); i++)
		{
			if(dbs.m_visited[ne_in[i].idx] == false)
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
////////////////////////////////////////////////////////////////////////////////////	
};



