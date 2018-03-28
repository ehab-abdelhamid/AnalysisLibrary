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

#include "geometric_partitioning.h"

namespace NWUClustering
{
	void get_extra_points(ClusteringAlgo& dbs)
	{
		#ifdef _DEBUG_GP
		MPI_Barrier(MPI_COMM_WORLD);
		double end, start = MPI_Wtime(), total_comp_time = 0, comp_time;
		#endif

		int k, rank, nproc, i, j;
        	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		interval* local_box = new interval[dbs.m_pts->m_i_dims];
		compute_local_bounding_box(dbs, local_box);

		float eps = sqrt(dbs.m_epsSquare);

		for(i = 0; i < dbs.m_pts->m_i_dims; i++)
		{
			local_box[i].upper += eps;
			local_box[i].lower -= eps; 
		}

		// all together all the extending bounding box
		interval* gather_local_box = new interval[dbs.m_pts->m_i_dims * nproc];

        	// gather the local bounding box first
        	MPI_Allgather(local_box, sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, gather_local_box,
        						sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, MPI_COMM_WORLD);

		bool if_inside, overlap;
		int count = 0, gcount;

		vector <float> empty;
		vector <vector <float> > send_buf;
		vector <vector <float> > recv_buf;
		send_buf.resize(nproc, empty);
		recv_buf.resize(nproc, empty);

	        vector <int> empty_i;
        	vector <vector <int> > send_buf_ind;
        	vector <vector <int> > recv_buf_ind;
        	send_buf_ind.resize(nproc, empty_i);
        	recv_buf_ind.resize(nproc, empty_i);

        	vector <vector <int> > send_buf_gid;
        	vector <vector <int> > recv_buf_gid;
        	send_buf_gid.resize(nproc, empty_i);
        	recv_buf_gid.resize(nproc, empty_i);

		#ifdef _DEBUG_GP
		MPI_Barrier(MPI_COMM_WORLD);
		end = MPI_Wtime();
		start = end;
		#endif

		for(k = 0; k < nproc; k++)
		{
			if(k == rank) // self
				continue;

			// check the two extended bounding box of proc rank and k. If they don't overlap, there must be no points 
			// SHOULD SUBTRACT EPS			

			overlap = true;
			for(j = 0; j < dbs.m_pts->m_i_dims; j++)
			{
				if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower < gather_local_box[k * dbs.m_pts->m_i_dims +j].lower)
				{
					//if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].upper < gather_local_box[k * dbs.m_pts->m_i_dims + j].lower)
					if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].upper - gather_local_box[k * dbs.m_pts->m_i_dims + j].lower < eps)	
					{
						overlap = false;
	                    			break;
					}
				}
				else
				{
                    			//if(gather_local_box[k * dbs.m_pts->m_i_dims + j].upper < gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower)
                    			if(gather_local_box[k * dbs.m_pts->m_i_dims + j].upper - gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower < eps)
					{
                        			overlap = false;
                        			break;
                    			}
				}
			}

			// the two bouding boxes are different, so continue to the next processors
			if(overlap == false)
				continue;

			// get the overlapping regions
			for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
			{
				if_inside = true;
				for(j = 0; j < dbs.m_pts->m_i_dims; j++)
				{
					if(dbs.m_pts->m_points[i][j] < gather_local_box[k * dbs.m_pts->m_i_dims + j].lower || 
						dbs.m_pts->m_points[i][j] > gather_local_box[k * dbs.m_pts->m_i_dims + j].upper)
					{
						if_inside = false;
						break;
					}
				}
				
				if(if_inside == true)
				{
					for(j = 0; j < dbs.m_pts->m_i_dims; j++)
						send_buf[k].push_back(dbs.m_pts->m_points[i][j]);
				
					send_buf_ind[k].push_back(dbs.m_pts->m_vec_local_ind[i]);
					send_buf_gid[k].push_back(dbs.m_pts->m_vec_global_point_ID[i]);

					count++;
				}
			}
		}

		//cout << "rank " << rank << " copied all points to send"<< endl;

        	#ifdef _DEBUG_GP
		MPI_Barrier(MPI_COMM_WORLD);
        	end = MPI_Wtime();
        	start = end;
        	#endif

		// now send buf have all the points. Send the size first to everyone
		vector <int> send_buf_size, recv_buf_size;
		send_buf_size.resize(nproc, 0);
		recv_buf_size.resize(nproc, 0);
		
		for(i = 0; i < nproc; i++)
			send_buf_size[i] = send_buf[i].size();
	
		MPI_Alltoall(&send_buf_size[0], 1, MPI_INT, &recv_buf_size[0], 1, MPI_INT, MPI_COMM_WORLD);

		int tag = 200, send_count, recv_count;
		MPI_Request req_send[3 * nproc], req_recv[3 * nproc];
		MPI_Status stat_send[3 * nproc], stat_recv;

		recv_count = 0;
		for(i = 0; i < nproc; i++)
		{
			if(recv_buf_size[i] > 0)
			{
				recv_buf[i].resize(recv_buf_size[i], 0);
				recv_buf_ind[i].resize(recv_buf_size[i] / dbs.m_pts->m_i_dims, -1);
				recv_buf_gid[i].resize(recv_buf_size[i] / dbs.m_pts->m_i_dims, -1);

				MPI_Irecv(&recv_buf[i][0], recv_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_recv[recv_count++]);
				MPI_Irecv(&recv_buf_ind[i][0], recv_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &req_recv[recv_count++]);
				MPI_Irecv(&recv_buf_gid[i][0], recv_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 2, MPI_COMM_WORLD, &req_recv[recv_count++]);
			}
		}

		send_count = 0;
        	for(i = 0; i < nproc; i++)
        	{
            		if(send_buf_size[i] > 0)
            		{
                		MPI_Isend(&send_buf[i][0], send_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_send[send_count++]);
				MPI_Isend(&send_buf_ind[i][0], send_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &req_send[send_count++]);
				MPI_Isend(&send_buf_gid[i][0], send_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 2, MPI_COMM_WORLD, &req_send[send_count++]);
			}
        	}

		int rtag, rsource, rpos;


		#if _GET_PARTION_STAT == 0
		//dbs.allocate_outer(dbs.m_pts->m_i_dims);
		#endif

		for(i = 0; i < recv_count; i++)
		{
			MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
		
			rtag = stat_recv.MPI_TAG;
			rsource = stat_recv.MPI_SOURCE;
		
			if(rtag == tag)
			{
				comp_time = MPI_Wtime();
				// process the request
				#if _GET_EXTRA_POINT_STAT == 0  // WHY THIS IS HERE??????????????????????????????????????????????????????????????
				//dbs.addPoints(rsource, recv_buf_size[rsource], dbs.m_pts->m_i_dims, recv_buf[rsource]);
				dbs.addPoints_NEW(rsource, recv_buf_size[rsource], dbs.m_pts->m_i_dims, recv_buf[rsource]);
				#endif
				total_comp_time += (MPI_Wtime() - comp_time);
				recv_buf[rsource].clear();
			}
			else if(rtag == tag + 1)
			{
				// postpond this computation and call update points later
				// processing immediately might lead to invalid computation
			}
			else if(rtag == tag + 2)
			{
				// postpond this computation and call update points later
				// processing immediately might lead to invalid computation
			}
		}	

		// MAY NOT NEED THIS
		if(send_count > 0)
			MPI_Waitall(send_count, &req_send[0], &stat_send[0]);

        	#ifdef _DEBUG_GP
		MPI_Barrier(MPI_COMM_WORLD);
        	end = MPI_Wtime();
        	start = end;
        	#endif

		// got all the points
		// now update the indices of the outer points

		#if _GET_EXTRA_POINT_STAT == 0 // WHY THIS IS HERE??????????????????????????????????????????????????????????????
		//dbs.updatePoints(recv_buf_ind);
		dbs.updatePoints_NEW(recv_buf_ind, recv_buf_gid);
		#endif

		MPI_Reduce(&count, &gcount, 1, MPI_INT, MPI_SUM, proc_of_interest, MPI_COMM_WORLD);
		
		//dbs.m_extra_point_per_processor = gcount/nproc; // save extra point per processor for future use if needed
		
		#ifdef _DEBUG
		if(rank == proc_of_interest)
		{
			//cout << "Total extra point " << gcount << endl;
			cout << "Gathered extra point per process (average) " << gcount/nproc << endl;
		}
		#endif

        	#ifdef _DEBUG_GP
		MPI_Barrier(MPI_COMM_WORLD);
        	end = MPI_Wtime();
        	start = end;
        	#endif

		empty.clear();
		send_buf.clear();
		recv_buf.clear();
		send_buf_size.clear();
		recv_buf_size.clear();
		send_buf_ind.clear();
		recv_buf_ind.clear();
		send_buf_gid.clear();
                recv_buf_gid.clear();

        	delete [] gather_local_box;
		delete [] local_box;
	}

	void start_partitioning(ClusteringAlgo& dbs)
	{
		int r_count, s_count, rank, nproc, i, j, k;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		// compute the local bouding box for each dimention
		interval* box = new interval[dbs.m_pts->m_i_dims];
		//compute_local_bounding_box(dbs, box);
		// we don't need to compute as we have stored this while reading
		
		for(i = 0; i < dbs.m_pts->m_i_dims; i++)
		{	
			box[i].upper = dbs.m_pts->m_box[i].upper;
			box[i].lower = dbs.m_pts->m_box[i].lower;
		}

		// compute the global bouding box for each dimention
		interval* gbox = new interval[dbs.m_pts->m_i_dims];
		compute_global_bounding_box(dbs, box, gbox, nproc);

		// find the loop count for nproc processors
		int internal_nodes, partner_rank, loops, b, color, sub_rank, d, max, sub_nprocs;
		//MPI_Comm new_comm;
		MPI_Status status;

		loops = 0;
		i = nproc;
		internal_nodes = 1;
		while((i = i >> 1) > 0)
		{
			loops++;
			internal_nodes = internal_nodes << 1;
		}
		
		internal_nodes = internal_nodes << 1;
		
		//gbox for each node in the tree [ONLY upto to reaching each processor]
		interval** nodes_gbox = new interval*[internal_nodes];
		for(i = 0; i < internal_nodes; i++)
			nodes_gbox[i] = new interval[dbs.m_pts->m_i_dims];
		
		copy_global_box_to_each_node(dbs, nodes_gbox, gbox, internal_nodes);
		// now each node in the tree has gbox

		vector <float> send_buf;
		vector <int>   invalid_pos_as;
		vector <float> recv_buf;
		vector <int> send_gid;
		vector <int> recv_gid;

	
		int pow2_i;
		float median;

		for(i = 0; i < loops; i++)
		{
			pow2_i = POW2(i);
			b  = nproc - (int) (nproc / pow2_i);
			color = (int)((rank & b) / POW2(loops - i ));
			partner_rank = rank ^ (int)(nproc/POW2(i + 1));

			MPI_Comm new_comm;
            		MPI_Comm_split(MPI_COMM_WORLD, color, rank, &new_comm);
           		MPI_Comm_rank(new_comm, &sub_rank);
	
			if(sub_rank == 0)
			{
				d = 0;
				for(j = 1; j < dbs.m_pts->m_i_dims; j++)
				{
					
					if(nodes_gbox[pow2_i + color][j].upper - nodes_gbox[pow2_i + color][j].lower > 
							nodes_gbox[pow2_i + color][d].upper - nodes_gbox[pow2_i + color][d].lower)
						d = j;
				}
			}	
	
			MPI_Bcast(&d, 1, MPI_INT, 0, new_comm);

			// compute the median in this dimension
			float median  = get_median(dbs, d, new_comm);		

			s_count = get_points_to_send(dbs, send_buf, send_gid, invalid_pos_as, median, d, rank, partner_rank);

			if (rank < partner_rank)
			{
				MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 4, &r_count, 1, MPI_INT, partner_rank, 5, MPI_COMM_WORLD, &status);
				recv_buf.resize(r_count * dbs.m_pts->m_i_dims, 0.0);
			    	MPI_Sendrecv(&send_buf[0], s_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 2,
							&recv_buf[0], r_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, &status);

				recv_gid.resize(r_count, -1);
				MPI_Sendrecv(&send_gid[0], s_count, MPI_INT, partner_rank, 6,
                                                        &recv_gid[0], r_count, MPI_INT, partner_rank, 7, MPI_COMM_WORLD, &status);

				send_buf.clear();
				send_gid.clear();
			}
			else
			{
				MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 5, &r_count, 1, MPI_INT, partner_rank, 4, MPI_COMM_WORLD, &status);
				recv_buf.resize(r_count * dbs.m_pts->m_i_dims, 0.0);
				MPI_Sendrecv(&send_buf[0], s_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 3, 
							&recv_buf[0], r_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 2, MPI_COMM_WORLD, &status);

				recv_gid.resize(r_count, -1);
				MPI_Sendrecv(&send_gid[0], s_count, MPI_INT, partner_rank, 7,
                                                        &recv_gid[0], r_count, MPI_INT, partner_rank, 6, MPI_COMM_WORLD, &status);

				send_buf.clear();
				send_gid.clear();
			}
	
			update_points(dbs, s_count, invalid_pos_as, recv_buf, recv_gid);
			recv_buf.clear();
			recv_gid.clear();
						
			copy_box(dbs, nodes_gbox[LOWER(pow2_i+color)], nodes_gbox[pow2_i+color]);
			nodes_gbox[LOWER(pow2_i+color)][d].upper =  median;
			copy_box(dbs, nodes_gbox[UPPER(pow2_i+color)], nodes_gbox[pow2_i+color]);
			nodes_gbox[UPPER(pow2_i+color)][d].lower =  median;	

			MPI_Comm_free(&new_comm);
		}
        
		// asign processor IDs to local points
		dbs.assign_prID_ind_local_points();

		// free the allocated memory
		for(i = 0; i < nproc; i++)
			delete [] nodes_gbox[i];

		delete [] nodes_gbox;
		delete [] gbox;
		delete [] box;
	}

	void update_points(ClusteringAlgo& dbs, int s_count, vector <int>& invalid_pos_as, vector <float>& recv_buf, vector <int>& recv_gid)
	{
		int i, j, k, l, r_count = recv_buf.size() / dbs.m_pts->m_i_dims, gid_pos;

		//cout << "r_count " << r_count << " s_count " << s_count << endl;
	
		if(r_count >= s_count)
		{
			//invalid_pos_as.reserve(dbs.m_pts->m_i_num_points + r_count - s_count);
			invalid_pos_as.resize(dbs.m_pts->m_i_num_points + r_count - s_count, 1);

			//dbs.m_pts->m_points.reserve(dbs.m_pts->m_i_num_points + r_count - s_count);
			//dbs.m_pts->m_points.resize(extents[dbs.m_pts->m_i_num_points + r_count - s_count][dbs.m_pts->m_i_dims]);

			dbs.m_pts->m_points.resize(dbs.m_pts->m_i_num_points + r_count - s_count);
                        for(int ll = 0; ll < dbs.m_pts->m_i_num_points + r_count - s_count; ll++)
                                dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

			dbs.m_pts->m_vec_global_point_ID.resize(dbs.m_pts->m_i_num_points + r_count - s_count, 1);

			j = 0;
			gid_pos = 0;
			for(i = 0; i < invalid_pos_as.size(); i++)
			{
				if(invalid_pos_as[i] == 1)
				{
					for(k = 0; k < dbs.m_pts->m_i_dims; k++)
						dbs.m_pts->m_points[i][k] = recv_buf[j++];

					dbs.m_pts->m_vec_global_point_ID[i] = recv_gid[gid_pos];
					gid_pos++;
				}
			}			

			dbs.m_pts->m_i_num_points = dbs.m_pts->m_i_num_points + r_count - s_count;
		}
		else
		{
			j = 0;
			i = 0;
			gid_pos = 0;
			if(recv_buf.size() > 0)
			{
				for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
				{
					if(invalid_pos_as[i] == 1)
					{
						for(k = 0; k < dbs.m_pts->m_i_dims; k++)
							dbs.m_pts->m_points[i][k] = recv_buf[j++];
					
						dbs.m_pts->m_vec_global_point_ID[i] = recv_gid[gid_pos];
						gid_pos++;	

						if(j == recv_buf.size())
						{
							i++;
							break;
						}
					}
				}
			}
			
			l = dbs.m_pts->m_i_num_points;
			for( ; i < invalid_pos_as.size(); i++)
			{
				if(invalid_pos_as[i] == 1)
				{
					while(l > i)
					{
						l--;
						if(invalid_pos_as[l] == 0)
							break;
					}

					if(invalid_pos_as[l] == 0)
					{
						for(k = 0; k < dbs.m_pts->m_i_dims; k++)
							dbs.m_pts->m_points[i][k] = dbs.m_pts->m_points[l][k];

						dbs.m_pts->m_vec_global_point_ID[i] = dbs.m_pts->m_vec_global_point_ID[l];
					}
				}
			}

			//dbs.m_pts->m_points.resize(extents[dbs.m_pts->m_i_num_points + r_count - s_count][dbs.m_pts->m_i_dims]);

			dbs.m_pts->m_points.resize(dbs.m_pts->m_i_num_points + r_count - s_count);
                        for(int ll = 0; ll < dbs.m_pts->m_i_num_points + r_count - s_count; ll++)
                                dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

			dbs.m_pts->m_vec_global_point_ID.resize(dbs.m_pts->m_i_num_points + r_count - s_count);
			dbs.m_pts->m_i_num_points = dbs.m_pts->m_i_num_points + r_count - s_count;
		}		
	}

	int get_points_to_send(ClusteringAlgo& dbs, vector <float>& send_buf, vector <int>& send_gid, vector <int>& invalid_pos_as, float median, int d, int rank, int partner_rank)
	{
		int i, count = 0, j;
		send_buf.reserve(dbs.m_pts->m_i_num_points * dbs.m_pts->m_i_dims);
		send_gid.reserve(dbs.m_pts->m_i_num_points);

		invalid_pos_as.clear();
		invalid_pos_as.resize(dbs.m_pts->m_i_num_points, 0);

		for(i = 0; i < dbs.m_pts->m_i_num_points; i++)
		{
			if (rank < partner_rank)
			{
				if(dbs.m_pts->m_points[i][d] > median)
				{
					invalid_pos_as[i] = 1;
					count++;
					for(j = 0; j < dbs.m_pts->m_i_dims; j++)
						send_buf.push_back(dbs.m_pts->m_points[i][j]);
					
					send_gid.push_back(dbs.m_pts->m_vec_global_point_ID[i]);
				}
			}
			else
			{
                		if(dbs.m_pts->m_points[i][d] <= median)
                		{
                    			invalid_pos_as[i] = 1;
                    			count++;
                    			for(j = 0; j < dbs.m_pts->m_i_dims; j++)
                        			send_buf.push_back(dbs.m_pts->m_points[i][j]);

					send_gid.push_back(dbs.m_pts->m_vec_global_point_ID[i]);
                		}
			}
		}

		return count;
	}	

	float get_median(ClusteringAlgo& dbs, int d, MPI_Comm& new_comm)
	{	
		/*
		double sum, g_sum;
		//float sum, g_sum;

		sum = 0.0;
		g_sum = 0.0;

    		for (int k = 0; k < dbs.m_pts->m_i_num_points; k++) 
      			sum += dbs.m_pts->m_points[k][d];
	
		//MPI_Allreduce(&sum, &g_sum, 1, MPI_FLOAT, MPI_SUM, new_comm);
		MPI_Allreduce(&sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, new_comm);

		int g_point_count = 0;
		MPI_Allreduce(&dbs.m_pts->m_i_num_points, &g_point_count, 1, MPI_INT, MPI_SUM, new_comm);

		//return g_sum / static_cast<float>(g_point_count);
		return g_sum / static_cast<double>(g_point_count);
		*/



		// ADDITIONAL CODE
		float median;
		
		vector <float> data;
		data.reserve(dbs.m_pts->m_i_num_points);
		data.resize(dbs.m_pts->m_i_num_points, 0);

		for (int k=0; k < dbs.m_pts->m_i_num_points; k++)
	      	//data.push_back(dbs.m_pts->m_points[k][d]);
	      		data[k] = dbs.m_pts->m_points[k][d];

		/*
       		#ifdef _DEBUG
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        	//MPI_Barrier(new_comm);
        	MPI_Barrier(MPI_COMM_WORLD);
       		//if(rank == proc_of_interest) 
       		cout << "AT median Pos 0 - rank " << rank << " data size " << data.size() << endl;
        	#endif
		*/
	
		median = findKMedian(data, data.size()/2);
		data.clear();
		
		#ifdef _DEBUG
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		//if(rank == proc_of_interest) cout << "Median VALUE " << median << endl;
		#endif
	
       		/*#ifdef _DEBUG
		//int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        	//MPI_Barrier(new_comm);
        	MPI_Barrier(MPI_COMM_WORLD);
       		//if(rank == proc_of_interest) 
		cout << "AT median Pos 1 - rank " << rank << " data size " << data.size() << endl;
        	#endif
		*/


		int proc_count;
		MPI_Comm_size(new_comm, &proc_count);

	    	//#ifdef _DEBUG
		//int rank;
		//MPI_Comm_rank(new_comm, &rank);
    		//if(rank == proc_of_interest) cout << "proc_count " << proc_count << endl;	
		//#endif

		//cout << "proc count " << proc_count << endl;

		vector <float> all_medians;
		all_medians.resize(proc_count, 0);

		MPI_Allgather(&median, sizeof(int), MPI_BYTE, &all_medians[0], sizeof(int), MPI_BYTE, new_comm);	

  		median = findKMedian(all_medians, all_medians.size()/2); 
		all_medians.clear();
		
		return median;	
	}

	void compute_local_bounding_box(ClusteringAlgo& dbs, interval* box)
	{
		int i, j;

		//we assume each processor has at least one point
		//if(dbs.m_pts->m_i_num_points > 0)
		{
			
			for(i = 0; i < dbs.m_pts->m_i_dims; i++)
			{
				box[i].upper = dbs.m_pts->m_points[0][i];
				box[i].lower = dbs.m_pts->m_points[0][i];
			}
		
			for(i = 0; i < dbs.m_pts->m_i_dims; i++)
			{
				for(j = 1; j < dbs.m_pts->m_i_num_points; j++)
				{
					if(box[i].lower > dbs.m_pts->m_points[j][i])
						box[i].lower = dbs.m_pts->m_points[j][i];
					else if(box[i].upper < dbs.m_pts->m_points[j][i])
                	    box[i].upper = dbs.m_pts->m_points[j][i];
				}
				//if(rank == 0)
    	    		//  cout << "proc " << rank << " upper " << box[i].upper << " lower " << box[i].lower << endl;
			}
		}
		/*else
		{
			// WHAT TO DO IF THERE IS NO POINTS
			// FOR THE TIME BEING set 0 
			for(i = 0; i < dbs.m_pts->m_i_dims; i++)
			{
				box[i].upper = 0;
				box[i].lower = 0;
			}
		}*/
	}

	void compute_global_bounding_box(ClusteringAlgo& dbs, interval* box, interval* gbox, int nproc)
	{
		int i, j, k;
	
		interval* gather_local_box = new interval[dbs.m_pts->m_i_dims * nproc];
	
		// gather the local bounding box first
		MPI_Allgather(box, sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, gather_local_box, 
					sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, MPI_COMM_WORLD);

		// compute the global bounding box
		for(i = 0; i < dbs.m_pts->m_i_dims; i++)
		{
			gbox[i].lower = gather_local_box[i].lower;
			gbox[i].upper = gather_local_box[i].upper;
			
			k = i;
			for(j = 0; j < nproc; j++, k += dbs.m_pts->m_i_dims)
			{
				if(gbox[i].lower > gather_local_box[k].lower)
					gbox[i].lower = gather_local_box[k].lower;
				
				if(gbox[i].upper < gather_local_box[k].upper)
                    			gbox[i].upper = gather_local_box[k].upper;
			}
		}
		
		delete [] gather_local_box;
	}

	void copy_global_box_to_each_node(ClusteringAlgo& dbs, interval** nodes_gbox, interval* gbox, int internal_nodes)
	{
        	int i, j;
        	for(i = 0; i < internal_nodes; i++)
        	{
            		for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            		{
                		nodes_gbox[i][j].upper = gbox[j].upper;
                		nodes_gbox[i][j].lower = gbox[j].lower;
            		}
        	}
	}
	
	void copy_box(ClusteringAlgo& dbs, interval* target_box, interval* source_box)
	{
			for(int j = 0; j < dbs.m_pts->m_i_dims; j++)
			{
				target_box[j].upper = source_box[j].upper;
				target_box[j].lower = source_box[j].lower;
			}
	}

	void print_points(ClusteringAlgo& dbs, int rank)
	{
		cout << "proc " << rank << " owned points: " << dbs.m_pts->m_points.size() << " verify " << dbs.m_pts->m_i_num_points << endl;
		int u, v;

		/*
		for(u = 0; u < dbs.m_pts->m_i_num_points; u++)
        {
        	for(v = 0; v < dbs.m_pts->m_i_dims; v++)
            	cout << dbs.m_pts->m_points[u][v] << " ";
            cout << endl;
        }
		*/

		cout << "proc " << rank << " outer points: " << dbs.m_pts_outer->m_points.size() << " verify " << dbs.m_pts_outer->m_i_num_points << endl;

		/*
		for(u = 0; u < dbs.m_pts_outer->m_i_num_points; u++)
        {
            for(v = 0; v < dbs.m_pts_outer->m_i_dims; v++)
                cout << dbs.m_pts_outer->m_points[u][v] << " ";
            cout << endl;
        }
		*/
	}

    void print_box(ClusteringAlgo& dbs, int rank, interval* box)
    {
        cout << "proc " << rank << " bbox: ";

		for(int v = 0; v < dbs.m_pts->m_i_dims; v++)
		{
        	if(v == dbs.m_pts->m_i_dims - 1)
				cout << "(" << box[v].upper << ", " << box[v].lower << ")";
			else
				cout << "(" << box[v].upper << ", " << box[v].lower << "), ";
		}

        cout << endl;
    }
};

