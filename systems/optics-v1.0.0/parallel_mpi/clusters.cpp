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

#include "clusters.h"

namespace NWUClustering
{
	Clusters::~Clusters()
	{
		if(m_pts)
		{
			m_pts->m_points.clear(); //resize(extents[0][0]); // there is no clear function
			delete m_pts;
			m_pts = NULL;
		}

		if(m_kdtree)
		{
			delete m_kdtree;
			m_kdtree = NULL;
		}

        if(m_pts_outer)
        {
			m_pts_outer->m_prIDs.clear();
			m_pts_outer->m_ind.clear();
            m_pts_outer->m_points.clear(); //resize(extents[0][0]); // there is no clear function
            delete m_pts_outer;
            m_pts_outer = NULL;
        }

        if(m_kdtree_outer)
        {
            delete m_kdtree_outer;
            m_kdtree_outer = NULL;
        }
	}

	bool Clusters::allocate_outer(int dims)
	{
		if(m_pts_outer == NULL)
        	{
	        	m_pts_outer = new Points_Outer;
        	    	m_pts_outer->m_prIDs.clear();
   	         	m_pts_outer->m_ind.clear();
        	    	m_pts_outer->m_i_dims = dims;
           	 	m_pts_outer->m_i_num_points = 0;
        	}
	}

	bool Clusters::addPoints_NEW(int source, int buf_size, int dims, vector<float>& raw_data)
	{
		int i, j, k, pos, num_points = buf_size / dims;

		//cout << "add points called" << endl; 

		// incorrect dimension
		if(m_pts->m_i_dims != dims)
			return false;

		pos = m_pts->m_i_num_points;
        	m_pts->m_i_num_points += num_points;
        	//m_pts->m_points.resize(extents[m_pts->m_i_num_points][dims]);

		m_pts->m_points.resize(m_pts->m_i_num_points);
                for(int ll = 0; ll < m_pts->m_i_num_points; ll++)
                	m_pts->m_points[ll].resize(dims);
		
		m_pts->m_vec_processor_ID.resize(m_pts->m_i_num_points, -1);
		
		k = 0;
		for(i = 0; i < num_points; i++)
		{
			for(j = 0; j < dims; j++)
				m_pts->m_points[pos][j] = raw_data[k++];
			
			m_pts->m_vec_processor_ID[pos] = source;
			
			pos++;
		}

		//cout << "outer " << m_pts_outer->m_i_num_points << endl;

		return true;
	}
	
	bool Clusters::updatePoints_NEW(vector< vector<int> >& raw_ind, vector< vector<int> >& raw_gid)
	{
		int i, source = -1, j = -1, prev_source = -1, pos;

		pos = m_pts->m_i_num_local_points;
		m_pts->m_vec_local_ind.resize(m_pts->m_i_num_points, -1);
		m_pts->m_vec_global_point_ID.resize(m_pts->m_i_num_points, -1);

		for(i = pos; i < m_pts->m_i_num_points; i++)
		{
			//source = m_pts_outer->m_prIDs[i];
			source = m_pts->m_vec_processor_ID[i];

			if(source != prev_source)
				j = 0;

			m_pts->m_vec_local_ind[i] = raw_ind[source][j];
			m_pts->m_vec_global_point_ID[i] = raw_gid[source][j];
			j++;	

			prev_source = source;
		}
		
		return true;
	}

	int Clusters::read_file(char* infilename, int isBinaryFile)
	{
		ssize_t numBytesRead;
		int     i, j, rank, nproc;
		int num_points, dims;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);

		if(isBinaryFile == 1)
        	{
			// each processor reads it own part ONLY
			ifstream file (infilename, ios::in|ios::binary);
			if(file.is_open())
  			{
				file.read((char*)&num_points, sizeof(int));
				file.read((char*)&dims, sizeof(int));
			
				if(rank == proc_of_interest) cout << "Number of points " << num_points << " Dimensions " << dims << endl;

				// compute the respective segments
				long long sch, lower, upper;

				if(num_points % nproc == 0)
					sch = num_points/nproc;
				else
					sch = num_points/nproc + 1;

				lower = sch * rank;
				upper = sch * (rank + 1);
				if(upper > num_points)
					upper = num_points;

				// allocate memory for points
				m_pts = new Points;
				m_pts->m_i_dims = dims;
				m_pts->m_i_num_points = upper - lower;//num_points;
				m_pts->m_total_point_count = num_points;

				m_pts->m_box = new interval[m_pts->m_i_dims];	
				m_pts->m_vec_global_point_ID.reserve(m_pts->m_i_num_points);

				//allocate memory for the points
				//m_pts->m_points.resize(extents[num_points][dims]);
				//m_pts->m_points.resize(extents[m_pts->m_i_num_points][dims]);
				m_pts->m_points.resize(m_pts->m_i_num_points);
				for(int ll = 0; ll < m_pts->m_i_num_points; ll++)
					m_pts->m_points[ll].resize(dims);

				point_coord_type* pt;					
			
				pt = (point_coord_type*) malloc(dims * sizeof(point_coord_type));
	
				// fseek to the respective position
				file.seekg(lower * dims * sizeof(point_coord_type), ios::cur);

				for (i = 0; i < upper - lower; i++)
				{
					file.read((char*)pt, dims*sizeof(point_coord_type));
					m_pts->m_vec_global_point_ID.push_back(lower + i);     		

					for (j = 0; j < dims; j++)
					{
						m_pts->m_points[i][j] = pt[j];
						
						if(i == 0)
						{	
							m_pts->m_box[j].upper = m_pts->m_points[i][j];
							m_pts->m_box[j].lower = m_pts->m_points[i][j];
						}
						else
						{
							if(m_pts->m_box[j].lower > m_pts->m_points[i][j])
								m_pts->m_box[j].lower = m_pts->m_points[i][j];
							else if(m_pts->m_box[j].upper < m_pts->m_points[i][j])
								m_pts->m_box[j].upper = m_pts->m_points[i][j];
						}
					}
				}
		
				free(pt);
				pt = NULL;

				file.close();
			}	
			else
			{
				if(rank == proc_of_interest) cout << "Error: no such file: " << infilename << endl;
				return -1;
			}
		}
		else
        	{
			if(rank == proc_of_interest) cout << "Invalid parameter" << endl;
			return -1;
		}

		return 0;
		
	}

	void Clusters::writePoints(ostream& o)
	{
		int i, j;
		for (i = 0; i < m_pts->m_i_num_points; i++)
		{	
			for(j = 0; j < m_pts->m_i_dims; j++)
                                o << " " << m_pts->m_points[i][j];
			o << endl;
		}
	}

	int Clusters::build_kdtree()
	{
		if(m_pts == NULL)
		{
			cout << "Point set is empty" << endl;
			return -1;
		}

		//m_kdtree = new kdtree2(m_pts->m_points, true);
		m_kdtree = new kdtree2(m_pts->m_points, false);
		
		if(m_kdtree == NULL)
		{
			cout << "Falied to allocate new kd tree for orginial points" << endl;
			return -1;
		}
		
		return 0;		
	} 

	int Clusters::build_kdtree_outer()
	{
		if(m_pts_outer == NULL)
		{
			cout << "Outer point set is empty" << endl;
			return -1;
		}

		//m_kdtree_outer = new kdtree2(m_pts->m_points_outer, true);
		m_kdtree_outer = new kdtree2(m_pts_outer->m_points, false);
		
		if(m_kdtree_outer == NULL)
		{
			cout << "Falied to allocate new kd tree for outer points" << endl;
			return -1;
		}
		
		return 0;		
	}

	void Clusters::assign_prID_ind_local_points()
	{
		int rank, i;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		m_pts->m_i_num_local_points = m_pts->m_i_num_points;

		m_pts->m_vec_processor_ID.reserve(m_pts->m_i_num_local_points);
		m_pts->m_vec_local_ind.reserve(m_pts->m_i_num_local_points);
			
		for(i = 0; i < m_pts->m_i_num_local_points; i++)
		{
			m_pts->m_vec_processor_ID.push_back(rank);
			m_pts->m_vec_local_ind.push_back(i);
		}


		// also set the upper and lower limit of global IDs
		
	}

	void Clusters::compute_upper_lower_limit_gids()
	{
		int i;
		//int rank;
		//MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		m_pts->m_owned_local_point_gID_upper = INT_MIN;
		m_pts->m_owned_local_point_gID_lower = INT_MAX;
		
		m_pts->m_not_owned_local_point_gID_upper = INT_MIN;
                m_pts->m_not_owned_local_point_gID_lower = INT_MAX;

		m_pts->m_all_local_point_gID_upper = INT_MIN;
                m_pts->m_all_local_point_gID_lower = INT_MAX;
		
		for(i = 0; i < m_pts->m_i_num_local_points; i++)
		{
			if(m_pts->m_vec_global_point_ID[i] > m_pts->m_owned_local_point_gID_upper)
				m_pts->m_owned_local_point_gID_upper = m_pts->m_vec_global_point_ID[i];

			if(m_pts->m_vec_global_point_ID[i] < m_pts->m_owned_local_point_gID_lower)
                                m_pts->m_owned_local_point_gID_lower = m_pts->m_vec_global_point_ID[i];
		}

		for(i = m_pts->m_i_num_local_points; i < m_pts->m_i_num_points; i++)
                {
                        if(m_pts->m_vec_global_point_ID[i] > m_pts->m_not_owned_local_point_gID_upper)
                                m_pts->m_not_owned_local_point_gID_upper = m_pts->m_vec_global_point_ID[i];

                        if(m_pts->m_vec_global_point_ID[i] < m_pts->m_not_owned_local_point_gID_lower)
                                m_pts->m_not_owned_local_point_gID_lower = m_pts->m_vec_global_point_ID[i];
                }
		
		if(m_pts->m_owned_local_point_gID_upper > m_pts->m_not_owned_local_point_gID_upper)
			m_pts->m_all_local_point_gID_upper = m_pts->m_owned_local_point_gID_upper;
		else
			m_pts->m_all_local_point_gID_upper = m_pts->m_not_owned_local_point_gID_upper;

		if(m_pts->m_owned_local_point_gID_lower < m_pts->m_not_owned_local_point_gID_lower)
                        m_pts->m_all_local_point_gID_lower = m_pts->m_owned_local_point_gID_lower;
                else
                        m_pts->m_all_local_point_gID_lower = m_pts->m_not_owned_local_point_gID_lower;

	} 
}
