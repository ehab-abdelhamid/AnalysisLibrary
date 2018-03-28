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

#ifndef _UTILS_
#define _UTILS_

#define _USE_ROUND_TIME // one of these two will be active at a time
//#define _USE_ROUND_COUNT
#define _CUT_OFF_TIME 5
#define _CUT_OFF_ROUND 2
#define _CUT_OFF_ACTIVE

// comment the following before production run
//#define _DEBUG
//#define _BORUVKA_FILTER_AND_KRUSKAL
//#define _KRUSKAL_SEQUENTIAL // Doesn't  work for large dataset as memory should overflow
#define _KRUSKAL //pairwise kruskal, Make sure to activate _GATHER_CORE_DIS_AT_A_TIME
//#define _KRUSKAL_SMART  // pairwise kruskal, Make sure to activate _GATHER_CORE_DIS_AT_A_TIME
#define _GATHER_CORE_DIS_AT_A_TIME
//#define _DEBUG_WRITE_INPUT
#define _DEBUG_GP // get extra point
#define _LOCK_BASED 1
#define proc_of_interest 0
#define _GET_EXTRA_POINT_STAT 0
#define _GET_LOCAL_POINT_COUNT 0
//#define _WRITE_CLUS_STAT_ 0
#define NOISE_ID -1

#include <mpi.h>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include "mutable_priority_queue.h"
#include <float.h>
#include <queue>
#include <climits>

#define POW2(x) (1 << (x))
#define ROOT 1
#define LOWER(i) (i<<1)
#define UPPER(i) ((i<<1)+1)
#define PARENT(i) (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)


using namespace std;


typedef float point_coord_type;
typedef vector <vector <point_coord_type> >         array2dfloat;

float findKMedian(vector<float>& A,int K);

#endif
