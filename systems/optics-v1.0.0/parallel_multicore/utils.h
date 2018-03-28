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

#ifndef _UTILS_
#define _UTILS_

// comment the following before production run
//#define _DEBUG
#define _KRUSKAL
//#define _BORUVKA_FILTER_AND_KRUSKAL // this requres _KRUSKAL to be activated
//#define _KRUSKAL_PARALLEL_MERGING
//#define _DEBUG_DETAILS
//#define _DEBUG_2
//#define _DEBUG_2
#define _LOCK_BASED
#define NOISE_ID -1
//#define CLUSTER_ID_STARTS 1

#include <omp.h>
#include <algorithm>
#include "mutable_priority_queue.h"

#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <climits>
#include <string>
#include <float.h>
#include <queue>

using namespace std;


typedef float point_coord_type;
typedef vector <vector <point_coord_type> >         array2dfloat;

float findKMedian(vector<float>& A,int K);

#endif
