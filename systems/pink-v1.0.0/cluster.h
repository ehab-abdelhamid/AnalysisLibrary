#ifndef __CLUSTER_H
#define __CLUSTER_H

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Macros controlling report output
#define PRINT_TIMING
//#define PRINT_STATS
#define PRINT_MST

//Define data types (precision) for input data and distances
typedef double data_t;
#define MPI_DATA MPI_DOUBLE
#define DATA_MAX DBL_MAX
#define PRINT_DATA "%lf"

typedef double dist_t;
#define DIST_MAX DBL_MAX
#define PRINT_UVDIST "(%d, %d, %lf)\n"

//Declares which distance metric to use for hierarchical clustering
#define DISTANCE_FUNCTION sq_euclid_dist

typedef struct
{
  dist_t distance;
  int left;
  int right;
} Node;

#ifdef PRINT_STATS
//Counters for collecting stats
extern long long distcount;
extern long long cpy;
#endif

//Before using MPI_NODE as a datatype, call init_type_MPI_NODE()
//Call free_type_MPI_NODE() when you no longer need MPI_NODE
#ifdef MPI_INCLUDED
MPI_Datatype MPI_NODE;
#define mpi_call(x) if ((x) != MPI_SUCCESS) {perror("MPI Error"); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#define init_type_MPI_NODE() { \
  mpi_call(MPI_Type_contiguous(sizeof(Node), MPI_BYTE, &MPI_NODE)); \
  mpi_call(MPI_Type_commit(&MPI_NODE)); \
}
#define free_type_MPI_NODE() mpi_call(MPI_Type_free(&MPI_NODE))
#endif

//Define comparator function
int nodecompare(const void* a, const void* b);

//Define distance functions
dist_t sq_euclid_dist(data_t* x, data_t* y, int len, dist_t max);
dist_t norm_cos_dist(data_t* x, data_t* y, int len, dist_t max);
dist_t sq_orb_dist(data_t* x, data_t* y, int len, dist_t max);

int prim_full(int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result);

int prim_full_shrink(int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result);

int prim_bipartite(int leftrows, int rightrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result);

int pslcluster(int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result);
#endif
