/**
 shrink.c
 
 Code for performing distributed single-linkage hierarchical agglomerative
 clustering.  This algorithm is described in "Parallel Hierarchical Clustering on Shared Memory Platforms" by Hendrix, Patwary, Agrawal, Liao, and Choudhary (HiPC 2012).
 
 This code is not designed for serial execution and will fail if using fewer
 than 3 OpenMP processes (see slink.c for serial execution of SLINK algorithm).
 
 The main work of SHRINK is divided into three phases:
 
 PHASE 0:  Data is distributed among the processes and data structures are
           initialized
 PHASE 1:  Each process computes partial dendrogram for local data using
           Prim's algorithm
 PHASE 2:  Partial dendrograms are merged, in binary fashion
 
 Author:  William Hendrix, Northwestern University, 2012
 */

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>

//cluster.h defines SLINK function and data structures
//Also defines macros controlling the printing of timing and stat reports
#include "cluster.h"
//unionfind.h defines UnionFind data structure and functions
#include "unionfind.h"

//Triangular numbers and inverse function
#define tri(x)  ((x) * ((x) - 1) >> 1)
#define itri(x) ((((int) sqrt(((x) << 3) + 1)) + 1) >> 1)

//Utility functions
#define chmalloc(x, type, nelem, msg) if ((x = (type*) malloc((nelem) * sizeof(type))) == NULL) {perror(msg); exit(EXIT_FAILURE);}

/**
 Driver routine for parallel SHRINK algorithm
 */
int main(int argc, char** argv)
{
  data_t* data;
  int* nedges;
  Node** minitree;
  int nthreads, part, pts, size, rem, rank;
  int dim[2]; // = {nrows, ncols}
  int i, j, start;
  FILE* fp;
  
  //PHASE 0:  Initialization and data distribution
  //Check command line argument(s)
  if (argc < 2)
  {
    printf("Usage:  %s <data file> <num_threads>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  //Read in and check number of threads
  nthreads = atoi(argv[2]);
  /*
  if (nthreads < 3)
  {
    printf("Must run %s with at least 3 threads\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  //*/

  #ifdef PRINT_TIMING
  double inittime = 0;
  double readtime = 0;
  double slinktime = 0;
  double mergetime = 0;
  double starttime = omp_get_wtime();
  double now;
  #endif
  
  #ifdef PRINT_STATS
  long long totaldcounts, totalcpy;
  fprintf(stderr, "Warning:  please re-compile with #define PRINT_STATS commented out\n");
  //Performance is bad otherwise
  #endif
  
  //Open file
  fp = fopen(argv[1], "r");
  if (fp == NULL)
  {
    perror("Error opening file");
    exit(EXIT_FAILURE);
  }
  
  //Read in metadata
  if (fread(dim, sizeof(int), 2, fp) != 2)
  {
    perror("Error reading dimensions from file");
    exit(EXIT_FAILURE);
  }
  
  //Allocate arrays
  #ifdef PRINT_TIMING
  now = omp_get_wtime();
  readtime += now - starttime;
  starttime = now;
  #endif
  //Calculate number of partitions and active processes
  part = itri(nthreads);
  nthreads = tri(part);
  omp_set_num_threads(nthreads);
  pts = dim[0] / part;
  rem = dim[0] % part;
  size = pts * dim[1];
  
  chmalloc(data, data_t, dim[0] * dim[1], "Error allocating array for data matrix");
  chmalloc(minitree, Node*, nthreads, "Error allocating array for MSTs");
  chmalloc(nedges, int, nthreads, "Error allocating array for MST sizes");
  //for (i = 0; i < nthreads; i++)
   //chmalloc(minitree[i], Node, (dim[0] - 1), "Error allocating array for MST mergetree");
  
  //Read in data
  #ifdef PRINT_TIMING
  now = omp_get_wtime();
  inittime += now - starttime;
  starttime = now;
  #endif
  if (fread(data, sizeof(data_t), dim[0] * dim[1], fp) != dim[0] * dim[1])
  {
    perror("Error reading in data from file");
    exit(EXIT_FAILURE);
  }
  fclose(fp);
  
  #ifdef PRINT_TIMING
  now = omp_get_wtime();
  readtime += now - starttime;
  starttime = now;
  #endif
  
  #pragma omp parallel private(i, j, rank, start) shared(data, dim, minitree, nthreads, pts, rem, nedges)
  {
    int* rowid;
    Node* mytree;
    chmalloc(rowid, int, (pts + 1) * 2, "Error allocating array for row ids");
    chmalloc(mytree, Node, (dim[0] - 1), "Error allocating array for MST mergetree");
    rank = omp_get_thread_num();
    
    //Calculate left and right partition IDs for local data
    int rpart = itri(rank);
    int lpart = rpart + tri(rpart) - rank - 1;
    
    //PHASE 1:  Run Prim's on local data partitions independently
    //Initialize array containing data point (row #) IDs, used for merging later
    if (lpart < rem)
      start = lpart * (pts + 1);
    else
      start = lpart * pts + rem;
    for (i = 0; i < pts; i++)
      rowid[i] = start + i;
    if (lpart < rem)
      rowid[i++] = start + i;
    if (rpart < rem)
      start = rpart * (pts + 1);
    else
      start = rpart * pts + rem;
    for (j = 0; j < pts; j++)
      rowid[i + j] = start + j;
    if (rpart < rem)
      rowid[i + j++] = start + j;
    
    //Calculate local dendrogram (Prim's)
    //if (pslcluster(i + j, dim[1], data, rowid, sq_euclid_dist, mytree))
    if (prim_full_shrink(i + j, dim[1], data, rowid, DISTANCE_FUNCTION, mytree))
    {
      perror("Error allocating local dendrogram");
      exit(EXIT_FAILURE);
    }
    minitree[rank] = mytree;
    nedges[rank] = i + j - 1;

    //Sort MST edges
    qsort(mytree, i + j - 1, sizeof(Node), nodecompare);

    //Declare and initialize variables
    #ifdef PRINT_TIMING
    if (rank == 0)
    {
      now = omp_get_wtime();
      slinktime += now - starttime;
      starttime = now;
    }
    #endif
    
    //PHASE 2:  Combine partial dendrograms to form full dendrogram (MST)
    Node* mergetree;
    int medges, delta;
    UnionFind* unionfind = init_unionfind(dim[0]);
    chmalloc(mergetree, Node, (dim[0] - 1), "Error allocating buffer for merging MSTs");
    
    //Code for receiving and merging partial dendrograms in binary fashion
    for (delta = 1; delta < nthreads; delta <<= 1)
    {
      #pragma omp barrier
      if ((rank & delta) != 0 | (rank + delta >= nthreads))
        continue;
      
      //"Receive" partial dendrogram (MST) from processor (rank ^ delta = rank + delta)
      Node* mytree = minitree[rank];
      Node* recvtree = minitree[rank + delta];
      int redges = nedges[rank + delta];
      
      dist_t t = minitree[rank][0].distance;
      int x;
      for (x = 1; x < nedges[rank]; x++)
      {
        if (mytree[x].distance > t)
          printf("ERROR!\n");
        t = mytree[x].distance;
      }
      
      //Reset variables
      reset_unionfind(unionfind);
      i = 0;
      j = 0;
      medges = 0;
      
      //Merge partial dendrograms (MSTs)
      while (1)
      {
        //Iteratively add min height merge from either tree
        if (mytree[i].distance < recvtree[j].distance)
        {
          //Local dendrogram has min height merge
          
          //Check whether merge vertices are already in same cluster
          if (unify(mytree[i].left, mytree[i].right, unionfind) == 0)
            mergetree[medges++] = mytree[i];
          
          //If all edges from local dendrogram have been added, iterate through received dendrogram
          if (++i >= nedges[rank])
          {
            while (j < redges)
            {
              if (unify(recvtree[j].left, recvtree[j].right, unionfind) == 0)
                mergetree[medges++] = recvtree[j];
              
              j++;
            }
            break;
          }
        }
        else
        {
          //Received dendrogram has min height merge
          
          //Check whether merge vertices are already in same cluster
          if (unify(recvtree[j].left, recvtree[j].right, unionfind) == 0)
            mergetree[medges++] = recvtree[j];
          
          //If all edges from received dendrogram have been added, iterate through local dendrogram
          if (++j >= redges)
          {
            while (i < nedges[rank])
            {
              if (unify(mytree[i].left, mytree[i].right, unionfind) == 0)
                mergetree[medges++] = mytree[i];
              
              i++;
            }
            break;
          }
        }
      }
      
      //Swap merged dendrogram with local dendrogram
      nedges[rank] = medges;
      minitree[rank] = mergetree;
      mergetree = mytree;
    }
    
    free(rowid);
    free(unionfind);
    free(mergetree);
  }
  
  #ifdef PRINT_TIMING
  mergetime += omp_get_wtime() - starttime;
  #endif
  //ALGORITHM COMPLETE
  
  //Output result, timing, and stats
  //Write out final MST
  #ifdef PRINT_MST
  for (i = 0; i < nedges[0]; i++)
    printf(PRINT_UVDIST, minitree[0][i].left, minitree[0][i].right, minitree[0][i].distance);
  #endif
  
  #ifdef PRINT_TIMING
  fprintf(stderr, "File:            %s\n", argv[1]);
  fprintf(stderr, "Processes:       %d\n", nthreads);
  fprintf(stderr, "Initialization:  %lf s\n", inittime);
  fprintf(stderr, "Reading:         %lf s\n", readtime);
  fprintf(stderr, "Prim's:          %lf s\n", slinktime);
  fprintf(stderr, "Merging:         %lf s\n", mergetime);
  fprintf(stderr, "Total:           %lf s\n", inittime + readtime + slinktime + mergetime);
  #endif
  
  #ifdef PRINT_STATS
  fprintf(stderr, "Distance_calcs:  %lld\n", totaldcounts);
  fprintf(stderr, "Memcopies:       %lld\n", totalcpy);
  #endif
  
  //Clean up
  free(minitree);
  free(data);
  
  exit(EXIT_SUCCESS);
}
