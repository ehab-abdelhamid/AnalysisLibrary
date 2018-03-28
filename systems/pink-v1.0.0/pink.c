/**
 pink.c
 
 Code for performing distributed single-linkage hierarchical agglomerative
 clustering.  This algorithm is described in "A Scalable Algorithm for 
 Single-Linkage Hierarchical Clustering on Distributed-Memory Architectures"
 by Hendrix, Palsetia, Patwary, Agrawal, Liao, and Choudhary (LDAV 2013).
 
 The main work of PINK is divided into three phases:
 
 PHASE 0:  Data is distributed among the processes and data structures are
           initialized
 PHASE 1:  Each process computes partial MST (dendrogram) for local data
           using Prim's algorithm
 PHASE 2:  Partial MSTs are merged in binary fashion
 
 Author:  William Hendrix, Northwestern University, 2013
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <mpi.h>

//cluster.h defines SLINK function and data structures
//Also defines macros controlling the printing of timing and stat reports
#include "cluster.h"
//unionfind.h defines UnionFind data structure and functions
#include "unionfind.h"

//Constants
#define EPSILON             1e-10
#define COMPLETE            1
#define BIPARTITE           0
#define NO_RIGHT_PARTITION -1

//Utility functions
MPI_Status status;
#define chmalloc(x, type, nelem, msg) if ((x = (type*) malloc((nelem) * sizeof(type))) == NULL) {perror(msg); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#define mpi_call(x) if ((x) != MPI_SUCCESS) {perror("MPI Error"); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
#define mpi_scall(x) mpi_call(x) if (status.MPI_ERROR) {perror("MPI error"); MPI_Abort(MPI_COMM_WORLD, status.MPI_ERROR);}

/**
 Helper function for reading in a binary file
 
 @param filename name of the file to read
 @param rank the rank of this process
 @param partitions the number of partitions into which to divide the file
 @param output_dims (output parameter) the dimensions (rows, cols) of the data
 @param output_left (output parameter) the number of rows in the left data partition for this process
 @param output_right (output parameter) the number of rows in the right data partition for this process
 @param output_complete (output parameter) the partition type for this process (COMPLETE, BIPARTITE, or NO_RIGHT_PARTITION)
 @param output_data (output parameter) the dataset for this process; the left-hand and right-hand partitions are concatenated
 @param output_rowid (output parameter) the IDs of the rows being read into the left- and right-hand data partitions
 */
void read_binary(char* filename, int rank, int partitions, int* output_dims, int* output_left, int* output_right, int* output_complete, data_t** output_data, int** output_rowid)
{
  int lpart, lpartsize, rpart, rpartsize, pts, rem, size, start;
  int block_length[2];
  int* rowid;
  data_t* data;
  MPI_Aint displacement[2];
  MPI_Datatype filetype;
  MPI_File fp;
  MPI_Offset size_check;
  
  //Calculate left and right partition IDs for local data
  rpart = (int) (sqrt(rank + rank + 1.5) + EPSILON);
  lpart = rank - (rpart * rpart >> 1);
  if (rpart == partitions)
  {
    lpart = partitions - 1;
    rpart = -1;
    *output_complete = NO_RIGHT_PARTITION;
  }
  else if (lpart == rpart)
  {
    lpart--;
    *output_complete = COMPLETE;
  }
  else
    *output_complete = BIPARTITE;

  mpi_call(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp));
  
  //Read in metadata
  mpi_scall(MPI_File_read_all(fp, output_dims, 2, MPI_INT, &status));
  pts = output_dims[0] / partitions;
  rem = output_dims[0] % partitions;
  size = pts * output_dims[1];
  
  //Verify file size
  mpi_call(MPI_File_seek(fp, 0, MPI_SEEK_END));
  mpi_call(MPI_File_get_position(fp, &size_check));
  if (rank == 0 && size_check != 2 * sizeof(int) + output_dims[0] * output_dims[1] * sizeof(data_t))
  {
    fprintf(stderr, "Error:  expected file size:  %d; found: %lld\n", 2 * sizeof(int) + output_dims[0] * output_dims[1] * sizeof(data_t), size_check);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  mpi_call(MPI_File_seek(fp, 8, MPI_SEEK_SET));
  
  //Allocate arrays
  chmalloc(data, data_t, (size + output_dims[1]) * 2, "Error allocating array for data matrix");
  chmalloc(rowid, int, (pts + 1) * 2, "Error allocating array for row ids");
  
  //Calculate offset, length, and IDs for left partition
  if (lpart < rem)
    start = lpart * (pts + 1);
  else
    start = lpart * pts + rem;
  displacement[0] = start * output_dims[1] * sizeof(data_t) + 2 * sizeof(int);
  for (lpartsize = 0; lpartsize < pts; lpartsize++)
    rowid[lpartsize] = start + lpartsize;
  if (lpart < rem)
  {
    rowid[lpartsize] = start + lpartsize;
    lpartsize++;
  }
  block_length[0] = lpartsize * output_dims[1];
  
  //Calculate offset, length, and IDs for right partition
  if (rpart >= 0)
  {
    if (rpart < rem)
      start = rpart * (pts + 1);
    else
      start = rpart * pts + rem;
    displacement[1] = start * output_dims[1] * sizeof(data_t) + 2 * sizeof(int);
    rowid += lpartsize;
    for (rpartsize = 0; rpartsize < pts; rpartsize++)
      rowid[rpartsize] = start + rpartsize;
    if (rpart < rem)
    {
      rowid[rpartsize] = start + rpartsize;
      rpartsize++;
    }
    rowid -= lpartsize;
    block_length[1] = rpartsize * output_dims[1];
    mpi_call(MPI_Type_create_hindexed(2, block_length, displacement, MPI_DATA, &filetype));
  }
  else
  {
    block_length[1] = rpartsize = 0;
    mpi_call(MPI_Type_create_hindexed(1, block_length, displacement, MPI_DATA, &filetype));
  }
  
  //Define MPI file view
  mpi_call(MPI_Type_commit(&filetype));
  mpi_call(MPI_File_set_view(fp, 0, MPI_BYTE, filetype, "native", MPI_INFO_NULL));
  
  //Read in left and right data partitions
  mpi_scall(MPI_File_read_all(fp, data, block_length[0] + block_length[1], MPI_DATA, &status));
  mpi_call(MPI_Type_free(&filetype));
  mpi_call(MPI_File_close(&fp));
  
  //Write data and row IDs to output arrays
  *output_left = lpartsize;
  *output_right = rpartsize;
  *output_data = data;
  *output_rowid = rowid;
}

/**
 * Driver routine for parallel PINK algorithm
 */
int main(int argc, char** argv)
{
  char* filename;
  data_t* data;
  int* rowid;
  Node* minitree;
  Node* recvtree;
  Node* mergetree;
  Node* temptree;
  int rank, nprocs, nparts, nleft, nright, nedges, redges, medges, isfull;
  int dim[2]; // = {nrows, ncols}
  int i, j, iter, delta;
  
  UnionFind* unionfind;
  MPI_Datatype MPI_NODE;
  
  #ifdef PRINT_STATS
  int recvcounter[14]; //14 receive phases: up to 16k processes
  long long edgecount, totaledges, totaldcounts, totalcpy;
  distcount = cpy = edgecount = totaledges = 0;
  #endif
  
  #ifdef PRINT_TIMING
  double starttime;
  double readtime = 0;
  double slinktime = 0;
  double mergetime = 0;
  double outputtime = 0;
  double maxreadtime, maxslinktime, maxmergetime;
  double avgreadtime, avgslinktime, avgmergetime;
  #endif
  
  //PHASE 0:  Initialization and data distribution
  //MPI initialization
  mpi_call(MPI_Init(&argc, &argv));
  #ifdef PRINT_TIMING
  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime();
  #endif
  mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
  
  //Check command line argument(s)
  if (argc < 2)
  {
    if (rank == 0)
      printf("Usage:  %s <data file>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  filename = argv[1];
  
  //Calculate number of partitions and active processes
  if (nprocs > 1)
  {
    nparts = (int) (sqrt(nprocs + nprocs) + EPSILON);
    nprocs = (nparts * nparts + 1) / 2;
  }
  else
    nparts = 1;
  
  //"Extra" processes go idle
  if (rank >= nprocs)
  {
    fprintf(stderr, "Warning:  Process %d not used\n", rank);
    
    //Needs to participate in MPI calls, even if not doing anything
    #ifdef PRINT_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    MPI_File fp;
    mpi_call(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp));
    mpi_scall(MPI_File_read(fp, dim, 2, MPI_INT, &status));
    mpi_call(MPI_File_set_view(fp, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL));
    mpi_scall(MPI_File_read_all(fp, data, 0, MPI_DATA, &status));
    mpi_call(MPI_File_close(&fp));
    
    #ifdef PRINT_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    readtime = 0;
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    #endif
    
    #ifdef PRINT_STATS
    edgecount = 0;
    mpi_call(MPI_Reduce(&edgecount, &totaledges, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&edgecount, &totaledges, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&edgecount, &totaledges, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    #endif
    
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }
  
  //Open file (MPI-I/O)
  read_binary(filename, rank, nparts, dim, &nleft, &nright, &isfull, &data, &rowid);
  
  //Allocate trees
  chmalloc(minitree, Node, dim[0], "Error allocating array for MST mergetree");
  chmalloc(recvtree, Node, dim[0], "Error allocating buffer for receiving MST");
  chmalloc(mergetree, Node, dim[0], "Error allocating buffer for merging MSTs");
  
  //PHASE 1:  Run Prim's algorithm on local data partitions independently
  #ifdef PRINT_TIMING
  readtime = MPI_Wtime() - starttime;
  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime();
  #endif
  
  nedges = nleft + nright - 1;
  if (isfull != BIPARTITE)
  {
    if (prim_full(nleft, dim[1], data, rowid, DISTANCE_FUNCTION, minitree))
    {
      perror("Error running Prim's on lefthand dataset");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (isfull != NO_RIGHT_PARTITION && prim_full(nright, dim[1], data + nleft * dim[1], rowid + nleft, DISTANCE_FUNCTION, minitree + nleft - 1))
    {
      perror("Error running Prim's on righthand dataset");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (isfull != NO_RIGHT_PARTITION)
      nedges--;
  }
  else if (prim_bipartite(nleft, nright, dim[1], data, rowid, DISTANCE_FUNCTION, minitree))
  {
    perror("Error running Prim's on bi-complete dataset");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  
  //Sort MST edges
  qsort(minitree, nedges, sizeof(Node), nodecompare);
  
  //PHASE 2:  Combine partial dendrograms to form full dendrogram (MST)
  #ifdef PRINT_TIMING
  slinktime = MPI_Wtime() - starttime;
  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime();
  #endif
  
  //Create MPI_Datatype for communicating partial dendrograms (MSTs)
  init_type_MPI_NODE();
  
  //Initialize Union-Find (Disjoint Set) data structure
  unionfind = init_unionfind(dim[0]);
  
  //Code for receiving and merging partial dendrograms in binary fashion
  for (iter = 0, delta = 1; (rank & delta) == 0 & delta < nprocs; iter++, delta <<= 1)
  {
    if (rank + delta >= nprocs)
    {
      //delta = 1 << (ffs(rank) - 1);
      while ((rank & delta) == 0)
        delta <<= 1;
      break;
    }
    
    //Receive partial dendrogram (MST) from processor (rank ^ delta = rank + delta)
    mpi_scall(MPI_Recv(&redges, 1, MPI_INT, rank + delta, delta, MPI_COMM_WORLD, &status));
    mpi_scall(MPI_Recv(recvtree, redges, MPI_NODE, rank + delta, delta, MPI_COMM_WORLD, &status));
    
    #ifdef PRINT_STATS
    recvcounter[iter] = redges;
    edgecount += redges;
    #endif
    
    //Reset variables
    reset_unionfind(unionfind);
    i = 0;
    j = 0;
    medges = 0;
    
    //Merge partial dendrograms (MSTs)
    while (1)
    {
      //Iteratively add min height merge from either tree
      if (minitree[i].distance < recvtree[j].distance)
      {
        //Local dendrogram has min height merge

        //Check whether merge vertices are already in same cluster
        if (unify(minitree[i].left, minitree[i].right, unionfind) == 0)
          mergetree[medges++] = minitree[i];
        
        //If all edges from local dendrogram have been added, iterate through received dendrogram
        if (++i >= nedges)
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
          while (i < nedges)
          {
            if (unify(minitree[i].left, minitree[i].right, unionfind) == 0)
              mergetree[medges++] = minitree[i];
            
            i++;
          }
          break;
        }
      }
    }
    
    //Swap merged dendrogram with local dendrogram
    nedges = medges;
    temptree = mergetree;
    mergetree = minitree;
    minitree = temptree;
  }
  
  //Send local partial dendrogram to another process
  if (rank != 0)
  {
    mpi_call(MPI_Send(&nedges, 1, MPI_INT, rank - delta, delta, MPI_COMM_WORLD));
    mpi_call(MPI_Send(minitree, nedges, MPI_NODE, rank - delta, delta, MPI_COMM_WORLD));
  }
  #ifdef PRINT_TIMING
  mergetime += MPI_Wtime() - starttime;
  starttime = MPI_Wtime();
  #endif
  //ALGORITHM COMPLETE
  
  //Output result, timing, and stats
  if (rank == 0)
  {
    #ifdef PRINT_MST
    //Write out final MST
    for (i = 0; i < nedges; i++)
      printf(PRINT_UVDIST, minitree[i].left, minitree[i].right, minitree[i].distance);
    #endif
    
  //Collect stats/timing
    #ifdef PRINT_TIMING
    outputtime = MPI_Wtime() - starttime;
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &avgreadtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgreadtime /= nprocs;
    mpi_call(MPI_Reduce(&slinktime, &maxslinktime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&slinktime, &avgslinktime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgslinktime /= nprocs;
    mpi_call(MPI_Reduce(&mergetime, &maxmergetime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&mergetime, &avgmergetime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgmergetime /= nprocs;
    fprintf(stderr, "File:                %s\n", argv[1]);
    fprintf(stderr, "Processes:           %d\n", nprocs);
    fprintf(stderr, "Max_reading:         %lf s\n", maxreadtime);
    fprintf(stderr, "Max_Prims:           %lf s\n", maxslinktime);
    fprintf(stderr, "Max_Merging:         %lf s\n", maxmergetime);
    fprintf(stderr, "Avg_reading:         %lf s\n", avgreadtime);
    fprintf(stderr, "Avg_Prims:           %lf s\n", avgslinktime);
    fprintf(stderr, "Avg_merging:         %lf s\n", avgmergetime);
    fprintf(stderr, "Output:              %lf s\n", outputtime);
    fprintf(stderr, "Total (p0):          %lf s\n", readtime + slinktime + mergetime + outputtime);
    #endif

    #ifdef PRINT_STATS
    mpi_call(MPI_Reduce(&edgecount, &totaledges, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Allreduce(&distcount, &totaldcounts, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&distcount, &totaldcounts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&cpy, &totalcpy, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    fprintf(stderr, "Edges_passed:        %lld\n", totaledges);
    //fprintf(stderr, "Distance_calcs:      %lld\n", totaldcounts);
    //fprintf(stderr, "Memcopies:           %lld\n", totalcpy);
    for (i = 0; i < iter; i++)
      fprintf(stderr, "Comm_%d:              %d\n", i, recvcounter[i]);
    #endif
  }
  else
  {
    //Collect stats/timing
    #ifdef PRINT_TIMING
    mpi_call(MPI_Reduce(&readtime, &maxreadtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&readtime, &avgreadtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgreadtime /= nprocs;
    mpi_call(MPI_Reduce(&slinktime, &maxslinktime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&slinktime, &avgslinktime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgslinktime /= nprocs;
    mpi_call(MPI_Reduce(&mergetime, &maxmergetime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    mpi_call(MPI_Reduce(&mergetime, &avgmergetime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    avgmergetime /= nprocs;
    #endif
    
    #ifdef PRINT_STATS
    mpi_call(MPI_Reduce(&edgecount, &totaledges, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Allreduce(&distcount, &totaldcounts, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&distcount, &totaldcounts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    //mpi_call(MPI_Reduce(&cpy, &totalcpy, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    #endif
  }
  
  //Clean up
  free_type_MPI_NODE();

  free(unionfind);
  free(minitree);
  free(recvtree);
  free(mergetree);
  
  free(rowid);
  free(data);
  
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
