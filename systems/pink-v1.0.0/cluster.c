#include "cluster.h"
#include <stdio.h>

#include <mpi.h>

#ifdef PRINT_STATS
long long distcount;
long long cpy;
#endif

/**
 Helper function for qsort to sort Node structs
 
 @param a pointer to first Node
 @param b pointer to second Node
 @return value > 0 if a > b, value < 0 if a < b, and 0 if a == b
 */
int nodecompare(const void* a, const void* b)
{
  const Node* node1 = (const Node*) a;
  const Node* node2 = (const Node*) b;
  const dist_t term1 = node1->distance;
  const dist_t term2 = node2->distance;
  if (term1 < term2)
    return -1;
  if (term1 > term2)
   return 1;
  return 0;
}

/**
 Function for calculating the square of the Euclidean distance between two
 points with len dimensions.  Aborts the distance computation early if it
 exceeds a given threshold.
 
 @param x    first vector
 @param y    second vector
 @param len  dimensionality for x and y
 @param max  max distance threshold; if d(x, y) > max, the value returned may be less than the true distance but > max
 
 @return squared Euclidean distance between x and y
 */
dist_t sq_euclid_dist(data_t* x, data_t* y, int len, dist_t max)
{
  int i;
  dist_t dist = (dist_t) 0;
  dist_t t;
  for (i = 0; i < len && dist < max; i++)
  {
    t = x[i] - y[i];
    dist += t * t;
  }
  
  return dist;
}

/**
 Function for calculating distance as 1 - cos(x, y), where x and y are
 assumed to be normalized vectors of length 1.  Aborts the distance 
 computation early if it exceeds a given threshold.
 
 @param x    first vector
 @param y    second vector
 @param len  dimensionality of x and y
 @param max  max distance threshold; if d(x, y) > max, the value returned may be less than the true distance but > max
 
 @return 1 - cos(x, y)
 */
dist_t norm_cos_dist(data_t* x, data_t* y, int len, dist_t max)
{
  int i;
  dist_t dist = (dist_t) 0;
  max = (dist_t) 1 - max;
  for (i = 0; i < len && dist < max; i++)
    dist += x[i] * y[i];
  
  //Check that 0 <= dist <= 1??
  
  return (dist_t) 1 - dist;
}

//Define constant values for orbital distance function (sq_orb_dist)
const dist_t ka = 1.25f;
const dist_t ke = 2.0f;
const dist_t ki = 2.0f;
const dist_t kw = 1e-6f;
const dist_t kf = 1e-6f;

/**
 Function for calculating dissimilarity between two orbits
 Data points are assumed to have the following 6 dimensions:
 
 n:  Mean daily motion
 a:  Semimajor axis
 e:  Eccentricity
 i:  Inclination
 w:  Argument of perihelion
 f:  Longitude of ascending node
 
 This function is described in more detail in:
 
 D. Nesvorny and D. Vokrouhlicky; New Candidates for Recent Asteroid
   Breakups; The Astronomical Journal, 132:1950-1958, November 2006.
 
 V. Zappala, A. Cellino, P. Farinella, and Z. Knezevic; Asteroid Families.
   I. Identification by Hierarchical Clustering and Reliability Assessment;
   The Astronomical Journal, 100:2030-2046, December 1990.
 
 @param x    first vector (length >= 6)
 @param y    second vector (length >= 6)
 @param len  (not used)
 @param max  (not used)
 @return the orbital distance between x and y
 */
dist_t sq_orb_dist(data_t* x, data_t* y, int len, dist_t max)
{
  //Define variables
  const dist_t n = (x[0] + y[0]) / 2;
  const dist_t a = (x[1] + y[1]) / 2;
  const dist_t da = x[1] - y[1];
  const dist_t de = x[2] - y[2];
  const dist_t di = sin(x[3]) - sin(y[3]);
  const dist_t dw = x[4] - y[4];
  const dist_t df = x[5] - y[5];
  
  //Calculate distance value
  dist_t retval = n * n * a * a * (ka * da / a * da / a + ke * de * de + ki * di * di + kw * dw * dw + kf * df * df);
    return retval;
}

/*
 Function for computing the MST of the complete weighted graph induced by the given data elements and distance function.  Builds the edge list of the MST using the vector of vertex IDs passed in, as these IDs will be critical when merging MSTs.  Difference from prim_full_shrink:  the data matrix is assumed to only contain the nrows X ncols data subselected from the full matrix, and row IDs are only used to determine the output MST.
 
 @param nrows     number of rows (records) in data matrix
 @param ncolumns  number of columns (dimensions) in data matrix
 @param data      (nrows X ncolumns) data matrix
 @param id        vector of row IDs
 @param distance  distance function for two records
 @param result    (output parameter) MST edge list, in order of increasing distance (assumed to be allocated already)
 
 @return 0 (success)
*/
int prim_full(int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result)
{
  int i, k, t;
  int size, minloc, minvert, minpar;
  dist_t min, dist;
  
  size = nrows - 1;
  for (i = 0; i < size; i++)
  {
    result[i].distance = DIST_MAX;
    result[i].left = i + 1;
  }
  
  minvert = 0;
  while (size > 0)
  {
    k = minvert;
    min = DIST_MAX;
    for (i = 0; i < size; i++)
    {
      t = result[i].left;
      dist = distance(&data[k * ncolumns], &data[t * ncolumns], ncolumns, result[i].distance);
      if (result[i].distance > dist)
      {
        result[i].right = k;
        result[i].distance = dist;
      }
      if (result[i].distance < min)
      {
        min = result[i].distance;
        minvert = t;
        minpar = result[i].right;
        minloc = i;
      }
    }
    size--;
    result[minloc] = result[size];
    result[size].distance = min;
    if (minvert < minpar) //Assumption:  id is monotonic
    {
      result[size].left = id[minvert];
      result[size].right = id[minpar];
    }
    else
    {
      result[size].left = id[minpar];
      result[size].right = id[minvert];
    }
  }
  
  return 0;
}

/*
 Function for computing the MST of the complete weighted bipartite graph induced by two data partitions and a given distance function.  The input matrix should be prepared so that the two data partitions are contiguous in memory, as a matrix with (left + right) rows and ncols columns.  Builds the edge list of the MST using the vector of vertex IDs passed in, as these IDs will be critical when merging MSTs.
 
 @param leftrows  number of rows (records) in left-hand data partition
 @param rightrows number of rows (records) in right-hand data partition
 @param ncolumns  number of columns (dimensions) in data matrix
 @param data      (leftrows + rightrows) X ncolumns data matrix
 @param id        vector of row IDs
 @param distance  distance function for two records
 @param result    (output parameter) MST edge list, in order of increasing distance (assumed to be allocated already)
 
 @return 0 (success)
*/
int prim_bipartite(int leftrows, int rightrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result)
{
  const int nrows = leftrows + rightrows;
  
  Node* left = result;
  Node* right = result + leftrows - 1;
  int lsize = leftrows - 1;
  int rsize = rightrows;
  dist_t lmin, rmin;
  int lidx, ridx, minvert, minpar;
  
  int i;
  dist_t d;
  data_t* sel;
  
  lmin = rmin = DIST_MAX;
  for (i = 0; i < lsize; i++)
  {
    left[i].distance = DIST_MAX;
    left[i].left = i + 1;
  }
  for (i = 0; i < rsize; i++)
  {
    d = distance(data, &data[(i + leftrows)*ncolumns], ncolumns, DIST_MAX);
    if (d < rmin)
    {
      rmin = d;
      ridx = i;
    }
    right[i].distance = d;
    right[i].left = i + leftrows;
    right[i].right = 0;
  }
  
  while (lsize > 0 && rsize > 0)
  {
    //Find min(left, right)
    if (lmin < rmin)
    {
      minvert = left[lidx].left;
      minpar = left[lidx].right;
      sel = &data[minvert * ncolumns];
      
      //Update distances for right side, calculate min
      rmin = DIST_MAX;
      for (i = 0; i < rsize; i++)
      {
        d = distance(sel, &data[right[i].left * ncolumns], ncolumns, right[i].distance);
        if (d < right[i].distance)
        {
          right[i].distance = d;
          right[i].right = minvert;
        }
        if (right[i].distance < rmin)
        {
          rmin = right[i].distance;
          ridx = i;
        }
      }
      
      //Swap min edge to end
      lsize--;
      left[lidx] = left[lsize];
      left[lsize].distance = lmin;
      if (minvert < minpar) //Assumption:  id is monotonic
      {
        left[lsize].left = id[minvert];
        left[lsize].right = id[minpar];
      }
      else
      {
        left[lsize].left = id[minpar];
        left[lsize].right = id[minvert];
      }
      
      //Recompute min for left side
      lmin = DIST_MAX;
      for (i = 0; i < lsize; i++)
        if (left[i].distance < lmin)
        {
          lmin = left[i].distance;
          lidx = i;
        }
    }
    else
    {
      minvert = right[ridx].left;
      minpar = right[ridx].right;
      sel = &data[right[ridx].left * ncolumns];
      
      //Update distances for left side, calculate min
      lmin = DIST_MAX;
      for (i = 0; i < lsize; i++)
      {
        d = distance(sel, &data[left[i].left * ncolumns], ncolumns, left[i].distance);
        if (d < left[i].distance)
        {
          left[i].distance = d;
          left[i].right = minvert;
        }
        if (left[i].distance < lmin)
        {
          lmin = left[i].distance;
          lidx = i;
        }
      }
      
      //Swap min edge to end
      rsize--;
      right[ridx] = right[rsize];
      right[rsize].distance = rmin;
      if (minvert < minpar) //Assumption:  id is monotonic
      {
        right[rsize].left = id[minvert];
        right[rsize].right = id[minpar];
      }
      else
      {
        right[rsize].left = id[minpar];
        right[rsize].right = id[minvert];
      }
      
      //Recompute min for right side
      rmin = DIST_MAX;
      for (i = 0; i < rsize; i++)
        if (right[i].distance < rmin)
        {
          rmin = right[i].distance;
          ridx = i;
        }
    }
  }
  
  //Correct remaining edges in left and right arrays
  for (i = 0; i < lsize; i++)
  {
    lidx = left[i].left;
    ridx = left[i].right;
    if (lidx < ridx) //Assumption:  id is monotonic
    {
      left[i].left = id[lidx];
      left[i].right = id[ridx];
    }
    else
    {
      left[i].left = id[ridx];
      left[i].right = id[lidx];
    }
  }
  for (i = 0; i < rsize; i++)
  {
    lidx = right[i].left;
    ridx = right[i].right;
    if (lidx < ridx) //Assumption:  id is monotonic
    {
      right[i].left = id[lidx];
      right[i].right = id[ridx];
    }
    else
    {
      right[i].left = id[ridx];
      right[i].right = id[lidx];
    }
  }
  
  return 0;
}

/*
 Function for computing the MST of the complete weighted graph induced by the data elements in the given rows using a given distance function.  Difference from prim_full:  the data matrix is assumed to contain the entire dataset, but prim_full_shrink only accesses those rows (data points) contained the id vector.
 
 @param nrows     number of rows (records) in data matrix
 @param ncolumns  number of columns (dimensions) in data matrix
 @param data      data matrix
 @param id        vector of row IDs
 @param distance  distance function for two records
 @param result    (output parameter) MST edge list, in order of increasing distance (assumed to be allocated already)
 
 @return 0 (success)
*/
int prim_full_shrink(int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result)
{
  int i, k, t;
  int size, minloc, minvert, minpar;
  dist_t min, dist;
  
  size = nrows - 1;
  for (i = 0; i < size; i++)
  {
    result[i].distance = DIST_MAX;
    result[i].left = id[i + 1];
  }
  
  minvert = id[0];
  while (size > 0)
  {
    k = minvert;
    min = DIST_MAX;
    for (i = 0; i < size; i++)
    {
      t = result[i].left;
      dist = distance(&data[k * ncolumns], &data[t * ncolumns], ncolumns, result[i].distance);
      if (result[i].distance > dist)
      {
        result[i].right = k;
        result[i].distance = dist;
      }
      if (result[i].distance < min)
      {
        min = result[i].distance;
        minvert = t;
        minpar = result[i].right;
        minloc = i;
      }
    }
    size--;
    result[minloc] = result[size];
    result[size].distance = min;
    if (minvert < minpar)
    {
      result[size].left = minvert;
      result[size].right = minpar;
    }
    else
    {
      result[size].left = minpar;
      result[size].right = minvert;
    }
  }
  
  return 0;
}


#define lt(x, y) (x.distance < y.distance || x.distance == y.distance && (x.left < y.left || x.left == y.left && x.right < y.right))
#define gt(x, y) (x.distance > y.distance || x.distance == y.distance && (x.left > y.left || x.left == y.left && x.right > y.right))

int pslcluster (int nrows, int ncolumns, data_t* data, int* id, dist_t (*distance)(data_t*, data_t*, int, dist_t), Node* result)
 /*******************************************************************************
  *
  * Purpose
  * =======
  *
  * The pslcluster routine performs single-linkage hierarchical clustering, using
  * either the distance matrix directly, if available, or by calculating the
  * distances from the data array. This implementation is based on the SLINK
  * algorithm, described in:
  *
  * Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link
  * cluster method. The Computer Journal, 16(1): 30-34.
  *
  * The output of this algorithm is identical to conventional single-linkage
  * hierarchical clustering, but is much more memory-efficient and faster. Hence,
  * it can be applied to large data sets, for which the conventional single-
  * linkage algorithm fails due to lack of memory.
  *
  * Authors
  * =======
  * Zhihua Du and Feng Lin ("A novel parallelization approach for hierarchical
  *   clustering") 
  * Modifications to interface with parallel algortihm added by William Hendrix
  * WARNING:  this algorithm is **not** guaranteed to produce the correct MST, 
  *           just a subgraph with an equal weight to the MST
  *
  * Arguments
  * =========
  *
  * nrows     (input) int
  * The number of rows in the data matrix, equal to the number of points
  *
  * ncolumns  (input) int
  * The number of columns in the data matrix, equal to the number of dim of
  * each point.
  *
  * data      (input) data_t[nrows * ncolumns]
  * The array containing the dataset.
  *
  * id        (input) int[nrows]
  * Array of ID values for each row of the data (modification for parallel SLINK)
  *
  * distance  (input) function:  (data_t*, data_t*, int) -> dist_t
  * Function for calculating distance between two data points (modification)
  *
  * result   (output) Node[>= nrows - 1]
  * Output array describing the Minimum Spanning Tree of the data.
  * Each Node will be of the form (ID_1, ID_2, dist), where dist is the 
  * distance between the points with IDs ID_1 and ID_2.  See cluster.h for a 
  * description of the Node structure. (modified for parallel SLINK)
  *
  *
  * Return value
  * ============
  *
  * pslcluster returns -1 if there was a memory allocation error and 0 otherwise
  */
{
  int i, j, m, n, k;
  const int nNodes = nrows - 1;
  int* vector;
  int* left;
  int* right;
  //Modification for parallel SLINK:  maintain endpoints as well as distance for merges
  Node* temp;
  int* index;
  int prank;
  int psize;

  temp = (Node*) malloc(nNodes * sizeof(Node));
  if (!temp)
    return -1;
  index = (int*) malloc(nrows * sizeof(int));
  if (index == NULL)
  {
    free(temp);
    return -1;
  }
  vector = (int*) malloc(nNodes * sizeof(int));
  if (vector == NULL)
  {
    free(index);
    free(temp);
    return -1;
  }
  
  for (i = 0; i < nNodes; i++)
    vector[i] = i;
  
  //for (i = 0; i < nrows; i++)
  //  result[i].distance = DATA_MAX;

  for (i = 0; i < nrows; i++)
  {
    //Modification for parallel SLINK:  Calculate distance rather than using distance matrix
    result[i].distance = DATA_MAX;
    for (j = 0; j < i; j++)
    {
      #ifdef PRINT_STATS
      distcount++;
      #endif
      temp[j].distance = distance(&data[i * ncolumns], &data[j * ncolumns], ncolumns, DIST_MAX);
      temp[j].left = id[j];
      temp[j].right = id[i];
    }
    
    for (j = 0; j < i; j++)
    {
      k = vector[j];
      if (gt(result[j], temp[j]))
      {
        if (lt(result[j], temp[k]))
        {
          temp[k] = result[j];
          #ifdef PRINT_STATS
          //Modification for parallel SLINK:  maintain a count for "memory copies"
          cpy++;
          #endif
        }
        result[j] = temp[j];
        #ifdef PRINT_STATS
        cpy++;
        #endif
        vector[j] = i;
      }
      else if (lt(temp[j], temp[k]))
      {
        temp[k] = temp[j];
        #ifdef PRINT_STATS
        cpy++;
        #endif
      }
    }
    for (j = 0; j < i; j++)
      if (gt(result[j], result[vector[j]]))
        vector[j] = i;
  }
 
  free(temp);
  free(index);
  free(vector);
  
  qsort(result, nNodes, sizeof(Node), nodecompare);
  
  //Modification for parallel SLINK:  terminate SLINK after sorting
  return 0;
  
  /*
  for (i = 0; i < nrows; i++)
    index[i] = i;
  for (i = 0; i < nNodes; i++)
  {   
    j = result[i].left;
    k = vector[j];
    result[i].left = index[j];
    result[i].right = index[k];
    index[k] = -i - 1;
  }
  if (prank == 0)
  {
    for (i = 0; i < nrows; i++) 
      printf("%d: index[%d]=%d, result[%d].left:%d, result[%d].right:%d\n", prank, i, index[i], i, result[i].left, i, result[i].right);
    printf("\n");
  }
  
  free(vector);
  free(index);

  result = realloc(result, nNodes * sizeof(Node));
  return(result);
  //*/
}
