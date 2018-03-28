/**
 unionfind.c
 
 Source code supporting Union-Find (Disjoint Set) data structure and relevant
 operations.
 
 Author: William Hendrix
 */

#include <stdlib.h>
#include <stdio.h>

#include "unionfind.h"

/**
 * Allocates and initializes new Union-Find (Disjoint Set) data structure
 * with n disjoint sets.
 */
UnionFind* init_unionfind(int n)
{
  int i;
  UnionFind* uf = (UnionFind*) malloc((n + 1) * sizeof(int));
  if (uf == NULL)
    return NULL;
  
  uf->length = n;
  return uf;
}

/**
 * Main utility function for Union-Find (Disjoint Set) data structure.
 * Performs a union on sets containing elements l and r in data structure uf.
 *
 * Returns 0 if sets were distinct, -1 if sets were already joined
 */
int unify(int l, int r, UnionFind* uf)
{
  int lc = l;
  int rc = r;
  int t;
  
  //lc = Find(l)
  while (lc < uf->array[lc])
    lc = uf->array[lc];
  
  //rc = Find(r)
  while (rc < uf->array[rc])
    rc = uf->array[rc];
  
  //If Find(l) != Find(r), perform union and possibly path compression
  if (lc < rc)
  {
    #ifdef PATH_COMPRESSION
    while (uf->array[l] != l)
    {
      t = l;
      l = uf->array[l];
      uf->array[t] = rc;
    }
    while (uf->array[r] != r)
    {
      t = r;
      r = uf->array[r];
      uf->array[t] = rc;
    }
    #endif
    
    uf->array[l] = rc;
    return 0;
  }
  else if (rc < lc)
  {
    #ifdef PATH_COMPRESSION
    while (uf->array[l] != l)
    {
      t = l;
      l = uf->array[l];
      uf->array[t] = lc;
    }
    while (uf->array[r] != r)
    {
      t = r;
      r = uf->array[r];
      uf->array[t] = lc;
    }
    #endif
    
    uf->array[r] = lc;
    return 0;
  }
  
  //If Find(l) == Find(r), sets were already joined
  #ifdef PATH_COMPRESSION
  while (uf->array[l] != l)
  {
    t = l;
    l = uf->array[l];
    uf->array[t] = lc;
  }
  while (uf->array[r] != r)
  {
    t = r;
    r = uf->array[r];
    uf->array[t] = rc;
  }
  #endif
  
  return -1;
}

/**
 * Initializes previously allocated Union-Find data structure to 
 * have n disjoint sets.
 */
void reset_unionfind(UnionFind* uf)
{
  int i;
  
  for (i = 0; i < uf->length; i++)
    uf->array[i] = i;
}

/**
 * Debug function:  prints contents of Union-Find data structure to stdout
 */
void print_unionfind(UnionFind* uf)
{
  int i;
  
  printf("%d", uf->array[0]);
  for (i = 1; i < uf->length; i++)
    printf(" %d", uf->array[i]);
  printf("\n");
}
