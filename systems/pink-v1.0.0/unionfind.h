#ifndef __UNIONFIND_H
#define __UNIONFIND_H

#define PATH_COMPRESSION

typedef struct
{
  int length;
  int array[];
} UnionFind;

//Allocate union-find data structure for n elements
UnionFind* init_unionfind(int n);

//Attempt to union components l and r, where l < r
//Returns -1 if l and r already belong to same component
int unify(int l, int r, UnionFind* uf);

//(Re-)initializes union-find data structure
void reset_unionfind(UnionFind* uf);

void print_unionfind(UnionFind* uf);
#endif
