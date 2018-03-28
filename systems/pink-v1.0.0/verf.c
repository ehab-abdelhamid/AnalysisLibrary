#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "cluster.h"
#include "unionfind.h"

int main(int argc, char** argv)
{
  int i, fp, nrow, ncol, left, right, lineno;
  size_t len;
  dist_t dist, tempdist;
  double total;
  const dist_t EPS = 1e-7;
  data_t* data;
  UnionFind* uf;
  FILE* in;
  char* linebuf;
  linebuf = (char*) malloc(128);
  
  if (argc != 3)
  {
    printf("Usage:  %s <PINK/SHRINK output> <data file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  fp = open(argv[2], O_RDONLY);
  if (fp == -1)
  {
    perror("Error opening data file");
    exit(EXIT_FAILURE);
  }
  if (read(fp, &nrow, sizeof(int)) != sizeof(int) || read(fp, &ncol, sizeof(int)) != sizeof(int))
  {
    perror("Error reading metadata from data file");
    exit(EXIT_FAILURE);
  }
  data = (data_t*) malloc(nrow * ncol * sizeof(data_t));
  if (data == NULL)
  {
    perror("Error allocating data array");
    exit(EXIT_FAILURE);
  }
  if (read(fp, data, nrow * ncol * sizeof(data_t)) != nrow * ncol * sizeof(data_t))
  {
    perror("Error reading data from data file");
    exit(EXIT_FAILURE);
  }
  close(fp);
  uf = init_unionfind(nrow);
  reset_unionfind(uf);
  
  in = fopen(argv[1], "r");
  lineno = 0;
  total = 0;
  while (getline(&linebuf, &len, in) != -1)
  {
    lineno++;
    if (linebuf[0] != '(')
      continue;
    
    if (sscanf(linebuf, PRINT_UVDIST, &left, &right, &dist) != 3)
    {
      fprintf(stderr, "Error reading values from line %d\n", lineno);
      perror("");
      exit(EXIT_FAILURE);
    }
    total += dist;
    
    //printf("Joining %d to %d at dist %lf\n", left, right, dist);
    if (unify(left, right, uf) != 0)
      printf("Mistake:  %d and %d are already joined (line %d)\n", left, right, lineno);
    
    if (abs((tempdist = sq_euclid_dist(&data[left * ncol], &data[right * ncol], ncol, DIST_MAX)) - dist) > EPS)
      printf("Mistake:  distance between %d and %d is %lf, not %lf (line %d)\n", left, right, tempdist, dist, lineno);
  }
  
  left = 0;
  while (left < uf->array[left])
    left = uf->array[left];

  for (i = 0; i < nrow; i++)
  {
    right = i;
    while (right < uf->array[right])
      right = uf->array[right];
    if (right != left)
      printf("Mistake:  %d is not connected to 0 (%d vs. %d)\n", i, right, left);
  }
  
  printf("Total weight:  %lf\n", total);
  fclose(in);
  free(data);
  free(uf);
  free(linebuf);
  
  exit(EXIT_SUCCESS);
}
