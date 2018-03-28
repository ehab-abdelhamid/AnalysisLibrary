#include <stdlib.h>
#include <stdio.h>

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    fprintf(stderr, "Usage:  %s [h (height) | k (clusters)] [threshold value] [input hierarchical clustering file] [output cluster file]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int k, kthresh, lines, readlines, u, v, numpoints, uc, vc, i, t, temp;
  double hthresh, h;
  FILE* infile, * outfile;
  int* id, * cluster, * uf, * size, * data;

  //Parse threshold type and value
  if ((argv[1][0] | 32) == 'h')
  {
    hthresh = atof(argv[2]);
    k = 0;
  }
  else if ((argv[1][0] | 32) == 'k')
  {
    kthresh = atoi(argv[2]);
    k = 1;
    if (kthresh < 1)
    {
      fprintf(stderr, "Argument 2:  Number of clusters must be at least one\n");
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "Argument 1:  threshold must be applied to dendrogram height (h) or or number of clusters (k)\n");
    exit(EXIT_FAILURE);
  }

  //Open input HC and output cluster files
  outfile = fopen(argv[4], "w");
  infile = fopen(argv[3], "r");
  if (infile == NULL)
  {
    perror("Error opening input HC file");
    exit(EXIT_FAILURE); 
  }
  if (outfile == NULL)
  {
    perror("Error opening output cluster file");
    exit(EXIT_FAILURE);
  }

  //Scan HC file to determine number of lines, points, and either clusters (k) or lines to read (h)
  readlines = numpoints = -1;
  lines = 0;
  while (fscanf(infile, "(%d, %d, %lf) ", &u, &v, &h) == 3)
  {
    if (u > numpoints) numpoints = u;
    if (v > numpoints) numpoints = v;
    if (k == 0 && readlines == -1 && h > hthresh)
      readlines = lines;
    lines++;
  }
  //Add vertex "0" to numpoints
  numpoints++;

  //Sanity check:  compare number of merges with number of data points
  if (numpoints > lines + 1)
    fprintf(stderr, "Warning:  dendrogram incomplete--%d lines in file for %d points\n", lines, numpoints);
  else if (numpoints < lines + 1)
    fprintf(stderr, "Warning:  dendrogram has too many merges--%d lines in file for %d points\n", lines, numpoints);
  
  //Calculate number of lines to read from file (k) or number of clusters (h)
  if (k != 0)
    readlines = numpoints - kthresh;
  else
    kthresh = numpoints - readlines;
  
  //Output something so user doesn't get nervous
  printf("%d point(s), %d line(s) in file, %d merge(s) to perform, %d cluster(s)\n", numpoints, lines, readlines, kthresh);
  
  //Sanity check:  check number of lines to read
  if (readlines < 0)
  {
    fprintf(stderr, "Error:  Number of clusters exceeds number of points in file\n");
    exit(EXIT_FAILURE);
  }
  else if (readlines > lines)
  {
    fprintf(stderr, "Error:  Too few merges in file to form %d clusters\n", kthresh);
    exit(EXIT_FAILURE);
  }
  
  //Allocate arrays
  cluster = (int*) malloc(kthresh * sizeof(int));
  id = (int*) malloc(numpoints * sizeof(int));
  uf = (int*) malloc(numpoints * sizeof(int));
  size = (int*) malloc(numpoints * sizeof(int));
  data = (int*) malloc(numpoints * sizeof(int));
  if (cluster == NULL || id == NULL || uf == NULL || size == NULL || data == NULL)
  {
    perror("Error allocating arrays for calculating clusters");
    exit(EXIT_FAILURE);
  }
  
  //Initialize disjoint set and set size arrays
  for (i = 0; i < numpoints; i++)
  {
    uf[i] = i;
    size[i] = 1;
  }
  
  //Perform readlines merges on disjoint set data structure
  rewind(infile);
  for (i = 0; i < readlines; i++)
  {
    fscanf(infile, "(%d, %d, %*lf) ", &u, &v);

    uc = u;
    vc = v;
    while (uc > uf[uc])
      uc = uf[uc];
    while (vc > uf[vc])
      vc = uf[vc];

    temp = u;
    while (temp > uc)
    {
      t = uf[temp];
      uf[temp] = uc;
      temp = t;
    }
    temp = v;
    while (temp > vc)
    {
      t = uf[temp];
      uf[temp] = vc;
      temp = t;
    }
    //printf("Merge (%d, %d) -> (%d, %d)\n", u, v, uc, vc);

    if (uc > vc)
    {
      uf[uc] = vc;
      size[vc] += size[uc];
    }
    else if (vc > uc)
    {
      uf[vc] = uc;
      size[uc] += size[vc];
    }
    else
    {
      fprintf(stderr, "Error in dendrogram (line %d):  Points %d and %d are already in same cluster!\n", i, u, v);
      exit(EXIT_FAILURE);
    }
  }
  fclose(infile);
  
  //Calculate cluster IDs
  for (i = 0, t = 0; i < numpoints; i++)
    if (uf[i] == i)
    {
      //printf("Cluster %d:  id %d, size %d\n", t, i, size[i]);
      id[i] = t++;
    }
  
  //Reorder cluster sizes by ID
  for (i = 0; i < numpoints; i++)
    if (uf[i] == i && id[i] != kthresh - 1)
      size[id[i] + 1] = size[i];
  
  //Calculate cluster offsets
  size[0] = 0;
  for (i = 1; i <= kthresh; i++)
    size[i] += size[i - 1];
  
  //Save clusters into data array
  for (i = 0; i < numpoints; i++)
  {
    t = uf[i];
    while (t < uf[t])
      t = uf[t];
    
    //printf("Point %d, cluster %d (id %d) -> pos %d\n", i, id[t], t, size[id[t]]);
    data[size[id[t]]++] = i;
  }
  
  //Output clusters to file
  fprintf(outfile, "Input file:     %s\n", argv[3]);
  fprintf(outfile, "# of clusters:  %d", kthresh);
  if (k == 0)
    fprintf(outfile, " (height threshold = %s)", argv[2]);
  fprintf(outfile, "\nCluster sizes:  %d", size[0]);
  for (i = 1; i < kthresh; i++)
    fprintf(outfile, " %d", size[i] - size[i - 1]);
  
  fprintf(outfile, "\nCluster membership:\n");
  for (i = 0, t = 0; i < numpoints; i++)
  {
    fprintf(outfile, "%d ", data[i]);
    if (i + 1 == size[t])
    {
      fprintf(outfile, "\n");
      t++;
    }
  }
  fclose(outfile);
  
  //Clean up
  free(cluster);
  free(id);
  free(uf);
  free(size);
  free(data);
  
  exit(EXIT_SUCCESS);
}
