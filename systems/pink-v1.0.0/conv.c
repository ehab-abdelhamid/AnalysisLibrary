#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "cluster.h"

int main(int argc, char** argv)
{
  int dim[2];
  data_t* data;
  int i, j;
  
  FILE* infile;
  MPI_File outfile;

  MPI_Init(&argc, &argv);
  if (argc < 3)
  {
    printf("Usage:  %s [ASCII input file] [binary output file]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  infile = fopen(argv[1], "r");
  if (infile == NULL)
  {
    perror("Error opening file");
    exit(EXIT_FAILURE);
  }
  if (fscanf(infile, "%d %d", &dim[0], &dim[1]) != 2)
  {
    printf("Error reading row and column counts from file\n");
    exit(EXIT_FAILURE);
  }
  printf("%d by %d matrix\n", dim[0], dim[1]);
  data = (double*) malloc(dim[0] * dim[1] * sizeof(data_t));
  if (data == NULL)
  {
    perror("Error allocating data matrix");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < dim[0] * dim[1]; i++)
  {
    if ((j = fscanf(infile, PRINT_DATA, &data[i])) != 1)
    {
      if (j == EOF)
        printf("File does not contain %d * %d data elements!\n", dim[0], dim[1]);
      else
        printf("Error reading data from file\n");
      exit(EXIT_FAILURE);
    }
    //printf(" %lf", data[i]);
    //if (i % dim[1] == dim[1] - 1)
    //  printf("\n");
  }
  fclose(infile);
  
  printf("Opening...\n");
  if (MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile) != MPI_SUCCESS)
  {
    perror("Error opening output file");
    exit(EXIT_FAILURE);
  }
  printf("Writing...\n");
  if (MPI_File_write(outfile, &dim[0], 2, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS)
  {
    perror("Error writing row and column count to output file");
    exit(EXIT_FAILURE);
  }
  if (MPI_File_write(outfile, data, dim[0] * dim[1], MPI_DATA, MPI_STATUS_IGNORE) != MPI_SUCCESS)
  {
    perror("Error writing data matrix to output file");
    exit(EXIT_FAILURE);
  }
  
  if (MPI_File_close(&outfile) != MPI_SUCCESS)
  {
    perror("Error closing output file (?)");
    exit(EXIT_FAILURE);
  }
  printf("Done\n");
  
  free(data);
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
