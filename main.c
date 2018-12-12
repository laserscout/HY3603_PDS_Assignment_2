#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "vp_helper.h"

MPI_Comm comm;
int processId;
int globalId;
int noProcesses;
int size;
int D;
int N;

int main (int argc, char **argv) {
  int reps,logp,l;
  FILE *f;
  struct points x;
  struct a_point vp;
  float *dist,median;

  D = 4;
  N = 128;
  MPI_Init (&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank (MPI_COMM_WORLD, &processId);
  MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);
  
  if (argc != 2) {
    printf("please enter a loop count, l");
    exit(1);
  }
  reps = atoi(argv[2]);

  
  /* f = fopen("data.bin", "r"); */
  /* if(f == NULL) */
  /*  { */
  /*     printf("Error opening file");    */
  /*     exit(1);              */
  /*  } */

  /* assert(read_data(f,&x) == 0); */
  //rand function for testing bellow
  rand_data(&x);

  logp=0;
  for (l=0;l<reps;l++)
    if (1<<l==noProcesses)
      logp=l;
  
  if (logp==0) {
    printf(" -np needs to be a power of two and larger than 2^l");
    exit(1);
  }
  
  for (l=0;l<logp;l++) {
    if (l!=0)
      assert(comm_split() == 0);

    set_vp(&x,&vp);
    find_dists(&x, &vp, dist);
    median = find_median(dist);
    points_exchange(&vp,&x,dist,median);

    free(dist);
    free(vp.data);
  }

  while(l<reps) {
    
    l++;
  }
  
  return 0;
}


