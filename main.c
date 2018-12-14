#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include "vp_helper.h"

MPI_Comm comm;
int processId;
int globalId;
int noProcesses;
int globalNo;
int size;
int D;
int N;

int main (int argc, char **argv) {
  int reps,logp,l;
  FILE *f;
  struct points x;
  struct a_point vp;
  float *dist;
  float median;

  D = 1;
  N = 8;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &globalId);
  MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);

  assert (MPI_Comm_dup(MPI_COMM_WORLD, &comm) == 0 );
  processId = globalId;
  globalNo = noProcesses;
  size = N*noProcesses; //sloppy
  //assert(MPI_Comm_split(MPI_COMM_WORLD,1,globalId,&comm) == 0);

  
  if (argc != 2) {
    printf("please enter a loop count, l\n");
    exit(1);
  }
  reps = atoi(argv[1]);

  logp=-1;
  for (l=0;l<reps;l++)
    if (1<<l==noProcesses)
      logp=l;

    if (logp<0) {
    printf("-np needs to be a power of two and larger than 2^l\n");
    exit(1);
  }

  /* f = fopen("data.bin", "r"); */
  /* if(f == NULL) */
  /*  { */
  /*     printf("Error opening file");    */
  /*     exit(1);              */
  /*  } */

  /* assert(read_data(f,&x) == 0); */
  //rand function for testing bellow
  rand_data(&x);

  dataprint(x.data,N);
    
  for (l=0;l<logp;l++) {
    printf("id:%d comm:%d np:%d size:%d l:%d logp:%d\n",processId,comm,noProcesses,size,reps,logp);
  
    set_vp(&x,&vp);
    //pointprint(vp.data);
    find_dists(&x, &vp, &dist);
    median = find_median(dist);
    points_exchange(&vp,&x,dist,median);
    //tree_grow(&vp, median, l);
    //free(dist);
    //free(vp.data);
    comm_split();
    printf("\ntest\n");
    dataprint(x.data,N);
  }

  while(l<reps) {
    
    l++;
  }

  MPI_Finalize();
  
  return 0;
}


