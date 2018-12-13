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

  D = 2;
  N = 8;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &globalId);
  MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);

  //  assert (MPI_Comm_dup(MPI_COMM_WORLD, &comm) == 0 );
  processId = globalId;
  globalNo = noProcesses;
  size = N*noProcesses; //sloppy
  assert(MPI_Comm_split(MPI_COMM_WORLD,1,globalId,&comm) == 0);

  
  if (argc != 2) {
    printf("please enter a loop count, l\n");
    exit(1);
  }
  reps = atoi(argv[1]);

  logp=-1;
  for (l=0;l<reps;l++)
    if (1<<l==noProcesses)
      logp=l;

  printf("id:%d\nnp:%d\nsize:%d\nl:%d\nlogp:%d\n",processId,noProcesses,size,reps,logp);
  
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
  printf("test\n");
    
  for (l=0;l<logp;l++) {
    set_vp(&x,&vp);
    pointprint(&vp);

    dist = (float *)malloc(N*sizeof(float));
    find_dists(&x, &vp, dist);

    /* MPI_Barrier; */
    /* for (int i=0;i<N;i++) */
    /*     printf("%d:dist[%d]=%f\n",processId,i,dist[i]); */
    /* MPI_Barrier; */

    /* double sum=0.0; */
    /* float *d_ptr; */
    /* d_ptr = (float *)malloc(N*sizeof(float)); */
    /* if (d_ptr==0) { */
    /*   printf("Couldn't allocate memory for the distance vector"); */
    /*   exit(1); */
    /* } */
    /* for (int i=0;i<N;i++){ */
    /*   for (int n=0;n<D;n++) */
    /* 	sum += exp((double)(x.data[i][n] - vp.data[n])); */
    /*   d_ptr[i]= (float)sqrt(sum); */
    /*   printf("%d:dist[%d]=%f\n",processId,i,d_ptr[i]); */
    /*   sum=0; */
    /* } */
    /* dist = d_ptr; */
    
    /* median = find_median(dist); */

    
    //points_exchange(&vp,&x,dist,median);
    //tree_grow(&vp, median, l);
    //free(dist);
    //free(vp.data);
    /* comm_split(); */
    dataprint(&x);
  }

  while(l<reps) {
    
    l++;
  }

  MPI_Finalize();
  
  return 0;
}


