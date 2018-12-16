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
int N;
int D;

int main (int argc, char **argv) {
  int reps,logp,l;
  int init_size;
  FILE *f;
  struct points x;
  struct a_point vp;
  struct vp_tree *tree;
  float *dist;
  float median;

  D = 3;
  N = 16;
  init_size=N;

  if (argc != 2) {
    printf("please enter a loop count, l\n");
    exit(1);
  }
  reps = atoi(argv[1]);

  logp=-1;
  for (l=1;l<=reps;l++) {
    if (1<<l==noProcesses){
      logp=l;
      break;
    }
  }
  if (logp<0) {
    printf("-np needs to be a power of two and larger than 2^l\n");
    exit(1);
  }

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &globalId);
  MPI_Comm_size (MPI_COMM_WORLD, &globalNo);

  assert (MPI_Comm_dup(MPI_COMM_WORLD, &comm) == 0 );
  processId = globalId;
  noProcesses = globalNo;
  size = N*noProcesses; //sloppy
  
  // if you want to actualy use data from a binary file, here is a
  // function that will read the data. Or, you can use random data
  // with the other function bellow.
  if (0) {
    f = fopen("data.bin", "r");
    if(f == NULL) {
      printf("Error opening file");
      exit(1);
    }
    assert(read_data(f,&x) == 0);
  }
  else
    rand_data(&x);

  
  dataprint(x.data,N);

  if (globalId==0)
    tree_init(&tree,reps);
  
  for (l=0;l<logp;l++) {
    set_vp(&x,&vp);
    //pointprint(vp.data);
    find_dists(&x, &vp, &dist);
    median = find_median(dist,x.data);
    //printf("id:%d comm:%d np:%d size:%d l:%d logp:%d, median:%f\n",processId,comm,noProcesses,size,reps,logp,median);
    points_exchange(&vp,&x,dist,median);
    tree_grow(&vp, median, l, 0, tree);
    //dataprint(x.data,N);

    free(dist);
    free(vp.data);
    comm_split();
  }

  int sl=0; 
  //continiue from where we left off
  while(l<reps) {
    int parts = 1<<sl;
    for (int i=0;i<parts;i++) {
      struct points x_part;
      x_part.data = x.data+i*size;
      
      //dataprint(x_part.data,size);
      set_vp(&x_part,&vp);
      //pointprint(vp.data);
      find_dists(&x_part,&vp,&dist);
      median = find_median(dist,x_part.data);
      //printf("id:%d size:%d, l:%d rep:%dof%d   and the median is:%f\n",globalId,size,l,i+1,parts,median);
      tree_grow(&vp, median, l, i, tree);
      
      free(vp.data);
      free(dist);
    }
    size = size/2;
    N = N/2;
    sl++;
    l++;
  }

  //Print a representation of the tree. Set printTree=1 to print
  int printTree = 1;
  if (globalId ==0 && printTree) {
    printf("=====The tree====\n");
    for (int i=0;i<(1<<reps)-1;i++) {
      printf("%d: ",i);
      for (int d=0;d<D;d++)
	printf("%f ",tree[i].center[d]);
      printf("R:%f, depth:%d\n",tree[i].radius,tree[i].depth);
    }
  }

  dataprint(x.data,init_size);

  MPI_Finalize();
  
  return 0;
}


