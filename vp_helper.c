#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include "vp_helper.h"

typedef double intype;

struct points {
  intype **data;
};

struct a_point {
  intype *data;
};


int read_data(FILE *bin, points *x) {
  intype **input_data;
  int length;

  length = N * sizeof(intype *) + N * D * sizeof(intype);
  input_data = (intype **)malloc(length);
  if (input_data == 0) {
    printf ("Could not allocate %d bytes for points struct creation", length);
    exit(1);
  }

  for(int i=0; i<N; i++) {
    fread(input_data[i], D * sizeof(intype), 1, bin);
  }

  fclose(bin);

  points->data = input_data;
  return 0;
}

int comm_split() {
  int color;
  MPI_Comm group_old;
  
  color = processId/(noProcesses/2);
  group_old = *comm;
  
  assert( MPI_Comm_split( group_old, color, globalId, comm ) == 0 );
  MPI_Comm_rank( comm, &processId );
  MPI_Comm_size( comm, &noProcesses);
  assert( MPI_Comm_free( &group_old ) == 0 ); // make sure that the first group isn't the global

  return 0;
}
 
int set_vp(a_point *vp) {

  int index, flag=1;
  MPI_Status stat;

  vp->data = (intype *)malloc(D * sizeof(intype));

  if (processId==0) {
    time_t t;
    int odds, rounds, parindex;
    int one=0;

    srand( (unsigned)time(&t) );
    index = rand() % size; // a number from 0 to size
    if (index>size) {
      printf{"Error,rand generated an index larger than the data size/n %d/n",index);
      exit(1);
    }
    rounds = size / noProcesses;
    odds = size % noProcesses;

    for (int i=0;i<noProcesses;i++) {
      parindex = index - rounds*i - one;
      if (odds>0){
	one++;
	odds--;
      }
      if (index < rounds*(i+1)+one) {
	MPI_Bcast(&i,1,MPI_INT,0,comm);
	if (i==0) {
	  for (int n=0; n<D; n++)
	    vp->data[n] = x->data[parindex][n];
	  MPI_Bcast(&vp->data,D,MPI_DOUBLE,0,comm); //used MPI_DOUBLE here!!
	}
	else {
	  MPI_Send(&parindex,1,MPI_INT,i,1,comm);
	  MPI_Bcast(&vp->data,D,MPI_DOUBLE,i,comm); //used MPI_DOUBLE here!!
	}
	flag=0;
	break;
      }
    }	
  }
  else {
    int chosen;
    MPI_Bcast(&chosen,1,MPI_INT,0,comm);
    if (chosen==processId){
      MPI_Recv(&index,1,MPI_INT,0,1,comm,&stat);
      for (int n=0; n<D; n++)
	vp->data[n] = x->data[index][n];
      MPI_Bcast(&vp->data,D,MPI_DOUBLE,chosen,comm); //used MPI_DOUBLE here!!
    }
    else
      MPI_Bcast(&vp->data,D,MPI_DOUBLE,chosen,comm); //used MPI_DOUBLE here!!
    flag=0;
  }
		  
  return flag;
}

int find_dists(points *x, a_point *vp, float *dist) {
  double sum=0.0;
  
  dist = (float *)malloc(N*sizeof(float));
  if (dist==0) {
    printf("Couldn't allocate memory for the distance vector");
    exit(1);
  }
  for (int i=0;i<N;i++){
    for (int n=0;n<D;n++)
      sum += exp((double)(x->data[i][n] - vp->data[n]));
    dist[i]= (float)sqrt(sum);
    sum=0;
  }
  return 0;
}

float find_median(float *dist) {
  float median;
  
  if (processId==0) {
    if (noProcesses==1){
      median=selection(dist,size);
      validationST(median,size,dist);
    }
    else
      median=masterPart(noProcesses,processId,size,N,dist,comm);
  }
  else
    slavePart(processId,N,dist,size,comm);
    
  return median;
}

 int points_exchange(a_point *vp, points *x, float *dist, float median) {

   
  return 0;
}

