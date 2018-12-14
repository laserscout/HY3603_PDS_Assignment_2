#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include "vp_helper.h"
#include "mpiFindMedian.h"


//extern MPI_Comm comm;

int read_data(FILE *bin, struct points *x) {
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

  x->data = input_data;
  return 0;
}

int rand_data(struct points *x) {
  intype **input_data, *ptr;
  int length;
  float a = 10;
  time_t t;

  //printf("%d:rand_data running...\n",processId);
  length = N * sizeof(intype *) + N * D * sizeof(intype);
  //printf("length=%d\n",length);
  input_data = (intype **)malloc(length);
  if (input_data == 0) {
    printf ("Could not allocate %d bytes for points struct creation", length);
    exit(1);
  }
  ptr =(intype *)input_data + N;
  for(int i=0;i<N;i++)
    input_data[i] = ptr + i*D; //https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/

  srand( (unsigned)time(&t)+processId );
  for(int i=0; i<N; i++) {
    for(int d=0;d<D;d++){
      input_data[i][d] = ((double)rand()/(double)(RAND_MAX))*a; //https://stackoverflow.com/a/13409133
      //printf("%f ",input_data[i][d]);
    }
  }

  x->data = input_data;
  return 0;
}

int comm_split() {
  int color;
  MPI_Comm group_old;
  
  color = processId/(noProcesses/2);
  group_old = comm;
  
  assert( MPI_Comm_split( group_old, color, globalId, &comm ) == 0 );
  MPI_Comm_rank( comm, &processId );
  MPI_Comm_size( comm, &noProcesses);
  assert( MPI_Comm_free( &group_old ) == 0 ); // make sure that the first group isn't the global
  size = size/2; //sloppy
  
  return 0;
}
 
int set_vp(struct points *x, struct a_point *vp) {

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
      printf("Error,rand generated an index larger than the data size/n %d/n",index);
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
	  MPI_Bcast(vp->data,D,MPI_DOUBLE,0,comm); //used MPI_DOUBLE here!!
	}
	else {
	  MPI_Send(&parindex,1,MPI_INT,i,1,comm);
	  MPI_Bcast(vp->data,D,MPI_DOUBLE,i,comm); //used MPI_DOUBLE here!!
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
      MPI_Bcast(vp->data,D,MPI_DOUBLE,chosen,comm); //used MPI_DOUBLE here!!
    }
    else
      MPI_Bcast(vp->data,D,MPI_DOUBLE,chosen,comm); //used MPI_DOUBLE here!!
    flag=0;
  }

  return flag;
}

int find_dists(struct points *x, struct  a_point *vp, float **dist) {
  double sum=0.0;
  float *fptr;

  fptr = (float *)malloc(N*sizeof(float));
  if (fptr==NULL) {
    printf("Couldn't allocate memory for the distance vector");
    exit(1);
  }
  for (int i=0;i<N;i++){
    for (int n=0;n<D;n++)
      sum += pow((double)(x->data[i][n] - vp->data[n]),2);
    fptr[i]= (float)sqrt(sum);
    /* printf("%d:dist[%d]=%f\n",processId,i,fptr[i]); */
    sum=0;
  }
  *dist = fptr;

  return 0;
}

float find_median(float *dist) {
  float median, *dist_copy;

  dist_copy = (float *)malloc(N*sizeof(float));
  if (dist_copy==NULL) {
    printf("Couldn't allocate memory for the distance vector");
    exit(1);
  }
  for (int i = 0;i<N;i++)
    dist_copy[i]=dist[i];

  /* for (int p=0;p<noProcesses;p++){ */
  /*   if (p==processId) { */
  /*     if (p==0) printf("\n\nBefore find_median\n"); */
  /*     for (int i=0;i<N;i++) */
  /* 	printf("%d:dist[%d]=%f\n",processId,i,dist[i]); */
  /*   } */
  /*   MPI_Barrier(comm); */
  /* }   */

  
  if (processId==0) {
    if (noProcesses==1){
      median=selection(dist_copy,size);
      validationST(median,size,dist);
    }
    else {
      median=masterPart(noProcesses,processId,size,N,dist_copy,comm);
      printf("median=%f\n",median);
    }
  }
  else
    slavePart(processId,N,dist_copy,size,comm);

  MPI_Bcast(&median,1,MPI_FLOAT,0,comm);

  /* for (int p=0;p<noProcesses;p++){ */
  /*   if (p==processId) { */
  /*     if (p==0) printf("\n\nAfter find_median\n"); */
  /*     for (int i=0;i<N;i++) */
  /* 	printf("%d:dist[%d]=%f\n",processId,i,dist[i]); */
  /*   } */
  /*   MPI_Barrier(comm); */
  /* }   */
    
  return median;
}

 int points_exchange(struct a_point *vp, struct  points *x, float *dist, float median) {
   int count=0;
   int noSwaps=0;
   int SRtable[noProcesses];
   int sw[N]; // indexes of elements to swap
   int hp = noProcesses/2;
   int priority=0;
   int countS=0;
   int swapid,send_start,send_finish,queue_start,length;
   intype **swap_data, *ptr;
   MPI_Status *status;
   MPI_Request *q;

   for (int p=0;p<noProcesses;p++){
     if (p==processId) {
       for (int i=0;i<N;i++) {
	 printf("%d:dist[%d]=%f",processId,i,dist[i]);
	 if (p<hp && dist[i]>median)
	   printf("<-");
	 if (p>=hp && dist[i]<=median)
	   printf("<-");
	 printf("\n");
       }
     }
     MPI_Barrier(MPI_COMM_WORLD);
   }

   //points toswap;

   //a_point temp;
   //temp->data = (*intype)malloc(8*sizeof(intype));

   if (processId<hp) {
     for (int n=0;n<N;n++) {
       //printf("%f<%f",dist[n],median);
       if (dist[n]>median) {
	 sw[noSwaps] = n;
	 noSwaps++;
	 //printf("y\n");
       }
       //else
	 //printf("n\n");
     }
   }
   else {
     for (int n=0;n<N;n++) {
       //printf("%f>=%f",dist[n],median);
       if (dist[n]<=median) {
	 sw[noSwaps] = n;
	 noSwaps++;
	 //printf("y\n");
       }
       //else
	 //printf("n\n");
     }
   }

   for (int p=0;p<noProcesses;p++){
     if (p==processId) {
       printf("%d:sw = ",processId);
       for (int i=0;i<noSwaps;i++)
	 printf("%d ",sw[i]);
       printf("\n");
     }
     MPI_Barrier(MPI_COMM_WORLD);
   }

   if (noSwaps==0) return 0;

   length = noSwaps * sizeof(intype *) + noSwaps * D * sizeof(intype);
   swap_data = (intype **)malloc(length);
   if (swap_data == NULL) {
     printf ("Could not allocate %d bytes for swap points struct creation", length);
     exit(1);
   }
   ptr = (intype *)swap_data + noSwaps;
   for (int i=0;i<noSwaps;i++)
     swap_data[i] = ptr + i*D;

   for (int n=0;n<noSwaps;n++)
     for (int i=0;i<D;i++)
       swap_data[n][i] = x->data[sw[n]][i];

   //dataprint(swap_data,noSwaps);   

   status = (MPI_Status *)malloc(noSwaps*2*sizeof(MPI_Status));
   q = (MPI_Request *)malloc(noSwaps*2*sizeof(MPI_Request));
   if (status==NULL || q==NULL) {
     printf ("Could not allocate memory for the MPI_Status and MPI_Request arrays");
     exit(1);
   }

   //   toswap->data = swap_data;

   //procceses 0 - np/2-1 keep the smaller elements and np/2 - np the larger
   MPI_Allgather(&noSwaps,1,MPI_INT,SRtable,1,MPI_INT,comm);

   for (int p=0;p<noProcesses;p++){
     if (p==processId) {
       printf("%d:",processId);
          for (int i=0;i<noProcesses;i++)
	    printf("%d ",SRtable[i]);
	  printf("\n");
     }
     MPI_Barrier(MPI_COMM_WORLD);
   }
   
   if (processId<hp){
     send_start = hp ;
     send_finish = noProcesses;
     queue_start = 0;
   }
   else {
     send_start = 0 ;
     send_finish = hp;
     queue_start = hp;
   }
   
   for (int p=queue_start;p<processId;p++)
     priority+=SRtable[p];
   swapid = priority;
     
   for (int p=send_start;p<send_finish;p++) {
     priority-=SRtable[p];
     printf("%d:priority=%d\n",processId,priority);
     while (priority<0 && countS<noSwaps) {
       MPI_Isend(swap_data[countS],D,MPI_DOUBLE,p,swapid,comm,q+countS);
       MPI_Irecv(x->data[sw[countS]],D,MPI_DOUBLE,p,swapid,comm,q+noSwaps+countS);
       priority++;
       swapid++;
       countS++;
       printf("%d:swap! %d of %d\n",processId,countS,noSwaps);
     }
     if(countS>=noSwaps){
       printf("%d:break",processId);
       break; // just to make sure it breaks, it's out of the while loop
     }
   }
   printf("%d: Before Wait_all\n",processId);
   assert (MPI_Waitall(noSwaps*2,q,status) == 0);
   printf("%d: I'm here\n",processId);
   free(q);
   free(status);
   free(swap_data);
  return 0;
}

struct vp_tree * newvp(struct a_point *center, float radius, int depth){
  struct vp_tree *node = (struct vp_tree *)malloc(sizeof(struct vp_tree));
  node->center = *center;
  node->radius = radius;
  node->bigger = NULL;
  node->smaller = NULL;
  node->depth = depth;
  //struct vp_tree *root;
  //int MPI_process_id_kati;
}

int tree_grow(struct vp_tree *root, struct a_point *center, float radius, int depth) {
  if (processId==0) {
    root->smaller = newvp(center, radius, depth);
  }
  return 0;
}

void dataprint(intype **x, int nom) {
  /* printf("dataprint running...\n"); */
  intype buff[D];
  MPI_Status stat;

  if (globalId==0) {
    printf("---0---\n");
    for (int n=0;n<nom;n++) {
      for (int d=0;d<D;d++) {
	printf("%f ", x[n][d]);
      }
      printf("\n");
    }

    for (int p=1;p<globalNo;p++) {
      printf("---%d--\n", p);
      for (int n=0;n<nom;n++) {
	MPI_Recv(buff,D,MPI_DOUBLE,p,1,MPI_COMM_WORLD,&stat);
	for (int d=0;d<D;d++) {
	  printf("%f ", buff[d]);
	}
	printf("\n");
      }
    }
  }
  else
    for (int i=0;i<nom;i++)
      MPI_Send(x[i],D,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
   
    MPI_Barrier(MPI_COMM_WORLD);
}

  
void pointprint(intype *x) {
  /* printf("dataprint running...\n"); */

  intype buff[D];
  MPI_Status stat;

  if (globalId==0) {
    printf("---0---\n");
    for (int d=0;d<D;d++) {
      printf("%f ", x[d]);
    }
    printf("\n");
    
    for (int p=1;p<globalNo;p++) {
      printf("---%d--\n", p);
      MPI_Recv(buff,D,MPI_DOUBLE,p,1,MPI_COMM_WORLD,&stat);
      for (int d=0;d<D;d++) {
	printf("%f ", buff[d]);
      }
      printf("\n");
    }
  }
  else
    MPI_Send(x,D,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
   
    MPI_Barrier(MPI_COMM_WORLD);
}
