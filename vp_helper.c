#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include "vp_helper.h"
#include "findMedian.h"

int read_data(FILE *bin, struct points *x) {
  intype **input_data, *ptr;
  int length;

  length = N * sizeof(intype *) + N * D * sizeof(intype);
  input_data = (intype **)malloc(length);
  if (input_data == 0) {
    printf ("Could not allocate %d bytes for points struct creation", length);
    exit(1);
  }

  ptr =(intype *)input_data + N;
  for(int i=0;i<N;i++)
    input_data[i] = ptr + i*D;
  //https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/

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

  length = N * sizeof(intype *) + N * D * sizeof(intype);
  input_data = (intype **)malloc(length);
  if (input_data == 0) {
    printf ("Could not allocate %d bytes for points struct creation", length);
    exit(1);
  }
  ptr =(intype *)input_data + N;
  for(int i=0;i<N;i++)
    input_data[i] = ptr + i*D;
  //https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/

  srand( (unsigned)time(&t)+processId );
  for(int i=0; i<N; i++) {
    for(int d=0;d<D;d++){
      input_data[i][d] = ((double)rand()/(double)(RAND_MAX))*a;
      //https://stackoverflow.com/a/13409133
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
  assert( MPI_Comm_free( &group_old ) == 0 );
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

    if (noProcesses==1) { //I'm on my own, do not broadcast
      for (int n=0; n<D; n++)
	vp->data[n] = x->data[index][n];
      return 0;
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
	  MPI_Bcast(vp->data,D,MPI_DOUBLE,0,comm); 
	}
	else {
	  MPI_Send(&parindex,1,MPI_INT,i,1,comm);
	  MPI_Bcast(vp->data,D,MPI_DOUBLE,i,comm); 
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
      MPI_Bcast(vp->data,D,MPI_DOUBLE,chosen,comm); 
    }
    else
      MPI_Bcast(vp->data,D,MPI_DOUBLE,chosen,comm);
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
    sum=0;
  }
  *dist = fptr;

  return 0;
}

float find_median(float *dist,intype **data) {
  float median;

  if (processId==0) {
    if (noProcesses==1){
      median=selection(dist,data,size);
      validationST(median, size, dist);
      return median;
    }
    else {
      median=masterPart(noProcesses,processId,size,N,dist,data,comm);
      validation(median, N, size,dist, processId, comm);
    }
  }
  else {
    slavePart(processId,N,dist,size,data,comm);
    validation(median, N, size,dist, processId, comm);
  }

  MPI_Bcast(&median,1,MPI_FLOAT,0,comm);

  return median;
}

 int points_exchange(struct  points *x, float *dist, float median) {
   int check;
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

   /* for (int p=0;p<noProcesses;p++){ */
   /*   if (p==processId) { */
   /*     for (int i=0;i<N;i++) { */
   /* 	 printf("%d:dist[%d]=%f",processId,i,dist[i]); */
   /* 	 if (p<hp && dist[i]>median) */
   /* 	   printf("<-"); */
   /* 	 if (p>=hp && dist[i]<=median) */
   /* 	   printf("<-"); */
   /* 	 printf("\n"); */
   /*     } */
   /*   } */
   /*   MPI_Barrier(MPI_COMM_WORLD); */
   /* } */
   if (processId<hp) {
     for (int n=0;n<N;n++) {
       if (dist[n]>median) {
	 sw[noSwaps] = n;
	 noSwaps++;
       }
     }
   }
   else {
     for (int n=0;n<N;n++) {
       if (dist[n]<=median) {
	 sw[noSwaps] = n;
	 noSwaps++;
       }
     }
   }
   
   /* for (int p=0;p<noProcesses;p++){ */
   /*   if (p==processId) { */
   /*     printf("%d:sw = ",processId); */
   /*     for (int i=0;i<noSwaps;i++) */
   /* 	 printf("%d ",sw[i]); */
   /*     printf("\n"); */
   /*   } */
   /*   MPI_Barrier(MPI_COMM_WORLD); */
   /* } */

   MPI_Allgather(&noSwaps,1,MPI_INT,SRtable,1,MPI_INT,comm);

   /* for (int p=0;p<noProcesses;p++){ */
   /*   if (p==processId) { */
   /*     printf("%d:",processId); */
   /*        for (int i=0;i<noProcesses;i++) */
   /* 	    printf("%d ",SRtable[i]); */
   /* 	  printf("\n"); */
   /*   } */
   /*   MPI_Barrier(MPI_COMM_WORLD); */
   /* } */

   if (processId==0){
     check=0;
     for (int i=0;i<hp;i++)
       check+=SRtable[i];
     for (int i=hp;i<noProcesses;i++)
       check-=SRtable[i];
     if (check!=0) {
       printf("Send and recive halves are not equal");
       exit(1);
     }
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
     while (priority<0 && countS<noSwaps) {
       MPI_Isend(swap_data[countS],D,MPI_DOUBLE,p,swapid,comm,q+countS);
       MPI_Irecv(x->data[sw[countS]],D,MPI_DOUBLE,p,swapid,comm,q+noSwaps+countS);
       priority++;
       swapid++;
       countS++;
       //printf("%d:swap! %d of %d\n",processId,countS,noSwaps);
     }
     if(countS>=noSwaps){
       //printf("%d:break",processId);
       break; // just to make sure it breaks, it's out of the while loop
     }
   }
   assert (MPI_Waitall(noSwaps*2,q,status) == 0);

   free(q);
   free(status);
   free(swap_data);
  return 0;
}

int tree_init(struct vp_tree **t, int depth){

  intype *c;
  int length;
  
  length = (1<<depth) - 1;
  *t = (struct vp_tree *)malloc(length*sizeof(struct vp_tree));
  if (*t == NULL) {
    printf("%d:Cound not allocate memory for the tree\n",globalId);
    exit(1);
  }

  c=(intype *)malloc(D*length*sizeof(intype));
  if (c==NULL) {
    printf("%d:Cound not allocate memory for the tree centers\n",globalId);
    exit(1);
  }
  for (int i=0;i<length;i++)
    (*t)[i].center = c+i*D;

  return 0;
}
	  
int tree_grow(struct a_point *center, float radius, int depth, int rep, struct vp_tree *t){
  int idx;
  float r;
  intype c[D];
  MPI_Status stat;

  idx = (1<<depth)-1+rep*globalNo;
  if (globalId == 0) {
    //printf("starting depth:%d idx:%d\n",depth,idx);
    t[idx].radius = radius;
    t[idx].depth = depth;
    for (int d=0;d<D;d++)
      t[idx].center[d] = center->data[d];

    int others=globalNo/noProcesses;
    for (int i=1; i<others;i++) {
      int p = i*noProcesses;
      MPI_Recv(c,D,MPI_DOUBLE,p,1,MPI_COMM_WORLD,&stat);
      MPI_Recv(&r,1,MPI_FLOAT,p,2,MPI_COMM_WORLD,&stat);
      //printf("recived depth:%d idx:%d\n",depth,idx+i);
      t[idx+i].radius = r;
      t[idx+i].depth = depth;
      for (int d=0;d<D;d++)
	t[idx+i].center[d] = c[d];
    }
  }
  else if (processId==0) {
    //printf("sending... depth:%d idx:%d\n",depth,idx+globalId);
    MPI_Send(center->data,D,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
    MPI_Send(&radius,1,MPI_FLOAT,0,2,MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
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
