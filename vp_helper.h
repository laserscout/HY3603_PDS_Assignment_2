#ifndef _VP_HELPER_H_
#define _VP_HELPER_H_

typedef double intype;

struct points {
  intype **data;
};

struct a_point {
  intype *data;
};

extern MPI_Comm comm;
extern int processId;
extern int globalId;
extern int noProcesses;
extern int size;
extern int D;
extern int N;


int read_data(FILE *bin, struct points *x);
int rand_data(struct points *x);
int comm_split();
int set_vp(struct points *x, struct a_point *vp);
int find_dists(struct points *x, struct a_point *vp, float *dist);
float find_median(float *dist);
int points_exchange(struct a_point *vp, struct points *x, float *dist, float median);

#endif
