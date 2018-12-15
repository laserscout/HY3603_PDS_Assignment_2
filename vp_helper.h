#ifndef _VP_HELPER_H_
#define _VP_HELPER_H_

typedef double intype;

struct points {
  intype **data;
};

struct a_point {
  intype *data;
};

struct vp_tree {
  intype *center;
  float radius;
  int owner;
  int depth;
};

extern MPI_Comm comm;
extern int processId;
extern int globalId;
extern int noProcesses;
extern int globalNo;
extern int size;
extern int N;
extern int D;


int read_data(FILE *bin, struct points *x);
int rand_data(struct points *x);
int comm_split();
int set_vp(struct points *x, struct a_point *vp);
int find_dists(struct points *x, struct a_point *vp, float **dist);
float find_median(float *dist,intype **data);
int points_exchange(struct points *x, float *dist, float median);
int tree_grow(struct a_point *center, float radius, int depth, int rep, struct vp_tree *t);
void dataprint(intype **x, int nom);
void pointprint(intype *x);
int tree_init(struct vp_tree **t, int depth);

#endif
