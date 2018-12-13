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
  struct a_point center;
  float radius;
  struct vp_tree *bigger;
  struct vp_tree *smaller;
  struct vp_tree *root;
  int MPI_process_id_kati;
  int depth;
};

extern MPI_Comm comm;
extern int processId;
extern int globalId;
extern int noProcesses;
extern int globalNo;
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
int tree_grow(struct vp_tree *root, struct a_point *center, float radius, int depth);
void dataprint(struct points *x);
void pointprint(struct a_point *x);
#endif
