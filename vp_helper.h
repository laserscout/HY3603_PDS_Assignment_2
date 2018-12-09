#ifndef _VP_HELPER_H_
#define _VP_HELPER_H_

struct points;
struct a_point;

extern MPI_Comm comm;
extern int processId;
extern int globalId;
extern int noProcesses;
extern int size;
extern int D;
extern int N;


int read_data(FILE *bin, points *x);
int comm_split();
int set_vp(a_point *vp);
int find_dists(points *x, a_point *vp, float *dist);
float find_median(float *dist);
int points_exchange(a_point *vp, points *x, float *dist, float median);

#endif
