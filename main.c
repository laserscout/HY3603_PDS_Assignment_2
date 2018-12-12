#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "vp_helper.h"

int main (int argc, char **argv) {
  int reps;
  FILE f;
  points x;
  a_point vp;
  float *dist,median;
  
  if (argc != 2) {
    printf("please enter a loop count, l");
    exit(1);
  }
  reps = atoi(argv[2]);

  assert(f = fopen("data.bin", "r") == 0);
  assert(read_data(f,&x) == 0);

  for (int l=0;l<reps;l++) {
    if (l!=0)
      assert(comm_split() == 0);

    set_vp(&vp);
    find_dists(&x, &vp, &dist);
    median = find_median(dist);
    points_exchange(&vp,&x,dist,median);

    free(dist);
    free(vp);
  }
  
  return 0;
}
