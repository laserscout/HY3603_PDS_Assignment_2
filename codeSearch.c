#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include "vp_helper.h"
#include "mpiFindMedian.h"

float distance(intype* pA, intype *pB);
int smallerRadius(intype *pA, struct vp_tree *vp);
struct vp_tree getSibling(struct vp_tree *input_node);
struct vp_tree *find_first(intype *pA, struct vp_tree *input_node);
intype *findNeighbours(intype *pA, int k, struct vp_tree input_node);

// k = number of neightbours with scanf
int k;
// depth of the tree
int depth;
// depth of the node
int nodeDepth;
// radius of the query point
float radiusQ;


struct vp_tree {
  //struct a_point center
  struct vp_tree *father; // father =NULL for root node
  intype *center; // vantage point
  float radius;
  struct vp_tree *bigger;
  struct vp_tree *smaller;
  //struct vp_tree *root;
  int MPI_process_id_kati;
  int depth;
};


intype *knn(intype *pA, struct vp_tree *treeStartNode, int k) {

  struct vp_tree *first_node;
  struct vp_tree *temp_node;
  intype neighbours[k];

  // find the first node
  first_node = find_first(pA, treeStartNode);
  neighbours = findNeighbours(pA, k, first_node, radiusQ);

  temp_node = getSibling(first_node);
  if((radiusQ + temp_node->radius) > distance(pA, temp_node->center)) {
    // search for neighbours in this leaf and calculate the new radius of the query point
  }

}


intype *findNeighbours(intype *pA, int k, struct vp_tree input_node) {
  // find the k nearest neighbours of pA point in the space with the vantage point of the input_node
  // return an array of these points
  // and save the radius of the 'neighbourhood' in searchRadius
}

struct vp_tree getSibling(struct vp_tree *input_node) {
  struct vp_tree *my_node;
  struct vp_tree *father_node;
  my_node = input_node;
  father_node = input_node->father;
  struct vp_tree *next_child;

  if(my_node == father_node->smaller)
    next_child = father_node->bigger;
  else
    next_child = father_node->smaller;

  return next_child;
}


// find the first vp node of pA
struct vp_tree *find_first(intype *pA, struct vp_tree *input_node) {

  int flag; // check for radius FLAG

// check if the node is a leaf
// if node is not a leaf, then..
  if(input_node->bigger =! NULL || input_node->smaller =! NULL) {
    // check the radius
    flag = smallerRadius(pA, input_node);
    // choose leaf
    if(flag == 0)
      find_first(pA, input_node->smaller);
    else
      find_first(pA, input_node->bigger);
  }
  else {
    // return the center of the bottom vantage point
    return input_node;
  }
}

// MIPWS THELEI KAI SYGKRISI KAI GIA ISO???
// check radius
int smallerRadius(intype *pA, struct vp_tree *vp) {

  float radiusA, radiusB;
  radiusA = distance(pA, vp->center);
  radiusB = vp->radius;

// if pA has smaller radius, return 0
  if(radiusA < radiusB)
    return 0;
// else return 1
  else
    return 1;
}


// find distance
float distance(intype* pA, intype *pB) {

  float dist;
  double temp;
  float sum = 0.0, dif = 0.0;
  squared = 0.0;

// D = dimensions
  for(int i = 0; i<D; i++) {
    dif = pA[i] - pB[i];
    dif = (float) abs((float) dif);
    squared = dif * dif;
    sum = sum + squared;
  }

  temp = sqrt((double) sum);
  dist = (float) temp;
  return dist;
}
