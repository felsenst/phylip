/* Version 4.0.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
*/


/*
    dist.h: included in fitch, kitsch, & neighbor
*/

#include "ml.h"

#define over            60


typedef long *intvector;

typedef node **pointptr;

/* debug: typedef struct ml_tree ml_tree; */

struct dist_tree {
  struct ml_tree ml_treepart;
} dist_tree;

/* debug:  ? typedef struct dist_tree dist_tree; */

typedef struct dist_node {                           /* subclass of ml_node */
  struct ml_node node;                        /* Base object, must be first */
  vector d, w;
  double t;
  boolean sametime;
  double dist;    /* debug: contml? (no) */   /* dist used in fitch, contml */
} dist_node;

#ifndef OLDC
/*function prototypes*/
void dist_node_init(dist_node* n, node_type type, long, long);
void dist_tree_init(struct dist_tree**, long, long);
void dist_tree_new(struct dist_tree**, long, long, int);
void dist_node_copy(node* src, node* dst);
void dist_node_free(node **np);
void alloctree(tree *, long);
void freetree(tree *, long);
void allocd(long, pointptr);
void freed(long, pointptr);
void allocw(long, pointptr);
void freew(long, pointptr);
void inputdata(boolean, boolean, boolean, boolean, vector *, intvector *);
void coordinates(node *, double, long *, double *, node *);
void drawline(long, double, node *, boolean);
void printree(node *, boolean, boolean);
void treeoutr(node *, long *, tree *);
void treeout(node *, long *, double, boolean, node *);
dist_node* dist_node_new(node_type type, long index);
/*function prototypes*/
#endif


/* End. */
