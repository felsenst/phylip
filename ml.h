/* Version 4.0. (c) Copyright 1993-2022 by the University of Washington.
   Written by Michal Palczewski and Joe Felsenstein
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* debug:  would it be better to split ml_tree / ml_node into two levels of
 * the hierarchy, so ml_node / ml_tree is for cases where the branch lengths
 * are iterated, and the codon_node / prot_node / dna_node are for cases
 * where it not only does that, it has molecular sequences
 * maybe call these  iterate.c / iterate.h  and  sequence.c / sequence.h ? */ 

#ifndef _ML_H_
#define _ML_H_
#endif


#include "bl.h"

extern boolean inserting, smoothit, smoothed, polishing;

struct allocx_t;
struct freex_t;

typedef void (*allocx_t)(struct node*, long, long);
typedef void (*freex_t)(node*);

typedef struct ml_tree {
  struct tree bl_tree;
} ml_tree;

typedef struct ml_node {                       /* subclass of generic node */
  struct node bl_node;                          /* Base object, must be first */
  allocx_t allocx;
  freex_t freex;    /* debug: stuff after here to be later moved to mldna.c ? */
  double* underflows;
  long endsite;
  long categs;
} ml_node;

typedef void (*makenewv_t)(tree *, node *);
typedef void (*nuview_t)(tree *, node *);
typedef void (*initialvtrav_t)(tree *, node *);

long endsite;

#ifndef OLDC /* prototypes */
void    ml_tree_new(struct tree **, long, long, long);
void    ml_tree_init(struct tree *, long, long);
void    ml_node_init(struct node*, node_type, long);
node*   ml_node_new(node_type, long, long);
void    ml_node_copy(node *, node *);
void    ml_node_free(node **);
void 	ml_node_reinit(node *);
void    ml_node_print(node *);
#endif

/* End.*/
