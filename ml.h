/* Version 4.0. (c) Copyright 1993-2025.
   Written by Michal Palczewski and Joe Felsenstein
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* debug:  would it be better to split ml_tree / ml_node into two levels of
 * the hierarchy, so ml_node / ml_tree is for cases where the branch lengths
 * are iterated, and the codon_node / prot_node / dna_node are for cases
 * where it not only does that, it has molecular sequences
 * maybe call these  iterate.c / iterate.h  and  sequence.c / sequence.h ? */ 

#ifndef ML_H
#define ML_H

#ifndef BL_H
#define BL_H
#include "bl.h"
#endif

/* debug: extern boolean inserting, smoothit, polishing; */

typedef struct ml_tree {
  struct bl_tree bl_tree;
} ml_tree;

typedef struct ml_node {                      /* subclass of generic bl_node */
  struct bl_node bl_node;                         /* base object, must be first */
  double* underflows;
  long endsite;
  long categs;
} ml_node;

typedef void (*allocx_t)(long, long, long, struct ml_node*);

#ifndef OLDC /* prototypes */
void    ml_tree_new(struct tree **, long, long, long);
void    ml_tree_init(struct tree *, long, long);
void    ml_node_init(struct node*, node_type, long);
struct  node* ml_node_new(node_type, long, long);
void    ml_node_copy(struct node *, struct node *);
void    ml_node_free(struct node **);
void    ml_node_print(struct node *);
#endif

#endif

/* end of conditional compilation if ML_H undefined */

/* End.*/
