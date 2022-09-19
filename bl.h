/* Version 4.0. (c) Copyright 1993-2022 by the University of Washington.
   Written by Michal Palczewski and Joe Felsenstein
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* debug:  would it be better to split ml_tree / ml_node into two levels of
 * the hierarchy, so ml_node / ml_tree is for cases where the branch lengths
 * are iterated, and the codon_node / prot_node / dna_node are for cases
 * where it not only does that, it has molecular sequences
 * maybe call these  iterate.c / iterate.h  and  sequence.c / sequence.h ? */ 

#ifndef _BL_H_
#define _BL_H_

#include "phylip.h"

extern boolean inserting, smoothit, smoothed, polishing;

struct allocx_t;
struct freex_t;

typedef void (*allocx_t)(struct node*, long, long);
typedef void (*freex_t)(node*);

typedef struct bl_tree {
  struct tree treepart;
} ml_tree;

typedef struct bl_node {                       /* subclass of generic node */
  struct node node;                          /* Base object, must be first */
} ml_node;

typedef void (*makenewv_t)(tree*, node*);
typedef void (*nuview_t)(tree*, node*);

typedef void (*initialvtrav_t)(tree*, node*);


long endsite;

#ifndef OLDC /* prototypes */
void    bl_tree_new(struct tree**, long, long, long);
void    bl_tree_init(struct tree*, long, long);
node*   bl_node_new(node_type, long, long);
void    bl_node_init(struct node *, node_type, long);
void    bl_node_copy(struct node *, struct  node *);
void    bl_node_free(node **);
void    bl_hookup(node *, node *);
void    bl_node_reinit(node *);
void    bl_node_print(node *);
void    bl_update(tree *, node *);
void    smooth_traverse(tree *, node *);
void    smooth(tree *, node %*);
void 	bl_tree_smoothall(tree *, node *);
void    bl_tree_do_branchl_on_insert(tree *, node *, node *);
void    bl_tree_insert_(tree *, node *, node *, boolean);
void    bl_tree_do_branchl_on_re_move(tree *, node *, node *);
void    bl_tree_re_move(tree *, node *, node **, boolean);
boolean bl_tree_try_insert_thorough(tree *, node *, node *, node *, 
                          double *, tree *, boolean, boolean, boolean);
boolean bl_tree_try_insert_(tree * , node * , node * , node * , double *,
                             tree *, boolean, boolean, boolean, double*);
void    blk_tree_insert_(tree *, node *, node *, boolean, boolean);
double  get_tyme(node *);
void    set_tyme (node *, double) ;
void    blk_tree_re_move(tree *, node *, node **, boolean);
double  min_child_tyme(node *);
double  parent_tyme(node *);
boolean valid_tyme(tree *, node *, double);
double  set_tyme_evaluate(tree *, node *, double);
void    blk_tree_makenewv(tree *, node *);
void    getthree(tree*, node *, double, double, double, double *, double *);
void    bl_treevaluate(tree *, boolean, boolean, boolean, boolean, tree *,
                        tree*, initialvtrav_t);
void    bl_initialvtrav(tree *, node *);
void    bl_treeoutrecurs(FILE *, tree *, node *, double, int *);
void    bl_treeout(FILE *, tree *, node *, double);
#endif

#endif /* _BL_H_ */


/* End. */
