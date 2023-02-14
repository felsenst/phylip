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
#endif

#ifndef PHYLIP_H
#include "phylip.h"
#endif

extern boolean inserting, smoothit, smoothed, polishing;

struct allocx_t;
struct freex_t;

typedef struct bl_tree {
  struct tree treepart;
} bl_tree;

typedef struct bl_node {                       /* subclass of generic node */
  struct node node;                          /* Base object, must be first */
} bl_node;

typedef void (*allocx_t)(struct bl_node*, long, long);
typedef void (*freex_t)(bl_node*);
typedef void (*makenewv_t)(tree*, bl_node*);
typedef void (*nuview_t)(tree*, bl_node*);

typedef void (*initialvtrav_t)(tree*, bl_node*);

long endsite;

#ifndef OLDC /* prototypes */
void    bl_tree_new(struct tree**, long, long, long);
bl_node* bl_node_new(node_type, long, long);
void    bl_node_copy(struct bl_node *, struct  bl_node *);
void    bl_node_init(struct node*, node_type, long);
void    bl_node_free(bl_node **);
void    bl_node_print(bl_node *);
void    bl_hookup(bl_node*, bl_node*);
void    allocx(long, long, long, bl_node**);
void    makevalues2(long, pointarray, long, long, sequence, steptr);
void    freex_notip(long, pointarray);
void    freex(long, pointarray);
void    bl_update(tree*, bl_node *);
void    smooth(tree*, bl_node *);
void    smooth_traverse(tree*, bl_node *);
void 	bl_tree_smoothall(tree*, bl_node*);
void 	bl_node_reinit(bl_node * n);
void    bl_tree_insert_(tree*, bl_node*, bl_node*, boolean);
void    bl_tree_re_move(tree*, bl_node*, bl_node**, boolean);
boolean bl_tree_try_insert_(tree* , bl_node* , bl_node* , bl_node* , double*, 
		tree*, boolean, boolean, boolean, double*);
boolean bl_tree_try_insert_thorough(tree*, bl_node*, bl_node*, bl_node*, 
                          double*, tree*, boolean, boolean, boolean);
void    bl_tree_do_branchl_on_insert(tree*, bl_node *, bl_node*);
void    bl_tree_do_branchl_on_re_move(tree*, bl_node*, bl_node*);
double  min_child_tyme(bl_node *);
double  parent_tyme(bl_node *);
boolean valid_tyme(tree *, bl_node *, double);
void    bl_treeoutrecurs(FILE*, tree*, bl_node*, double, int*);
void    bl_treeout(FILE*, tree*, bl_node*, double);
void    bl_treevaluate(tree*, boolean, boolean, boolean, boolean, tree*,
                        tree*, initialvtrav_t);
void    bl_initialvtrav(tree*, bl_node *);
#endif

/* End. */
