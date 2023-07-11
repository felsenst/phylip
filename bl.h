/* Version 4.0. (c) Copyright 1993-2023 by the University of Washington.
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

typedef void (*tree_save_lr_nodes_t)(tree*,node*,node*);
typedef void (*tree_restore_lr_nodes_t)(tree*,node*,node*);

typedef struct bl_tree {
  struct tree treepart;
  long endsite;
  tree_save_lr_nodes_t save_lr_nodes;
  tree_restore_lr_nodes_t restore_lr_nodes;
  tree_save_traverses_t save_traverses;
} bl_tree;

typedef struct bl_node {                       /* subclass of generic node */
  struct node node;                          /* Base object, must be first */
  double v, tyme, deltav, oldlen, ssq;        /* ssq used only in contrast */
  boolean iter;                       /* iter used in dnaml, fitch & restml */
} bl_node;

typedef void (*makenewv_t)(struct bl_tree*, struct bl_node*);
typedef void (*nuview_t)(struct bl_tree*, struct bl_node*);
typedef void (*initialvtrav_t)(struct bl_tree*, struct bl_node*);

extern boolean smoothed;

#ifndef OLDC /* prototypes */
void    bl_tree_new(struct bl_tree**, long, long, long);
void    bl_tree_init(struct bl_tree*, long, long);
struct  bl_node* bl_node_new(node_type, long, long);
boolean bl_node_good(struct tree*, struct node*);
void    bl_node_copy(struct bl_node *, struct  bl_node *);
void    bl_node_init(struct bl_node*, node_type, long);
void    bl_node_free(struct bl_node **);
void    bl_node_print(struct bl_node *);
void    bl_hookup(struct bl_node*, struct bl_node*);
void    allocx(long, long, long, struct bl_node**);
double  get_tyme(struct bl_node *);
void    set_tyme(struct bl_node*, double);
void    bl_update(struct bl_tree*, struct bl_node *);
void    smooth(struct bl_tree*, struct bl_node *);
void    smooth_traverse(struct bl_tree*, bl_node *);
void 	bl_tree_smoothall(struct bl_tree*, bl_node*);
void 	bl_node_reinit(struct bl_node *);
void    bl_tree_insert_(struct bl_tree*, struct bl_node*, 
                          struct bl_node*, boolean);
void    unrooted_tree_save_lr_nodes(struct tree*, struct node*, struct node*);
void    unrooted_tree_restore_lr_nodes(struct tree*, struct node*, 
		                         struct node*);
void    blk_tree_makenewv(struct tree*, struct node*);
void    bl_tree_re_move(struct bl_tree*, struct bl_node*, 
		          struct bl_node**, boolean);
void    blk_tree_insert_(struct bl_tree*, struct bl_node*, struct bl_node*, 
                           boolean , boolean );
void    blk_tree_re_move(struct bl_tree*, struct bl_node *, struct bl_node**,
                           boolean);
boolean bl_tree_try_insert_(struct bl_tree*, struct bl_node*, struct bl_node*,
                              struct bl_node*, double*, struct bl_tree*, 
                              boolean, boolean, boolean, double*);
boolean bl_tree_try_insert_thorough(struct bl_tree*, struct bl_node*, 
                                      struct bl_node*, struct bl_node*, 
                                      double*, struct bl_tree*, 
                                      boolean, boolean, boolean);
void    bl_tree_save_traverses(struct tree*, struct node*, struct node*);
void    bl_tree_restore_traverses(struct tree*, struct node*, struct node*);
void    bl_tree_do_branchl_on_insert(struct bl_tree*, bl_node *, bl_node*);
void    bl_tree_do_branchl_on_re_move(struct bl_tree*, bl_node*, bl_node*);
double  min_child_tyme(struct bl_node *);
double  parent_tyme(struct bl_node *);
boolean valid_tyme(struct bl_tree *, struct bl_node *, double);
double  set_tyme_evaluate(struct bl_tree*, struct bl_node*, double);
void    addelement2(struct tree*, struct node*, Char*, long*, FILE*, boolean, 
	  double*, boolean*, long*, long*, long, boolean*, boolean, long);
void    treeread2 (struct tree*, FILE*, struct node**, boolean, 
          double*, boolean*, boolean*, long*, boolean, long);
void    bl_treeoutrecurs(FILE*, struct tree*, bl_node*, double, int*);
void    bl_treeout(FILE*, struct tree*, struct bl_node*, double);
void    getthree(struct tree*, struct node*, double, 
                   double, double, double*, double*);
void    bl_treevaluate(struct tree*, boolean, boolean, boolean, 
		         boolean, struct tree*, 
                         struct tree*, initialvtrav_t);
void    bl_initialvtrav(struct tree*, bl_node *);
#endif

#endif

/* end of  #ifndef  that conditions on this header file not already used */

/* End. */

