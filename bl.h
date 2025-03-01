/* Version 4.0. (c) Copyright 1993-2023 by the University of Washington.
   Written by Michal Palczewski and Joe Felsenstein
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* debug:  would it be better to split ml_tree / ml_node into two levels of
 * the hierarchy, so ml_node / ml_tree is for cases where the branch lengths
 * are iterated, and the codon_node / prot_node / dna_node are for cases
 * where it not only does that, it has molecular sequences
 * maybe call these  iterate.c / iterate.h  and  sequence.c / sequence.h ? */ 

#ifndef BL_H
#define BL_H

#ifndef PHYLIP_H
#include "phylip.h"
#endif

extern boolean inserting, smoothit, polishing;
extern FILE *infile, *outfile, *intree, *intree2, *outtree;
typedef void (*tree_restore_lr_nodes_t)(struct tree*, struct node*);
/* debug:   declare here? where?  extern boolean smoothed;   debug */

typedef struct bl_tree {
  struct tree treepart;
  tree_save_lr_nodes_t save_lr_nodes;
  tree_restore_lr_nodes_t restore_lr_nodes;
  tree_save_traverses_t save_traverses;
} bl_tree;

typedef struct bl_node {                        /* subclass of generic node */
  struct node node;                           /* base object, must be first */
  double v, tyme, deltav, oldlen, ssq;         /* ssq used only in contrast */
  boolean iter;                       /* iter used in dnaml, fitch & restml */
} bl_node;

typedef void (*makenewv_t)(struct tree*, struct node*);  /* debug: should these three be just tree*, node* ? */
typedef void (*nuview_t)(struct tree*, struct node*);
typedef void (*initialvtrav_t)(struct tree*, struct node*);

#ifndef OLDC /* prototypes */
void    bl_tree_new(struct tree**, long, long, long);
void    bl_tree_init(struct tree*, long, long);
struct  node* bl_node_new(node_type, long, long);
boolean bl_node_good(struct tree*, struct node*);
void    bl_node_copy(struct node*, struct node*);
void    bl_node_init(struct node*, node_type, long);
void    bl_node_free(struct node **);
void    bl_node_print(struct node *);
void    bl_hookup(struct node*, struct node*);
void    allocx(long, long, long, struct node**);
void    bl_update(struct tree*, struct node *);
void    smooth(struct tree*, struct node *);
void    smooth_traverse(struct tree*, node *);
void 	bl_tree_smoothall(struct tree*, node*);
void 	bl_node_reinit(struct node *);
void    bl_tree_insert_(struct tree*, struct node*, 
                          struct node*, boolean);
void    unrooted_tree_save_lr_nodes(struct tree*, struct node*);
void    unrooted_tree_restore_lr_nodes(struct tree*, struct node*);
void    blk_tree_makenewv(struct tree*, struct node*);
void    bl_tree_re_move(struct tree*, struct node*, struct node**, boolean);
void    blk_tree_insert_(struct tree*, struct node*, struct node*, 
                           boolean, boolean);
void    blk_tree_re_move(struct tree*, struct node *, struct node**,
                           boolean);
boolean bl_tree_try_insert_(struct tree*, struct node*, struct node*,
                              struct node*, double*, struct tree*, 
                              boolean, boolean, boolean, double*);
boolean bl_tree_try_insert_thorough(struct tree*, struct node*, 
                                      struct node*, struct node**, 
                                      double*, struct tree*, 
                                      boolean, boolean, boolean);
void    bl_tree_save_traverses(struct tree*, struct node*);
void    bl_tree_restore_traverses(struct tree*, struct node*);
void    bl_tree_do_branchl_on_insert(struct tree*, node *, node*);
void    bl_tree_do_branchl_on_re_move(struct tree*, node*, node*);
double  get_tyme(struct node *);
void    set_tyme(struct node*, double);
double  min_child_tyme(struct node *);
double  parent_tyme(struct node *);
boolean valid_tyme(struct tree *, struct node *, double);
double  set_tyme_evaluate(struct tree*, struct node*, double);
void    addelement2(struct tree*, struct node*, Char*, long*, FILE*, boolean, 
	  double*, boolean*, long*, long*, long, boolean*, boolean, long);
void    treeread2 (struct tree*, FILE*, struct node**, boolean, 
          double*, boolean*, boolean*, long*, boolean, long);
void    bl_treeoutrecurs(FILE*, struct tree*, node*, double, int*);
void    bl_treeout(FILE*, struct tree*, struct node*, double);
void    getthree(struct tree*, struct node*, double, 
                   double, double, double*, double*);
void    bl_treevaluate(struct tree*, boolean, boolean, boolean, 
		         boolean, struct tree*, 
                         struct tree*, initialvtrav_t);
void    bl_initialvtrav(struct tree*, bl_node *);
#endif

#endif

/* end of conditional compilation if BL_H initially undefined */

/* End. */

