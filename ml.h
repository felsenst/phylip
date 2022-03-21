/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Michal Palczewski
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* debug:  would it be better to split ml_tree / ml_node into two levels of
 * the hierarchy, so ml_node / ml_tree is for cases where the branch lengths
 * are iterated, and the codon_node / prot_node / dna_node are for cases
 * where it not only does that, it has molecular sequences
 * maybe call these  iterate.c / iterate.h  and  sequence.c / sequence.h ? */ 

#ifndef _ML_H_
#define _ML_H_

#include "phylip.h"

extern boolean inserting, smoothit, smoothed, polishing;

struct allocx_t;
struct freex_t;

typedef void (*allocx_t)(node*, long, long);
typedef void (*freex_t)(node*);

typedef struct ml_tree {
  struct tree treepart;
} ml_tree;

typedef struct ml_node {                       /* subclass of generic node */
  struct node node;                          /* Base object, must be first */
  allocx_t allocx;
  freex_t freex;
  double* underflows;
  long endsite;
  long categs;
} ml_node;

typedef void (*makenewv_t)(tree*, node*);
typedef void (*nuview_t)(tree*, node*);

typedef void (*initialvtrav_t)(tree*, node*);

#ifndef OLDC /* prototypes */
void    ml_tree_new(struct tree**, long, long, long);
void    ml_tree_init(struct tree*, long, long);
node*   ml_node_new(node_type, long, long);
void    ml_node_init(struct node*, node_type, long);
void    ml_node_free(node **);
void    ml_node_print(node *);
void    ml_node_copy(node*, node*);
void    ml_hookup(node*, node*);
void    allocx(long, long, long, ml_node**);
void    makevalues2(long, pointarray, long, long, sequence, steptr);
void    freex_notip(long, pointarray);
void    freex(long, pointarray);
void    ml_update(tree*, node *);
void    smooth(tree*, node *);
void    smooth_traverse(tree*, node *);
void 	ml_tree_smoothall(tree*, node*);
void 	ml_node_reinit(node * n);
void    ml_tree_insert_(tree*, node*, node*, boolean);
void    ml_tree_re_move(tree*, node*, node**, boolean);
boolean ml_tree_try_insert_(tree* , node* , node* , node* , double*, tree*,
                            boolean, boolean, boolean, double*);
boolean ml_tree_try_insert_thorough(tree*, node*, node*, node*, 
                          double*, tree*, boolean, boolean, boolean);
void    ml_tree_do_branchl_on_insert(tree*, node *, node*);
void    ml_tree_do_branchl_on_re_move(tree*, node*, node*);
void    mlk_tree_insert_(tree*, node *, node *, boolean, boolean);
double  get_tyme(node *);
void    set_tyme (node*, double) ;
void    mlk_tree_re_move(tree* t, node *item, node** where, boolean recompute);
void    getthree(tree*, node *, double, double, double, double *, double *);
double  min_child_tyme(node *);
double  parent_tyme(node *);
boolean valid_tyme(tree *, node *, double);
double  set_tyme_evaluate(tree *, node *, double);
void    mlk_tree_makenewv(tree*, node *);
void    empiricalfreqs(double *, double *, double *, double *, steptr, pointarray);
void    ml_treevaluate(tree*, boolean, boolean, boolean, boolean, tree*,
                        tree*, initialvtrav_t);
void    ml_initialvtrav(tree*, node *);
#endif

#endif /* _ML_H_ */


/* End.*/
