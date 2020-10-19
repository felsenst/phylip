/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Michal Palczewski
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifndef _ML_H_
#define _ML_H_

#include "phylip.h"

extern boolean inserting;

typedef struct ml_node ml_node;

typedef void (*allocx_t)(ml_node*, long, long);
typedef void (*freex_t)(ml_node*);
typedef void (*makenewv_t)(tree*, node*);
typedef void (*nuview_t)(tree*, node*);

typedef struct ml_tree {
  tree tree;
  nuview_t nuview;
  makenewv_t makenewv;
} ml_tree;

struct ml_node {
  node node;            /* Base object, must be first */
  allocx_t allocx;
  freex_t freex;
  double* underflows;
  long endsite;
  long categs;
};

typedef struct codon_node {
  ml_node ml_node;      /* Base object, must be first */
  cphenotype codonx;
} codon_node;

typedef struct prot_node {
  ml_node ml_node;      /* Base object, must be first */
  pphenotype x;
} prot_node;

typedef struct dna_node{
  ml_node ml_node;      /* Base object, must be first */
  phenotype x;
} dna_node;

typedef void (*initialvtrav_t)(tree*, node*);

#ifndef OLDC /* prototypes */
node *  dna_node_new(node_type, long);
void    dna_node_init(node *, node_type, long);
void    ml_node_init(node *, node_type, long);
void    ml_node_free(node **);
void    ml_node_print(node *);
node *  prot_node_new(node_type, long);
void    prot_node_init(node *, node_type, long);
void    prot_node_allocx(ml_node*, long, long);
node *  codon_node_new(node_type, long);
void    codon_node_init(node *, node_type, long);
void    codon_node_allocx(ml_node*, long, long);
void    dna_node_allocx(ml_node*, long, long);
void    ml_tree_init(tree*, long, long);
void    ml_node_copy(node*, node*);
void    ml_hookup(node*, node*);
void    allocx(long, long, long, ml_node**);
void    fix_x(dna_node*, long, double, long);
void    fix_protx(prot_node*, long, double, long);
void    prot_node_copy(node*, node*);
void    codon_node_copy(node*, node*);
void    dna_node_copy(node*, node*);
void    codon_node_freex(ml_node*);
void    prot_node_freex(ml_node*);
void    dna_node_freex(ml_node*);
void    makevalues2(long, pointarray, long, long, sequence, steptr);
void    codon_freex_notip(long, pointarray);
void    prot_freex_notip(long, pointarray);
void    freex_notip(long, pointarray);
void    freex(long, pointarray);
void    ml_update(tree*, node *);
void    smooth(tree*, node *);
void    smooth_traverse(tree*, node *);
void    ml_tree_insert_(tree*, node*, node*, boolean);
void    ml_tree_re_move(tree*, node*, node**, boolean);
boolean ml_tree_try_insert_(tree* , node* , node* , node* , double*, tree*,
                            boolean, boolean, boolean, double*);
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


// End.
