/* version 3.6. (c) Copyright 1993-2010 by the University of Washington.
   Written by Michal Palczewski
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include "phylip.h"

extern boolean inserting;

typedef struct ml_node ml_node;

typedef void (*allocx_t)(ml_node*, long, long);
typedef void (*freex_t)(ml_node*);
typedef void (*makenewv_t)(tree*, node*); 
/* EWFIX.REMOVE
typedef void (*do_branchl_on_insert_t)(tree*, node *, node*);
typedef void (*do_branchl_on_re_move_t)(tree* ,node*, node*);
*/

typedef struct ml_tree {
  tree tree;
  makenewv_t makenewv;
  /* EWFIX.REMOVE
  do_branchl_on_insert_t do_branchl_on_insert;
  do_branchl_on_re_move_t do_branchl_on_re_move;
  */
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


typedef void (*inittravtree_t)(tree*, node*);

#ifndef OLDC /* prototypes */
node  * ml_node_new(boolean,long);
node *  dna_node_new(node_type type, long index);
void    dna_node_init(node *n, node_type type, long index);
void    ml_node_init(node*, node_type, long);
void    ml_node_allocx(ml_node* n,long);
void    ml_node_free(node **np);
void    ml_node_print(node *n);
void    prot_node_init(node *n, node_type type, long index);
node *  prot_node_new(node_type type, long index);
void    prot_node_allocx(ml_node* n,long endsite,long rcategs);
void    codon_node_init(node *n, node_type type, long index);
node *  codon_node_new(node_type type, long index);
void    codon_node_allocx(ml_node* n,long endsite,long rcategs);
void    dna_node_allocx(ml_node* n, long endsite, long rcategs);
void    ml_tree_allocx(long nonodes, long endsite, long param, 
    ml_node** treenode, boolean usertree);
ml_tree* ml_tree_new(long nonodes, long spp);
void    ml_tree_init(tree* t, long nonodes, long spp);
void    ml_node_copy(node* src, node* dest);
void    allocx(long nonodes, long endsite, long param, ml_node** treenode);
void    ml_unrooted_tree_locrearrange(tree*, node*,boolean,tree* , tree* );
void    fix_x(dna_node* p,long site, double maxx, long rcategs);
void    fix_protx(prot_node* p,long site,double maxx, long rcategs);
void    fix_codonx(codon_node* p,long site,double maxx, long rcategs);
void    prot_node_copy(node* src, node* dest);
void    codon_node_copy(node* src, node* dest);
void    dna_node_copy(node* src, node* dest);
void    codon_node_freex(ml_node* n);
void    prot_node_freex(ml_node* n);
void    dna_node_freex(ml_node* n);
void    makevalues2(long, pointarray, long, long, sequence, steptr);
void    codon_freex_notip(long nonodes, pointarray treenode);
void    codon_freex(long nonodes, pointarray treenode);
void    prot_freex_notip(long nonodes, pointarray treenode);
void    prot_freex(long nonodes, pointarray treenode);
void    freex_notip(long, pointarray);
void    freex(long, pointarray);
void    update(tree*, node *p);
void    smooth(tree*, node *p);
void    ml_tree_insert_(tree* t, node*p, node*q,boolean dooinit,boolean multf);
void    ml_tree_re_move(tree* t, node*p, node** q,boolean);
boolean ml_tree_try_insert_(tree* ,node* ,node* ,node** ,double*, tree *,tree*,
    boolean, boolean* );
void    ml_tree_do_branchl_on_insert(tree*t, node *p, node* q);
void    ml_tree_do_branchl_on_re_move(tree* t, node* p, node*q);
void    mlk_tree_insert_(tree *t,node *newtip, node *below, boolean,boolean);
double  get_tyme(node *p);
void    set_tyme (node* p,double tyme) ;
void    mlk_tree_re_move(tree* t,node *item, node** where, boolean recompute);
void    getthree(tree* t,node *p, double thigh, double tlow,double tdelta,
    double *x,double *lnl);
double  min_child_tyme(node *p);
double  parent_tyme(node *p);
boolean valid_tyme(tree *t, node *p, double tyme);
double  set_tyme_evaluate(tree *t, node *p, double tyme);
void    mlk_tree_makenewv(tree* t, node *p);
void    empiricalfreqs(double *,double *,double *,double *,steptr,pointarray);
void    ml_treevaluate(tree* curtree, boolean improve, boolean reusertree, boolean global,
    boolean progress, tree* priortree, tree* bestree, inittravtree_t inittravtree);
void ml_inittravtree(tree* t,node *p);
#endif
