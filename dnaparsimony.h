/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Michal Palczewski.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#include "parsimony.h"


typedef struct dnapars_node {
  pars_node pars_node;
  nucarray cumlengths;           /* bookkeeps cumulative minimum lengths    */
  nucarray numreconst;           /* bookkeeps number of  reconstructions    */
  nucarray *numnuc;              /* bookkeeps number of nucleotides         */
  baseptr base;                  /* the sequence in dnapars/comp/penny      */
} dnapars_node;


typedef struct dnapars_tree {
  pars_tree pars_tree;
} dnapars_tree;

struct LOC_hyptrav {
  boolean bottom;
  node *r;
  long *hypset;
  boolean maybe, nonzero;
  long tempset, anc;
} ;


double  dnapars_tree_evaluate(tree* t, node *r, boolean dummy);
void    dna_initmin(dnapars_node *p, long sitei, boolean internal);
void    dna_compmin(node *, node *);
void    dna_treelength(node *, long, pointarray);
void    dna_initbase(node *, long);
long    dna_getlargest(long *);
void    dna_hyprint(tree*, long, long, struct LOC_hyptrav *, Char *);
void    dna_hyptrav(tree* t, node *, long *, long, long, boolean, Char *);
void    dna_hypstates(tree*, long, Char *);
node*   dnapars_node_new(node_type type, long index);
void    dnapars_node_init(node* p, node_type type, long index);
void    dnapars_node_reinit(node* p);
void    dnapars_node_free(node **pp);
void    dnapars_node_copy(node* src, node* dst);
void    dna_makevalues(tree* t, boolean usertree);
boolean dna_branchcollapsible(tree* t, node* n);
void    dnapars_tree_nuview(tree* t, node* p);
tree*   dnapars_tree_new(long nonodes, long spp);
void    dnapars_tree_init(tree*, long, long);


// End.
