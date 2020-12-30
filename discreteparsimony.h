/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/*
  discrete.h: included in pars
*/


#include "parsimony.h"

typedef struct discretepars_node {
  pars_node pars_node;
  discnucarray discnumreconst;   /* bookkeeps number of  reconstructions   */
  discnucarray disccumlengths;   /* bookkeeps cummulative minimum lengths  */
  discbaseptr discbase;          /* the sequence in pars                   */
  discbaseptr olddiscbase;       /* record previous sequence               */
  discnucarray *discnumnuc;      /* bookkeeps number of nucleotides        */
} discretepars_node;

typedef struct discretepars_tree {
  pars_tree pars_tree;
} discretepars_tree;

struct LOC_hyptrav {
  boolean bottom;
  node *r;
  discbaseptr hypset;
  boolean maybe, nonzero;
  unsigned char tempset, anc;
};

extern long nonodes, endsite, outgrno, which;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff;
extern steptr weight, category, alias, location, ally;
extern sequence inputSequences;
extern sequence convtab;

#ifndef OLDC
/*function prototypes*/
tree*  discretepars_tree_new(long nonodes, long spp);
void   discretepars_tree_init(tree*, long nonodes, long spp);
void   discretepars_tree_nuview(tree* t, node*p);
double discretepars_tree_evaluate(tree* t, node *n, boolean saveit);
boolean discretepars_tree_branchcollapsible(tree* t, node* n);

node*  discretepars_node_new(node_type type, long index);
void   discretepars_node_init(node *n, node_type type, long index);
void   discretepars_node_free(node **np);

void   inputdata(long);
void   sitesort(long, steptr);
void   sitecombine(long);
void   sitescrunch(long);
void   makevalues(tree*, boolean);
long   discgetlargest(long *);
void   findoutgroup(node *, boolean *);
void   dischyptrav(tree* t, node *r_, discbaseptr hypset_, long b1, long b2, boolean bottom_);
void   disc_hypstates(tree*, long);
void   discinitbranchlen(node *);
void   discinitbase(node *, long);
void   discinittreetrav(node *, long);
void   discminpostorder(node *, pointarray);
void   discbranchlength(node *, node *, double *, pointarray);
void   discbranchlentrav(node *, node *, long, long, double *, pointarray);
void   drawline(long, double, node *);
void   treeout(node *, long, long *, node *);
void   standev(long, long, long, double, double *, long **, longer);
void   initdiscmin(discretepars_node *p, long sitei, boolean internal);
void   discmultifillin(node *p, node *q, long dnumdesc);
void   disccompmin(node *p, node *desc);
void   discinitmin(discretepars_node *p, long sitei, boolean internal);
void   dischyprint(tree* t, long b1, long b2, struct LOC_hyptrav *htrav);
/*function prototypes*/
#endif


// End.
