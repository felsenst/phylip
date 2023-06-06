/* Version 4.0.  Copyright, 2022.
   Written by Michal Palczewski and Joe Felsenstein */

/* These are versions of functions, and support functions, for programs
 * computing likelihoods. Functions that infer branch lengths on the tree
 * but are not specific to likelihood inference are instead in bl.c */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <assert.h>

#include "bl.h"
#include "ml.h"

#define DEBUG
#define MAKENEWV_DEBUG
/* #define USE_NEW_MAKENEWV */

const double MIN_BRANCH_LENGTH = epsilon/4.0;
const double MIN_ROOT_TYME = -10;

/* TODO check to see which of these are needed here */
long endsite;                           // RSGdebug: Check this.
extern long nextree, which;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
extern steptr weight, category, alias, location, ally;
extern sequence inputSequences;
extern struct node** lrsaves;
extern long rcategs;
extern boolean usertree, lngths, smoothit, smoothed, polishing;
boolean inserting;


void ml_tree_new(struct ml_tree **tp, long nonodes, long spp, long treesize)
{ /* make a new ml_tree.  Calls to generic_tree_new,
   * casting ml_tree** to tree** as we call it 
   * then call  ml_tree_init */
  struct bl_tree **bltp;

  bltp = (struct bl_tree**)tp;
  bl_tree_new(bltp, nonodes, spp, treesize);      /* next up tree hierarchy */
  ml_tree_init(*tp, nonodes, spp);
} /* ml_tree_new */


void ml_tree_init(struct ml_tree* t, long nonodes, long spp)
{ /* set up function variables in ml_tree.  Currently these are actually
   * attributes of the generic tree that need ml function versions */
  struct bl_tree *blt;

  blt = (struct bl_tree*)t;
  bl_tree_init(blt, nonodes, spp);                 /* go up class hierarchy */
/* debug: need here?   ((ml_tree*)t)->nuview = ml_tree_nuview;
  (t.tree)->makenewv_t = ml_tree->makenewv_t;
 * */
} /* ml_tree_init */


struct ml_node* ml_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct ml_node* n;

  n = (struct ml_node*)bl_node_new(type, index, nodesize);
  ml_node_init(n, type, index);
  return n;
} /* ml_node_new */


void ml_node_init(struct ml_node *n, node_type type, long index)
{
  /* initialize a node for ml trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */
  struct node* nn;
  struct bl_node* bn;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  bn = (struct bl_node*)n;
  nn = (struct node*)n;
  generic_node_init(nn, type, index);                /* go up node hierarchy */
  nn->node_print_f = (node_print_t)bl_node_print;
  bn->tyme = 0;
} /* ml_node_init */


void ml_node_copy(ml_node* src, ml_node* dest)
{ /* copy contents of an ml_node but not its pointers */
  bl_node *srcb = (bl_node *)src;
  bl_node *destb = (bl_node *)dest;
  bl_node_copy(srcb, destb);                              /* go up hierarchy */
  dest->categs = src->categs;
  dest->endsite = src->endsite;
  set_tyme((bl_node*)dest, ((bl_node*)src)->tyme);

  if(dest->underflows)                  // RSGbugfix
    memcpy(dest->underflows, src->underflows, src->endsite * sizeof(double));
  else
    assert(src->underflows == NULL);    // RSGdebug
} /* ml_node_copy */


void ml_node_free(struct ml_node **np)
{
  /* free a node for ml trees */
  struct ml_node *n = (ml_node*)*np;

  n->freex((ml_node*)n);
  generic_node_free((struct node**)np);
} /* ml_node_free */


void ml_node_reinit(struct ml_node *n)
{
  /* reset things for an ml tree node */
  bl_node * bln = (bl_node*)n;

  bln->tyme = 0;
  // BUG.970 -- does freex need refreshing ?
  // BUG.970 -- leave for dna_node and prot_node ?
  bl_node_reinit(bln);                               /* go up the hierarchy */
} /* ml_node_reinit */


void ml_node_print(struct ml_node * n)
{
  /* for debugging */
  struct bl_node * bn = (bl_node*)n;

  bl_node_print(bn);
  /* debug:  ?? printf(" ml(bn.endsite:%ld tyme:%lf)", ((struct bl_tree*)mn)->endsite, mn->bl_node.tyme); */
} /* ml_node_print */

/* End. */

