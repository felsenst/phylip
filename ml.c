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
#include "phylip.h"

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
extern node** lrsaves;
extern long rcategs;
extern boolean usertree, lngths, smoothit, smoothed, polishing;
boolean inserting;


void ml_tree_new(struct tree **tp, long nonodes, long spp, long treesize)
{ /* make a new ml_tree.  Calls to generic_tree_new,
   * casting ml_tree** to tree** as we call it 
   * then call  ml_tree_init */

  bl_tree_new(tp, nonodes, spp, treesize);        /* next up tree hierarchy */
} /* ml_tree_new */


void ml_tree_init(struct tree* t, long nonodes, long spp)
{ /* set up function variables in ml_tree.  Currently these are actually
   * attributes of the generic tree that need ml function versions */

  bl_tree_init(t, nonodes, spp);                   /* go up class hierarchy */
  t->smoothall = ml_tree_smoothall;
  t->insert_ = (tree_insert_t)ml_tree_insert_;
  t->re_move = ml_tree_re_move;
  t->try_insert_ = (tree_try_insert_t)ml_tree_try_insert_;
  t->do_branchl_on_insert_f = ml_tree_do_branchl_on_insert;
  t->do_branchl_on_re_move_f = ml_tree_do_branchl_on_re_move;
/* debug: need here?   ((ml_tree*)t)->nuview = ml_tree_nuview;
  (t.tree)->makenewv_t = ml_tree->makenewv_t;
 * */
} /* ml_tree_init */


void ml_node_init(struct node *n, node_type type, long index)
{
  /* initialize a node for ml trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */
  ml_node* nn;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  generic_node_init(n, type, index);                /* go up node hierarchy */
  n->copy = ml_node_copy;
  n->free = ml_node_free;
  n->reinit = ml_node_reinit;
  n->node_print_f = ml_node_print;
  nn = (ml_node*)n;
  nn->freex = NULL;         /* x is only defined for dna_node and prot_node */
  nn->node.tyme = 0;
} /* ml_node_init */


struct node* ml_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct node* nn;

  nn = bl_node_new(type, index, nodesize);
  ml_node_init(nn, type, index);
  return nn;
} /* ml_node_new */


void ml_node_copy(node* srcn, node* destn) // RSGbugfix
{ /* copy an ml_node */
  ml_node *src = (ml_node *)srcn;
  ml_node *dest = (ml_node *)destn;
  assert(srcn);                         // RSGdebug
  assert(destn);                        // RSGdebug
  bl_node_copy(srcn, destn);
  dest->categs = src->categs;
  dest->endsite = src->endsite;
  set_tyme((node*)dest, src->node.tyme);

  if(dest->underflows)                  // RSGbugfix
    memcpy(dest->underflows, src->underflows, src->endsite * sizeof(double));
  else
    assert(src->underflows == NULL);    // RSGdebug
} /* ml_node_copy */


void ml_node_free(node **np)
{
  /* free a node for ml trees */
  ml_node *n = (ml_node*)*np;
  n->freex((node*)n);
  generic_node_free(np);
} /* ml_node_free */


void ml_node_reinit(node * n)
{
  /* reset things for an ml tree node */
  ml_node * mln = (ml_node*)n;

  mln->node.tyme = 0;
  // BUG.970 -- does freex need refreshing ?
  // BUG.970 -- leave for dna_node and prot_node ?
  bl_node_reinit(n);                                 /* go up the hierarchy */
} /* ml_node_reinit */


void ml_node_print(node * n)
{
  /* for debugging */
  bl_node_print(n);
  ml_node * mn = (ml_node*)n;
  printf(" ml(endsite:%ld tyme:%lf)", mn->endsite, mn->node.tyme);
} /* ml_node_print */

/* End. */

