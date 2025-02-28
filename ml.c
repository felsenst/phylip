/* Version 4.0.  Copyright, 2022.
   Written by Michal Palczewski and Joe Felsenstein */

/* These are versions of functions, and support functions, for programs
 * computing likelihoods. Functions that infer branch lengths on the tree
 * but are not specific to likelihood inference are instead in bl.c */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>

#ifndef BL_H
#include "bl.h"
#define BL_H
#endif

#ifndef ML_H
#include "ml.h"
#define ML_H
#endif

#define DEBUG
#define MAKENEWV_DEBUG
/* #define USE_NEW_MAKENEWV */

extern const double MIN_BRANCH_LENGTH;
extern const double MIN_ROOT_TYME;

/* TODO check to see which of these are needed here */
extern long nextree, which;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
extern steptr weight, category, alias, location, ally;
extern sequence inputSequences;
extern struct node** lrsaves;
long rcategs;                    /* number of rate categories, default is 1 */


extern allocx_t allocx_f;


void ml_tree_new(struct tree **tp, long nonodes, long spp, long treesize)
{ /* make a new ml_tree.  Calls to generic_tree_new,
   * casting ml_tree** to tree** as we call it 
   * then call  ml_tree_init */

  bl_tree_new(tp, nonodes, spp, treesize);      /* next up tree hierarchy */
  ml_tree_init(*tp, nonodes, spp);
} /* ml_tree_new */


void ml_tree_init(struct tree* t, long nonodes, long spp)
{ /* set up function variables in ml_tree.  Currently these are actually
   * attributes of the generic tree that need ml function versions */

  bl_tree_init(t, nonodes, spp);                 /* go up class hierarchy */
/* debug: need here?   ((ml_tree*)t)->nuview = ml_tree_nuview;
  (t.tree)->makenewv_t = ml_tree->makenewv_t;
 * */
} /* ml_tree_init */


struct node* ml_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct node *n;

  n = bl_node_new(type, index, nodesize);
  return n;
} /* ml_node_new */


void ml_node_init(struct node *n, node_type type, long index)
{
  /* initialize a node for ml trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */
  long i;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  bl_node_init(n, type, index);
  n->node_print_f = (node_print_t)ml_node_print;
  for (i = 0; i < ((struct ml_node *)n)->endsite; i++)
    ((struct ml_node*)n)->underflows[i] = 0.0;
} /* ml_node_init */


void ml_node_copy(struct node* src, struct node* dest)
{ /* copy contents of an ml_node but not its pointers */
  ml_node *srcmln = (ml_node *)src;
  ml_node *destmln = (ml_node *)dest;

  bl_node_copy(src, dest);                              /* go up hierarchy */
  destmln->categs = srcmln->categs;
  destmln->endsite = srcmln->endsite;
  if(srcmln->underflows)                  // RSGbugfix
    memcpy(&(destmln->underflows), &(srcmln->underflows), 
		     sizeof(srcmln->underflows));
  else
    assert(destmln->underflows == NULL);    // RSGdebug
} /* ml_node_copy */


void ml_node_free(struct node **np)
{ /* free a node for ml trees */

/* debug:  something to free the data goes here */
  generic_node_free((struct node**)np);
} /* ml_node_free */


void ml_node_print(struct node * n)
{
  /* for debugging */

  bl_node_print(n);
  /* debug:  ?? printf(" ml(bn.endsite:%ld tyme:%lf)", ((struct bl_tree*)mn)->endsite, mn->bl_node.tyme); */
} /* ml_node_print */

/* End. */

