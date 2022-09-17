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


struct node* ml_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct node* nn;

  nn = generic_node_new(type, index, nodesize);
  ml_node_init(nn, type, index);
  return nn;
} /* ml_node_new */


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


void ml_node_copy(node* srcn, node* destn) // RSGbugfix
{ /* copy an ml_node */
  ml_node *src = (ml_node *)srcn;
  ml_node *dest = (ml_node *)destn;
  assert(srcn);                         // RSGdebug
  assert(destn);                        // RSGdebug
  generic_node_copy(srcn, destn);
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
  generic_node_print(n);
  ml_node * mn = (ml_node*)n;
  printf(" ml(endsite:%ld tyme:%lf)", mn->endsite, mn->node.tyme);
} /* ml_node_print */


void getthree(tree* t, node *p, double thigh, double tlow, double tdelta, double *x, double *lnl)
{
  /* compute score at a new triple of points */
  int i;
  double tt = ((ml_node*)p)->node.tyme;
  double td = fabs(tdelta);

  x[0] = tt - td;
  x[1] = tt;
  x[2] = tt + td;

  if ( x[0] < tlow + epsilon )
  {
    x[0] = tlow + epsilon;
    x[1] = ( x[0] + x[2] ) / 2;
  }

  if ( x[2] > thigh - epsilon )
  {
    x[2] = thigh - epsilon;
    x[1] = ( x[0] + x[2] ) / 2;
  }

  for ( i = 0 ; i < 3 ; i++ )
  {
    set_tyme(p, x[i]);
    t->nuview(t, p);
    lnl[i] = t->evaluate(t, p, 0);
  }
}  /* getthree */


void ml_treevaluate(tree* curtree, boolean improve, boolean reusertree,
                    boolean global, boolean progress, tree* priortree,
                    tree* bestree, initialvtrav_t initialvtrav)
{
  double bestyet;
  /* evaluate a user tree */

  smoothit = improve;
  if (reusertree)
  {
    arbitrary_resolve(curtree);
    curtree->smoothall(curtree, curtree->root);
    if (global)
      curtree->globrearrange(curtree, bestree, progress, smoothit, &bestyet);
    else
      curtree->locrearrange(curtree, curtree->root->back, smoothit, &bestyet,
                             bestree, priortree, false, &bestyet);
    polishing = true;
    smoothit = true;
    curtree->smoothall(curtree, curtree->root);
    polishing = false;
  }
  else
  {
    if (!lngths) {
      inittrav(curtree, curtree->root);
      inittrav(curtree, curtree->root->back);
    }
    polishing = true;
    smoothit = true;
    curtree->evaluate(curtree, curtree->root, 0);     /* get current value */
    if (!lngths)
      curtree->smoothall(curtree, curtree->root);
    smoothit = improve;
    polishing = false;
  }
  curtree->evaluate(curtree, curtree->root, true);
}  /* ml_treevaluate */


void ml_initialvtrav(tree* t, node *p)
{
  /* traverse tree to set branch lengths  v  to initial values
   * must be called twice the first time, at both ends of
   * a branch such as the root branch.  Is separate from the
   * task of setting initialized booleans for views to false   */
  node* q;

  if ((!lngths) || p->iter)
  {
    p->v = initialv;
    p->back->v = initialv;
  }

  if (!p->tip)
  {
    q = p->next;
    while ( q != p )
    {
      ml_initialvtrav(t, q->back);
      q = q->next;
    }
  }
}  /* ml_initialvtrav */


void ml_treeoutrecurs(FILE* outtreefile, tree* t, node* p, double bl_scale, int* col)
{ 
  /* write out to output file a subtree, recursively.  This is the version 
   * with branch lengths and a scale factor,  bl_scale  */
  long i, n, w;
  Char c;
  double x;
  node *q, *qfirst;
  boolean inloop;

  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {       /* find out how long the name is */
      if (nayme[p->index-1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {                      /* ... then write it out */
      c = nayme[p->index-1][i];
      if (c == ' ')                          /* convert blank to underscore */
        c = '_';
      putc(c, outtree);
    }
    (*col) += n;                  /* ... and update where on is in the line */
  }
  else {                                           /* if this is a fork ... */
    qfirst = p;                       /* save node where you entered circle */
    if ((t->root == p) && (p->back != NULL))   /* if root has non-null back */
      q = p;
    else
      q = p->next;                           /* if null back or not at root */
    putc('(', outtree);                     /* open the paren for this fork */
    (*col)++;
    inloop = false;
    do {
      if (inloop) {                                 /* if not in first furc */
        putc(',', outtree);
        (*col)++;
        if (*col > 45) {                                  /* if got too far */
          putc('\n', outtree);
          *col = 0;
        }
      }
      if (q->back != NULL) {                /* just making sure is not null */
        ml_treeoutrecurs(outtreefile, t, q->back, bl_scale, col); /* go out */
        inloop = true;                  /* will need comma before next furc */
      }
      q = q->next;                                  /* continue around fork */
    } while (q != qfirst); /* until you get to where you entered the circle */
    putc(')', outtree);             /* then close the paren for this circle */
    (*col)++;
  }
  x = p->v * bl_scale;                   /* now write out the branch length */
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else {
    if (x == 0.0)
      w = 0;
    else
      w = (long)(0.43429448222 * log(-x)) + 1;
  }
  if (w < 0)
    w = 0;
  if (p == t->root)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
} /* ml_treeoutrecurse */


void ml_treeout(FILE* outtreefile, tree* t, node* p, double bl_scale)
{
  /* write out file with representation of final tree2 */
  int col;
  boolean found;
  node *q;

  assert(p->index > 0);                 // RSGdebug

  q = findrootmostandroot(t, p, &found);
  if (found)
    p = q;
  col = 0;
  ml_treeoutrecurs(outtreefile, t, p, bl_scale, &col);
}  /* ml_treeout */

/* End. */

