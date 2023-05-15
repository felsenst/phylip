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

#define DEBUG
#define MAKENEWV_DEBUG
/* #define USE_NEW_MAKENEWV */

const double MIN_BRANCH_LENGTH = epsilon/4.0;
const double MIN_ROOT_TYME = -10;

long endsite;
extern sequence inputSequences;
extern boolean lngths, smoothit, polishing;
boolean inserting;


void bl_tree_new(struct tree **tp, long nonodes, long spp, long treesize)
{ /* make a new bl_tree.  Calls to generic_tree_new,
   * casting bl_tree** to tree** as we call it 
   * then call  bl_tree_init */

  generic_tree_new(tp, nonodes, spp, treesize);   /* next up tree hierarchy */
} /* bl_tree_new */


void bl_tree_init(struct tree* t, long nonodes, long spp)
{ /* 
   * attributes of the generic tree that need ml function versions */

  generic_tree_init(t, nonodes, spp);              /* go up class hierarchy */
} /* bl_tree_init */


struct bl_node* bl_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct node* n;
  struct bl_node* nn;

  n = generic_node_new(type, index, nodesize);
  nn = (struct bl_node*)n;
  bl_node_init(nn, type, index);
  return (bl_node*)nn;
} /* bl_node_new */


void bl_node_init(struct bl_node *n, node_type type, long index)
{
  /* initialize a node for ml trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */
  struct bl_node* nn;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  generic_node_init((struct node*)n, type, index);                /* go up node hierarchy */
  nn = (bl_node*)n;
  nn->tyme = 0;
} /* bl_node_init */


void bl_node_copy(struct bl_node* srcn, struct bl_node* destn)
{ /* copy a bl_node */
/* debug: shouldn't length of node be involved? */
  struct node *src = (struct node *)srcn;
  struct node *dest = (struct node *)destn;
  assert(srcn);                         // RSGdebug
  assert(destn);                        // RSGdebug
  generic_node_copy(src, dest);
  set_tyme(destn, srcn->tyme);
} /* bl_node_copy */


void bl_node_free(struct bl_node **np)
{
  /* free a node for bl trees */
  struct node *n;
	 
  n = (struct node*)np;
  generic_node_free(&n);
} /* bl_node_free */


void bl_hookup(struct bl_node* p, struct bl_node* q){
/* hook up two nodes, set branch length to initial value
   (one of the nodes may be in a fork circle) */

  hookup((struct node*)p, (struct node*)q);
  p->v = initialv;
  q->v = initialv;
} /* bl_hookup */


void bl_node_reinit(struct bl_node * bln)
{
  /* reset things for an ml tree node */
  struct node * n;

  n = (struct node*)bln;
  bln->tyme = 0.0;
  generic_node_reinit(n);
} /* bl_node_reinit */


void bl_node_print(struct bl_node * bln)
{
  /* for debugging only */
  struct node * n;
 
  n = (struct node*)bln;

  generic_node_print(n);
  printf(" bl(tyme:%lf)", bln->tyme);
} /* bl_node_print */


void bl_update(struct tree *t, struct bl_node *pp)
{ /* calls nuview to make views at both ends of a branch.  Each is
   * made by recursive calls outward from there, as needed,
   * indicated by boolean initialized
   * the nuviews for the specific program are in turn called from
   * generic_tree_nuview  in phylip.c  */
/* debug:   I think redundant with calls in phylip.c  */
  struct node *p;

  p = (struct node*)pp;
  if (p != NULL) {                                /* if not a NULL node ... */
    if (!p->tip)
      generic_tree_nuview(t, p);                    /* recurse from one end */
    if ((p->back != NULL) && (!p->back->tip)) {
      if (!p->back->tip)
        generic_tree_nuview(t, p->back);          /* recurse from the other */
    }
  }
}  /* bl_update */


void smooth_traverse(tree* t, bl_node *pp)
{ /* start traversal, smoothing branch lengths, in both directions from
   * this branch */
 /* debug: in which file should this be defined? bl.c? ml.c? */
  struct node *p;

  p = (struct node*)pp;
  smooth(t, pp);
  smooth(t, ((struct bl_node*)(p->back)));
} /* smooth_traverse */


void smooth(tree* t, bl_node *pp)
{  /* repeatedly and recursively do one step of smoothing on a branch */
 /* debug: in which file should this be defined? bl.c? ml.c? */
 struct node *p, *sib_ptr;

  p = (struct node*)pp;
  if ( p == NULL )
    return;
/* debug:    if ( p->back == NULL )
    return;
debug */
  smoothed = false;

  bl_update(t, pp);      /* get views at both ends updated, maybe recursing */
  t->makenewv (t, p);                         /* new value of branch length */
  inittrav (t, p);                 /* set inward-looking pointers false ... */
  inittrav (t, p->back);               /* ... from both ends of this branch */

  if ( p->tip )
    return;
  if ( (smoothed && !polishing) || !smoothit )
    return;

  for ( sib_ptr = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next )
  {       /* recursion out one end, the  p  end, to do this on all branches */
    if ( sib_ptr->back )
    {
      smooth(t, (struct bl_node*)(sib_ptr->back));     /* go out from there */
      sib_ptr->initialized = false;  /* inward-looking views need adjusting */
    }
  }
  bl_update(t, ((struct bl_node*)(p->back)));          /* update ends views */
}  /* smooth */


void bl_tree_smoothall(tree* t, bl_node* pp)
{
  /* go through the tree multiple times re-estimating branch lengths
   * using makenewv, with "initialized" reset and views updated
   * as needed.  It may seem like we are doing too many smooths, but sometimes
   * branch near p may already be completely smoothed from an
   * insert, this insures we spread out in the tree */
 /* debug: in which file should this be defined? bl.c? ml.c? */
  boolean save;
  int i;
  struct node* p, *q;

  save = smoothit;
  smoothit = true;
  p = (struct node*)pp;
  if ( p->tip )
    p = p->back;

/* debug:   editing mistake near here when removing debugging prints? */
  for ( i = 0 ; i < smoothings ; i++ )
  {
    smooth(t, ((struct bl_node*)(p->back)));
    if ( p->tip )
      return;
    for ( q = p->next ; q != p ; q = q->next)
      smooth(t, ((struct bl_node*)(q->back)));
  }
  smoothit = save;
} /* bl_tree_smoothall */


void bl_tree_do_branchl_on_insert(tree* t, struct bl_node* forknode, 
                                    struct bl_node* qq)
{ /* split original  qq->v  branch length evenly beween forknode->next and 
   * forknode->next->next.  forknode->back  must be subtree or tip.
   * This assumes interior node is a bifurcation.  Not able to cope if we 
   * insert at a rootmost branch one of whose ends is NULL
   * forknode should be where tip was hooked to.
   * that connection set to initial v for *both* directions.  */
  double newv;
  struct node *forkn, *q;

  forkn = (struct node*)forknode;
  if (forkn->back != NULL) {              /* condition should never fail */
    forknode->v = initialv; 
    ((struct bl_node*)(forkn->back))->v = initialv;
  }

  q = (struct node*)qq;
  if (q->back != NULL)
    newv = qq->v * 0.5;
  else
    newv = initialv;

  if (forkn->next->back != NULL) { /* forknode->next for both directions */
    ((struct bl_node*)(forkn->next))->v = newv ;
    ((struct bl_node*)(forkn->next->back))->v = newv ;
  }

  if (forkn->next->back != NULL) {     /* next->next for both directions */
    ((struct bl_node*)(forkn->next->next))->v = newv;
    ((struct bl_node*)(forkn->next->next->back))->v = newv;
  }

  /* debug:  BUG.970 -- might consider invalidating views here or in generic */
  /* debug:  do values of ->v get set earlier anyway?  */
  inittrav(t, forkn);         /* some of this block of code unnexessary? */
  inittrav(t, forkn->back);
  inittrav(t, forkn->next);
  if (forkn->next->back != NULL)
    inittrav(t, forkn->next->back);
  inittrav(t, forkn->next->next);
  if (forkn-> next->next->back != NULL)
    inittrav(t, forkn->next->next->back);
} /* bl_tree_do_branchl_on_insert */


void bl_tree_insert_(tree *t, struct bl_node *pp, struct bl_node *qq, 
                       boolean multif)
{
 /* 
  * After inserting via generic_tree_insert, branch length gets initialv. If
  * t->do_newbl is true, all branches optimized.
  *
  * Insert p near q 
  * p is the interior fork connected to the inserted subtree or tip */
  long i;
  struct node *p, *q;

  p = (struct node*)pp;
  q = (struct node*)qq;
  generic_tree_insert_(t, p, q, multif);  /* debug:  maybe "multif"? */

  if ( !t->do_newbl )
  {
    invalidate_traverse(p);        /* set initialized false on views ... */
    invalidate_traverse(p->next);  /* ... looking in towards this fork */
    invalidate_traverse(p->next->next);
    p->initialized = false;        /* set initialized false on views ... */
    p->next->initialized = false;  /* ... out from the interior node */
    p->next->next->initialized = false;  
    inserting = true;
    bl_update(t, pp);               /* update the views outward */
    bl_update(t, (struct bl_node*)(p->next));
    bl_update(t, (struct bl_node*)(p->next->next));
    inserting = false;
  }
  else    /* this is the case where we recurse outwards, smoothing */
  {
    inittrav(t, p);        /* set inward-looking pointers false */
    inittrav(t, p->back);
    bl_update(t, (struct bl_node*)p);
    for ( i = 0 ; i < smoothings ; i++)
    {
      smooth_traverse(t, (struct bl_node*)p);        /* go around fork, out */
    }
  }
} /* bl_tree_insert */


void bl_tree_do_branchl_on_re_move(tree* t, struct bl_node* pp,
            	                              struct bl_node* qq)
{
  /* sum up branch lengths on the two other neighbors of  p
   * that are connected to fork  q */
 /* debug:  only works for bifurcations */
  /*
   * debug: BUG.970 -- add this when moved into re_move
   * assert(q->next->next->next == q);
   *
   * also, should we call generic_do_branchl_on_re_move(t, p, q); ??
   */
  double combinedEdgeWeight;
  struct node *q;

  q = (struct node*)qq;
  combinedEdgeWeight = qq->v;
  if (q->back != NULL)
    combinedEdgeWeight = qq->v + ((struct bl_node*)(q->back))->v;
  qq->v = combinedEdgeWeight;
  if (q->back != NULL)
    ((struct bl_node*)(q->back))->v = combinedEdgeWeight;
} /* bl_tree_do_branchl_on_re_move */


void bl_tree_re_move(tree *t, struct bl_node *pp, struct bl_node **qq, 
                       boolean do_newbl)
{
  /* remove  p  and record in  q  where it was
   * assumes bifurcations
   * do_newbl is boolean which tells whether branch lengths get redone   */
  long i;
  struct node *p, **q;

  p = (struct node*)pp;
  q = (struct node**)qq;
  generic_tree_re_move(t, p, q, do_newbl);

  if ( do_newbl )
  {
    for (i = 0 ; i < smoothings ; i++ )
    {
      if ( smoothit ) {
        smooth_traverse(t, *qq);
        smooth_traverse(t, (struct bl_node*)((*q)->back));
      }
    }
  }
  else {   /* update views at both ends of branch connected to q */
    if (!((*q)->tip))
      bl_update(t, ((struct bl_node*)(*q)));
    if (!((*q)->back->tip))
      bl_update(t, ((struct bl_node*)((*q)->back)));
  }
} /* bl_tree_re_move */


boolean bl_tree_try_insert_thorough(tree *t, struct bl_node *pp, 
                                     struct bl_node *qq, 
                                     struct bl_node *qqwherein,
                                     double *bestyet, tree *bestree, boolean 
                                     thorough, boolean storing, 
                                     boolean atstart)
{
 /* Temporarily inserts  p  at  q  and evaluates. If the rearrangement is
  * better than bestyet, updates bestyet and returns true.  If this is the
  * first place to insert, set  bestyet  to the current likelihood and set
  * qwhere  to the current place  q  */
  double like;
  boolean succeeded, bettertree;
  struct node* whereRemoved;
  struct node *p, *q;

  succeeded = false;
  t->save_traverses(t, p, q);
  t->insert_(t, p, q, false);
  t->smoothall(t, t->root);
  like = t->evaluate(t, p, false);                    /* get score for tree */
printf("t->score, bestyet, like are now  %14.8f, %14.8f, %14.8f\n", t->score, *bestyet, like);   /* debug */

  if (atstart) {          /* save the tree if it is the first one or better */
    bettertree = true;
    *bestyet = like;
printf("set *bestyet to  %14.8f\n", like);   /* debug */
  } else {
    bettertree = (like > *bestyet);
printf("*bestyet, like are %14.8f, %14.8f\n", *bestyet, like);   /* debug */
if(bettertree) printf("found better tree, t->score = %14.8f\n", t->score); /* debug */
    succeeded = bettertree;
    }
  if (bettertree) {                    /* set variables for return, and ...*/
    *bestyet = like;
printf("set *bestyet to  %14.8f\n", like);   /* debug */
    qqwherein = qq;
    t->copy(t, bestree);              /* save the tree in bestree, and ... */
printf("bestree->score is now  %14.8f\n", bestree->score);   /* debug */
  }
  t->re_move(t, p, &whereRemoved, false);    /* then remove inserted stuff */

/* debug: not sure what whereRemoved is doing for us:  assert(whereRemoved == q);  */
/* debug:  probably redundant: */   t->restore_traverses(t, p, q); /*  debug */

  /* Update t->score of tree on which placements are being tested */
  like = t->evaluate(t, q, 0);   /* evaluate restored tree to update views */

  return succeeded;
} /* bl_tree_try_insert_thorough */


boolean bl_tree_try_insert_(tree* t, struct bl_node* pp, struct bl_node* qq, 
                          struct bl_node* qwherein, double* bestyet, 
                          tree* bestree, boolean thorough,
                          boolean storing, boolean atstart, double* bestfound)
{
 /* Passes to bl_tree_try_insert_thorough or bl_tree_try_insert_notthorough
  * depending on the value of thorough. If multf is given, sets to
  * false.  */
  boolean succeeded = false;
  struct node *p, *q;

  p = (struct node*)pp;
  q = (struct node*)qq;
  if ( thorough )
    succeeded = bl_tree_try_insert_thorough(t, pp, qq, qwherein, bestyet,
                                           bestree, thorough, false, atstart);
  else  /* debug:  need to have a _notthorough function here instead? */
    generic_tree_insert_(t, p, q, false);

  return succeeded;
} /* bl_tree_try_insert_ */


void blk_tree_insert_(tree *t, struct bl_node *nnewtip, 
		        struct bl_node *bbelow, 
                        boolean dummy, boolean dummy2)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  boolean done;
  struct node *p, *newfork, *newtip, *below;
  struct bl_node* nnewfork;
  /* debug GOT TO HERE */

  newtip = (struct node*)nnewtip;
  below = (struct node*)bbelow;
  /* first stick it in the right place */
  rooted_tree_insert_(t, newtip, below, dummy);

  below = t->nodep[below->index - 1];
  newfork = t->nodep[newtip->back->index - 1];
  nnewfork = (struct bl_node*)newfork;
  newtip = t->nodep[newtip->index-1];
  /* now for the tyme stuff */
  if (((struct bl_node*)newtip)->tyme 
         < ((struct bl_node*)below)->tyme)
    p = newtip;
  else p = below;

  set_tyme(nnewfork, ((struct bl_node*)p)->tyme);
  if (newfork->back != NULL)
  {
    /* here we rescale the tree to fit the subtree being added      *
     * note that if we are sticking a new node into the tree and    *
     * the branches are only epsilon appart, allow the branch       *
     * lengths to be 1/2 epsilon, so that we interfere with the     *
     * tree minimally                                               */
    if (((struct bl_node*)p)->tyme 
           > ((struct bl_node*)newfork->back)->tyme)
      set_tyme(nnewfork, (((struct bl_node*)p)->tyme +
               ((bl_node*)newfork->back)->tyme) / 2.0);
    else
      set_tyme(nnewfork, ((bl_node*)p)->tyme - (epsilon/2));
    do {
      p = t->nodep[p->back->index - 1];
      done = (p == t->root);
      if (!done) {
        done = (((struct bl_node*)t->nodep[p->back->index - 1])->tyme 
                   < ((struct bl_node*)p)->tyme);
        set_tyme((struct bl_node*)(p->back), 
		  ((struct bl_node*)p)->tyme - epsilon/2);
      }
    } while (!done);
  }
  else
    set_tyme(nnewfork, ((struct bl_node*)newfork)->tyme - initialv);

  if ( !smoothit ) {
    smooth(t, nnewfork);
    smooth(t, ((struct bl_node*)(newfork->back)));
  }
  else {
    inittrav(t, newtip);
    inittrav(t, newtip->back);
    for (i = 0 ; i < smoothings ; i++) {
      smooth(t, nnewfork);
      smooth(t, ((struct bl_node*)(newfork->back)));
    }
  }
}  /* blk_tree_insert_ */


double get_tyme(node *p)
{ /* return the tyme of a bl_node. p must point to struct bl_node. */
  return ((struct bl_node *)p)->node.tyme;
} /* get_tyme */


void set_tyme (struct bl_node* p, double tyme)
{ /* Set the tyme of a node and its sibs. p must point to struct bl_node. */
  struct bl_node *sib_ptr;

  sib_ptr = p;
  if ( p->next )
    do {
      ((bl_node*)sib_ptr)->node.tyme = tyme;
      /* added because changing tymes usually invalidates data likelihood.
       * This set seems to fix a failure to find the best tree in some
       * cases, but if the flags are being properly maintained it shouldn't...
       * apparent fix to bug#296, JY and MK 2015/05/18 */
      ((bl_node*)sib_ptr)->node.initialized = false;
      sib_ptr = sib_ptr->next;
    } while (sib_ptr != p );
  else
    ((bl_node*)p)->node.tyme = tyme;
} /* set_tyme */


void blk_tree_re_move(tree* t, struct bl_node *item, struct bl_node** where, 
                        boolean do_newbl) {
  /* Removes nodes item and its ancestor, where, from the tree.
     The new descendant of where's ancestor is made to be where's second descendant (other than item).
     Also returns pointers to the deleted nodes, item and where, and records where they were deleted from. */
  long i;
  struct bl_node* whereloc;

  rooted_tree_re_move(t, item, &whereloc, do_newbl);
  if ( where )  *where = whereloc;

  if ( do_newbl ) {
    inittrav(t, whereloc);
    inittrav(t, whereloc->back);
    for ( i = 0 ;  i < smoothings ; i++) {
      smooth(t, whereloc);
      smooth(t, whereloc->back);
    }
  }
  else smooth(t, whereloc->back);
}  /* blk_tree_re_move */


#ifdef USE_NEW_MAKENEWV

/******* PROPAGATED FROM 3.6 ************/

double min_child_tyme(node *p)
{
  /* Return the minimum tyme of all children. p must be a parent nodelet */
  double min;
  struct bl_node *q;

  min = 1.0;                                /* Tymes are always nonpositive */
  for ( q = p->next; q != p; q = q->next ) {
    if ( get_tyme(q->back) < min )
      min = get_tyme(q->back);
  }
  return min;
} /* min_child_tyme */


double parent_tyme(node *p)
{
  /* Return the tyme of the parent of node p.  p must be a parent node. */
  if (p->back)
    return get_tyme(p->back);
  /* else */
  return MIN_ROOT_TYME;
} /* parent_tyme */


boolean valid_tyme(tree *t, struct bl_node *p, double tyme) {
  /* Return true if tyme is a valid tyme to assign to node p. tyme must be
   * finite, not greater than any of p's children, and not less than p's
   * parent. Also, tip nodes can only be assigned 0. Otherwise false is
   * returned. */

  p = t->nodep[p->index - 1];

#ifdef __USE_C99        /* debug: TODO Find a way to check without this. */
  if ( !isfinite(tyme) ) return false;
#endif
  if ( p->tip == true && tyme != 0.0 ) return false;
  if ( tyme > min_child_tyme(p) ) return false;
  if ( tyme < parent_tyme(p) ) return false;
  return true;
} /* valid_tyme */


double set_tyme_evaluate(tree *t, struct bl_node *p, double tyme)
{
  /* Change tyme of node p and return likelihood
   * Views should be invalidated and regenerated before calling
   * evaluate() anywhere else in the tree. */

  /* node *sib_ptr;
     long num_sibs, i;  debug */

  assert( valid_tyme(t, p, tyme) );

  set_tyme(p, tyme);
  t->nuview(t, p);

  return t->evaluate(t, p, false);
} /* set_tyme_evaluate */


void blk_tree_makenewv(tree* t, struct bl_node *p)
{
  /* Improve a node tyme using Newton-Raphson
   *
   * Slope and curvature are estimated at the current point and used to
   * interpolate a new point. If the curvature is positive, the next point
   * the estimations are pushed uphill by a fraction of the total range.
   * If any interpolation fails to produce a better result, the result is
   * retracted by a given factor toward the original point and tested again.
   *
   * The function terminates when max_iterations have been performed, or
   * when the likelihood improvement for any iteration is less than epsilon,
   * or when a retraction fails. If the iterations are exhausted,
   * 'smoothed' is set false, indicating that further improvement may be
   * possible by additional calls to makenewv(). Otherwise 'smoothed' is left
   * untouched.
   *
   * Define MAKENEWV_DEBUG to get lots of junk printed to stdout.
   * Each loop prints a character, as follows:
   *   '/' ->   Positive curvature, positive slope, moved +
   *   '\' ->   Positive curvature, negative slope, moved -
   *   ')' ->   Negative curvature, moved +
   *   '(' ->   Negative curvature, moved -
   *   '<' ->   Retracting back by retract_factor
   *   'X' ->   Retraction failed, keeping current point
   */
  /* Tuning constants */
  const double likelihood_epsilon = epsilon/1000.0;
                     /* Any change in score less than this, and we're done. */
  const double tyme_delta = epsilon;  /* Small tyme difference used
                                            to compute slope and curvature. */
  const double min_tyme_delta = tyme_delta / 10.0;     /* absolute smallest */
  const double uphill_step_factor = 0.05; /* Fraction of the current branch 
                                             length to move uphill in positive
                                             curvature regions. */
  const double retract_factor = 0.5; /* This defines how far back we go if the
                                        interpolated point is lower than the 
                                        original value */
  const double min_tdelta = epsilon; /* Minimum to which we will retract 
                                        before giving up. */
  const long max_iterations = 100; /* Maximum iterations - typically we stop 
                                      much sooner */
  double min_tyme, max_tyme;
  double current_tyme, new_tyme;
  double current_likelihood, new_likelihood;
  double x[3], lnl[3];       /* tyme (x) and score (lnl) points below, at, and
                                above the current tyme */
/* debug: rename lnl? change "likelihood" to "score" in comments? */
  double s21, s10, slope;
  double curv;
  double uphill_step;
  double tdelta;                  /* interpolated point minus current point */
  long iteration;
  boolean done;
  long num_sibs, i;
  struct bl_node *sib_ptr;

  if ( p->tip )                                               /* skip tips. */
    return;

  struct node *s = t->nodep[p->index - 1];

#ifdef MAKENEWV_DEBUG
  double start_tyme = get_tyme(s);
  double start_likelihood = t->score;
  long uphill_steps = 0;
#endif /* MAKENEWV_DEBUG */

  if (s == t->root)                  /* Tyme cannot be less than parent ... */
    min_tyme = MIN_ROOT_TYME;
  else
    min_tyme = get_tyme(s) + MIN_BRANCH_LENGTH;

  max_tyme = min_child_tyme(s) - MIN_BRANCH_LENGTH;    /* or > any children */

  /* Nothing to do if we can't move */
  if ( max_tyme < min_tyme + 2.0*min_tyme_delta ) {  /* done if can't move! */
    done = true;
    return;
  }

  current_tyme = get_tyme(s);
  current_likelihood = t->evaluate(t, s, false);

  uphill_step = (max_tyme - min_tyme) * uphill_step_factor;

  done = false;
  for ( iteration = 0; iteration < max_iterations; iteration++) {
    x[0] = current_tyme - tyme_delta;     /* three points for interpolation */
    if ( x[0] < min_tyme )
      x[0] = min_tyme;
    x[2] = current_tyme + tyme_delta;
    if ( x[2] > max_tyme )
      x[2] = max_tyme;
    x[1] = (x[0] + x[2]) / 2.0;

    lnl[0] = set_tyme_evaluate(t, s, x[0]);   /* scores of the three points */
    lnl[1] = set_tyme_evaluate(t, s, x[1]);
    lnl[2] = set_tyme_evaluate(t, s, x[2]);

    s21 = (lnl[2] - lnl[1]) / (x[2] - x[1]);              /* compute slopes */
    s10 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    slope = s21 + s10 / 2.0;

    curv = (s21 - s10) / ((x[2] - x[0]) / 2);          /* compute curvature */
    if (curv >= 0.0) { /* In negative curvature regions, just move uphill by a
                        * fraction of the current length */
      tdelta = copysign(uphill_step, slope);
#ifdef MAKENEWV_DEBUG
      uphill_steps++;
      if (tdelta > 0) putchar('/');
      else putchar('\\');
#endif /* MAKENEWV_DEBUG */
    }
    else {                              /* otherwise guess where slope is 0 */
      tdelta = -(slope / curv);
#ifdef MAKENEWV_DEBUG
      if (tdelta > 0) putchar(')');
      else putchar('(');
#endif /* MAKENEWV_DEBUG */
    }

    new_tyme = current_tyme + tdelta;              /* propose to move there */
    if (new_tyme <= min_tyme) {               /* but don't go past min_tyme */
      new_tyme = min_tyme;
      tdelta = new_tyme - current_tyme;
    }
    else if (new_tyme >= max_tyme) {                /* ... or past max_tyme */
        new_tyme = max_tyme;
        tdelta = new_tyme - current_tyme;
      }
    new_likelihood = set_tyme_evaluate(t, s, new_tyme);
    while ( new_likelihood < current_likelihood ) {
            /* if our estimate is worse, retract until we find a better one */
#ifdef MAKENEWV_DEBUG
      putchar('<');
#endif /* MAKENEWV_DEBUG */
      tdelta *= retract_factor;
      uphill_step *= retract_factor;
      if (fabs(tdelta) < min_tdelta) {      /* if can't retract far enough ...
                                         /* keep the current point and quit */
        new_likelihood = set_tyme_evaluate(t, s, current_tyme);
        done = true;
#ifdef MAKENEWV_DEBUG
        putchar('X');
#endif /* MAKENEWV_DEBUG */
        break;
      }
      new_tyme = current_tyme + tdelta;
      new_likelihood = set_tyme_evaluate(t, s, new_tyme);
    }
    if ( new_likelihood - current_likelihood < likelihood_epsilon ) {
      done = true;
    }
    current_likelihood = new_likelihood;
    if (done) break;
  }
  num_sibs = count_sibs(s);              /* invalidate and regenerate views */
  sib_ptr = p;
  for (i = 0 ; i < num_sibs; i++ ) {
    sib_ptr = sib_ptr->next;
    inittrav (t, sib_ptr);
  }
  if (!done) smoothed = false;
#ifdef MAKENEWV_DEBUG
  fprintf(stdout, "\nmakenewv(): node %ld: %ld iterations (%f,%f) => (%f,%f)\n", p->index, iteration+1, start_tyme, start_likelihood, current_tyme, current_likelihood);
#endif
}  /* blk_tree_makenewv */


/******* END PROPAGATED FROM 3.6 ************/

#else /* ifndef USE_NEW_MAKENEWV */

void blk_tree_makenewv(tree* t, struct node *p) {
  /* improve a node time */
  long it, imin, imax, i;
  double tt, tfactor, tlow, thigh, oldlike, oldx, ymin, ymax, s32, s21, yold;
  boolean done, already;
  struct node *s, *sib_ptr, *sib_back_ptr;
  double tdelta, curv, slope, lnlike;
  double  x[3], lnl[3];

  if ( p->tip )                                            /* don't do tips */
    return;

  s = t->nodep[p->index - 1];
  oldx = ((bl_node*)s)->node.tyme;                        /* store old tyme */
  lnlike = oldlike = t->evaluate(t, p, 0);  /* evaluate and store old score */
  if (s == t->root)
    tlow = -10.0;                           /* default minimum tyme at root */
  else
    tlow = ((bl_node*)(s->back))->node.tyme;    /* otherwise tyme >= parent */

  sib_ptr = s;                   /* set maximum tyme to smallest child tyme */
  thigh = ((bl_node*)s->next->back)->node.tyme;
  for (sib_ptr = s->next ; sib_ptr != s ; sib_ptr = sib_ptr->next) {
    sib_back_ptr = sib_ptr->back;
    if (((bl_node*)sib_back_ptr)->node.tyme < thigh)
      thigh = ((bl_node*)sib_back_ptr)->node.tyme;
  }
  if (thigh - tlow < 4.0*epsilon)   /* if thigh and tlow are close to equal */
    return;
  if (s != t->root)
    tdelta = (thigh - tlow) / 10.0;
  else {
    tdelta = (thigh - ((bl_node*)s)->node.tyme) / 5.0;
    if (tdelta  < 2 * epsilon ) tdelta = 2 * epsilon;
  }
  getthree(t, s, thigh, tlow, tdelta, x, lnl);          /* get three points */
  it = 0;
  tfactor = 1.0;
  done = false;
  while (it < iterations && !done) {
    ymax = lnl[0];
    imax = 0;
    for (i = 1; i <= 2; i++) {  /* figure out which point has largest score */
      if (lnl[i] > ymax) {
        ymax = lnl[i];
        imax = i;
      }
    }
    if (imax != 1) {             /* swap points so that x[1] scores highest */
      /* debug: TODO Explain why we are doing this */
      ymax = x[1];                                /* ymax is temporary only */
      x[1] = x[imax];
      x[imax] = ymax;
      ymax = lnl[1];
      lnl[1] = lnl[imax];
      lnl[imax] = ymax;
    }
    tt = x[1];
    yold = tt;
    s32 = (lnl[2] - lnl[1]) / (x[2] - x[1]);               /* average slope */
    /* avg slope near (x[1]+x[0])/2 */
    s21 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    if (fabs(x[2] - x[0]) > epsilon)                   /* average curvature */
      curv = (s32 - s21) / ((x[2] - x[0]) / 2);
    else
      curv = 0.0;
    slope = (s32 + s21) / 2 - curv * (x[2] - 2 * x[1] + x[0]) / 4;
                                               /* interpolate slope at x[1] */
    if (curv >= 0.0) {
      if (slope < 0)
        tdelta = -fabs(tdelta);
      else
        tdelta = fabs(tdelta);
    }
    else
      tdelta = -(tfactor * slope / curv);     /* approximate Newton-Raphson */
    if (tt + tdelta <= tlow + epsilon)           /* don't let it go too far */
      tdelta = tlow + epsilon - tt;
    if (tt + tdelta >= thigh - epsilon)
      tdelta = thigh - epsilon - tt;
    tt += tdelta;
    done = (fabs(yold - tt) < epsilon || fabs(tdelta) < epsilon);
    set_tyme(s, tt);
    t->nuview(t, s);
    lnlike = t->evaluate(t, s, false);
    ymin = lnl[0];
    imin = 0;
    for (i = 1; i <= 2; i++) {    /* figure out which of the three original */
      if (lnl[i] < ymin) {                /* points has the lowest ln score */
        ymin = lnl[i];
        imin = i;
      }
    }
    already = (tt == x[0]) || (tt == x[1]) || (tt == x[2]);
    if (!already && ymin < lnlike) {  /* if the minimum point is lower than */
      x[imin] = tt;                 /* our new interpolated point then take */
      lnl[imin] = lnlike;                /* that point and put it where the */
    }                                              /* interpolated point is */
    if (already || lnlike < oldlike) {
      tt = oldx;                    /* if either our interpolated point has */
      set_tyme(s, oldx);            /* a lower score or is equivalent to    */
      tfactor /= 2;                     /* our original, reinterpolate this */
      tdelta /= 2;                              /* time go only half as far */
      t->score = oldlike;
      lnlike = oldlike;
    }
    else {
      tfactor = 1.0;
      oldlike = lnlike;
      oldx = tt;
    }
    if (!done)                          /* apply it to the sibs */
    {
      set_tyme(p, tt);
      t->nuview(t, p);
      for (sib_ptr = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next)
        t->nuview(t, sib_ptr);
    }
    it++;
  }
  if ( smoothit )
    inittrav(t, p);
  p->initialized = false;
  for (sib_ptr = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next) {
    sib_ptr->initialized = false;
    if ( smoothit )
      inittrav(t, sib_ptr);
  }
  t->score = lnlike;
  smoothed = smoothed && done;
}  /* blk_tree_makenewv */

#endif /* USE_NEW_MAKENEWV */


void getthree(tree* t, struct node *p, double thigh, double tlow, 
               double tdelta, double *x, double *lnl)
{
  /* compute scores at a new triple of points */
  int i;
  double tt = ((bl_node*)p)->node.tyme;
  double td = fabs(tdelta);

  x[0] = tt - td;                        /* points are on either side of tt */
  x[1] = tt;
  x[2] = tt + td;
  if (x[0] < tlow + epsilon) {                /* make sure not to go to low */
    x[0] = tlow + epsilon;
    x[1] = ( x[0] + x[2] ) / 2;
  }
  if (x[2] > thigh - epsilon) {                          /* ... or too high */
    x[2] = thigh - epsilon;
    x[1] = ( x[0] + x[2] ) / 2;
  }
  for ( i = 0 ; i < 3 ; i++ ) {                 /* get scores for all three */
    set_tyme(p, x[i]);
    t->nuview(t, p);
    lnl[i] = t->evaluate(t, p, 0);  /* debug: make sure score of tree is not reset */
  }
}  /* getthree */


void bl_treevaluate(tree* curtree, boolean improve, boolean reusertree,
                    boolean global, boolean progress, tree* priortree,
                    tree* bestree, initialvtrav_t initialvtrav)
{
  /* evaluate a user tree */
  double bestyet;

  smoothit = improve;
  if (reusertree)                             /* if rearrange on user trees */
  {
    arbitrary_resolve(curtree);              /* put in zero-length branches */
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
  else {
    if (!lngths) {                  /* put initial lengths on all branches */
      inittrav(curtree, curtree->root);
      inittrav(curtree, curtree->root->back);
    }
    polishing = true;
    smoothit = true;
    curtree->evaluate(curtree, curtree->root, 0);     /* get current value */
    if (!lngths)
      curtree->smoothall(curtree, curtree->root); /* improve branch lengths */
    smoothit = improve;
    polishing = false;
  }
  curtree->evaluate(curtree, curtree->root, true);             /* get score */
}  /* bl_treevaluate */


void bl_initialvtrav(tree* t, struct bl_node *p)
{
  /* traverse tree to set branch lengths  v  to initial values
   * must be called twice the first time, at both ends of
   * a branch such as the root branch.  Is separate from the
   * task of setting initialized booleans for views to false   */
  struct bl_node* q;

  if (p == NULL)                       /* if this is a NULL branch bail out */
    return;
  if ((!lngths) || p->iter) {     /* set length of this branch to  initialv */
    p->v = initialv;
    p->back->v = initialv;
  }
  if (!p->tip) {     /* go around circle, calling initialvtrav on all backs */
    q = p->next;
    while ( q != p ) {
      bl_initialvtrav(t, q->back);
      q = q->next;
    }
  }
}  /* bl_initialvtrav */


void bl_treeoutrecurs(FILE* outtreefile, tree* t, struct bl_node* p, 
                        double bl_scale, int* col)
{ 
  /* write out to output file a subtree, recursively.  This is the version 
   * with branch lengths and a scale factor,  bl_scale  */
  long i, n, w;
  Char c;
  double x;
  struct bl_node *q, *qfirst;
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
        bl_treeoutrecurs(outtreefile, t, q->back, bl_scale, col); /* go out */
        inloop = true;                  /* will need comma before next furc */
      }
      q = q->next;                                  /* continue around fork */
    } while (q != qfirst); /* until you get to where you entered the circle */
    putc(')', outtree);             /* then close the paren for this circle */
    (*col)++;
  }
  x = p->v * bl_scale;                   /* now write out the branch length */
  if (x > 0.0)           /* widths in decimal places hence divide by ln(10) */
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
} /* bl_treeoutrecurse */


void bl_treeout(FILE* outtreefile, tree* t, struct bl_node* p, 
                  double bl_scale)
{
  /* write out file with representation of final tree2 */
  int col;
  boolean found;
  struct bl_node *q;

  assert(p->index > 0);                 // RSGdebug

  q = findrootmostandroot(t, p, &found);
  if (found)
    p = q;
  col = 0;
  bl_treeoutrecurs(outtreefile, t, p, bl_scale, &col);
}  /* bl_treeout */

/* End. */

