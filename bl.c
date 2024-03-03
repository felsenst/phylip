/* Version 4.0.  Copyright, 2023.
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

extern long endsite;
extern sequence inputSequences;
extern boolean lngths;
boolean inserting;


void bl_tree_new(struct tree **tp, long nonodes, long spp, long treesize)
{ /* make a new bl_tree.  Calls to generic_tree_new, casting bl_tree** to 
   * tree** as we call it, then call  bl_tree_init */
  struct tree *t;

  generic_tree_new(tp, nonodes, spp, treesize);                  /* next up */
  t = *tp;
  bl_tree_init(t, nonodes, spp);         /* initialize tree at this level */
} /* bl_tree_new */


void bl_tree_init(struct tree* t, long nonodes, long spp)
{ /* attributes of the generic tree that need bl function versions */

/* debug: if anything to initialize, do this here, but none right now */
  t->insert_ = bl_tree_insert_;
  t->save_lr_nodes = unrooted_tree_save_lr_nodes;
  t->restore_lr_nodes = unrooted_tree_restore_lr_nodes;
  t->save_traverses = bl_tree_save_traverses;
  t->restore_traverses = bl_tree_restore_traverses;
} /* bl_tree_init */


void bl_tree_save_traverses(struct tree* t, struct node* p, struct node* q) {
  /* bl_tree version ofn saving traverses */
	/* debug: more stuff here! */
} /* bl_tree_save_traverses */


struct node* bl_node_new(node_type type, long index, long nodesize) {
  /* go up hierarchy creating a node, initializing it */
  struct node* n;

  n = generic_node_new(type, index, nodesize);
  return n;
} /* bl_node_new */


void bl_node_init(struct node *n, node_type type, long index)
{
  /* initialize a node for bl trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  struct bl_node *bln;

  assert(index >= 0);

  generic_node_init(n, type, index);               /*  go up node hierarchy */
  bln = (struct bl_node*)n;
  bln->tyme = 0.0;
  bln->v = initialv;     /* debug: should be demoted to bl.h/bl.c ? */
  n->reinit = bl_node_reinit;
} /* bl_node_init */


boolean bl_node_good(struct tree *t, struct node *n)
{ /* check whether node exists and has forward and back lengths equal */
  boolean ok;
  struct bl_node *bln, *blnb;

  ok = generic_node_good(t, n);        /* first check if back branch exists */
  if (ok) {
    bln = (struct bl_node *)n;
    blnb = (struct bl_node *)(n->back);
    ok = (blnb->v == bln->v);    /* is sensible to check equality of reals! */
  }
  return ok;
} /* bl_node_good */


void bl_node_copy(struct node* srcn, struct node* destn)
{ /* copy a bl_node */
  struct bl_node *srcbln, *destbln;

  srcbln = (struct bl_node *)srcn;
  destbln = (struct bl_node *)destn;

  assert(srcn);                         // RSGdebug
  assert(destn);                        // RSGdebug
  generic_node_copy(srcn, destn);                        /* go up hierarchy */
  set_tyme(destn, srcbln->tyme);
  destbln->v = srcbln->v;
  destbln->deltav = srcbln->deltav;
  destbln->iter = srcbln->iter;
} /* bl_node_copy */


void bl_node_free(struct node **p)
{
  /* free a node for bl trees */
  struct node *n;
	 
  n = *p;
  generic_node_free(&n);
} /* bl_node_free */


void bl_hookup(struct node* p, struct node* q){
/* hook up two nodes, set branch length to initial value
   (one of the nodes may be in a fork circle) */

  hookup((struct node*)p, (struct node*)q);
  ((struct bl_node*)p)->v = initialv;
  ((struct bl_node*)q)->v = initialv;
} /* bl_hookup */


void bl_node_reinit(struct node* n)
{
  /* reset things for a bl tree node */
  struct bl_node *bln;

  bln = (struct bl_node*)n;
  generic_node_reinit(n);
  bln->tyme = 0.0;
  bln->v = initialv;
  bln->iter = true;   /* debug: do we know we need true? */
} /* bl_node_reinit */


void bl_node_print(struct node * n)
{
  /* for debugging only */
 
  generic_node_print(n);
  printf(" bl(tyme:%lf)", ((struct bl_node *)n)->tyme);
} /* bl_node_print */


void bl_update(struct tree *t, struct node *p)
{ /* calls nuview to make views at both ends of a branch.  Each is
   * made by recursive calls outward from there, as needed,
   * indicated by boolean initialized
   * the nuviews for the specific program are in turn called from
   * generic_tree_nuview  in phylip.c  */
/* debug:   I think redundant with calls in phylip.c  */

  if (p != NULL) {                                /* if not a NULL node ... */
    if (!p->tip)
      generic_tree_nuview(t, p);                    /* recurse from one end */
    if (p->back != NULL) {
      if (!p->back->tip)
        generic_tree_nuview(t, p->back);          /* recurse from the other */
    }
  }
}  /* bl_update */


void smooth(struct tree* t, node *p)
{  /* recursively do one step of smoothing on a branch, where
      smoothing includes getting views at both ends and using the 
      appropriate function to get a new branch length */
 /* debug: in which file should this be defined? bl.c? ml.c? */
 struct node *sib_ptr;

  if ( p == NULL )
    return;
  smoothed = false;

  bl_update(t, p);       /* get views at both ends updated, maybe recursing */
  if (p != NULL) {
    if (p->back != NULL) {
      t->makenewv (t, p);   /* new branch length using appropriate function */
      inittrav (t, p);    /* and thus set inward-looking pointers false ... */
      inittrav (t, p->back);                /* ... from both ends of branch */

      if ( p->tip )
        return;
      if ( (smoothed && !polishing) || !smoothit )
        return;
      }
      bl_update(t, p->back);                           /* update ends views */
    }
    for ( sib_ptr = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next ) {
          /* recursion out one end, the  p  end, to do this on all branches */
      if ( sib_ptr->back )
      {
        smooth(t, sib_ptr->back);                      /* go out from there */
        sib_ptr->initialized = false;    /* adjust its inward-looking views */
      }
    } 
}  /* smooth */


void bl_tree_smoothall(struct tree* t, node* p)
{
  /* go through the tree multiple times re-estimating branch lengths
   * using makenewv, with "initialized" reset and views updated
   * as needed.  It may seem like we are doing too many smooths, but sometimes
   * branch near p may already be completely smoothed from an
   * insert, this insures we spread out in the tree */
 /* debug: in which file should this be defined? bl.c? ml.c? */
  struct node *q;
  boolean save;
  int i;

  save = smoothit;
  smoothit = true;
  if (p == NULL) {        /* set outward-looking views uninitialized */
    inittrav(t, p);
    inittrav(t, p->back);
  }



  if ( p->tip )
    p = p->back;

  for ( i = 0 ; i < smoothings ; i++ )
  {
    smooth(t, p->back);
    if ( p->tip )
      return;
    for ( q = p->next ; q != p ; q = q->next)
      smooth(t, q->back);
  }
  smoothit = save;
} /* bl_tree_smoothall */


void bl_tree_do_branchl_on_insert(struct tree* t, struct node* forknode,
                                    struct node* q)
{ /* split original  qq->v  branch length evenly beween forknode->next and 
   * forknode->next->next.  forknode->back  must be subtree or tip.
   * This assumes interior node is a bifurcation.  Not able to cope if we 
   * insert at a rootmost branch one of whose ends is NULL
   * forknode should be where tip was hooked to.
   * that connection set to initial v for *both* directions.  */
  double newv;
  struct bl_node *qq;

  if (forknode->back != NULL) {              /* condition should never fail */
    ((struct bl_node*)forknode)->v = initialv; 
    ((struct bl_node*)(forknode->back))->v = initialv;
  }

  qq = (struct bl_node*)q;
  if (q->back != NULL)
    newv = qq->v * 0.5;
  else
    newv = initialv;

  if (forknode->next->back != NULL) { /* forknode->next for both directions */
    ((struct bl_node*)(forknode->next))->v = newv ;
    ((struct bl_node*)(forknode->next->back))->v = newv ;
  }

  if (forknode->next->next->back != NULL) {   /* next->next both directions */
    ((struct bl_node*)(forknode->next->next))->v = newv;
    ((struct bl_node*)(forknode->next->next->back))->v = newv;
  }

  /* debug:  BUG.970 -- might consider invalidating views here or in generic */
  /* debug:  do values of ->v get set earlier anyway?  */
  inittrav(t, forknode);                  /* some of this code unnecessary? */
  inittrav(t, forknode->back);
  inittrav(t, forknode->next);
  if (forknode->next->back != NULL)
    inittrav(t, forknode->next->back);
  inittrav(t, forknode->next->next);
  if (forknode->next->next->back != NULL)
    inittrav(t, forknode->next->next->back);
} /* bl_tree_do_branchl_on_insert */


void bl_tree_insert_(struct tree *t, struct node *p, 
                       struct node *q, boolean multif)
{
 /* 
  * After inserting via generic_tree_insert, branch length gets initialv. If
  * t->do_newbl is true, all branches optimized.
  *
  * Insert p near q 
  * p is the interior fork connected to the inserted subtree or tip */
  long i;

  generic_tree_insert_((struct tree*)t, p, q, multif);  /* debug:  maybe "multif"? */

  if ( !(((tree*)t)->do_newbl) )
  {
    invalidate_traverse(p);           /* set initialized false on views ... */
    invalidate_traverse(p->next);       /* ... looking in towards this fork */
    invalidate_traverse(p->next->next);
    p->initialized = false;           /* set initialized false on views ... */
    p->next->initialized = false;         /* ... out from the interior node */
    p->next->next->initialized = false;  
    inserting = true;
    bl_update(t, p);                            /* update the views outward */
    bl_update(t, p->next);
    bl_update(t, p->next->next);
    inserting = false;
  }
  else             /* this is the case where we recurse outwards, smoothing */
  {
    inittrav(t, p);                    /* set inward-looking pointers false */
    inittrav(t, p->back);
    bl_update(t, p);
    for ( i = 0 ; i < smoothings ; i++)
    {
      smooth(t, p);                      /* go around fork, out each branch */
    }
  }
} /* bl_tree_insert_ */


void generic_tree_save_lr_nodes(tree* t, node* p, node* r) {
  /* null operations if not replaced by polymorphic variant */
} /* generic_tree_save_lr_nodes */


void generic_tree_restore_lr_nodes(tree* t, node* p, node* r) {
  /* null operations if not replaced by polymorphic variant */
} /* generic_tree_restore_lr_nodes */


void unrooted_tree_save_lr_nodes(tree* t, node* p, node* r)
{
  /* save views and branch lengths around fork that is removed. */

  r->node_copy(r, t->lrsaves[0]);
  r->next->node_copy(r->next->back, t->lrsaves[1]);
  r->next->next->node_copy(r->next->next->back, t->lrsaves[2]);
  p->next->node_copy(p->next, t->lrsaves[3]);
  p->next->next->node_copy(p->next->next, t->lrsaves[4]);
  t->rb = r;                       /* pointers to the nodes of the fork ... */
  t->rnb = r->next;                                /* ... that contains  r  */
  t->rnnb = r->next->next;          /* (the "b" in their names is in error) */
} /* unrooted_tree_save_lr_nodes */


void unrooted_tree_restore_lr_nodes(tree* t, node* p, node* r)
{
    /* restore  r  fork nodes and inward views at  p  in unrooted tree case */
  struct bl_node *tr, *trb, *trnb, *trnnb, *trn, *trnn, *pnb, *pnnb;

  t->lrsaves[0]->node_copy(t->lrsaves[0], r->back);       /* these restore views */
  t->lrsaves[1]->node_copy(t->lrsaves[1], r->next->back);
  t->lrsaves[2]->node_copy(t->lrsaves[2], r->next->next->back);
  t->lrsaves[3]->node_copy(t->lrsaves[3], r->next);      /* inward-looking views */
  t->lrsaves[4]->node_copy(t->lrsaves[4], r->next->next);

  tr = (struct bl_node*)r;                               /* r  as a bl_node */
  trb = (struct bl_node*)(r->back);
  trnb = (struct bl_node*)(r->next->back);   /* to two of neighboring nodes */
  trnnb = (struct bl_node*)(r->next->next->back);
  trn = (struct bl_node*)(r->next);                 /* for view back in ... */
  trnn = (struct bl_node*)(r->next->next);         /* ... and the other one */
  trb->v = tr->v                         ;     /* branch lengths around  r  */
  trn->v = trnb->v;
  trnn->v = trnnb->v;
  pnb = (struct bl_node*)(p->next->back);
  pnnb = (struct bl_node*)(p->next->next->back);
  pnb->v = ((struct bl_node*)(p->next))->v; /* ... and on branches beyond p */
  pnnb->v = ((struct bl_node*)(p->next->next))->v;

  inittrav(t, r->back);    /*  to make sure initialized booleans are OK ... */
  inittrav(t, r->next->back);                /* ... in the neighbors of  r  */
  inittrav(t, r->next->next->back);
} /* unrooted_tree_restore_lr_nodes */


void bl_tree_do_branchl_on_re_move(struct tree* t, struct node* p,
            	                              struct node* q)
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
  struct bl_node *qq;

  qq = (struct bl_node*)q;
  combinedEdgeWeight = qq->v;
  if (q->back != NULL)
    combinedEdgeWeight = qq->v + ((struct bl_node*)(q->back))->v;
  qq->v = combinedEdgeWeight;
  if (q->back != NULL)
    ((struct bl_node*)(q->back))->v = combinedEdgeWeight;
} /* bl_tree_do_branchl_on_re_move */


void bl_tree_re_move(struct tree *t, struct node *p, 
                       struct node **q, boolean do_newbl)
{
  /* remove  p  and record in  q  where it was
   * assumes bifurcations
   * do_newbl is boolean which tells whether branch lengths get redone   */
  long i;

  generic_tree_re_move(t, p, q, do_newbl);

  if ( do_newbl )
  {
    for (i = 0 ; i < smoothings ; i++ )
    {
      if ( smoothit ) {
        smooth(t, *q);
        smooth(t, (*q)->back);
      }
    }
  }
  else {   /* update views at both ends of branch connected to q */
    if (!((*q)->tip))
      bl_update(t, *q);
    if (!((*q)->back->tip))
      bl_update(t, (*q)->back);
  }
} /* bl_tree_re_move */


void bl_tree_restore_traverses(struct tree *t, struct node *p, 
		                                 struct node* q) {
  /* restore branch lengths and mark views (un?)initialized */
  struct bl_node *pp, *qq, *ppb, *qqb;

  generic_tree_restore_traverses(t, p, q);
  pp = (struct bl_node*)p;
  ppb = (struct bl_node*)(p->back);
  if ( p->back )
  {
    ppb->v = pp->v;
    inittrav(t, p->back);      /* ... and similarly for other end if branch */
  }
  qq = (struct bl_node*)q;
  qqb = (struct bl_node*)(q->back);
  if ( q->back )
  {
    qqb->v = qq->v;
    inittrav(t, q->back);
  }
} /* bl_tree_restore_traverses */


boolean bl_tree_try_insert_thorough(struct tree *t, struct node *pp, 
                                     struct node *qq, 
                                     struct node *qqwherein,
                                     double *bestyet, struct tree *bestree,
                                     boolean thorough, boolean storing, 
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
  struct tree *tt;

  p = (struct node*)pp;
  q = (struct node*)qq;
  tt = (struct tree*)t;
  succeeded = false;
  tt->save_traverses(tt, p, q);
  tt->insert_(tt, p, q, false);
  tt->smoothall(tt, tt->root);
  like = tt->evaluate(tt, p, false);                  /* get score for tree */
printf("t->score, bestyet, like are now  %14.8f, %14.8f, %14.8f\n", tt->score, *bestyet, like);   /* debug */

  if (atstart) {          /* save the tree if it is the first one or better */
    bettertree = true;
    *bestyet = like;
printf("set *bestyet to  %14.8f\n", like);   /* debug */
  } else {
    bettertree = (like > *bestyet);
printf("*bestyet, like are %14.8f, %14.8f\n", *bestyet, like);   /* debug */
if(bettertree) printf("found better tree, tt->score = %14.8f\n", tt->score); /* debug */
    succeeded = bettertree;
    }
  if (bettertree) {                    /* set variables for return, and ...*/
    *bestyet = like;
printf("set *bestyet to  %14.8f\n", like);   /* debug */
    qqwherein = qq;
    tt->copy(tt, (struct tree*)bestree);  /* save tree in bestree, and ... */
printf("bestree->score is now  %14.8f\n", ((struct tree*)bestree)->score);   /* debug */
  }
  tt->re_move(tt, p, &whereRemoved, false);  /* then remove inserted stuff */

/* debug: not sure what whereRemoved is doing for us:  assert(whereRemoved == q);  */
/* debug:  probably redundant: */   tt->restore_traverses(tt, p, q); /*  debug */

  /* Update t->score of tree on which placements are being tested */
  generic_update(tt, q);                             /* update views on restored tree */

  return succeeded;
} /* bl_tree_try_insert_thorough */


boolean bl_tree_try_insert_(struct tree* tt, struct node* pp, 
                          struct node* qq, struct node* qwherein, 
                          double* bestyet, struct tree* bestree, 
                          boolean thorough, boolean storing, boolean atstart, 
                          double* bestfound)
{
 /* Passes to bl_tree_try_insert_thorough or bl_tree_try_insert_notthorough
  * depending on the value of thorough. If multf is given, sets to
  * false.  */
  boolean succeeded = false;
  struct node *p, *q;
  struct tree* t;

  p = (struct node*)pp;
  q = (struct node*)qq;
  t = (struct tree*)tt;
  if ( thorough )
    succeeded = bl_tree_try_insert_thorough(tt, pp, qq, qwherein, bestyet,
                                           bestree, thorough, false, atstart);
  else  /* debug:  need to have a _notthorough function here instead? */
    generic_tree_insert_(t, p, q, false);

  return succeeded;
} /* bl_tree_try_insert_ */


void blk_tree_insert_(struct tree *t, struct node *newtip, 
		        struct node *below, boolean dummy, boolean dummy2)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  boolean done;
  struct node *p, *newfork;

  /* first stick it in the right place */
  rooted_tree_insert_(t, newtip, below, dummy);

  below = t->nodep[below->index - 1];
  newfork = t->nodep[newtip->back->index - 1];
  newtip = t->nodep[newtip->index-1];
  /* now for the tyme stuff */
  if (((struct bl_node*)newtip)->tyme 
         < ((struct bl_node*)below)->tyme)
    p = newtip;
  else p = below;

  set_tyme(newfork, ((struct bl_node*)p)->tyme);
  if (newfork->back != NULL)
  {
    /* here we rescale the tree to fit the subtree being added      *
     * note that if we are sticking a new node into the tree and    *
     * the branches are only epsilon appart, allow the branch       *
     * lengths to be 1/2 epsilon, so that we interfere with the     *
     * tree minimally                                               */
    if (((struct bl_node*)p)->tyme 
           > ((struct bl_node*)newfork->back)->tyme)
      set_tyme(newfork, (((struct bl_node*)p)->tyme +
               ((struct bl_node*)newfork->back)->tyme) / 2.0);
    else
      set_tyme(newfork, ((bl_node*)p)->tyme - (epsilon/2));
    do {
      p = t->nodep[p->back->index - 1];
      done = (p == t->root);
      if (!done) {
        done = ((struct bl_node*)(t->nodep[p->back->index - 1]))->tyme 
                   < ((struct bl_node*)p)->tyme;
        set_tyme(p->back,  ((struct bl_node*)p)->tyme - epsilon/2);
      }
    } while (!done);
  }
  else
    set_tyme(newfork, ((struct bl_node*)newfork)->tyme - initialv);

  if ( !smoothit ) {
    smooth(t, newfork);
    smooth(t, newfork->back);
  }
  else {
    inittrav(t, newtip);
    inittrav(t, newtip->back);
    for (i = 0 ; i < smoothings ; i++) {
      smooth(t, newfork);
      smooth(t, newfork->back);
    }
  }
}  /* blk_tree_insert_ */


double get_tyme(struct node *p)
{ /* return the tyme of a bl_node. p must point to struct bl_node. */
  return (((struct bl_node*)p)->tyme);
} /* get_tyme */


void set_tyme (struct node* p, double tyme)
{ /* Set the tyme of a node and its sibs. */
  struct node *q, *sib_ptr;

  q = p;
  sib_ptr = q;
  if (p->next)
    do {
      ((struct bl_node*)sib_ptr)->tyme = tyme;
      /* added because changing tymes usually invalidates data likelihood.
       * This set seems to fix a failure to find the best tree in some
       * cases, but if the flags are being properly maintained it shouldn't...
       * apparent fix to bug#296, JY and MK 2015/05/18 */
      sib_ptr->initialized = false;
      sib_ptr = sib_ptr->next;
    } while (sib_ptr != q);
  ((struct bl_node*)p)->tyme = tyme;
} /* set_tyme */


void blk_tree_re_move(struct tree* t, struct node *item, 
                        struct node** where, boolean do_newbl) {
  /* Removes nodes item and its ancestor, where, from the tree.
     The new descendant of where's ancestor is made to be where's second
     descendant (other than item).  Also returns pointers to the deleted 
     nodes, item and where, and records where they were deleted from. */
  long i;
  struct node *whereloc;

  rooted_tree_re_move(t, item, &whereloc, do_newbl);
  if ( where )  where = &whereloc;

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


double min_child_tyme(struct node *p)
{
  /* Return the minimum tyme of all children. p must be a parent nodelet */
  double min;
  struct node *q, *qb;

  min = 1.0;                                /* tymes are always nonpositive */
  for ( q = p->next; q != p; q = q->next ) {
    qb = q->back;
    if ( get_tyme(qb) < min )
      min = get_tyme(qb);
  }
  return min;
} /* min_child_tyme */


double parent_tyme(struct node *p)
{
  /* Return the tyme of the parent of node p.  p must be a parent node. */

  if (p->back)
    return get_tyme(p->back);
  else
    return MIN_ROOT_TYME;
} /* parent_tyme */


boolean valid_tyme(struct tree *t, struct node *p, double tyme) {
  /* Return true if tyme is a valid tyme to assign to node p. tyme must be
   * finite, not greater than any of p's children, and not less than p's
   * parent. Also, tip nodes can only be assigned 0. Otherwise false is
   * returned. */

  p = t->nodep[p->index - 1];

#ifdef __USE_C99        /* debug: TODO Find a way to check without this. */
  if ( !isfinite(tyme) ) return false;
#endif
  if ((p->tip == true) && (tyme != 0.0) ) return false;
  if ( tyme > min_child_tyme(p) ) return false;
  if ( tyme < parent_tyme(p) ) return false;
  return true;
} /* valid_tyme */


double set_tyme_evaluate(struct tree *t, struct node *p, double tyme)
{
  /* Change tyme of node p and return likelihood
   * Views should be invalidated and regenerated before calling
   * evaluate() anywhere else in the tree. */

  assert( valid_tyme(t, p, tyme) );

  set_tyme(p, tyme);
  t->nuview(t, p);

  return t->evaluate(t, p, false);
} /* set_tyme_evaluate */


void blk_tree_makenewv(struct tree* t, struct node *p)
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
  struct node *s, *sib_ptr;

  if ( p->tip )                                               /* skip tips. */
    return;

  s = t->nodep[p->index - 1];

#ifdef MAKENEWV_DEBUG
  double start_tyme = get_tyme(s);
  double start_likelihood = t->score;
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
                                            keep the current point and quit */
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

void getthree(struct tree* t, struct node *p, double thigh, double tlow, 
               double tdelta, double *x, double *lnl)
{
  /* compute scores at a new triple of points */
  int i;
  double tt = ((bl_node*)p)->tyme;
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


void bl_treevaluate(struct tree* curtree, boolean improve, boolean reusertree,
                    boolean global, boolean progress, struct tree* priortree,
                    struct tree* bestree, initialvtrav_t initialvtrav)
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


void bl_initialvtrav(struct tree* t, struct bl_node *p)
{
  /* traverse tree to set branch lengths  v  to initial values
   * must be called twice the first time, at both ends of
   * a branch such as the root branch.  Is separate from the
   * task of setting initialized booleans for views to false   */
  struct node *pp, *qq;

  pp = (struct node*)p;
  if (pp == NULL)                      /* if this is a NULL branch bail out */
    return;
  if ((!lngths) || p->iter) {     /* set length of this branch to  initialv */
    p->v = initialv;
    ((struct bl_node*)(pp->back))->v = initialv;
  }
  if (!pp->tip) {     /* go around circle, calling initialvtrav on all backs */
    qq = pp->next;
    while ( qq != pp ) {
      bl_initialvtrav(t, ((struct bl_node*)(qq->back)));
      qq = qq->next;
    }
  }
}  /* bl_initialvtrav */


void addelement2(tree* t, struct node *qq, Char *ch, long *parens,
                 FILE *treefile, boolean lngths, double *trweight,
                 boolean *goteof, long *nextnode, long *ntips, 
                 long no_species, boolean *haslengths, boolean unifok,
                 long maxnodes)
{ /* recursive procedure adds nodes to user-defined tree
   * -- old-style bifurcating-only version used only by treeread2
   * which is used only in Contml, Fitch, Kitsch, and Restml.  */
  struct node *pfirst = NULL, *p;
  struct bl_node *q;
  long i, current_loop_index;
  boolean notlast, minusread;
  Char str[MAXNCH];
  double valyew, divisor;
  long furcs = 0;

  if ((*ch) == '(') {

    current_loop_index = (*nextnode) + t->spp;
    (*nextnode)++;

    if ( (maxnodes != -1) && (current_loop_index > maxnodes)) {
      sprintf(progbuf,
            "ERROR in intree file: Attempting to allocate too many nodes.\n");
      print_progress(progbuf);
      sprintf(progbuf,
                  "This is usually caused by a unifurcation.  To use this\n");
      print_progress(progbuf);
      sprintf(progbuf,
                  "intree with this program, use Retree to read and write\n");
      print_progress(progbuf);
      sprintf(progbuf, "this tree.\n");
      print_progress(progbuf);
      exxit(-1);
    }
    /* This is an assignment of an interior node */
    p = t->nodep[current_loop_index];
    pfirst = p;
    notlast = true;
    while (notlast) {      /* This while loop goes through a circle (triad for
                                         the case of bifurcations) of nodes */
      furcs++;
      p = p->next;
      /* added to ensure that non base nodes in loops have indices */
      p->index = current_loop_index + 1;

      getch(ch, parens, treefile);

      addelement2(t, p, ch, parens, treefile, lngths, trweight, goteof,
                   nextnode, ntips, no_species, haslengths, unifok, maxnodes);
      /* recursive call for subtrees */

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
                 (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    if ( furcs <= 1 && !unifok ) {
      sprintf(progbuf,
               "ERROR in intree file: A Unifurcation was detected.\n");
      print_progress(progbuf);
      sprintf(progbuf,
            "To use this intree with this program, use Retree to read and\n");
      print_progress(progbuf);
      sprintf(progbuf, " write this tree.\n");
      print_progress(progbuf);
      exxit(-1);
    }

  } else if ((*ch) != ')') {                       /* read the species name */
    for (i = 0; i < MAXNCH; i++)
      str[i] = '\0';
    take_name_from_tree (ch, str, treefile);
    match_names_to_data (str, t->nodep, &p, spp);
    pfirst = p;
    if ((*ch) == ')')
      (*parens)--;
    (*ntips)++;
    strncpy (p->nayme, str, MAXNCH);
  } else
    getch(ch, parens, treefile);

  if ((*ch) == '[')          /* getting tree weight from last comment field */
  {
    if (!eoln(treefile))
    {
      if(fscanf(treefile, "%lf", trweight) < 1)
      {
        printf("\n\nERROR reading tree file./n/n");
        exxit(-1);
      }
      getch(ch, parens, treefile);
      if (*ch != ']')
      {
        sprintf(progbuf, "\n\nERROR:  Missing right square bracket.\n\n");
        print_progress(progbuf);
        exxit(-1);
      }
      else
      {
        getch(ch, parens, treefile);
        if (*ch != ';') {
          sprintf(progbuf,
                  "\n\nERROR:  Missing semicolon after square brackets.\n\n");
          print_progress(progbuf);
          exxit(-1);
        }
      }
    }
  }
  else if ((*ch) == ';') {
    (*trweight) = 1.0 ;
    if (!eoln(treefile))
      sprintf(progbuf, "WARNING:  Tree weight set to 1.0\n");
    print_progress(progbuf);
  }
  else
    (*haslengths) = (*haslengths) && (qq == NULL);
  q = (struct bl_node*)qq;
  if (qq != NULL)
    hookup(qq, pfirst);
  /* debug:   if (q != NULL) {
    if (q->branchnum < pfirst->branchnum)
    pfirst->branchnum = q->branchnum;
    else
    q->branchnum = pfirst->branchnum;
    }  FIXME check if we need this for restml */

  if ((*ch) == ':') {                               /* read a branch length */
    processlength(&valyew, &divisor, ch, &minusread, treefile, parens);
    if (qq != NULL) {
      if (!minusread)
        q->oldlen = valyew / divisor;
      else
        q->oldlen = initialv;
      if (lngths) {
        q->v = valyew / divisor;
        ((struct bl_node*)(qq->back))->v = q->v;
        q->iter = false;
        ((struct bl_node*)(qq->back))->iter = false;
      }
    }
  }
}  /* addelement2 */


void treeread2 (tree* t, FILE *treefile, node **root, boolean lngths,
                 double *trweight, boolean *goteof, boolean *haslengths,
                 long *no_species, boolean unifok, long maxnodes)
{
  /* read in user-defined tree and set it up
     -- old-style bifurcating-only version used only in Fitch, Kitsch,
     Contml, and Restml.  Needs to be replaced by generic treeread */
  char  ch;
  long parens = 0;
  long ntips = 0;
  long nextnode;

  (*goteof) = false;
  nextnode = 0;

  /* Eats all blank lines at start of file */
  while (eoln(treefile) && !eoff(treefile))
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  }

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
     * encounter an open-paren */
    getch(&ch, &parens, treefile);
  }

  addelement2(t, NULL, &ch, &parens, treefile, lngths, trweight, goteof,
              &nextnode, &ntips, (*no_species), haslengths, unifok, maxnodes);
  (*root) = t->nodep[*no_species];

  /*eat blank lines */
  while (eoln(treefile) && !eoff(treefile))
    scan_eoln(treefile);

  (*root)->oldlen = 0.0;

  if (parens != 0) {
    sprintf(progbuf, "\n\nERROR in tree file:  unmatched parentheses.\n\n");
    print_progress(progbuf);
    exxit(-1);
  }
}  /* treeread2 */


void bl_treeoutrecurs(FILE* outtreefile, struct tree* t, struct node* p, 
                        double bl_scale, int* col)
{ 
  /* write out to output file a subtree, recursively.  This is the version 
   * with branch lengths and a scale factor,  bl_scale  */
  long i, n, w;
  Char c;
  double x;
  struct node *q, *qfirst; 
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
        bl_treeoutrecurs(outtreefile, t, q->back, 
			  bl_scale, col);  /* go out the branch recursively */
        inloop = true;                  /* will need comma before next furc */
      }
      q = q->next;                                  /* continue around fork */
    } while (q != qfirst); /* until you get to where you entered the circle */
    putc(')', outtree);             /* then close the paren for this circle */
    (*col)++;
  }
  x = ((struct bl_node*)p)->v * bl_scale;    /* now write out branch length */
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


void bl_treeout(FILE* outtreefile, struct tree* t, struct node* p, 
                  double bl_scale)
{
  /* write out file with representation of final tree2 */
  int col;
  boolean found;
  struct node *q;

  assert(p->index > 0);                 // RSGdebug

  q = findrootmostandroot(t, p, &found);
  if (found)
    p = q;
  col = 0;
  bl_treeoutrecurs(outtreefile, t, p, bl_scale, &col);
}  /* bl_treeout */

/* End. */

