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

extern FILE *outtree;

const double MIN_BRANCH_LENGTH = epsilon/4.0;
const double MIN_ROOT_TYME = -10;

extern long endsite;
extern sequence inputSequences;
extern boolean lngths, smoothit, polishing;
boolean inserting;


void bl_tree_new(struct bl_tree **tp, long nonodes, long spp, long treesize)
{ /* make a new bl_tree.  Calls to generic_tree_new, casting bl_tree** to 
   * tree** as we call it, then call  bl_tree_init */
  struct tree **tt;
  struct bl_tree *bltt;

  tt = (struct tree**)tp;
  generic_tree_new(tt, nonodes, spp, treesize);                  /* next up */
  bltt = (struct bl_tree *)(*tt);
  bl_tree_init(bltt, nonodes, spp);         /* initialize tree at this level */
} /* bl_tree_new */


void bl_tree_init(struct bl_tree* t, long nonodes, long spp)
{ /* attributes of the generic tree that need bl function versions */

/* debug: if anything to initialize, do this here, but none right now */
  ((struct tree*)t)->save_lr_nodes = unrooted_tree_save_lr_nodes;
  ((struct tree*)t)->restore_lr_nodes = unrooted_tree_restore_lr_nodes;
  ((struct tree*)t)->save_traverses = bl_tree_save_traverses;
  ((struct tree*)t)->restore_traverses = bl_tree_restore_traverses;
} /* bl_tree_init */


void bl_tree_save_traverses(struct tree* t, struct node* p, struct node* q) {
  /* bl_tree version ofn saving traverses */
	/* debug: more stuff here! */
} /* bl_tree_save_traverses */


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
  /* initialize a node for bl trees */
/* debug: not needed for dist_node creation but needed for sequence types.  Needs nodesize argument? probably not */

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  generic_node_init((struct node*)n, type, index);                /* go up node hierarchy */
  n->tyme = 0;
  n->v = initialv;     /* debug: should be demoted to bl.h/bl.c ? */
/* debug: initialize branch lengths here too? */
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
  struct bl_node *srcbln;
  struct bl_node *destbln;

  srcbln = (struct bl_node *)srcn;
  destbln = (struct bl_node *)destn;

  assert(srcn);                         // RSGdebug
  assert(destn);                        // RSGdebug
  generic_node_copy(srcn, destn);                /* first call generic copy */
  set_tyme(destbln, srcbln->tyme);
  destbln->v = srcbln->v;
  destbln->deltav = srcbln->deltav;
  destbln->iter = srcbln->iter;
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


void bl_node_reinit(struct bl_node* bln)
{
  /* reset things for an ml tree node */
  struct node *n;

  n = (struct node*)bln;
  generic_node_reinit(n);
  bln->tyme = 0.0;
  bln->v = initialv;
  bln->iter = true;   /* debug: do we know we need true? */
} /* bl_node_reinit */


void bl_node_print(struct bl_node * bln)
{
  /* for debugging only */
  struct node * n;
 
  n = (struct node*)bln;

  generic_node_print(n);
  printf(" bl(tyme:%lf)", bln->tyme);
} /* bl_node_print */


void bl_update(struct bl_tree *blt, struct bl_node *pp)
{ /* calls nuview to make views at both ends of a branch.  Each is
   * made by recursive calls outward from there, as needed,
   * indicated by boolean initialized
   * the nuviews for the specific program are in turn called from
   * generic_tree_nuview  in phylip.c  */
/* debug:   I think redundant with calls in phylip.c  */
  struct node *p;
  struct tree *t;

  t = (struct tree*)blt;
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


void smooth_traverse(struct bl_tree* t, bl_node *pp)
{ /* start traversal, smoothing branch lengths, in both directions from
   * this branch */
 /* debug: in which file should this be defined? bl.c? ml.c? */
  struct node *p;

  p = (struct node*)pp;
  smooth(t, pp);
  smooth(t, ((struct bl_node*)(p->back)));
} /* smooth_traverse */


void smooth(struct bl_tree* t, bl_node *pp)
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
  ((struct tree*)t)->makenewv ((struct tree *)t, p);   /* new branch length */
  inittrav (((struct tree*)t), p); /* set inward-looking pointers false ... */
  inittrav (((struct tree*)t), p->back);    /* ... from both ends of branch */

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


void bl_tree_smoothall(struct bl_tree* t, bl_node* pp)
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


void bl_tree_do_branchl_on_insert(struct bl_tree* t, struct bl_node* forknode,
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
  inittrav((struct tree*)t, forkn);    /* some of this code unnecessary? */
  inittrav((struct tree*)t, forkn->back);
  inittrav((struct tree*)t, forkn->next);
  if (forkn->next->back != NULL)
    inittrav((struct tree*)t, forkn->next->back);
  inittrav((struct tree*)t, forkn->next->next);
  if (forkn-> next->next->back != NULL)
    inittrav((struct tree*)t, forkn->next->next->back);
} /* bl_tree_do_branchl_on_insert */


void bl_tree_insert_(struct bl_tree *t, struct bl_node *pp, 
                       struct bl_node *qq, boolean multif)
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
    bl_update(t, pp);                           /* update the views outward */
    bl_update(t, (struct bl_node*)(p->next));
    bl_update(t, (struct bl_node*)(p->next->next));
    inserting = false;
  }
  else             /* this is the case where we recurse outwards, smoothing */
  {
    inittrav((struct tree*)t, p);      /* set inward-looking pointers false */
    inittrav((struct tree*)t, p->back);
    bl_update(t, (struct bl_node*)p);
    for ( i = 0 ; i < smoothings ; i++)
    {
      smooth_traverse(t, (struct bl_node*)p);        /* go around fork, out */
    }
  }
} /* bl_tree_insert */


void generic_tree_save_lr_nodes(tree* t, node* p, node* r) {
  /* null operations if not replaced by polymorphic variant */
} /* generic_tree_save_lr_nodes */


void generic_tree_restore_lr_nodes(tree* t, node* p, node* r) {
  /* null operations if not replaced by polymorphic variant */
} /* generic_tree_restore_lr_nodes */


void unrooted_tree_save_lr_nodes(tree* t, node* p, node* r)
{
  /* save views and branch lengths around fork that is removed. */

  r->copy(r, t->lrsaves[0]);
  r->next->copy(r->next->back, t->lrsaves[1]);
  r->next->next->copy(r->next->next->back, t->lrsaves[2]);
  p->next->copy(p->next, t->lrsaves[3]);
  p->next->next->copy(p->next->next, t->lrsaves[4]);
  t->rb = r;                       /* pointers to the nodes of the fork ... */
  t->rnb = r->next;                                /* ... that contains  r  */
  t->rnnb = r->next->next;          /* (the "b" in their names is in error) */
} /* unrooted_tree_save_lr_nodes */


void unrooted_tree_restore_lr_nodes(tree* t, node* p, node* r)
{
    /* restore  r  fork nodes and inward views at  p  in unrooted tree case */
  struct bl_node *tr, *trb, *trnb, *trnnb, *trn, *trnn, *pnb, *pnnb;

  t->lrsaves[0]->copy(t->lrsaves[0], r->back);       /* these restore views */
  t->lrsaves[1]->copy(t->lrsaves[1], r->next->back);
  t->lrsaves[2]->copy(t->lrsaves[2], r->next->next->back);
  t->lrsaves[3]->copy(t->lrsaves[3], r->next);      /* inward-looking views */
  t->lrsaves[4]->copy(t->lrsaves[4], r->next->next);

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


void bl_tree_do_branchl_on_re_move(struct bl_tree* t, struct bl_node* pp,
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


void bl_tree_re_move(struct bl_tree *t, struct bl_node *pp, 
                       struct bl_node **qq, boolean do_newbl)
{
  /* remove  p  and record in  q  where it was
   * assumes bifurcations
   * do_newbl is boolean which tells whether branch lengths get redone   */
  long i;
  struct node *p, **q;

  p = (struct node*)pp;
  q = (struct node**)qq;
  generic_tree_re_move((struct tree*)t, p, q, do_newbl);

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


boolean bl_tree_try_insert_thorough(struct bl_tree *t, struct bl_node *pp, 
                                     struct bl_node *qq, 
                                     struct bl_node *qqwherein,
                                     double *bestyet, struct bl_tree *bestree,
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
  like = tt->evaluate(tt, q, 0); /* evaluate restored tree to update views */

  return succeeded;
} /* bl_tree_try_insert_thorough */


boolean bl_tree_try_insert_(struct bl_tree* tt, struct bl_node* pp, 
                          struct bl_node* qq, struct bl_node* qwherein, 
                          double* bestyet, struct bl_tree* bestree, 
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


void blk_tree_insert_(struct bl_tree *tt, struct bl_node *nnewtip, 
		        struct bl_node *bbelow, 
                        boolean dummy, boolean dummy2)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  boolean done;
  struct node *p, *newfork, *newtip, *below;
  struct bl_node* nnewfork;
  struct tree* t;

  newtip = (struct node*)nnewtip;
  below = (struct node*)bbelow;
  t = (struct tree*)tt;
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
    smooth(tt, nnewfork);
    smooth(tt, ((struct bl_node*)(newfork->back)));
  }
  else {
    inittrav(t, newtip);
    inittrav(t, newtip->back);
    for (i = 0 ; i < smoothings ; i++) {
      smooth(tt, nnewfork);
      smooth(tt, ((struct bl_node*)(newfork->back)));
    }
  }
}  /* blk_tree_insert_ */


double get_tyme(struct bl_node *p)
{ /* return the tyme of a bl_node. p must point to struct bl_node. */
  return (p->tyme);
} /* get_tyme */


void set_tyme (struct bl_node* p, double tyme)
{ /* Set the tyme of a node and its sibs. p must point to struct bl_node. */
  struct bl_node *sib_ptr;
  struct node *pp, *ssib_ptr;

  sib_ptr = p;
  pp = (struct node*)p;
  ssib_ptr = pp;
  if (pp->next)
    do {
      sib_ptr->tyme = tyme;
      /* added because changing tymes usually invalidates data likelihood.
       * This set seems to fix a failure to find the best tree in some
       * cases, but if the flags are being properly maintained it shouldn't...
       * apparent fix to bug#296, JY and MK 2015/05/18 */
      ssib_ptr->initialized = false;
      ssib_ptr = ssib_ptr->next;
      sib_ptr = (struct bl_node*)ssib_ptr;
    } while (ssib_ptr != pp);
  else
    p->tyme = tyme;
} /* set_tyme */


void blk_tree_re_move(struct bl_tree* t, struct bl_node *item, 
                        struct bl_node** where, boolean do_newbl) {
  /* Removes nodes item and its ancestor, where, from the tree.
     The new descendant of where's ancestor is made to be where's second
     descendant (other than item).  Also returns pointers to the deleted 
     nodes, item and where, and records where they were deleted from. */
  long i;
  struct bl_node *whereloc;
  struct node *wwhereloc;
  struct tree* tt;

  tt = (struct tree*)t;
  rooted_tree_re_move(tt, (struct node*)item, &wwhereloc, do_newbl);
  whereloc = (struct bl_node*)wwhereloc;
  if ( where )  where = &whereloc;

  if ( do_newbl ) {
    whereloc = (struct bl_node *)wwhereloc;
    inittrav(tt, wwhereloc);
    inittrav(tt, wwhereloc->back);
    for ( i = 0 ;  i < smoothings ; i++) {
      smooth(t, whereloc);
      smooth(t, (struct bl_node*)(wwhereloc->back));
    }
  }
  else smooth(t, (struct bl_node*)(wwhereloc->back));
}  /* blk_tree_re_move */


double min_child_tyme(struct bl_node *pp)
{
  /* Return the minimum tyme of all children. p must be a parent nodelet */
  double min;
  struct node *p, *q;
  struct bl_node* blnqb;

  p = (struct node*)pp;
  min = 1.0;                                /* tymes are always nonpositive */
  for ( q = p->next; q != p; q = q->next ) {
    blnqb = (struct bl_node*)(q->back);
    if ( get_tyme(blnqb) < min )
      min = get_tyme(blnqb);
  }
  return min;
} /* min_child_tyme */


double parent_tyme(struct bl_node *pp)
{
  /* Return the tyme of the parent of node p.  p must be a parent node. */
  struct node *p;

  p = (struct node*)pp;
  if (p->back)
    return get_tyme((struct bl_node*)(p->back));
  else
    return MIN_ROOT_TYME;
} /* parent_tyme */


boolean valid_tyme(struct bl_tree *blt, struct bl_node *pp, double tyme) {
  /* Return true if tyme is a valid tyme to assign to node p. tyme must be
   * finite, not greater than any of p's children, and not less than p's
   * parent. Also, tip nodes can only be assigned 0. Otherwise false is
   * returned. */
  struct tree* t;
  struct node* p;

  p = (struct node*)pp;
  t = (struct tree*)blt;
  p = t->nodep[p->index - 1];

#ifdef __USE_C99        /* debug: TODO Find a way to check without this. */
  if ( !isfinite(tyme) ) return false;
#endif
  if ( (p->tip == true) && (tyme != 0.0) ) return false;
  if ( tyme > min_child_tyme(pp) ) return false;
  if ( tyme < parent_tyme(pp) ) return false;
  return true;
} /* valid_tyme */


double set_tyme_evaluate(struct bl_tree *blt, struct bl_node *pp, double tyme)
{
  /* Change tyme of node p and return likelihood
   * Views should be invalidated and regenerated before calling
   * evaluate() anywhere else in the tree. */
  struct tree* t;
  struct node* p;

  assert( valid_tyme(blt, pp, tyme) );
  t = (struct tree*)blt;
  p = (struct node*)pp;

  set_tyme(pp, tyme);
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
  struct bl_tree *blt;
  struct bl_node *bls;

  if ( p->tip )                                               /* skip tips. */
    return;

  s = t->nodep[p->index - 1];
  bls = (struct bl_node*)s;
  blt = (struct bl_tree*)t;

#ifdef MAKENEWV_DEBUG
  double start_tyme = get_tyme(bls);
  double start_likelihood = t->score;
#endif /* MAKENEWV_DEBUG */

  if (s == t->root)                  /* Tyme cannot be less than parent ... */
    min_tyme = MIN_ROOT_TYME;
  else
    min_tyme = get_tyme(bls) + MIN_BRANCH_LENGTH;

  max_tyme = min_child_tyme(bls) - MIN_BRANCH_LENGTH;  /* or > any children */

  /* Nothing to do if we can't move */
  if ( max_tyme < min_tyme + 2.0*min_tyme_delta ) {  /* done if can't move! */
    done = true;
    return;
  }

  current_tyme = get_tyme(bls);
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

    lnl[0] = set_tyme_evaluate(blt, bls, x[0]);   /* scores of the three points */
    lnl[1] = set_tyme_evaluate(blt, bls, x[1]);
    lnl[2] = set_tyme_evaluate(blt, bls, x[2]);

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
    new_likelihood = set_tyme_evaluate(blt, bls, new_tyme);
    while ( new_likelihood < current_likelihood ) {
            /* if our estimate is worse, retract until we find a better one */
#ifdef MAKENEWV_DEBUG
      putchar('<');
#endif /* MAKENEWV_DEBUG */
      tdelta *= retract_factor;
      uphill_step *= retract_factor;
      if (fabs(tdelta) < min_tdelta) {      /* if can't retract far enough ...
                                            keep the current point and quit */
        new_likelihood = set_tyme_evaluate(blt, bls, current_tyme);
        done = true;
#ifdef MAKENEWV_DEBUG
        putchar('X');
#endif /* MAKENEWV_DEBUG */
        break;
      }
      new_tyme = current_tyme + tdelta;
      new_likelihood = set_tyme_evaluate(blt, bls, new_tyme);
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
    set_tyme((struct bl_node*)p, x[i]);
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


void bl_treeoutrecurs(FILE* outtreefile, struct tree* t, struct bl_node* pp, 
                        double bl_scale, int* col)
{ 
  /* write out to output file a subtree, recursively.  This is the version 
   * with branch lengths and a scale factor,  bl_scale  */
  long i, n, w;
  Char c;
  double x;
  struct node *p, *q, *qfirst; 
  boolean inloop;

  p = (struct node*)pp;
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
        bl_treeoutrecurs(outtreefile, t, ((struct bl_node*)(q->back)), 
			  bl_scale, col); /* go out the branch recursiovely */
        inloop = true;                  /* will need comma before next furc */
      }
      q = q->next;                                  /* continue around fork */
    } while (q != qfirst); /* until you get to where you entered the circle */
    putc(')', outtree);             /* then close the paren for this circle */
    (*col)++;
  }
  x = pp->v * bl_scale;                  /* now write out the branch length */
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


void bl_treeout(FILE* outtreefile, struct tree* t, struct bl_node* pp, 
                  double bl_scale)
{
  /* write out file with representation of final tree2 */
  int col;
  boolean found;
  struct node *p, *q;
  struct bl_node *blp;

  p = (struct node*)pp;
  assert(p->index > 0);                 // RSGdebug

  q = findrootmostandroot(t, p, &found);
  if (found)
    p = q;
  blp = (struct bl_node*)p;
  col = 0;
  bl_treeoutrecurs(outtreefile, t, blp, bl_scale, &col);
}  /* bl_treeout */

/* End. */

