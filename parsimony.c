/* Version 4.0.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Michal Palczewski.
   */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "parsimony.h"

/* globals best defined here */
extern long nonodes, endsite, outgrno, which;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff;
extern steptr weight, category, alias, location, ally;
extern node** lrsaves;
extern tree* curtree, bestree;

long maxtrees;
double *threshwt;
boolean usertree;
boolean reusertree;
node *temp, *temp1, *temp2, *tempsum, *temprm, *tempadd, *tempf, *tmp, *tmp1, *tmp2, *tmp3, *tmprm, *tmpadd;
boolean lastrearr, recompute;
double bestyet, bestlike, bstlike2, rebestyet;
bestelm *bestrees, **rebestrees;
long *place;
boolean mulf;
sequence convtab;
long renextree, nextree;

void setbottomtraverse(node *p);
static void oldsavetree(tree* t, long *place);
static void bintomulti(tree *t, node **root, node **binroot);
static void reroot3(tree* t, node *outgroup, node *root, node *root2, node *lastdesc);
static void backtobinary(tree* t, node **root, node *binroot);
static boolean outgrin(node *root, node *outgrnode);
static long get_numdesc(node* root, node* p);
static void moveleft(node *root, node *outgrnode, node **flipback);
static void flipnodes(node *nodea, node *nodeb);


node* root_tree(tree* t, node* here)
{
  /* Setup a root in a tree, a 3-node circle between  here  and  here->back.
   * This is useful for functions that expect a
   * rooted tree such as oldsavetree.  This does NOT reorient the tree so
   * that all nodep pointers point to the rootmost node of that interior
   * node.  */
  long k;
  node *nuroot1, *nuroot2, *nuroot3, *there;  /* the three nodes in the new circle
                                               * and a saving of  here->back */

  k = generic_tree_findemptyfork(t);
  there = here->back;
  nuroot1 = t->get_forknode(t, k+1);
  hookup(here, nuroot1);
  nuroot2 = t->get_forknode(t, k+1);
  nuroot1->next = nuroot2;
  hookup(nuroot2, there);
  nuroot3 = t->get_forknode(t, k+1);
  nuroot2->next = nuroot3;
  nuroot3->next = nuroot1;
  nuroot3->back = NULL;
  t->nodep[k] = nuroot3;
  return nuroot3;
} /* root_tree */


void reroot_tree(tree* t)
{
  /* Removes a root from a tree; useful after a return from functions that
   * expect a rooted tree (e.g. oldsavetree()). Then reroots the tree before
   * releasing a forknode to the freelist, to avoid tree components
   * pointing (even temporarily) into garbage. */
  long i;
  node *p, *fakeroot;

  fakeroot = t->root;     /* make sure it becomes the bottom node in circle */
  while (fakeroot->back != NULL)
    fakeroot = fakeroot->next;
  if ( count_sibs(fakeroot) > 2 )
  {                         /* debug:  doing what?  for multifurcation case */
    for (p = fakeroot ; p->next != fakeroot ; p = p->next)
      p->next = fakeroot->next;           /* bypass node fakeroot points to */
    if ( t->nodep[fakeroot->index - 1 ] == fakeroot)
      t->nodep[fakeroot->index - 1 ] = p;    /* have fakeroot point to fork */
  }
  else
  {     /* are only 2 sibs.  Be careful to hook two real sibs to each other */
    hookup(fakeroot->next->back, fakeroot->next->next->back); /* debug: always OK? */
  }
  if ( t->root == fakeroot) /* set root of tree if was pointing to fakeroot */
  {
    if (t->nodep[outgrno-1]->back != NULL)    /* if that tip is on the tree */
      t->root = t->nodep[outgrno - 1]->back;
    else {            /* find the tip of lowest number actually on the tree */
      for (i = 0; t->nodep[i]->back == NULL; i++) { }
      t->root = t->nodep[i]->back;                        /* set it to root */
      }
  }
  t->release_fork(t, fakeroot);
  t->root = t->nodep[outgrno - 1]->back;
} /* reroot_tree */


boolean pars_tree_try_insert_(tree * t, node * item, node * p, node * there,
                          double* bestyet, tree* bestree, boolean thorough,
                          boolean storing, boolean atstart, double* bestfound)
{
  /* insert item at p, if resulting tree has a better score, update bestyet
   * and there
   * This version actually does the hookups which are quickly dissolved,
   * however none of the changes are propegated in the tree and it is as
   * if it never got inserted. If we are on the last rearrangement save a
   * bestscoring insert to the bestrees array
   * item  should be an interior fork hooked to a tip or subtree which
   * is item->back */
  double like;
  boolean succeeded = false;
  node* dummy;
  boolean found = false;
  long pos = 0;

  t->save_traverses(t, item, p);     /* need to restore to leave tree same  */
  t->insert_(t, item, p->back, false);
  initializetrav(t, t->root);
  initializetrav(t, t->root->back);
  like = t->evaluate(t, p, false);
  t->score = like;
printf(" score = %lf, bestyet = %lf, bestfound = %lf", like, *bestyet, *bestfound); /* debug */
  if (like > *bestyet) {
    generic_tree_copy(t, bestree);
printf(" (new bestyet)");  /* debug */
    *bestyet = like;
    there = p;
  }
/* debug */ printf("\n");
  if (storing) {
    savetree(t, place);           /* make coded representation of this tree */
    if (atstart) {                    /* case when this is first tree tried */
      pos = 0;                       /* put it at the beginning of bestrees */
      found = false;
      if (nextree == 0) {
        *bestfound = like;
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (Initial)\n", like, *bestyet, *bestfound);  /* debug */
        addbestever(pos, &nextree, maxtrees, false, place, bestrees, like);
printf("Added an initial tree to bestrees, now %ld of them\n", nextree);
      }
      *bestyet = like;
      succeeded = true;
    } 
    else {
      if ( like == *bestfound )                 /* deciding on a later tree */
      {                /* find where it goes in numerical order in bestrees */
        findtree(&found, &pos, nextree, place, bestrees);
        succeeded = true;
        if (!found) {                  /* if found same tree, do not add it */
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (Tied)\n", like, *bestyet, *bestfound);  /* debug */
          addtiedtree(&pos, &nextree, maxtrees, false, place, bestrees, like);
printf("Added another tied tree to bestrees, now %ld of them\n", nextree);
        }
      } else {
        if (like > *bestfound) {       /* replacing all or adding one more? */
          *bestfound = like;
          *bestyet = like;
          pos = 0;                   /* put it at the beginning of bestrees */
          found = false;
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (Better)\n", like, *bestyet, *bestfound);  /* debug */
          addbestever(pos, &nextree, maxtrees, false, place, bestrees, like);
printf("Added new best tree to bestrees, score = %lf, now %ld of them\n", like, nextree);
          succeeded = true;
          *bestyet = like;
        }
      }
    }
  }
  if (succeeded)
  {
    there = p;
/* debug:    *multf = false;   */
  }
  t->re_move(t, item, t->root, true);       /* pull the branch back off the tree */
/* debug:  is preceding statement correct?  &dummy?  */
  t->restore_traverses(t, item, p);           /* debug: what is tis doing? */
  t->evaluate(t, p, 0);   /* debug:   as in dnaml, but may not be needed */


  found = false;                /* debug: why this? May not have any effect */
  pos = 0;

  /* debug:  Uncommenting the following code will allow for a multifurcating
   * search, However you will generally find every multifurcation as a
   * separate tree including ambiguous ones.  */
#if 0
  if ( p->tip == false )
  {
    t->insert_(t, item, p, false);
    like = t->evaluate(t, p, false);
    if (like >= *bestyet || *bestyet == UNDEFINED)
    {
      *multf = true;
      *there = p;
      if ( like > *bestyet )
        succeeded  = true;
      if ( lastrearr )
      {
        savetree(t, place);
        findtree(&found, &pos, nextree, place, bestrees);
        if ( ((like > *bestyet) && !found) || nextree == 1)
          addbestever(pos, &nextree, maxtrees, false, place, bestrees);
        else if ( !found )
          addtiedtree(pos, &nextree, maxtrees, false, place, bestrees);
      }
      *bestyet = like;
    }
    t->re_move(t, item, p, true);
  }
#endif

  return succeeded;
} /* pars_tree_try_insert */


tree* pars_tree_new(long nonodes, long spp)
{
  /* make a new pars_tree */
  tree* t = Malloc(sizeof(pars_tree));

  generic_tree_init(t, nonodes, spp);
  pars_tree_init(t, nonodes, spp);
  return t;
} /* pars_tree_new */


void pars_tree_init(tree* t, long nonodes, long spp)
{
  /* setup of a new tree  with  spp  tips */

  t->globrearrange = pars_globrearrange;
  t->try_insert_ = (tree_try_insert_t)pars_tree_try_insert_;
  t->evaluate = pars_tree_evaluate;
} /* pars_tree_init */


void pars_node_init(node* p, node_type type, long index)
{
  /* set up a new node for a pars_tree */
  pars_node *pn = (pars_node *)p;

  generic_node_init(p, type, index);
  p->copy = pars_node_copy;
  p->init = pars_node_init;
  p->reinit = pars_node_reinit;
  p->free = pars_node_free;
  p->node_print_f = pars_node_print;

  if (pn->numsteps)
    free(pn->numsteps);
  pn->numsteps = Malloc(endsite * sizeof(long));
} /* pars_node_init */


void pars_node_reinit(node * n)
{
  /* re-setup a pars_tree node */
  generic_node_reinit(n);
  pars_node *pn = (pars_node *)n;
  if (pn->numsteps)
    free(pn->numsteps);
  pn->numsteps = Malloc(endsite * sizeof(long));
} /* pars_node_reinit */


void pars_node_print(node * n)
{
  /* print out steps for a pars_tree node
   * is this just for debugging? */
  generic_node_print(n);
  pars_node * pn = (pars_node*)n;
  if(pn->numsteps == NULL) printf(" numsteps:<empty>");
  else
  {
    long i;
    printf(" numsteps:");
    for(i = 0; i < endsite; i++)
    {
      printf(" %ld", pn->numsteps[i]);
    }
  }
} /* pars_node_print */


void pars_node_free(node **pp)
{ 
  /* free a node from a pars_tree */
  pars_node *pn = (pars_node *)*pp;
  free(pn->numsteps);
  generic_node_free(pp);
} /* pars_node_free */


void pars_node_copy(node* srcn, node* dstn)
{ 
  /* copy a pars_tree node */
  pars_node *src = (pars_node *)srcn;
  pars_node *dst = (pars_node *)dstn;

  generic_node_copy(srcn, dstn);
  if (dst->numsteps == NULL )
    dst->numsteps = Malloc(endsite * sizeof(long));
  memcpy(dst->numsteps, src->numsteps, endsite * sizeof(long));
} /* pars_node_copy */


void collapsebestrees(tree *t, bestelm *bestrees, long *place, long chars,
                       boolean progress, long *finalTotal)
{
  /* Goes through all best trees, collapsing trees where possible,
   * and deleting trees that are not unique.    */
  long i, j, k, pos ;
  boolean found;
  long treeLimit = nextree < maxtrees ? nextree : maxtrees;

  for(i = 0 ; i < treeLimit ; i++)
  {
    bestrees[i].collapse = true;
  }
  if(progress)
  {
    sprintf(progbuf, "\nCollapsing best trees\n   ");
    print_progress(progbuf);
  }
  k = 0;
  for(i = 0 ; i < treeLimit ; i++) {
    if(progress)
    {
      if(i % ((treeLimit / 72) + 1) == 0)
      {
        sprintf(progbuf, ".");
        print_progress(progbuf);
      }
    }
    while(!bestrees[k].collapse)
      k++;
    load_tree(t, k, bestrees);                         /* Reconstruct tree. */
    while ( treecollapsible(t, t->nodep[0]))
      collapsetree(t, t->nodep[0]);
    savetree(t, place);          /* set aside collapsed tree in place array */
    if ( k != (treeLimit-1) ) {          /* if not at the last tree already */
      for (j = k ; j < (treeLimit - 1) ; j++) /* move rest of trees forward */
      {
        memcpy(bestrees[j].btree, bestrees[j + 1].btree, spp * sizeof(long));
        bestrees[j].gloreange = bestrees[j + 1].gloreange;
        bestrees[j + 1].gloreange = false;
        bestrees[j].locreange = bestrees[j + 1].locreange;
        bestrees[j + 1].locreange = false;
        bestrees[j].collapse = bestrees[j + 1].collapse;
      }
    }
    treeLimit--;         /* because there is now one fewer tree in bestrees */
    pos = 0;
    findtree(&found, &pos, treeLimit, place, bestrees);

    if (!found)      /* put the new tree in the the list if it wasn't found */
    {
      addtree(pos, &treeLimit, false, place, bestrees);
    }
  }
  if (progress)
  {
    sprintf(progbuf, "\n");
    print_progress(progbuf);
    phyFillScreenColor();
  }
  *finalTotal = treeLimit;
} /* collapsebesttrees */


static long get_numdesc(node* root, node* p)
{
  /* we used to bookkeep a numdesc variable,  this is no longer
   * necessary, however some older functions still like it. */
  if ( p->tip )
    return 0;
  if ( root->index == p->index && root != p)
    return count_sibs(p) - 1;
  if ( (p != root && p->back == NULL) ||
       (p->next->back == NULL && p->next != root))
    return 0;
  return count_sibs(p);
} /* get_numdesc */


void reroot(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. used if
   * the root is binary.
   * used in dnacomp & dnapars */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* reroot */


void reroot2(node *outgroup, node *root)
{
  /* reorients tree, putting outgroup in desired position. */
  /* used in dnacomp & dnapars */
  node *p;

  p = outgroup->back->next;
  while (p->next != outgroup->back)
    p = p->next;
  root->next = outgroup->back;
  p->next = root;
}  /* reroot2 */


void reroot3(tree* t, node *outgroup, node *root, node *root2, node *lastdesc)
{
  /* reorients tree, putting back outgroup in original position.
   * used in dnacomp & dnapars */
  node *p;

  p = root->next;
  while (p->next != root)
    p = p->next;
  t->release_forknode(t, root);
  p->next = outgroup->back;
  root2->next = lastdesc->next;
  lastdesc->next = root2;
}  /* reroot3 */


static void bintomulti(tree *t, node **root, node **binroot)
{
  /* Make a binary tree multifurcating:
   * attaches root's left child to its right child and makes the right
   * child the new root. */
  node *left, *right, *newnode, *temp;

  right = (*root)->next->next->back;
  left = (*root)->next->back;
  if (right->tip)
  {
    (*root)->next = right->back;
    (*root)->next->next = left->back;
    temp = left;
    left = right;
    right = temp;
    right->back->next = *root;
  }
  newnode = t->get_forknode(t, right->index);
  newnode->next = right->next;
  newnode->back = left;
  left->back = newnode;
  right->next = newnode;
  (*root)->next->back = (*root)->next->next->back = NULL;
  *binroot = *root;
  *root = right;
  (*root)->back = NULL;
} /* bintomulti */


void oldsavetree(tree* t, long *place)
{
   /* record in array  place  where each species has to be
    * added to reconstruct this tree. This code assumes a root
    * this is the older function,  a new function roots the tree and calls this
    * function to save the tree */
  long i, j, nextnode, nvisited, newforknum;
  long* lineagenumber;
  node *p, *q, *r = NULL, *root2, *lastdesc, *outgrnode, *binroot, *flipback;
  boolean done, newfork;

  binroot = NULL;
  lastdesc = NULL;
  root2 = NULL;
  flipback = NULL;
  outgrnode = t->nodep[outgrno - 1];
  lineagenumber = (long *)Malloc(nonodes*sizeof(long));
  setbottomtraverse(t->root);  /* set booleans indicating which way is down */
  nextnode = spp + 1;
  for (i = 0; i < spp; i++)                     /* initialize "place" array */
    place[i] = 0;
  for (i = 0; i < nonodes; i++)          /* which lineage each tree node is */
    lineagenumber[i] = 0;                          /* ... starts out zeroed */
  place[0] = 1;                                     /* this one is always 1 */
  lineagenumber[0] = 1;                        /* first lineage is number 1 */
  for (i = 1; i <= spp; i++)                           /* for each tip, ... */
  {
    p = t->nodep[i - 1];                           /* start with species  i */
    if (p->back != NULL) {                    /* if this tip is in the tree */
      p = p->back;               /* go to the interior node connected to it */
      while (lineagenumber[p->index - 1] == 0)    /* if no number yet there */
      {
        lineagenumber[p->index - 1] = i;  /* set value to index of that tip */
        while (!p->bottom)             /* go around circle to find way down */
          p = p->next;
        p = p->back;                             /* go down to earlier fork */
        if (p == NULL)                   /* if we went past bottom fork ... */
          break;                                 /* blast out of while loop */
      }
      if (p != NULL) {              /* we ran into a nonzero lineage number */
        place[i-1] = lineagenumber[p->index -1];   /* record in place array */
        newforknum = spp + i - 1;       /* number of new fork when attaches */
        while (lineagenumber[p->index - 1] == place[i-1])  /* going on down */
        {
          lineagenumber[p->index - 1] = newforknum;       /* from here down */
          while (!p->bottom)           /* go around circle to find way down */
            p = p->next;
          if (p->back == NULL)         /* blast out of loop if reached root */
            break;
          else
            p = p->back;                         /* go down to earlier fork */
        }
      }
    }
    if (i > 1)         /* this is for dealing with multifurcations, somehow */
    {
#if 0                  /* debug: commenting out for now.  Not sure how this */
      q = t->nodep[i - 1];  /* ... works, how  r  is initialized */
      newfork = true;
      nvisited = sibsvisited(q, place);
      if (nvisited == 0)
      {
        if (parentinmulti(r, t->root))
        {
          nvisited = sibsvisited(r, place);
          if (nvisited == 0)
            place[i - 1] = place[p->index - 1];
          else if (nvisited == 1)
            place[i - 1] = smallest(r, place);
          else
          {
            place[i - 1] = -smallest(r, place);
            newfork = false;
          }
        }
        else
          place[i - 1] = place[p->index - 1];
      }
      else if (nvisited == 1)
      {
        place[i - 1] = place[p->index - 1];
      }
      else
      {
        place[i - 1] = -smallest(q, place);
        newfork = false;
      }
      if (newfork)
      {
        j = place[p->index - 1];
        done = false;
        while (!done)
        {
          place[p->index - 1] = nextnode;
          while (!p->bottom) {
            p = p->next;
          }
          p = p->back;
          done = (p == NULL);
          if (!done)
            done = (place[p->index - 1] != j);
          if (done)
          {
            nextnode++;
          }
        }
      }
#endif     /* debug: end of commented-out code */
    }
    if (flipback)
      flipnodes(outgrnode, flipback->back);
    else
    {
      if (root2)
      {
        reroot3(t, outgrnode, t->root, root2, lastdesc);
        t->root = root2;
      }
    }
    if (binroot)
      backtobinary(t, &t->root, binroot);
  }
}  /* oldsavetree */


void savetree(tree* t, long *place)
{
  /* Record in  place  where each species has to be added to reconstruct
   * this tree. This code finds out whether there is already a root fork
   * that has a NULL ancestor, if not it  roots the tree.  Then it calls
   * and calls oldsavetree to save it.  Then, if necessary, it removes
   * the new temporary root fork */
  boolean wasrooted;
  node *oldroot, *p, *q, *outgrnode;

  wasrooted = false;
  outgrnode = t->nodep[outgrno - 1];
  p = outgrnode->back;
  if (p != NULL) {
    wasrooted = (p->back == NULL);  /* Check: was it already rooted by ... */
    q = p;                    /* going around circle to see whether it was */
    if (!wasrooted) {
      q = q->next;
      while (q != p) {
        if (q->back == NULL) {
          wasrooted = true;
          oldroot = q;
          break;
        }
        q = q->next;
      }
    }
  }
  if (!wasrooted) {
    oldroot = t->nodep[outgrno-1];
    t->root = root_tree(t, oldroot);           /* put in a "fake" root fork */
  }
  oldsavetree(t, place);        /* now save this rooted tree in array place */
  if (!wasrooted) {                /* remove the fake root if one was added */
    reroot_tree(t); 
    t->root = oldroot;
  }
}  /* savetree */


void addbestever(long pos, long *nextree, long maxtrees, boolean collapse,
                  long *place, bestelm *bestrees, double score)
{
  /* adds first best tree. If we are rearranging on usertrees, 
   * add it to the second array of trees if the score is good enough
   * pos is the position where it will be added which is 0   */
  long repos;
  boolean found;

  pos = 0;
  *nextree = 0;

  addtree(pos, nextree, collapse, place, bestrees);
  if ( reusertree )
  {
    if ( score == UNDEFINED ) return;
    if ( score != UNDEFINED && (score > rebestyet || rebestyet == UNDEFINED))
    {
      renextree = 1;
      rebestyet = score;
      addtree(1, &renextree, collapse, place, rebestrees[1]);  /* debug: correct? */
      renextree = 1;
    }
    else if ( score != UNDEFINED && score == rebestyet )
    {
      findtree(&found, &repos, renextree, place, rebestrees[1]);
      if ( !found && renextree <= maxtrees )
        addtree(repos, &renextree, collapse, place, rebestrees[1]);
    }
  }
} /* addbestever */


void addtiedtree(long* pos, long *nextree, long maxtrees, boolean collapse, long *place, bestelm *bestrees, double score)
{
  /* add a tied tree.   pos is the position in the range  0 .. (nextree-1) */
  boolean found;
  long repos;

  if (*nextree <= maxtrees)
    addtree(*pos, nextree, collapse, place, bestrees);
  if ( reusertree )    /* debug:  this part needs more debugging */
  {
    if ( rebestyet == score )
    {
      findtree(&found, &repos, renextree, place, rebestrees[1]);
      if ( !found && renextree <= maxtrees )
        addtree(repos, &renextree, collapse, place, rebestrees[1]);
    }
  }
} /* addtiedtree */


void add_to_besttrees(tree* t, long score, bestelm* bestrees,
                       double* bestfound)
{
  /* take the tree we have found and try to add it to the array bestrees:
   * if none are already there, make it the first one, if it is better than
   * the ones that are there then toss them and start over with just this
   * one, if tied with them add it in too */
/* debug:  may not need in view of pars_tree_try_insert  */

  boolean found = false;
  long *pos = 0;
  
  if ( (nextree = 1) || !(score < *bestfound)) {    /* if save this one ... */
    savetree(t, place);
    if (score > *bestfound) {        /* if it will be the lone new best one */
      addbestever(*pos, &nextree, maxtrees, false, place, bestrees, score);
    } else {                            /* it is another tree tied for best */
      findtree(&found, pos, nextree-1, place, bestrees);  /* already there? */
      if (!found)                      /* save it only if not already there */
        addtiedtree(pos, &nextree, maxtrees, false, place, bestrees, score);
    }
  }
} /* add_to_besttrees */


boolean pars_addtraverse(tree* t, node* p, node* q, boolean contin,
                         node* qwherein, double* bestyet, bestelm* bestrees,
                         boolean thorough, boolean storing, boolean atstart,
                         double* bestfound)
{
  /* wrapper for addraverse, calling generic addtraverse
   * and then taking the best tree found and adding it to the array
   * of tied best trees found. Function like this works for parsimony-like
   * criteria where there are exact ties, not for likelihood or distance
   * criteria */
/* debug:  not yet called from anywhere */
   boolean success;

/* debug:  does this make any sense?  Already saving best tree yet in generic version
   success = generic_tree_addtraverse(t, p, q, contin, qwherein,
                   bestyet, &bestree, thorough, storing, atstart, bestfound);
   add_to_besttrees(t, t->score, bestrees, bestfound);     debug */
   return success;
} /* pars_addtraverse */


void flipnodes(node *nodea, node *nodeb)
{
  /* flip nodes */
  node *backa, *backb;

  backa = nodea->back;
  backb = nodeb->back;
  backa->back = nodeb;
  backb->back = nodea;
  nodea->back = backb;
  nodeb->back = backa;
} /* flipnodes */


void moveleft(node *root, node *outgrnode, node **flipback)
{
  /* makes outgroup node to leftmost child of root */
  node *p;
  boolean done;

  p = root->next;
  done = false;
  while (p != root && !done)
  {
    if (p->back == outgrnode)
    {
      *flipback = p;
      flipnodes(root->next->back, p->back);
      done = true;
    }
    p = p->next;
  }
} /* moveleft */


void printbranchlengths(node *p)
{ 
  /* print branch lengths */
  node *q;
  long i;

  if (p->tip)
    return;
  q = p->next;
  do {
    fprintf(outfile, "%6ld      ", q->index - spp);
    if (q->back->tip)
    {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->back->index - 1][i], outfile);
    }
    else
      fprintf(outfile, "%6ld    ", q->back->index - spp);
    fprintf(outfile, "   %f\n", q->v);
    if (q->back)
      printbranchlengths(q->back);
    q = q->next;
  } while (q != p);
} /* printbranchlengths */


void initbranchlen(node *p)
{
  /* initial values of branch lengths */
  node *q;

  p->v = 0.0;
  if (p->back)
    p->back->v = 0.0;
  if (p->tip)
    return;
  q = p->next;
  while (q != p)
  {
    initbranchlen(q->back);
    q = q->next;
  }
  q = p->next;
  while (q != p)
  {
    q->v = 0.0;
    q = q->next;
  }
} /* initbranchlen */


boolean alltips(node *forknode, node *p)
{
  /* returns true if all descendants of forknode except p are tips;
     false otherwise.  */
  node *q, *r;
  boolean tips;

  tips = true;
  r = forknode;
  q = forknode->next;
  do {
    if (q->back && q->back != p && !q->back->tip)
      tips = false;
    q = q->next;
  } while (tips && q != r);
  return tips;
} /* alltips */


void flipindexes(long nextnode, pointarray treenode)
{
  /* flips index of nodes between nextnode and last node.
   * This is intended to move any empty fork to be numerically last  */
  long last;
  node *temp;

  last = nonodes;
  while (treenode[last - 1]->back == NULL)
    last--;                     /* go earlier in forks to find nonempty one */
  if (last > nextnode)       /* swap it with place where next fork is to be */
  {
    temp = treenode[nextnode - 1];
    treenode[nextnode - 1] = treenode[last - 1];
    treenode[last - 1] = temp;
    newindex(nextnode, treenode[nextnode - 1]);
    newindex(last, treenode[last - 1]);
  }
} /* flipindexes */


long sibsvisited(node *anode, long *place)
{
  /* computes the number of nodes which are visited after anode among
     the siblings of its fork circle */
  node *p;
  long nvisited;

  while (!anode->bottom)                   /* go around circle, find bottom */
    anode = anode->next;
  p = anode->back->next; /* go down to ancestor, and on to its next sibling */
  nvisited = 0;
  do {       /* how many of those aunts/uncles are tips already encountered */
    if (!p->bottom && place[p->back->index - 1] != 0)
      nvisited++;
    p = p->next;
  } while (p != anode->back);
  return nvisited;
}  /* sibsvisited */


boolean parentinmulti(node *anode, node* root)
{
  /* sees if anode's parent has more than 2 children */
  node *p;

  while (!anode->bottom)
    anode = anode->next;
  p = anode->back;
  while (!p->bottom)
    p = p->next;
  return (get_numdesc(root, p) > 2);
} /* parentinmulti */


long smallest(node *anode, long *place)
{
  /* finds the smallest index of sibling of anode */
  node *p;
  long min;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  if (p->bottom) p = p->next;
  min = nonodes;
  do {
    if (p->back && place[p->back->index - 1] != 0)
    {
      if (p->back->index <= spp)
      {
        if (p->back->index < min)
          min = p->back->index;
      }
      else
      {
        if (place[p->back->index - 1] < min)
          min = place[p->back->index - 1];
      }
    }
    p = p->next;
    if (p->bottom) p = p->next;
  } while (p != anode->back);
  return min;
}  /* smallest */


void backtobinary(tree* t, node **root, node *binroot)
{ /* restores binary root */
  node *p;

  binroot->next->back = (*root)->next->back;
  (*root)->next->back->back = binroot->next;
  p = (*root)->next;
  (*root)->next = p->next;
  binroot->next->next->back = *root;
  (*root)->back = binroot->next->next;
  t->release_forknode(t, p);
  *root = binroot;
} /* backtobinary */


void newindex(long i, node *p)
{
  /* assigns index i to fork that  p is in */

  while (p->index != i)
  {
    p->index = i;
    p = p->next;   /* ... and move on around circle. */
  }
} /* newindex */


void load_tree(tree* t, long treei, bestelm* bestrees)
{
  /* restores tree  treei  from array bestrees (treei is the index
   * of the array, so tree 5 has treei = 4).  Add all the tips to a tree one
   * by one in order.  The array element  bestree[treei].btree[k]  indicates
   * that the k-th tip is to be connected to a new fork that is below tip or
   * fork  btree[k].  If negative, it indicates that it is to be added as an
   * extra furc to the fork that is already at the bottom of the branch that
   * is below fork (or tip) abs(btree[k])  */
  long i, j, nsibs;
  boolean foundit = false;
  node *p, *q, *below, *bback, *forknode, *newtip, *beforewhere, *afterwhere;

  release_all_forks(t);              /* to make sure all interior nodes
                                        are on the list at free_fork_nodes  */
                                     /* then make tree of first two species */
  forknode = t->get_fork(t, spp);      /* was put on nodep, index is  spp+1 */
  hookup(t->nodep[1], forknode->next);
  hookup(t->nodep[0], forknode->next->next);

  for ( j = 3; j <= spp ; j++ )     /* adding one by one species, 3, 4, ... */
  {
    newtip = t->nodep[j-1];

    if ( bestrees[treei].btree[j-1] > 0 )    /* j-th entry in "place" array */
    {                                                    /*  if bifurcation */
      forknode = t->get_fork(t, spp+j-2);       /* put a new fork circle in */
      hookup(newtip, forknode);
      below = t->nodep[bestrees[treei].btree[j - 1] - 1];
      bback = below->back;
      hookup(forknode->next, below);
      if ( bback )
        hookup(forknode->next->next, bback);
      t->nodep[spp+j-2] = forknode->next->next;   /* know which way is down */
    }
    else
    {          /*  if goes into a multifurcation put a new node into circle */
      bback= t->nodep[t->nodep[-bestrees[treei].btree[j-1]-1]->back->index-1];
      beforewhere = bback->next;            /* where will that node be put? */ 
      afterwhere = bback->next->next;   /* ... the nodes before, after that */
      do {   /* move around fork circle until just before the downward link */
        beforewhere = afterwhere;
        afterwhere = afterwhere->next;
      } while ((afterwhere->next) != bback);
      forknode = t->get_forknode(t, below->index);        /* get a new node */
      hookup(newtip, forknode);             /* hook the tip to the new node */
      beforewhere->next = forknode;               /* put it the right place */
      forknode->next = afterwhere;
    }
  }

  forknode = NULL;
  for (i = spp; i < nonodes; i++) {      /* check all interior node circles */
    p = t->nodep[i];
    if (p != NULL) {
      q = p;
      do {
        if (q->back == NULL) {
          forknode = q;            /* find a node that has nothing below it */
          foundit = true;
          }
        q = q->next;
      } while (q != p); 
    }
  }
  if (foundit) {    /* remove the interior node which has an empty neighbor */
    nsibs = count_sibs(forknode); 
    if ( nsibs > 2 )
    {                            /* find the circle member that precedes it */
      for ( q = forknode ; q->next != forknode ; q = q->next);
      q->next = forknode->next;                      /* and connect past it */
      t->nodep[q->index - 1] = q;           /* and have nodep point to that */
      t->release_forknode(t, forknode);                      /* and toss it */
    }
    else
    {                   /* if the interior node has only two real neighbors */
      hookup(forknode->next->back, forknode->next->next->back);
      t->release_fork(t, forknode);        /* release the whole fork circle */
    }
  }

  t->root = t->nodep[outgrno - 1]->back;
} /* load_tree */


void setbottomtraverse(node *p)
{ 
  /* set boolean "bottom" to true on one of the nodes in a fork circle
   * at each interior fork to show which way is down.  Go around
   * fork circle and set others to false.  Traverse out through tree
   * to do this on all nodes.  */
  node *q;

  p->bottom = true;                       /* set the one you arrive at true */
  if (p->tip)
    return;
  q = p->next;
  while (q != p)                  /* go around circle, set all others false */
  {
    q->bottom = false;
    setbottomtraverse(q->back);
    q = q->next;
  }
}  /* setbottomtraverse */


boolean outgrin(node *root, node *outgrnode)
{
  /* checks if outgroup node is a child of root */
  node *p;

  p = root->next;
  while (p != root)
  {
    if (p->back == outgrnode)
      return true;
    p = p->next;
  }
  return false;
} /* outgrin */


void pars_globrearrange(tree* curtree, tree* bestree, boolean progress,
                         boolean thorough, double* bestfound)
{ /* does global rearrangements */
  /* The more general generic_unrooted_locrearrange also works but this is
   * much faster because it gets to take advantage of some of the speedups
   * available in the parsimony programs */
  int i, j, k, num_sibs, num_sibs2;
  node *where, *sib_ptr, *sib_ptr2, *qwhere;
  double bestyet;
  boolean mulf;
  node* removed;

  bestyet = curtree->evaluate(bestree, bestree->root, 0);

  if (progress)
  {
    sprintf(progbuf, "   ");
    print_progress(progbuf);
  }

  for ( i = 0 ; i < curtree->nonodes ; i++ )
  {
    sib_ptr  = curtree->nodep[i];
    if (sib_ptr == NULL)
      continue;          /* skip this case if no interior node circle there */
    if ( sib_ptr->tip )
      num_sibs = 0;
    else
      num_sibs = count_sibs(sib_ptr);

    if ( progress && ((i - spp) % (( curtree->nonodes / 72 ) + 1 ) == 0 ))
    {
      sprintf(progbuf, ".");
      print_progress(progbuf);
    }

    for ( j = 0 ; j <= num_sibs ; j++ )
    {
      sib_ptr = curtree->nodep[i];
      for ( k = 0 ; k < j ; k++ ) {
        sib_ptr = sib_ptr->next;
        if ( sib_ptr->back == NULL || sib_ptr->back->tip )
          continue;                     /* skip over rest of loop this time */

        removed = sib_ptr;   /* debug: or is it  sib_ptr->back ? */
        mulf = 2 != count_sibs(removed->back);
        curtree->re_move(curtree, removed, &where, true);
/* printf(" remove %ld:%ld\n", removed->index, removed->back->index);   debug */
        qwhere = where;

        if ( where->tip)
        {
          num_sibs2 = 0;
          sib_ptr2 = where->back;
        }
        else
        {
          num_sibs2 = count_sibs(where);
          sib_ptr2 = where;
        }
        for ( k = 0 ; k <= num_sibs2 ; k++ )
        {
          curtree->addtraverse(curtree, removed, sib_ptr2->back, true,
                     qwhere, &bestyet, bestree, true, true, false, bestfound);
          sib_ptr2 = sib_ptr2->next;
        }
        curtree->insert_(curtree, removed, qwhere, false);
        curtree->root = curtree->nodep[0]->back;
/* debug: why?        bestyet = curtree->evaluate(curtree, curtree->root, 0);   debug */
      }
    }
  }
  if (progress)
  {
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
} /* pars_globrearrange */


boolean treecollapsible(tree* t, node* n)
{
 /* find out whether there is any collapsible branch on the tree */
  node *sib;
  boolean collapsible = false;

  if ( ((pars_tree*)t)->branchcollapsible(t, n) )
    return true;

  if ( n->back->tip == true )
    return false;

  for ( sib = n->back->next ; sib != n->back ; sib = sib->next )
  {
    collapsible =  treecollapsible(t, sib) || collapsible;
  }
  return collapsible;
} /* treecollapsible */


void collapsebranch(tree* t, node* n)
{ /* remove a branch and merge the forks at both ends */
  node* m, *sib, *newfork;
  long nsibs;

  m = n->back;

  nsibs = count_sibs(m);
  (void)nsibs;               // RSGnote: Variable set but never referenced.

  for ( sib = m->next ; sib != m ; sib = sib->next )
  {
    if ( sib == m->next )
      newfork = n;
    else
    {
      newfork = t->get_forknode(t, n->index);
      newfork->next = n->next;
      n->next = newfork;
    }
    hookup(sib->back, newfork);
  }
  t->release_fork(t, m); 

  inittrav(t, n);   /* now make initialized pointer false (ones looking in? */
  inittrav(t, n->back);
  n->initialized = false;
} /* collapsebranch */


void collapsetree(tree* t, node* n)
{  /* collapse all branches that are designated as collapsible */
  node *sib;
  if ( ((pars_tree*)t)->branchcollapsible(t, n) )
    collapsebranch(t, n);
  if ( n->back->tip == true)
    return;
  for ( sib = n->back->next ; sib != n->back ; sib = sib->next )
  {
    collapsetree(t, sib);
  }
} /* collapsetree */


void printree(tree* t)
{ /* prints out diagram of the tree */
  long tipy;
  double scale, tipmax;
  long i;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(t, t->root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline3(i, scale, t->root);
  fprintf(outfile, "\n  remember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
  putc('\n', outfile);
}  /* printree */


void coordinates(tree* t, node *p, double lengthsum, long *tipy,
                 double *tipmax)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

  if (p == NULL)
    return;
  if (p->tip)
  {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    xx = q->v;
    if (xx > 100.0)
      xx = 100.0;
    coordinates(t, q->back, lengthsum + xx, tipy, tipmax);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if ((p == t->root) || count_sibs(p) > 2)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline3(long i, double scale, node *start)
{
  /* draws one row of the tree diagram by moving up tree */
  /* used in dnapars */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = start;
  q = start;
  extra = false;
  if (i == (long)p->ycoord)
  {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
  }
  else
    fprintf(outfile, "  ");
  do {
    if (!p->tip)
    {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax)
        {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (r == p)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done)
    {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip)
      {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      }
      else
      {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    }
    else if (!p->tip)
    {
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord || p == start))
      {
        putc('|', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      }
      else
      {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    }
    else
    {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip)
  {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline3 */


void writesteps(tree* t, long chars, boolean weights, steptr oldweight)
{
  /* used in dnacomp, dnapars, & dnapenny */
  long i, j, k, l;
  k=0;

  /*calculate the steps */
  if (t->root->initialized == false ) t->nuview(t, t->root);

  /* print them */
  putc('\n', outfile);
  if (weights)
    fprintf(outfile, "weighted ");
  fprintf(outfile, "steps in each site:\n");
  fprintf(outfile, "      ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%4ld", i);
  fprintf(outfile, "\n     r------------------------------------");
  fprintf(outfile, "-----\n");
  for (i = 0; i <= (chars / 10); i++)
  {
    fprintf(outfile, "%5ld", i * 10);
    putc('|', outfile);
    for (j = 0; j <= 9; j++)
    {
      k = i * 10 + j;
      if (k == 0 || k > chars)
        fprintf(outfile, "    ");
      else
      {
        l = location[ally[k - 1] - 1];
        if (oldweight[k - 1] > 0)
          fprintf(outfile, "%4ld",
                   oldweight[k - 1]
                    * (((pars_node*)t->root)->numsteps[l - 1] / weight[l - 1]));
        else
          fprintf(outfile, "%4ld", (((pars_node*)t->root)->numsteps[k - 1] ));
      }
    }
    putc('\n', outfile);
  }
} /* writesteps */


void grandrearr(tree* t, tree* bestree, boolean progress,
                 boolean rearrfirst, double* bestfound)
{
  /* calls "global" (SPR) rearrangement on best trees */
  long treei;
  long i, oldbestyet;
  boolean done = false;

  lastrearr = true;
/* debug:  whay do this, best trees are in array bestrees and  t  is not necessarily one of them
  savetree(t, place);
  addbestever(pos, &nextree, maxtrees, false, place, bestrees, UNDEFINED);
debug:   */

  for ( i = 0 ; i <= nextree-1 ; i++)    /* set all saved trees as as-yet */
    bestrees[i].gloreange = false;       /* globally un-rearranged */

  oldbestyet = bestree->score;
  while (!done)
  {
    treei = findunrearranged(bestrees, nextree, true); /* find one not done */
    if (treei < 0)
      done = true;
    else
      bestrees[treei].gloreange = true;

    if (!done)
    {
      load_tree(t, treei, bestrees);                /* reconstruct the tree */
      bestyet = t->evaluate(t, t->root, 0);                /* get its score */
      t->globrearrange(t, bestree, progress, true, bestfound); /* rearrange */
      done = rearrfirst || (oldbestyet == bestyet);    /* if not any better */
    }
  }
} /* grandrearr */


void treeout3(node *p, long nextree, long *col, long indent, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnapars -- writes branch lengths */
  node *q;
  long i, n, w;
  double x;
  Char c;

  if (p == root)
    indent = 0;
  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)
    {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)
    {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  }
  else
  {
    putc('(', outtree);
    (*col)++;
    indent++;                         /* increment amount of line indent */
    q = p->next;
    while (q != p)
    {
      treeout3(q->back, nextree, col, indent, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60)
      {
        putc('\n', outtree);
        *col = 0;
        for (i = 1; i <= indent; i++) /* write indent at beginning of line */
          putc(' ', outtree);
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p != root)
  {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    *col += w + 8;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout3 */


void initparsnode(tree *treep, node **p, long len, long nodei, long *ntips, long *parens, initops whichinit, pointarray treenode, Char *str, Char *ch, FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit)
  {
    case bottom:
      *p = treep->get_forknode(treep, nodei);
      treenode[nodei - 1] = *p;
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      break;
    case tip:
      match_names_to_data (str, treenode, p, spp);
      break;
    case length:         /* if there is a length, read it and discard value */
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      break;
    default:            /*cases hslength, hsnolength, treewt, unittrwt, iter, */
      break;
  }
} /* initparsnode */


double pars_tree_evaluate(tree* t, node*p, boolean dummy)
{
  /* not used but could be, would work for all parsimony programs
   * slower than a native version */
  node *root  = NULL;
  node *left  = NULL;
  node *right = NULL;
  long i;
  double sum = 0, steps;
  boolean newfork = false;

  generic_tree_evaluate(t, p, dummy);

  if (p->back )
  {
    newfork = true;
    /* make a new ring of nodes */
    root = t->get_forknode(t, t->nonodes);
    left = t->get_forknode(t, t->nonodes);
    right = t->get_forknode(t, t->nonodes);

    root->next = left;
    left->next = right;
    right->next = root;

    /* graft onto the tree */
    left->back = p;
    right->back = p->back;
  }
  else
  {
    root = p; /* in case this tree is already rooted */
  }

  /* get an updated view and count the steps */
  t->nuview(t, root);
  for ( i = 0 ; i < endsite ; i++ )
  {
    steps = ((pars_node*)root)->numsteps[i];
    if (((pars_tree*)t)->supplement)
      steps += ((pars_tree*)t)->supplement(t, i);
    if ( steps > threshwt[i] )
      steps  = threshwt[i];
    sum += steps;
  }

  if ( newfork )
  {
    t->release_forknode(t, root);
    t->release_forknode(t, left);
    t->release_forknode(t, right);
  }
  t->score = -sum;
  return -sum;
} /* pars_tree_evaluate */


bestelm* allocbestree(void)
{ 
  /* Malloc space for bestelm and bestrees */
  long i;
  bestelm* bestrees;

  bestrees = (bestelm *)Malloc(maxtrees * sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes * sizeof(long));
  return bestrees;
} /* allocbestree */


bestelm** allocbestrees(void)
{
  /* alloc space for array of bestrees */
  long i;
  bestelm **rebestrees = Malloc(2 * sizeof(bestelm*));

  for ( i = 0 ; i < 2 ; i++ )
    rebestrees[i] = allocbestree();
  return rebestrees;
} /* allocbestrees */


/* End. */
