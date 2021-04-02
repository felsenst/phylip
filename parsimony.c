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
  node *nuroot1, *nuroot2, *nuroot3, *there;  /* the three nodes in new circle
                                               * and a saving of here->back */
  k = generic_tree_findemptyfork(t);
  there = here->back;
  nuroot1 = t->get_forknode(t, k+1);
  hookup(here, nuroot1);
  nuroot2 = t->get_forknode(t, k+1);
  nuroot1->next = nuroot2;
  if (there != NULL)
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
  {                       /* pulling off bottom node in multifurcation case */
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


boolean pars_tree_try_insert_(tree* t, node* item, node* p, node* there,
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
  long i, pos = 0;

  t->save_traverses(t, item, p);     /* need to restore to leave tree same  */
/* printf("try inserting %ld on %ld:%ld\n", item->index, p->back->index, p->index); debug */
  t->insert_(t, item, p->back, false);
  initializetrav(t, t->root);                /* make sure updates all views */
  initializetrav(t, t->root->back);
/* printf("\n"); debug */
  like = t->evaluate(t, p, false);
/* printf("\n"); debug */
  t->score = like;
printf(" score = %lf, bestyet = %lf, bestfound = %lf\n", like, *bestyet, *bestfound); /* debug */
  if (like > *bestyet) {
    generic_tree_copy(t, bestree);
printf(" (new bestyet)"); /*  debug */
    *bestyet = like;
    there = p;
printf("\n"); /* debug */
    if (storing) {
      savetree(t, place);     /* storable coded representation of this tree */
      if (atstart) {                       /* when this is first tree tried */
        pos = 0;                     /* put it at the beginning of bestrees */
        found = false;
        if (nextree == 0) {
          *bestfound = like;                   /* score of the stored trees */
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (Initial, as Tree #1))\n", like, *bestyet, *bestfound);  /* debug */
          addbestever(pos, &nextree, maxtrees, false, place, bestrees, like);
printf("Added an initial tree to bestrees, now %ld of them\n", nextree); /*    debug */
        }
        *bestyet = like;          /* same value as *bestfound.  Why needed? */
        succeeded = true;         /* to be updated when "tryinsert" returns */
      } 
      else {
        if ( like == *bestfound )               /* deciding on a later tree */
        {              /* find where it goes in numerical order in bestrees */
          findtree(&found, &pos, nextree, place, bestrees);
          succeeded = true;
          if (!found) {                /* if found same tree, do not add it */
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (tied, as Tree %ld)\n", like, *bestyet, *bestfound, pos+1);  /* debug */
            addtiedtree(&pos, &nextree, maxtrees, false, place, bestrees, like);
printf("TREE %ld: ", pos+1);for (i = 0; i < spp; i++) printf("%ld ", place[i]);printf("\n");  /* debug */
printf("Added another tied tree to bestrees, now %ld of them\n", nextree); /*   debug */
          }
        } else {          /* since  like  is not the same as the best score */
          if (like > *bestfound) {                        /* replacing all? */
            *bestfound = like;
            *bestyet = like;
            pos = 0;                 /* put it at the beginning of bestrees */
            found = false;
printf(" score = %lf, bestyet = %lf, bestfound = %lf  (Better, as Tree #1)\n", like, *bestyet, *bestfound); /* debug: */
            addbestever(pos, &nextree, maxtrees, false, place, bestrees, like);
printf("TREE %ld: ", pos+1);for (i = 0; i < spp; i++) printf("%ld ", place[i]);printf("\n");  /* debug */
printf("Added new best tree to bestrees, score = %lf, now %ld of them\n", like, nextree);  /* debug */
            succeeded = true;
            *bestyet = like;
          }
        }
      }
    }
  }
  if (succeeded)
  {
    there = p;
/* debug:    *multf = false;   */
  }
  t->re_move(t, item, &dummy, true);   /* pull the branch back off the tree */
/* printf("then remove %ld from %ld:%ld\n", item->index, p->back->index, p->index); debug */
/* debug:  is preceding statement correct?  &dummy?  */
  t->restore_traverses(t, item, p);           /* debug: what is this doing? */
  t->evaluate(t, p, 0);   /* debug:   as in dnaml, but may not be needed */

  found = false;                /* debug: why this? May not have any effect */
  pos = 0;

  /* debug:  Uncommenting the following code will allow for a multifurcating
   * search, However doing a search only of resolved trees will also find
   * multifurcations when collapsible branches are collapsed */
#if 0
  if ( p->tip == false )
  {
    t->insert_(t, item, p, false);
    like = t->evaluate(t, p, false);
    t->score = like;
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


void parsimony_tree_init(tree* t, long nonodes, long spp)
{
  /* setup some functions of a new tree,  t,   with  spp  tips */
/* debug: are nonodes and spp arguments needed in this function? */

  t->globrearrange = pars_globrearrange;
  t->try_insert_ = (tree_try_insert_t)pars_tree_try_insert_;
  t->evaluate = pars_tree_evaluate;
} /* parsimony_tree_init */


void pars_node_init(node* p, node_type type, long index)
{
  /* set up an already-allocated new node for a parsimony_tree */
  pars_node *pn = (pars_node *)p;

  generic_node_init(p, type, index);
  p->copy = pars_node_copy;
  p->init = pars_node_init;
  p->reinit = pars_node_reinit;
  p->free = pars_node_free;
  p->node_print_f = pars_node_print;

  if (pn->numsteps)
    free(pn->numsteps);
  pn->numsteps = Malloc(endsite * sizeof(double));
} /* pars_node_init */


void pars_node_reinit(node * n)
{
  /* re-setup a pars_tree node */
  generic_node_reinit(n);
  pars_node *pn = (pars_node *)n;
  if (pn->numsteps)
    free(pn->numsteps);
  pn->numsteps = Malloc(endsite * sizeof(double));
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
   * and deleting trees that can be further collapsed. Continues this
   * until there are no further changes.   */
  long i, j, k, pos;
  double like;
  boolean found, collapsible, collapsed, somecollapsed;
  long treeLimit;
  node* p;

  treeLimit = nextree < maxtrees ? nextree : maxtrees;
  for (i = 0 ; i < treeLimit ; i++)           /* mark all trees collapsible */
  {
    bestrees[i].collapse = true;
  }
  if (progress)
  {
    sprintf(progbuf, "\nCollapsing best trees\n   ");
    print_progress(progbuf);
  }
/* debug:  comment out outer loop to see whether it is really necessary  */
    k = 0;
    do {           /* go to end of saved trees, trying to collapse branches */
      if (progress)
      {
        if ( (i % ((treeLimit / 72) + 1) ) == 0)  /* (% = mod) progress as  */
        { /* debug: fix formula for 72 */  /* ... row of dots up to 71 long */
          sprintf(progbuf, ".");
          print_progress(progbuf);
        }
      }
      while (!bestrees[k].collapse) /* go along bestrees until find one ... */
        k++;                               /* ... that might be collapsible */
      if (k >= treeLimit)          /* bail if all trees have been looked at */
        break;
      load_tree(t, k, bestrees);                       /* Reconstruct tree. */
      initializetrav(t, t->root);            /* make sure updates all views */
      initializetrav(t, t->root->back);
      like = t->evaluate(t, t->root, false);
      t->score = like;
      collapsible = false; 
      p = NULL;                /* for recording where tree can be collapsed */
  printf("(nextree before checking collapsibility is: %ld)\n", nextree);
  printf("STARTING treecollapsible on tree  %ld\n", k+1); /* debug */
  printf("LOOKING AT TREE #%ld: ", k+1);for (i = 0; i < spp; i++) printf("%ld ", place[i]);printf("\n");  /* debug */
      while ( treecollapsible(t, t->nodep[outgrno-1], &p, collapsible) ) {
printf("\nOuter call of treecollapsible on %ld:%ld\n", outgrno, t->nodep[outgrno-1]->back->index);  /* debug */
        collapsed = false;
        collapsetree(t, p->back, &collapsed);     /* collapse branch if can */
      }
      if (collapsed) {                        /* if something was collapsed */
        somecollapsed = somecollapsed || collapsible;
  printf("TREE #%ld collapsed\n", k+1);
        savetree(t, place);         /* record collapsed tree in place array */
  printf("COLLAPSED TREE: ");for (i = 0; i < spp; i++) printf("%ld ", place[i]);printf("\n");  /* debug */
        if ( k < (treeLimit-1) ) {       /* if not at the last tree already */
  printf("Shifting trees %ld through %ld down one\n", k+2, treeLimit);  /* debug */
          for (j = k ; j < (treeLimit - 1) ; j++)    /* shift rest of trees */
          {                   /* (in the process, overwriting the k-th tree */
            memcpy(bestrees[j].btree, bestrees[j+1].btree, spp * sizeof(long));
            bestrees[j].gloreange = bestrees[j + 1].gloreange;
            bestrees[j].locreange = bestrees[j + 1].locreange;
            bestrees[j + 1].gloreange = false;
            bestrees[j + 1].locreange = false;
            bestrees[j].collapse = bestrees[j + 1].collapse;
          }
        }
        treeLimit--;      /* because there is now one fewer tree in bestrees */
        nextree--;
printf("(nextree now %ld)\n", nextree); /* debug */
        pos = 0;
        findtree(&found, &pos, treeLimit, place, bestrees); /* find where   */
               /* ... the collapsed tree is to go, or whether already there */
        if (!found)  /* put the new tree in the the list if it wasn't found */
        {                       /* (note: treeLimit is increased as needed) */
          addtree(pos, &treeLimit, false, place, bestrees);
          nextree++;
          bestrees[pos].collapse = false;
  printf("ADDING NEW TREE: as number %ld: ", pos+1);
  for (i = 0; i < spp; i++) printf("%ld ", place[i]);printf("\n"); /* debug */
if (pos < treeLimit) {  /* debug */
printf(" and shifting trees %ld through %ld up one\n", pos+2, treeLimit); }/* debug*/
          if (pos <= k)        /* keep  k  pointing at next tree to examine */
            k++;
printf("(nextree now %ld)\n", nextree); /* debug */
printf("Next tree to examine is tree #%ld\n", k+1);  /* debug */
        }
    else printf("ALREADY THERE at %ld\n", pos); /* debug */
    }
  } while (k < treeLimit);
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
    * added to reconstruct this tree. This code assumes a root.
    * Trees get saved by storing that array in a 2D array.
    * This is the older function, a new function roots the tree
    * and calls this function, then stores "place" to save the tree.
    * The code imagines us adding tips to a tree and giving numbers
    * to the new interior forks, those numbers are the "lineage
    * numbers", which are not same as the current node index. */
  long i, newforknum, hitlineage, hitforknum;
  long* lineagenumber;
  node *p, *q, *rootnode;
  boolean topfork, donelineage;

  rootnode = t->nodep[outgrno - 1]->back; /* find fork attached to outgroup */
  while ( !(rootnode->back == NULL) ) {  /* which node in fork is at bottom */
    rootnode = rootnode->next;
  }                                         /* then call recursively to ... */
  setbottomtraverse(rootnode); /* set booleans indicating which way is down */
  lineagenumber = (long *)Malloc(nonodes*sizeof(long));
  /* find out where lineage down from each species connects */
  for (i = 0; i < nonodes; i++)      /* which lineage each tree node is ... */
    lineagenumber[i] = 0;                          /* ... starts out zeroed */
  lineagenumber[0] = 1;                        /* first lineage is number 1 */
  place[0] = 1;
  newforknum = spp;  /* one less than  number of next new fork to be put in */
  topfork = true;
  for (i = 1; i <= spp; i++)                            /* for each tip ... */
  {
    p = t->nodep[i - 1];                           /* start with species  i */
    lineagenumber[p->index - 1] = i;           /* set to number of that tip */
    if (p == NULL)                     /* if species  i  is not in the tree */
      break;                          /* ... then go on to the next species */
    if (p->back != NULL) {           /* if its back node is in the tree ... */
      p = p->back;                     /* ... go down to that interior node */
      while (lineagenumber[p->index - 1] == 0)    /* if no number yet there */
      {
        lineagenumber[p->index - 1] = i;       /* set to number of that tip */
/* printf("set species %ld lineagenumber to %ld\n", p->index, i); debug */
        while (!p->bottom)             /* go around circle to find way down */
          p = p->next;
        p = p->back;                             /* go down to earlier fork */
        if (p == NULL)                   /* if we went past bottom fork ... */
          break;                              /* ... bail out of while loop */
      }
      /* start a new lineage from where it connects, unless already a fork --
           in that case, set the lineage number negative for multifurcation */
      if (p != NULL) {           /* if we ran into a nonzero lineage number */
        hitlineage = lineagenumber[p->index - 1];         /* which lineage? */
        topfork = true;           /* will use to see if at top of a lineage */
        q = p;                     /* checking all descendants of this fork */
        do {     /* go around circle seeing if forks above are same lineage */
          if ( !(q->bottom) ) {                  /* node doesn't point down */
            if (q->back != NULL) {  /* node descendant not same lineage ... */
              topfork = topfork &&
                       ( !(lineagenumber[q->back->index - 1] == hitlineage) );
            }
          }
          q = q->next;
        } while (q != p);     /*  topfork  is true if none are same lineage */
        if (topfork) {        /* if this fork is the top one in its lineage */
          hitforknum = hitlineage;
          lineagenumber[p->index - 1] = -abs(hitforknum);  /* make negative */
          place[i-1] = -abs(hitforknum);
/* printf("set fork %ld lineagenumber to %ld\n", p->index, lineagenumber[p->index-1]); debug */
        } else {                          /* going on down that lineage ... */
          place[i-1] = hitlineage;                       /* set place value */
          newforknum++;
          do {  /* ... while still on that branch and no other new fork yet */
            lineagenumber[p->index - 1] = newforknum;    /* ... renumbering */
/* printf("set species %ld lineagenumber to %ld\n", p->index, newforknum); debug */
            while (!p->bottom)           /* go around circle to find way down */
              p = p->next;
            if (p->back == NULL)         /* blast out of loop if reached root */
              break;
            else
              p = p->back;                         /* go down to earlier fork */
            donelineage = topfork;
            if (!donelineage) 
              donelineage = (lineagenumber[p->index - 1] != hitlineage);
          } while (!donelineage);
        }
      }
    }
  }
/* printf("got place values for that tree: ");  debug */
/* for (i=0; i<spp; i++) printf(" %ld", place[i]); printf("\n");  debug */
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
  if (!wasrooted) {  /* debug:  change here?  */
    oldroot = t->nodep[outgrno-1];
    t->root = root_tree(t, oldroot);           /* put in a "fake" root fork */
  }
  moveroottooutgroup(t);
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


void addtiedtree(long* pos, long *nextree, long maxtrees, boolean collapse,
                  long *place, bestelm *bestrees, double score)
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
printf("Adding as new best tree (Tree 1), %ld\n", score); /*  debug */
    } else {                            /* it is another tree tied for best */
      findtree(&found, pos, nextree-1, place, bestrees);  /* already there? */
      if (!found) {                    /* save it only if not already there */
        addtiedtree(pos, &nextree, maxtrees, false, place, bestrees, score);
printf("Adding as tied tree (Tree %ld), %ld\n", *pos, score);  /* debug */
      } else {
/* printf("found that tree already there\n");  debug */
      }
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

/* debug:  does this make any sense?  Already saving best tree yet in generic version ...
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
  /* print the table of branch lengths */
  node *q, *r;
  long i;

  if (p->tip)
    return;
  r = p;                      /* save the node you start with for this fork */
  if (p->back == NULL) {   /* don't try to go out a branch to the null root */
    q = p->next;
  } else {       /* ... but do if the root branch leads to the outgroup ... */
    if (p == curtree->root)
      q = p;
    else                    /* ... but don't if you arrived via that branch */
      q = p->next;
  }
  do {                    /* go around fork circle, recursing out as needed */
    if (q->back != NULL) {  /* unless are at the bottom fork of rooted tree */
      fprintf(outfile, "%6ld      ", q->index - spp); /* print fork number, */
      if (q->back->tip)
      {
        for (i = 0; i < nmlngth; i++)         /* ... then name if a tip ... */
          putc(nayme[q->back->index - 1][i], outfile);
      }
      else
        fprintf(outfile, "%6ld    ", q->back->index - spp);  /* else number */
      fprintf(outfile, "   %f\n", q->v);
      printbranchlengths(q->back);       /* on our way out through the tree */
    }
    q = q->next;
  } while (q != r);
} /* printbranchlengths */


void initbranchlen(node *p)
{
  /* initial values of branch lengths */
  node *q;

  p->v = 0.0;
  if (p->back != NULL)             /* set length at start of branch to zero */
    p->back->v = 0.0;
  if (p->tip)
    return;
  q = p->next;
  while (q != p)
  {
    if (q->back != NULL)                        /* recurse out through tree */
      initbranchlen(q->back);
    q = q->next;
  }
  q = p->next;
  while (q != p)
  {
    q->v = 0.0;      /* ... then set branch length at end of branch to zero */
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


void moveroottooutgroup(tree* t) {
  /* if root is on a multifurcating node, remove it from that fork circle,
     toss the circle node, make a new fork and use it to put the root
     on the outgroup lineage. */
  boolean foundit;
  node *forknode, *preforknode, *p, *q;

  if (count_sibs(t->root) > 2) { /* if the root is on a multifurcating fork */
    foundit = false;
    p = t->nodep[outgrno - 1]->back;   /* for the fork nearest the outgroup */
    q = p;
    do {
      if (q->back == NULL) {
        forknode = q;           /* find a node that has nothing below it */
        foundit = true;
      }
      q = q->next;
    } while (q != p); 
    if (foundit) {    /* remove the interior node which has an empty neighbor */
        preforknode = precursor(forknode);     /* find fork node preceding it */
        preforknode->next = forknode->next;            /* and connect past it */
        t->nodep[q->index - 1] = preforknode; /* and have nodep point to that */
        t->release_forknode(t, forknode);                      /* and toss it */
        t->root = root_tree(t, t->nodep[outgrno - 1]);    /* put root fork in */
    }
  }
} /* moveroottooutgroup */


void load_tree(tree* t, long treei, bestelm* bestrees)
{
  /* restores tree  treei  from array bestrees (treei is the index
   * of the array, so tree 5 has treei = 4).  Add all the tips to a tree one
   * by one in order.  The array element  bestree[treei].btree[k], if
   * positive, indicates that tip  k+1  is to be connected to a new fork
   * that is just below tip or fork  btree[k].  The numbers here are not the
   * values of "index" but tips and forks are numbered in order of their
   * consideration, so tips from 1 .. spp and forks are numbered from  spp+1
   * on as they are added to the tree.  If negative, btree[k] indicates that
   * species  k+1  is to be added as an extra furc to the fork has number
   * -btree[k], */
  long i, j, numofnewfork, belowindex;
  long  nsibs;  /* debug */
  boolean foundit = false; /* debug */
  node *p, *q, *below, *bback, *forknode, *newtip, *bbot, *afterwhere;

  release_all_forks(t);              /* to make sure all interior nodes
                                        are on the list at free_fork_nodes  */
                                     /* then make tree of first two species */
  forknode = t->get_fork(t, spp);     /* fork put on nodep, index is  spp+1 */
  hookup(t->nodep[1], forknode->next);
  hookup(t->nodep[0], forknode->next->next);
  forknode->back = NULL;
  numofnewfork = spp;                           /* index-1 of next new fork */
  t->root = t->nodep[0]->back;

  for ( j = 3; j <= spp ; j++ )     /* adding one by one species, 3, 4, ... */
  {
    newtip = t->nodep[j-1];
    belowindex = bestrees[treei].btree[j - 1];             /* where it goes */
    if ( belowindex > 0 ) {          /*  if the fork is to be a bifurcation */
      numofnewfork++;
      forknode = t->get_fork(t, numofnewfork);  /* put a new fork circle in */
      below = t->nodep[belowindex - 1]; 
      bback = below->back;
      hookup(forknode->next, below);         /* hook to it node it is below */
      hookup(newtip, forknode->next->next);               /* hook tip to it */
      if ( bback != NULL )      /* if below the new fork is not a NULL node */
        hookup(forknode, bback);
      t->nodep[numofnewfork] = forknode;        /* points to bottom of fork */
    }
    else
    {          /*  if goes into a multifurcation put a new node into circle */
      bbot = t->nodep[-belowindex-1];      /* the fork circle's bottom node */
      afterwhere = precursor(bbot);     /* find fork node that points to it */
      forknode = t->get_forknode(t, bbot->index);         /* get a new node */
      hookup(newtip, forknode);             /* hook the tip to the new node */
      afterwhere->next = forknode;                /* put it the right place */
      forknode->next = bbot;              /* namely, the last in the circle */
    }
  }
  /* if root is at multifurcation, move it to be next to outgroup instead */
  moveroottooutgroup(t);
} /* load_tree */


void setbottomtraverse(node *p)
{ 
  /* set boolean "bottom" to true on one of the nodes in a fork circle
   * at each interior fork to show which way is down.  Go around
   * fork circle and set others to false.  Traverse out through tree
   * to do this on all nodes.  */
  node *q;

  if (p == NULL)                  /* return if this node is not on the tree */
    return;
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


void pars_globrearrange(tree* curtree, tree* bestree, boolean progress,
                         boolean thorough, double* bestfound)
{ /* does global rearrangements.  The more general
   * generic_unrooted_locrearrange also works but this is much faster because
   * it gets to take advantage of some of the speedups available in the
   * parsimony programs */
  int i;
  node *p, *where, *there, *qwhere;
  double bestyet;
  boolean success, successaftertraverse, dontremove, donttrythere;
  node* removed;

/*  bestyet = curtree->evaluate(bestree, bestree->root, 0); debug */
  bestyet = *bestfound;

  success = false;        /* to be  true  when a better tree has been found */
  do {           /* loop until a global rearrangement pass does not improve */
    if (progress)
    {
      sprintf(progbuf, "   ");
      print_progress(progbuf);
    }
    for ( i = 0 ; i < curtree->nonodes ; i++ )
    {
      if ( (curtree->root->index == (i+1)) ||
           (curtree->root->back->index == (i+1)) )
        continue;    /* skip this case if this branch connects to root fork */
      p = curtree->nodep[i]->back;
      if (p == NULL)
        continue;        /* skip this case if no interior node circle there */

      if ( progress &&
             (((i-spp)+1) % (( (curtree->nonodes - 1) / 72 ) + 1 ) == 0 ))
      {                     /* print progress characters, up to 72 per line */
        sprintf(progbuf, ".");
        print_progress(progbuf);
      }

      if ( p->back == NULL )                               /* skip to next */
        continue;                      /* ... if nothing connected here ... */
      dontremove = (p->index == curtree->root->index) ||
                    (p->back->index == curtree->root->index);
      if (!dontremove) {
        removed = p; 
/*  printf(" remove %ld:%ld\n", removed->index, removed->back->index); debug */
        curtree->re_move(curtree, removed, &where, true);
        qwhere = where;                  /* to hold best place to insert it */
        success = pars_tree_try_insert_(curtree, removed, where, there,
                        &bestyet, bestree, true, true, true, bestfound);
        donttrythere = (where->tip ||
                            where->back == curtree->root->index);
        if (!donttrythere) {
          successaftertraverse = generic_tree_addtraverse_1way(curtree,
                                   removed, where, true, qwhere, &bestyet,
                                   bestree, true, true, false, bestfound);
          success = success || successaftertraverse;
        }
        donttrythere = (where->back->index == curtree->root->index) ||
                         (where->back->tip);     /* a tip or rootmost fork? */
        if ( !donttrythere ) {    /* if not, do traversal also at other end */
          successaftertraverse = generic_tree_addtraverse_1way(curtree,
                                   removed, where->back, true, qwhere,
                                   &bestyet, bestree, true, true, false,
                                   bestfound);
        }
        success = success || successaftertraverse;
/* printf("inserting at %ld:%ld\n", qwhere->index, qwhere->back->index); debug */
        curtree->insert_(curtree, removed, qwhere, false); /* in best place */
        curtree->root = curtree->nodep[outgrno-1]->back;        /* set root */
/* printf("setting root as: %ld\n", curtree->root->index); debug */
/* debug: why?        bestyet = curtree->evaluate(curtree, curtree->root, 0);   debug */
      }
    }
  } while (success);
  if (progress)
  {
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
} /* pars_globrearrange */


boolean treecollapsible(tree* t, node* n, node** p, boolean collapsible)
{
 /* find out whether there is any collapsible branch on the tree.
  * In initial call of the recursion,  n  should be a node that is
  * not a tip.  p is initially NULL but is set to n if that branch
  * can be collapsed */

  node *sib;

printf("\ncalled treecollapsible with %ld:%ld\n", n->index, n->back->index); /*  debug */
  if ( n == NULL )                /* in case it is called on branch at root */
    return false;
  if ( n->back == NULL )                /* in case it is called on branch at root */
    return false;

printf("calling branchcollapsible with branch %ld-%ld\n", n->index, n->back->index); /*  debug */
  if ( ((pars_tree*)t)->branchcollapsible(t, n) ) {    /* check this branch */
printf(" (collapsible) \n"); /*  debug */
    *p = n;                     /* record the node where it can be collapsed */
    return true;             /* then bail out and do not recurse further in */
  }
  else   /* debug */
/*     printf(" (not collapsible) \n");   debug */

  if ( n->back->tip == true )         /* in case we've reached a tip branch */
    return false;
printf("going around circle for fork %ld\n", n->back->index); /* debug */
  for ( sib = n->back->next ; sib != n->back ; sib = sib->next )
  {                                                  /* recurse further out */
    if (sib->back != NULL) {
printf("collapsible was %ld, now do recursive call on %ld-%ld\n", (long)collapsible, sib->index, sib->back->index);  /* debug */
    collapsible = treecollapsible(t, sib, p, collapsible) || collapsible;
    }
   else
      printf("root branch skipped\n"); /* debug */
  }
printf("collapsible: %d\n", collapsible);  /* debug */
  return collapsible;
} /* treecollapsible */


void collapsebranch(tree* t, node* n)
{ /* remove a branch and merge the forks at both ends
   * node  n  must have its back pointer point to a fork circle */
  node *m, *prem, *pren, *p, *q; 
  long i, j, k;
  node *sib, *newfork;  /* debug: can delete this once code block deleted */
  long nsibs;   /* debug: this too */

  inittrav(t, n);    /* make initialized pointers looking in to fork  false */
  if (n == NULL)      /* make sure both ends of branch exist and are forks */
    return;
  if (n->tip)
    return;
  m = n->back;              /* get other end of branch too, do same checks */
  inittrav(t, m->back);        /* (recurse out from fork circle doing that) */
    /* debug:  the preceding may be unnecessary, they may already be OK */
/* printf("COLLAPSING branch %ld:%ld\n",n->index,m->index); debug */
  if (m == NULL)
    return;
  if (m->tip)
    return;
  prem = precursor(m);            /* node in fork circle that precedes  m  */
  pren = precursor(n);                            /* ... and ditto for  n  */
  i = n->index;                          /* get indexes of the two circles */
  j = m->index;
  prem->next = n->next;                      /* merge the two fork circles */
  pren->next = m->next;
  if (i > j) {                 /* find the smaller of the fork indexes ... */
    k = j;
    j = i;
    i = k;                                           /* ... it will be  i  */
  }
  if ( (t->nodep[i-1] == n) )
    t->nodep[i-1] = n->next;
  if ( (t->nodep[i-1] == m) )
    t->nodep[i-1] = m->next;
  t->nodep[j-1] = NULL;     /* debug: necessary? Done by release_forknode? */
  m->next = NULL;         /* set these disconnected nodes to point nowhere */
  n->next = NULL;
  p = t->nodep[i-1];                               /* start of merged fork */
  p->index = i; 
  for (q = p->next; q != p; q = q->next) { /* renumber nodes in one of ... */
    q->index = i;                         /* ... the original fork circles */
  }
  t->release_forknode(t, m); /* now recycle  m, n  as are no longer needed */
  t->release_forknode(t, n);
  t->score = t->evaluate(t, t->nodep[outgrno-1], false); 
/* debug:  need to replace all this code ... */
#if 0
  for ( sib = m->next ; sib != m ; sib = sib->next )    /* go around circle */
  {
    if ( sib == m->next )   /* debug:   huh?  check! */
      newfork = n;
    else
    {
      newfork = t->get_forknode(t, n->index);            /* get a new node */
      newfork->next = n->next;                /* put it in the fork circle */
      n->next = newfork;
    }
    hookup(sib->back, newfork);       /* hook stuff in back to the new node */
  }
  t->release_fork(t, m);               /* toss  m  back onto free node list */
#endif
} /* collapsebranch */


void collapsetree(tree* t, node* n, boolean* collapsed)
{  /* collapse all branches that are designated as collapsible
    * node  n  should be an outgroup tip or a fork circle node from
    * which we are proceeding out the back subtree of that node */
  node *sib;

  if (n == NULL)
    return;
  if (n->back == NULL)
    return;
  if ( n->back->tip == true)           /* if back is a tip don't go further */
    return;
  if ( ((pars_tree*)t)->branchcollapsible(t, n) ) {
    collapsebranch(t, n);                 /* collapse this branch if we can */
    *collapsed = true;
    return;
  }
  else /* go around circle, for all but initial node, collapse back subtree */
    for ( sib = n->back->next ; sib != n->back ; sib = sib->next ) {
      collapsetree(t, sib, collapsed);                    /* try further on */
    }
} /* collapsetree */


void printree(tree* t)
{ /* prints out diagram of the tree */
  long tipy;
  double scale, tipmax;
  long i;
  boolean *found;
  node *p, *q;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;            /* this will be the number for each line printed out */
  tipmax = 0.0;     /* will keep track of how far to right the tree extends */
  p = findroot(t, t->root, found);      /* find node with null back pointer */
  q = p;
  coordinates(t, p, 0.0, &tipy, &tipmax);       /* get coordinates of nodes */
  scale = 1.0 / (long)(tipmax + 1.000);      /* horizontal scaling for tree */
  for (i = 1; i <= (tipy - down); i++)       /* draw one line at a time ... */
    drawline3(i, scale, p);      /* ... each time going up from root to tip */
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
  if (p->tip)                                            /* if at a tip ... */
  {
    p->xcoord = (long)(over * lengthsum + 0.5);      /* set coordinates ... */
    p->ycoord = (*tipy);
    p->ymin = (*tipy);        /* ... as well as min, max  y  coordinate ... */
    p->ymax = (*tipy);    /* ... which for a tip is just its  y cooprdinate */
    (*tipy) += down;      /* increment the  tipy  value, ready for next tip */
    if (lengthsum > (*tipmax))  /* update how far out to right tree extends */
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;      /* for interior fork, go around the circle of nodes ... */
  do {
    xx = q->v;        /* get the value of the branch length leading out ... */
    if (xx > 100.0)                /* ... from that node in the fork circle */
      xx = 100.0;      /*   debug:  why this?  I think it works */
    coordinates(t, q->back, lengthsum + xx, tipy, tipmax);   /* recurse out */
    q = q->next;                /* continue around the fork circle of nodes */
  } while (p != q);
  first = p->next->back;    /* around circle again, getting pointers to ... */
  q = p;       /* ... the descendants of first and last nodes in the circle */
  while (q->next != p)           /* ... not counting the one you arrived to */
    q = q->next;
  last = q->back;   /* here's the pointer to the descendant of the last one */
  p->xcoord = (long)(over * lengthsum + 0.5);
  if ((p == t->root) || count_sibs(p) > 2)
    p->ycoord = p->next->next->back->ycoord;     /* where 2nd descendant is */
  else                              /* (that works OK for 3, 4 descendants) */
    p->ycoord = (first->ycoord + last->ycoord) / 2;  /* mean, if bifurcates */
  p->ymin = first->ymin;  /* on way back down tree propagate ymin, ymax ... */
  p->ymax = last->ymax;            /* ... to come from first, last subtrees */
}  /* coordinates */


void drawline3(long i, double scale, node *start)
{
  /* draws one row of the tree diagram by moving up tree
   * used in pars and dnapars */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first =NULL, *last =NULL;
  boolean done;

  p = start;
  q = start;
  extra = false;
  if (i == (long)p->ycoord)                /* print number of interior fork */
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
      first = p->next->back;         /* find first immediate descendant ,,, */
      r = p;
      while (r->next != p) {
        r = r->next;
        if (r->back != NULL)                            /* ... and last one */
          last = r->back;
      }
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);    /* how many chars */
    if (n < 3 && !q->tip)      /* branch must be at least 3 characters long */
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done)  /* are we where a corner is needed? */
    {
      if ((long)p->ycoord != (long)q->ycoord)     /* if at corner of branch */
        putc('+', outfile);
      else                                  /* otherwise print another dash */
        putc('-', outfile);
      if (!q->tip)                         /* if we're not yet at a tip ... */
      {
        for (j = 1; j <= n - 2; j++)  /* print enough dashes to get to fork */
          putc('-', outfile);
        if (q->index - spp >= 10)                   /* print number at fork */
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      }
      else
      {      /* if it's to lead to a tip, print enough dashes but no number */
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    }
    else if (!p->tip)
    {                     /* if a branch crosses here, print a vertical bar */
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          (i != (long)p->ycoord))
      {
        putc('|', outfile);      /* print the vertical bar where it crosses */
        for (j = 1; j < n; j++)             /* ... then print enough spaces */
          putc(' ', outfile);
      }
      else
      {                                  /* ... otherwise just print spaces */
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    }
    else
    {             /*  debug: (correct?)  if  p  is a tip, print spaces */
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip)      /* now we print the species name */
  {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index-1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline3 */


void writesteps(tree* t, long chars, boolean weights, steptr oldweight)
{
  /* used in Dnacomp, Dnapars, Dnapenny, Pars, and Penny */
  long i, j, k, l;
  boolean *found;
  node *p;
  k=0;

  /*calculate the steps */
  p = findroot(t, t->root, found);
  if (p->initialized == false ) t->nuview(t, p);

  /* print them */
  putc('\n', outfile);
  if (weights)
    fprintf(outfile, "weighted ");
  fprintf(outfile, "steps in each site:\n");
  fprintf(outfile, "      ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%4ld", i);
  fprintf(outfile, "\n     ,-------------------------------------");
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
        l = location[ally[k-1]-1];
        if (oldweight[k] > 0)
          fprintf(outfile, "%4ld",
                   oldweight[k]
                    * (((pars_node*)p)->numsteps[l-1])/weight[l-1]);
        else
          fprintf(outfile, "%4ld", (((pars_node*)p)->numsteps[l - 1] ));
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
  /* write out file with representation of final tree
   * used in pars and dnapars -- writes branch lengths
   * This version also indents additional lines of trees appropriately */
  node *q;
  long i, n, w;
  double x;
  boolean atstart = false;
  Char c;

  if (p == NULL)                                /* bail out if no node here */
    return;
  if (p == root) {
    indent = 0;
    atstart = true;
  }
  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)
    {                                           /* how long is species name */
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)                       /* write out species name */
    {
      c = nayme[p->index - 1][i];
      if (c == ' ')                    /* with underscores replacing blanks */
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  }
  else
  {
    putc('(', outtree);                   /* left paren if encounter a fork */
    (*col)++;
    indent++;                            /* increment amount of line indent */
    q = p;
    do                                           /* loop over nodes in fork */
    {
      if ((q != p) || atstart) {        /* ... unless we are done with fork */
        treeout3(q->back, nextree, col, indent, root);       /* recurse out */
        atstart = false;
        if ((q->back != NULL) && (q->next != p))
          putc(',', outtree);         /* ... printing a comma for next furc */
        (*col)++;
        if (*col > 60)  /* move on to new line if get too far right on line */
        {
          putc('\n', outtree);
          *col = 0;
          for (i = 1; i <= indent; i++) /* write indent at beginning of line */
            putc(' ', outtree);
        }
      }
      q = q->next;
    } while (q != p);
    putc(')', outtree);        /* finish with fork by printing right-paren */
    (*col)++;
  }
  x = p->v;                                      /* printing branch length */
  if (x > 0.0)       /* hard part is figuring out how wide number is to be */
    w = (long)(0.43429448222 * log(x));
  else { if (x == 0.0)
      w = 0;
    else
      w = (long)(0.43429448222 * log(-x)) + 1;
  }
  if (w < 0)                                           /* gets bizarre ... */
    w = 0;
  if (p != root)
  {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    *col += w + 8;
  }
  if (p != root)                                       /* end of a recurse */
    return;
  if (nextree > 2)              /* weights of trees if there are tied ones */
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
