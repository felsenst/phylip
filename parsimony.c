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

static void savetraverse(node *p);
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
  node *nuroot1, *nuroot2, *nuroot3, *there;  /* the three nodes in the new circle
                                               * and a saving of  here->back */
  there = here->back;
  nuroot1 = t->get_forknode(t, nonodes);
  hookup(here, nuroot1);
  nuroot2 = t->get_forknode(t, nonodes);
  nuroot1->next = nuroot2;
  hookup(nuroot2, there);
  nuroot3 = t->get_forknode(t, nonodes);
  nuroot2->next = nuroot3;
  nuroot3->next = nuroot1;
  nuroot3->back = NULL;
  t->nodep[nonodes-1] = nuroot3;
  return nuroot3;
} /* root_tree */


void reroot_tree(tree* t, node* fakeroot) // RSGbugfix: Name change.
{
  /* Removes a root from a tree; useful after a return from functions that
   * expect a rooted tree (e.g. oldsavetree()). Then reroots the tree before
   * releasing a FORKNODE or FORKRING to a FREELIST, to avoid tree components
   * pointing (even temporarily) into garbage. */
  node *p;

  if ( count_sibs(fakeroot) > 2 )
  {
    for (p = fakeroot ; p->next != fakeroot ; p = p->next);
    p->next = fakeroot->next;
    if ( t->nodep[fakeroot->index - 1 ] == fakeroot)
      t->nodep[fakeroot->index - 1 ] = p;
    if ( t->root == fakeroot)   // RSGbugfix: Reroot before sending FORKNODE to freelist.
      t->root = t->nodep[outgrno - 1]->back;
    t->release_forknode(t, fakeroot);
  }
  else
  {
    hookup(fakeroot->next->back, fakeroot->next->next->back);
    if ( t->root == fakeroot)           // RSGbugfix: Reroot before sending FORKRING to freelist.
      t->root = t->nodep[outgrno - 1]->back;
    t->release_fork(t, fakeroot);
  }
} /* reroot_tree */


boolean pars_tree_try_insert_(tree * t, node * item, node * p, node * there,
  double * bestyet, tree * bestree, tree * priortree, boolean thorough,
  boolean *multf )
{
  /* insert item at p, if the resulting tree has a better score, update bestyet
   * and there
   * This version actually does the hookups which are quickly dissolved,
   * however none of the changes are propegated in the tree and it is like as
   * if it never got inserted. If we are on the last rearrangement save a
   * bestscoring insert to the bestrees array
   * item  should be an interior fork hooked to a tip or subtree which
   * is item->back */
  double like;
  boolean succeeded = false;
  node* dummy;
  boolean found = false;
  long pos = 0;

/* debug:    t->save_traverses(t, item, p);      */
  t->insert_(t, item, p, false);
  like = t->evaluate(t, p, false);

  if (like > *bestyet || *bestyet == UNDEFINED)
  {
    there = p;
    succeeded = true;
    *multf = false;
  }

  if ( lastrearr && (like >= *bestyet || *bestyet == UNDEFINED))
  {
    savetree(t, place);
    findtree(&found, &pos, nextree, place, bestrees);
    if ( !found )
    {
      if (*bestyet < like || nextree == 1 )
        addbestever(&pos, &nextree, maxtrees, false, place, bestrees, like);
      else
        addtiedtree(pos, &nextree, maxtrees, false, place, bestrees, like);
    }
  }
  if (like > *bestyet || *bestyet == UNDEFINED)
    *bestyet = like;
  t->re_move(t, item, &dummy, true);
/* debug:   t->restore_traverses(t, item, p);   */
  t->evaluate(t, p, 0);   // as in dnaml, but may not be needed

  found = false;
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
        if ( (like > *bestyet  && !found) || nextree == 1)
          addbestever(&pos, &nextree, maxtrees, false, place, bestrees);
        else if ( !found )
          addtiedtree(pos, &nextree, maxtrees, false, place, bestrees);
      }
      *bestyet = like;
    }
    t->re_move(t, item, &dummy, true);
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
                       boolean progress, long * finalTotal)
{
  /* Goes through all best trees, collapsing trees where possible,
   * and deleting trees that are not unique.    */
  long i, j, k, pos ;
  boolean found;
  long treeLimit = nextree - 1 < maxtrees ? nextree - 1 : maxtrees;

  (void)chars;                          // RSGnote: Parameter never used.

  for(i = 0 ; i < treeLimit ; i++)
  {
    bestrees[i].collapse = true;
  }
  if(progress)
  {
    sprintf(progbuf, "Collapsing best trees\n   ");
    print_progress(progbuf);
  }
  k = 0;
  for(i = 0 ; i < treeLimit ; i++)
  {
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
    /* Reconstruct tree. */
    load_tree(t, k, bestrees);
    while ( treecollapsible(t, t->nodep[0]))
      collapsetree(t, t->nodep[0]);
    savetree(t, place);
    /* move everything down in the bestree list */
    for(j = k ; j < (treeLimit - 1) ; j++)
    {
      memcpy(bestrees[j].btree, bestrees[j + 1].btree, spp * sizeof(long));
      bestrees[j].gloreange = bestrees[j + 1].gloreange;
      bestrees[j + 1].gloreange = false;
      bestrees[j].locreange = bestrees[j + 1].locreange;
      bestrees[j + 1].locreange = false;
      bestrees[j].collapse = bestrees[j + 1].collapse;
    }
    treeLimit--;

    pos=0;
    findtree(&found, &pos, treeLimit, place, bestrees);

    /* put the new tree in the the list if it wasn't found */
    if(!found)
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
  long i, j, nextnode, nvisited;
  node *p, *q, *r = NULL, *root2, *lastdesc, *outgrnode, *binroot, *flipback;
  boolean done, newfork;

  binroot = NULL;
  lastdesc = NULL;
  root2 = NULL;
  flipback = NULL;
  outgrnode = t->nodep[outgrno - 1];
#if 0                        /* debug: don't need this?  */
  if (get_numdesc(t->root, t->root) == 2)
    bintomulti(t, &t->root, &binroot);
#endif
  if (outgrin(t->root, outgrnode))
  {
    if (outgrnode != t->root->next->back)
      moveleft(t->root, outgrnode, &flipback);
  }
  else
  {
    root2 = t->root;
    lastdesc = t->root->next;
    while (lastdesc->next != t->root)
      lastdesc = lastdesc->next;
    lastdesc->next = t->root->next;
    t->root = t->get_forknode(t, outgrnode->back->index);
    reroot2(outgrnode, t->root);
  }
  savetraverse(t->root);
  nextnode = spp + 1;
  for (i = nextnode; i <= nonodes; i++)
    if (get_numdesc(t->root, t->nodep[i - 1]) == 0)
      flipindexes(i, t->nodep);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[t->root->index - 1] = 1;
  for (i = 1; i <= spp; i++)
  {
    p = t->nodep[i - 1];
    while (place[p->index - 1] == 0)
    {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      r = p;
      p = p->back;
    }
    if (i > 1)
    {
      q = t->nodep[i - 1];
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
          while (!p->bottom)
            p = p->next;
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
    }
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
}  /* oldsavetree */


void savetree(tree* t, long *place)
{
  /* Record in  place  where each species has to be added to reconstruct
   * this tree. This code roots the tree and calls oldsavetree to save it. */
  node *oldroot, *p, *outgrnode;

  outgrnode = t->nodep[outgrno - 1];
  p = outgrnode->back;
  oldroot = p;
  t->root = root_tree(t, p);
  oldsavetree(t, place);
  reroot_tree(t, t->root);              // RSGbugfix: Name change.
  t->root = oldroot;
}  /* savetree */


void addbestever(long *pos, long *nextree, long maxtrees, boolean collapse,
                  long *place, bestelm *bestrees, double score)
{
  /* adds first best tree. If we are rearranging on usertrees, 
   * add it to the second array of trees if the score is good enough */
  long repos;
  boolean found;

  *pos = 1;
  *nextree = 1;

  addtree(*pos, nextree, collapse, place, bestrees);
  if ( reusertree )
  {
    if ( score == UNDEFINED ) return;
    if ( score != UNDEFINED && (score > rebestyet || rebestyet == UNDEFINED))
    {
      renextree = 1;
      rebestyet = score;
      addtree(1, &renextree, collapse, place, rebestrees[1]);
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


void addtiedtree(long pos, long *nextree, long maxtrees, boolean collapse, long *place, bestelm *bestrees, double score)
{
  /* add a tied tree */
  boolean found;
  long repos;

  if (*nextree <= maxtrees)
    addtree(pos, nextree, collapse, place, bestrees);
  if ( reusertree )
  {
    if ( rebestyet == score )
    {
      findtree(&found, &repos, renextree, place, rebestrees[1]);
      if ( !found && renextree <= maxtrees )
        addtree(repos, &renextree, collapse, place, rebestrees[1]);
    }
  }
} /* addtiedtree */

/* debug:  may not need in view of pars_try_insert  */
#if 0
void add_to_besttrees(tree* t, long score, bestelm* bestrees)
{
  /* take the tree we have found and try to add it to the array bestrees:
   * if none are already there, make it the first one, if it is better than
   * the ones that are there then toss them and start over with just this
   * one, if tied with them add it in too */

  boolean better, worse, tied;
  long bestscoreyet;
  
  bestscoreyet = bestyet;
  if (!(score > bestscoreyet)) {           /* if going to save this one ... */
    savetree(t, place);
    if (score < bestscoreyet) {      /* if it will be the lone new best one */
      addbestever(&pos, nextree, maxtrees, collapse, place, bestrees, score);
    } else {                            /* it is another tree tied for best */
      addtiedtree(&pos, *nextree, maxtrees, collapse, place, bestrees, score);
    }
  }
    if ( !found )
    {
      if (*bestyet < like || nextree == 1 )
        addbestever(&pos, &nextree, maxtrees, false, place, bestrees, like);
      else
        addtiedtree(&pos, &nextree, maxtrees, false, place, bestrees, like);
    }
} /* add_to_besttrees */
#endif     /* debug: end of commented-out stuff */


boolean pars_addtraverse(tree* t, node* p, node* q, boolean contin,
                          node* qwherein, double* bestyet, bestelm* bestrees,
                          boolean thorough)
{
  /* wrapper for addraverse, calling generic addtraverse
   * and then taking the best tree found and adding it to the array
   * of tied best trees found. Function like this works for parsimony-like
   * criteria where there are exact ties, not for likelihood or distance
   * criteria */
   boolean success;

   success = generic_tree_addtraverse(t, p, q, contin, qwherein,
                                       bestyet, &bestree, thorough);
/* debug:    add_to_besttrees(tree* t, long score, bestelm* bestrees);  */
   return success;
} /* pars_addtraverse */


static void flipnodes(node *nodea, node *nodeb)
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


static void moveleft(node *root, node *outgrnode, node **flipback)
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
  /* flips index of nodes between nextnode and last node.  */
  long last;
  node *temp;

  last = nonodes;
  while (treenode[last - 1]->back == NULL)
    last--;
  if (last > nextnode)
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
  /* computes the number of nodes which are visited earlier than anode among
     its siblings */
  node *p;
  long nvisited;

  while (!anode->bottom) anode = anode->next;
  p = anode->back->next;
  nvisited = 0;
  do {
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

  while (!anode->bottom) anode = anode->next;
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
  /* assigns index i to node p */

  while (p->index != i)
  {
    p->index = i;
    p = p->next;
  }
} /* newindex */


void load_tree(tree* t, long treei, bestelm* bestrees)
{
  /* restores a tree from bestrees */
  long i, j, nsibs;
  boolean foundit = false;
  node *p, *q, *below, *bback, *forknode, *newtip;

  destruct_tree(t);       /* to make sure all interior nodes are on 
                           * the list at free_fork_nodes  */
  /* restore the tree */
  forknode = t->get_fork(t, spp);      /* was put on nodep, index is  spp+1 */
  hookup(t->nodep[1], forknode->next);
  hookup(t->nodep[0], forknode->next->next);

  for ( j = 3; j <= spp ; j++ )     /* adding one by one species, 3, 4, ... */
  {
    newtip = t->nodep[j-1];

    if ( bestrees[treei].btree[j-1] > 0 )    /* j-th entry in "place" array */
    {              /*  if bifurcation */
      forknode = t->get_fork(t, spp+j-2);       /* put a new fork circle in */
      hookup(newtip, forknode);
      below = t->nodep[bestrees[treei].btree[j - 1] - 1];
      bback = below->back;
      hookup(forknode->next, below);
      if ( bback )
        hookup(forknode->next->next, bback);
    }
    else
    {          /*  if goes into a multifurcation put a new node into circle */
      below = t->nodep[t->nodep[-bestrees[treei].btree[j-1]-1]->back->index-1];
      forknode = t->get_forknode(t, below->index);
      forknode->next = below->next;
      below->next = forknode;
    }
  }

  forknode = NULL;
  for (i = spp; i <= nonodes; i++) {     /* check all interior node circles */
    p = t->nodep[i];
    for (q = p; q != p; q = q->next) {                 /* go all way around */
      if (q->back == NULL) {
        forknode = q;              /* find a node that has nothing below it */
        foundit = true;
      }
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
  t->score = bestyet;
} /* load_tree */


static void savetraverse(node *p)
{ 
  /* set boolean "bottom" on each interior node to show which way is down */
  node *q;

  p->bottom = true;
  if (p->tip)
    return;
  q = p->next;
  while (q != p)
  {
    q->bottom = false;
    savetraverse(q->back);
    q = q->next;
  }

}  /* savetraverse */


static boolean outgrin(node *root, node *outgrnode)
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


void pars_globrearrange(tree* curtree, boolean progress, boolean thorough)
{ /* does global rearrangements */
  /* The more general generic_unrooted_locrearrange also works but this is
   * much faster because it gets to take advantage of some of the speedups
   * available in the parsimony programs */
  int i, j, k, num_sibs, num_sibs2;
  node *where, *sib_ptr, *sib_ptr2, *qwhere;
  double oldbestyet;
  double bestyet;
  boolean multf, mulf;
  node* removed;

  bestyet = oldbestyet = curtree->evaluate(curtree, curtree->root, 0);

  if (progress)
  {
    sprintf(progbuf, "   ");
    print_progress(progbuf);
  }

  for ( i = 0 ; i < curtree->nonodes ; i++ )
  {
    sib_ptr  = curtree->nodep[i];
    if (sib_ptr == NULL)
      continue;         /* skip this case if no interior node circle there */
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

        removed = sib_ptr->back;
        mulf = 2 != count_sibs(removed->back);
        curtree->re_move(curtree, removed, &where, true);
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
                                qwhere, &bestyet, &bestree, multf);
          sib_ptr2 = sib_ptr2->next;
        }
        curtree->insert_(curtree, removed, where, mulf);
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


void grandrearr(tree* t, boolean progress, boolean rearrfirst)
{
  /* calls global rearrangement on best trees */
  long treei;
  long i, pos;
  boolean done = false;

  lastrearr = true;
  savetree(t, place);
  addbestever(&pos, &nextree, maxtrees, false, place, bestrees, UNDEFINED);

  for ( i = 0 ; i < nextree ; i++)
    bestrees[i].gloreange = false;

  while (!done)
  {
    treei = findunrearranged(bestrees, nextree, true);
    if (treei < 0)
      done = true;
    else
      bestrees[treei].gloreange = true;

    if (!done)
    {
      load_tree(t, treei, bestrees);
      t->evaluate(t, t->root, 0);
      t->globrearrange(t, progress, true);
      done = rearrfirst;
    }
  }
} /* grandrearr */


void treeout3(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dnapars -- writes branch lengths */
  node *q;
  long i, n, w;
  double x;
  Char c;

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
    q = p->next;
    while (q != p)
    {
      treeout3(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60)
      {
        putc('\n', outtree);
        *col = 0;
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

  (void)len;                            // RSGnote: Parameter never used.
  (void)ntips;                          // RSGnote: Parameter never used.

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
} /* bestelm */


/* End. */
