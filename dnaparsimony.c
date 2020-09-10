/* Version 4.0.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Michal Palczewski.
   */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "dnaparsimony.h"

extern long nonodes, endsite, outgrno, nextree, which;
extern boolean transvp;
extern steptr weight, category, alias, location, ally;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff, usertree;
extern double *threshwt;
extern sequence inputSequences;

double nsteps[maxuser], minsteps;
long **fsteps;
long minwhich;


tree* dnapars_tree_new(long nonodes, long spp)
{
  /* make new dnapars_tree */ 

  tree *t;
  
  t = generic_tree_new(nonodes, spp);
  pars_tree_init(t, nonodes, spp);
  dnapars_tree_init(t, nonodes, spp);
  return t;
} /* dnapars_tree_new */


void dnapars_tree_init(tree* t, long nonodes, long spp)
{
  /* set up functions for dnapars_tree */

  t->evaluate = dnapars_tree_evaluate;
  t->nuview = dnapars_tree_nuview;
  ((pars_tree*)t)->branchcollapsible = dna_branchcollapsible;
} /* dnapars_tree_init */


node* dnapars_node_new(node_type type, long index) // RSGbugfix
{
  /* make a new dnapars_node */

  node* n = Malloc(sizeof(dnapars_node));
  dnapars_node_init(n, type, index);
  return n;
} /* dnapars_node_new */


void dnapars_node_init(node* n, node_type type, long index)
{
  /* mostly, set up local functiona for a dnapars_node */
  dnapars_node *dn = (dnapars_node *)n;

  pars_node_init(n, type, index);

  n->copy = dnapars_node_copy;
  n->init = dnapars_node_init;
  n->reinit = dnapars_node_reinit;
  n->free = dnapars_node_free;

  if (dn->base)
    free(dn->base);
  dn->base = Malloc(endsite * sizeof(long));

  if (dn->numnuc) free(dn->numnuc);
  dn->numnuc = Malloc(endsite * sizeof(nucarray));
} /* dnapars_node_init */


void dnapars_node_reinit(node * n)
{
  /* reinitialize a dnapars_node */

  pars_node_reinit(n);

  dnapars_node *dn = (dnapars_node *)n;

  if (dn->base)
    free(dn->base);
  dn->base = Malloc(endsite * sizeof(long));

  if (dn->numnuc)
    free(dn->numnuc);
  dn->numnuc = Malloc(endsite * sizeof(nucarray));
} /* dnapars_node_reinit */


void dnapars_node_free(node **pp)
{
  /* free a dnapars_node */
  dnapars_node * dp = (dnapars_node *)*pp;

  free(dp->base);
  free(dp->numnuc);
  pars_node_free(pp);
} /* dnapars_node_free */


void dnapars_node_copy(node* srcn, node* dstn)
{
  /* copy a dnapars_node */
  dnapars_node *src = (dnapars_node *)srcn;
  dnapars_node *dst = (dnapars_node *)dstn;

  pars_node_copy(srcn, dstn);
  if (dst->base == NULL )
  {
    dst->base = Malloc(endsite * sizeof(long));
    dst->numnuc = Malloc(endsite * sizeof(nucarray));
  }
  memcpy(dst->base, src->base, endsite * sizeof(long));
  memcpy(dst->numnuc, src->numnuc, endsite * sizeof(nucarray));
} /* dnapars_node_copy */


void dna_initmin(dnapars_node *p, long sitei, boolean internal)
{
  /* some bookkeeping involving inferring branch lengths */
  long i;

  if (internal)
  {
    for (i = (long)A; i <= (long)O; i++)
    {
      p->cumlengths[i] = 0;
      p->numreconst[i] = 1;
    }
  }
  else
  {
    for (i = (long)A; i <= (long)O; i++)
    {
      if (p->base[sitei - 1] & (1 << i))
      {
        p->cumlengths[i] = 0;
        p->numreconst[i] = 1;
      }
      else
      {
        p->cumlengths[i] = -1;
        p->numreconst[i] = 0;
      }
    }
  }
} /* initmin */


void dna_compmin(node *p, node *desc)
{
  /* computes minimum lengths up to p */
  long i, j, minn, cost, desclen, descrecon=0, maxx;

  maxx = 10 * spp;
  for (i = (long)A; i <= (long)O; i++)
  {
    minn = maxx;
    for (j = (long)A; j <= (long)O; j++)
    {
      if (transvp)
      {
        if ((((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G)))
            || (((j == (long)C) || (j == (long)T))
                && ((i == (long)C) || (i == (long)T))))
          cost = 0;
        else
          cost = 1;
      }
      else
      {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (((dnapars_node*)desc)->cumlengths[j] == -1)
      {
        desclen = maxx;
      }
      else
      {
        desclen = ((dnapars_node*)desc)->cumlengths[j];
      }
      if (minn > cost + desclen)
      {
        minn = cost + desclen;
        descrecon = 0;
      }
      if (minn == cost + desclen)
      {
        descrecon += ((dnapars_node*)desc)->numreconst[j];
      }
    }
    ((dnapars_node*)p)->cumlengths[i] += minn;
    ((dnapars_node*)p)->numreconst[i] *= descrecon;
  }
  p->initialized = true;
} /* compmin */


void branchlength(node *subtr1, node *subtr2, double *brlen, pointarray treenode)
{
  /* computes a branch length between two subtrees for a given site */
  long i, j, minn, cost, nom, denom;
  node *temp;

  if (subtr1->tip)
  {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  if (subtr1->index == outgrno)
  {
    temp = subtr1;
    subtr1 = subtr2;
    subtr2 = temp;
  }
  minpostorder(subtr1, treenode);
  minpostorder(subtr2, treenode);
  minn = 10 * spp;
  nom = 0;
  denom = 0;
  for (i = (long)A; i <= (long)O; i++)
  {
    for (j = (long)A; j <= (long)O; j++)
    {
      if (transvp)
      {
        if ((((i == (long)A) || (i == (long)G))
             && ((j == (long)A) || (j == (long)G)))
            || (((j == (long)C) || (j == (long)T))
                && ((i == (long)C) || (i == (long)T))))
          cost = 0;
        else
          cost = 1;
      }
      else
      {
        if (i == j)
          cost = 0;
        else
          cost = 1;
      }
      if (((dnapars_node*)subtr1)->cumlengths[i] != -1 &&
          (((dnapars_node*)subtr2)->cumlengths[j] != -1))
      {
        if (((dnapars_node*)subtr1)->cumlengths[i] +
            cost + ((dnapars_node*)subtr2)->cumlengths[j] < minn)
        {
          minn = ((dnapars_node*)subtr1)->cumlengths[i] + cost +
            ((dnapars_node*)subtr2)->cumlengths[j];
          nom = 0;
          denom = 0;
        }
        if (((dnapars_node*)subtr1)->cumlengths[i] + cost +
            ((dnapars_node*)subtr2)->cumlengths[j] == minn)
        {
          nom += ((dnapars_node*)subtr1)->numreconst[i] *
            ((dnapars_node*)subtr2)->numreconst[j] * cost;
          denom += ((dnapars_node*)subtr1)->numreconst[i] *
            ((dnapars_node*)subtr2)->numreconst[j];
        }
      }
    }
  }
  *brlen = (double)nom/(double)denom;
} /* branchlength */


void inittreetrav(node *p, long sitei)
{
  /* traverse tree to clear boolean initialized and set up base */
  node *q;

  if (p->tip)
  {
    dna_initmin((dnapars_node*)p, sitei, false);
    p->initialized = true;
    return;
  }
  q = p->next;
  while (q != p)
  {
    inittreetrav(q->back, sitei);
    q = q->next;
  }
  dna_initmin((dnapars_node*)p, sitei, true);
  p->initialized = false;
  q = p->next;
  while (q != p)
  {
    dna_initmin((dnapars_node*)q, sitei, true);
    q->initialized = false;
    q = q->next;
  }
} /* inittreetrav */


void minpostorder(node *p, pointarray treenode)
{
  /* traverses an n-ary tree, computing minimum steps at each node */
  node *q;

  if (p->tip)
  {
    return;
  }
  q = p->next;
  while (q != p)
  {
    if (q->back)
      minpostorder(q->back, treenode);
    q = q->next;
  }
  if (!p->initialized)
  {
    q = p->next;
    while (q != p)
    {
      if (q->back)
        dna_compmin(p, q->back);
      q = q->next;
    }
  }
}  /* minpostorder */


void branchlentrav(node *p, node *root, long sitei, long chars, double *brlen, pointarray treenode)
{
  /*  traverses the tree computing tree length at each branch */
  node *q;

  if (p->tip)
    return;
  if (p->index == outgrno)
    p = p->back;
  q = p->next;
  do {
    if (q->back)
    {
      branchlength(q, q->back, brlen, treenode);
      q->v += (weight[sitei - 1] * (*brlen)/chars);
      q->back->v += (weight[sitei - 1] * (*brlen)/chars);
      if (!q->back->tip)
        branchlentrav(q->back, root, sitei, chars, brlen, treenode);
    }
    q = q->next;
  } while (q != p);
}  /* branchlentrav */


void dna_treelength(node *root, long chars, pointarray treenode)
{ 
  /*  calls branchlentrav at each site */
  long sitei;
  double trlen;

  initbranchlen(root);
  for (sitei = 1; sitei <= endsite; sitei++)
  {
    trlen = 0.0;
    dna_initbase(root, sitei);
    inittreetrav(root, sitei);
    branchlentrav(root, root, sitei, chars, &trlen, treenode);
  }
} /* treelength */


void dna_initbase(node *p, long sitei)
{
  /* traverse tree to initialize base at internal nodes */
  node *q;
  long i, largest;

  if (p->tip)
    return;
  q = p->next;
  while (q != p)
  {
    if (q->back)
    {
      memcpy(((dnapars_node*)q)->numnuc, ((dnapars_node*)p)->numnuc, endsite * sizeof(nucarray));
      for (i = (long)A; i <= (long)O; i++)
      {
        if (((dnapars_node*)q->back)->base[sitei - 1] & (1 << i))
          ((dnapars_node*)q)->numnuc[sitei - 1][i]--;
      }
      if (p->back)
      {
        for (i = (long)A; i <= (long)O; i++)
        {
          if (((dnapars_node*)p->back)->base[sitei - 1] & (1 << i))
            ((dnapars_node*)q)->numnuc[sitei - 1][i]++;
        }
      }
      largest = dna_getlargest(((dnapars_node*)q)->numnuc[sitei - 1]);
      ((dnapars_node*)q)->base[sitei - 1] = 0;
      for (i = (long)A; i <= (long)O; i++)
      {
        if (((dnapars_node*)q)->numnuc[sitei - 1][i] == largest)
          ((dnapars_node*)q)->base[sitei - 1] |= (1 << i);
      }
    }
    q = q->next;
  }
  q = p->next;
  while (q != p)
  {
    dna_initbase(q->back, sitei);
    q = q->next;
  }
} /* initbase */


long dna_getlargest(long *numnuc)
{
  /* find the largest in array numnuc */
  long i, largest;

  largest = 0;
  for (i = (long)A; i <= (long)O; i++)
    if (numnuc[i] > largest)
      largest = numnuc[i];
  return largest;
} /* dna_getlargest */


void dna_hyptrav(tree* t, node *r_, long *hypset_, long b1, long b2, boolean bottom_, Char *basechar)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  baseptr  ancset;
  nucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  ancset = Malloc(endsite * sizeof(long));

  tempnuc = (nucarray *)Malloc(endsite * sizeof(nucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    memset(((dnapars_node*)Vars.r)->numnuc, 0, endsite * sizeof(nucarray));
  for (i = b1 - 1; i < b2; i++)
  {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip)
    {
      p = Vars.r->next;
      for (k = (long)A; k <= (long)O; k++)
        if (Vars.anc & (1 << k))
          ((dnapars_node*)Vars.r)->numnuc[j - 1][k]++;
      do {
        for (k = (long)A; k <= (long)O; k++)
          if (((dnapars_node*)p->back)->base[j - 1] & (1 << k))
            ((dnapars_node*)Vars.r)->numnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest = dna_getlargest(((dnapars_node*)Vars.r)->numnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)A; k <= (long)O; k++)
      {
        if (((dnapars_node*)Vars.r)->numnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      ((dnapars_node*)Vars.r)->base[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = ((dnapars_node*)t->nodep[Vars.r->back->index - 1])
        ->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (((dnapars_node*)Vars.r)->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || ((dnapars_node*)Vars.r)->base[j - 1] != Vars.anc);
  }
  dna_hyprint(t, b1, b2, &Vars, basechar);
  Vars.bottom = false;
  if (!Vars.r->tip)
  {
    memcpy(tempnuc, ((dnapars_node*)Vars.r)->numnuc, endsite * sizeof(nucarray));
    q = Vars.r->next;
    do
    {
      memcpy(((dnapars_node*)Vars.r)->numnuc, tempnuc, endsite * sizeof(nucarray));
      for (i = b1 - 1; i < b2; i++)
      {
        j = location[ally[i] - 1];
        for (k = (long)A; k <= (long)O; k++)
          if (((dnapars_node*)q->back)->base[j - 1] & (1 << k))
            ((dnapars_node*)Vars.r)->numnuc[j - 1][k]--;
        largest = dna_getlargest(((dnapars_node*)Vars.r)->numnuc[j - 1]);
        ancset[j - 1] = 0;
        for (k = (long)A; k <= (long)O; k++)
          if (((dnapars_node*)Vars.r)->numnuc[j - 1][k] == largest)
            ancset[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset[j - 1];
      }
      dna_hyptrav(t, q->back, ancset, b1, b2, Vars.bottom, basechar);
      q = q->next;
    } while (q != Vars.r);
  }
  free(tempnuc);
  free(ancset);
}  /* hyptrav */


void dna_hypstates(tree* t, long chars, Char *basechar)
{
  /* fill in and describe states at interior nodes */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, n;
  baseptr nothing;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  if (dotdiff)
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (baseptr)Malloc(endsite * sizeof(long));
  memset(nothing, 0, endsite * sizeof(long));
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++)
  {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    dna_hyptrav(t, t->root, nothing, i * 40 - 39, n, true, basechar);
  }
  free(nothing);
}  /* dna_hypstates */


void dna_hyprint(tree* t, long b1, long b2, struct LOC_hyptrav *htrav, Char *basechar)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k, n;
  boolean dot;
  bases b;

  if (htrav->bottom)
  {
    if (!outgropt)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  }
  else
    fprintf(outfile, "%4ld   ", htrav->r->back->index - spp);
  if (htrav->r->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[htrav->r->index - 1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", htrav->r->index - spp);
  if (htrav->bottom)
    fprintf(outfile, "          ");
  else if (htrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (htrav->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++)
  {
    j = location[ally[i - 1] - 1];
    htrav->tempset = ((dnapars_node*)htrav->r)->base[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc = ((dnapars_node*)t->nodep[htrav->r->back->index - 1])->base[j - 1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)
      putc('.', outfile);
    else if (htrav->tempset == (1 << A))
      putc('A', outfile);
    else if (htrav->tempset == (1 << C))
      putc('C', outfile);
    else if (htrav->tempset == (1 << G))
      putc('G', outfile);
    else if (htrav->tempset == (1 << T))
      putc('T', outfile);
    else if (htrav->tempset == (1 << O))
      putc('-', outfile);
    else
    {
      k = 1;
      n = 0;
      for (b = A; b <= O; b = b + 1)
      {
        if (((1 << b) & htrav->tempset) != 0)
          n += k;
        k += k;
      }
      putc(basechar[n - 1], outfile);
    }
    if (i % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


double dnapars_tree_evaluate(tree* t, node *n, boolean saveit)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, steps;
  double term;
  double sum;
  dnapars_node* p = (dnapars_node*)n;
  dnapars_node* q = (dnapars_node*)n->back;
  long base1, base2;
  sum = 0.0;

  generic_tree_evaluate(t, n, saveit);

  for (i = 0; i < endsite; i++)
  {
    steps = ((pars_node*)p)->numsteps[i] + ((pars_node*)q)->numsteps[i];
    base1 = p->base[i];
    base2 = q->base[i];
    if ( transvp )
    {
      if (base1 & purset) base1 = purset;
      if (base1 & pyrset) base1 = pyrset;
    }
    if ( (base1 & base2) == 0 )
      steps += weight[i];
    if ( ((pars_tree*)t)->supplement)
      steps += ((pars_tree*)t)->supplement(t, i);

    if (steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser)
  {
    nsteps[which - 1] = sum;
    if (which == 1)
    {
      minwhich = 1;
      minsteps = sum;
    }
    else if (sum < minsteps)
    {
      minwhich = which;
      minsteps = sum;
    }
  }
  t->score = -sum;
  return t->score;
}  /* dnapars_tree_evaluate */


void dnapars_tree_nuview(tree* t, node* p)
{ 
  /* calculate the view for all endsite sites 
   * includes possibility of multifurcations */
  node *q;
  dnapars_node *qback;
  long i, j, base1, newbase, steps, largest;
  long numnuc[5] = {0, 0, 0, 0, 0};
  boolean bif;
  long root = 0;

/* debug:   generic_tree_nuview(t, p);      not needed, generic will call appropriate one */
  bif = (count_sibs(p) == 2);

  for ( i = 0 ; i < endsite ; i++ )
  {
    newbase = 0xff;
    steps = 0;
    for ( q = p->next ; q != p ; q = q->next )
    {
      qback = (dnapars_node*)q->back;
      if ( qback == NULL )
      {
        root = 1;
        continue;
      }
      /* if we root the tree we can safely ignore the root*/
      base1 = qback->base[i];
      if ( transvp )
      {
        if (base1 & purset) base1 = purset;
        if (base1 & pyrset) base1 = pyrset;
      }
      newbase &= base1;
      steps += ((pars_node*)qback)->numsteps[i];
    }
    if ( newbase == 0 )
    {
      if ( !bif )
      {
        memset(numnuc, 0, sizeof(numnuc));
        for (j = (long)A; j <= (long)O; j++)
        {
          for ( q = p->next ; q != p ; q = q->next )
          {
            qback = (dnapars_node*) q->back;
            if ( qback == NULL )
              continue;
            if ( qback->base[i] & (1 << j) )
              numnuc[j]++;
          }
        }
        largest = dna_getlargest(numnuc);
        for (j = (long)A; j <= (long)O; j++)
        {
          if (numnuc[j] == largest )
          {
            newbase |= 1 << j;
          }
        }
        steps += (weight[i]) * (count_sibs(p) - largest - root);
      }
      else
      { /* optimization for bifurcation, code above still works though*/
        newbase = ((dnapars_node*)(p->next->back))->base[i] | ((dnapars_node*)(p->next->next->back))->base[i];
        steps += weight[i];
      }
    }

    ((dnapars_node*)p)->base[i] = newbase;
    ((pars_node*)p)->numsteps[i] = steps;
  }

  p->initialized = true;
} /* dnapars_tree_nuview */


boolean dna_branchcollapsible(tree* t, node* n)
{
  /* dnapars version of checking whether a branch is collapsible */
  boolean collapsible = true;
  long i;
  node* q;

  if ( n->tip == true || n->back->tip == true )
    return false;

  q = n->back;
  if ( q->initialized == false ) t->nuview(t, q);
  if ( n->initialized == false ) t->nuview(t, n);

  for ( i = 0 ; i < endsite ; i++ )
  {
    if ( (((dnapars_node*)q)->base[i] & ((dnapars_node*)n)->base[i] ) == 0)
      return false;
  }
  return collapsible;
} /* dna_branchcollapse */


void dna_makevalues(tree* t, boolean usertree)
{
  /* set up fractional likelihoods at tips */
  /* used by dnacomp, dnapars, & dnapenny */
  long i, j;
  char ns = 0;

  for (j = 0; j < endsite; j++)
  {
    for (i = 0; i < spp; i++)
    {
      switch (inputSequences[i][alias[j] - 1])
      {
        case 'A':
          ns = 1 << A;
          break;

        case 'C':
          ns = 1 << C;
          break;

        case 'G':
          ns = 1 << G;
          break;

        case 'U':
          ns = 1 << T;
          break;

        case 'T':
          ns = 1 << T;
          break;

        case 'M':
          ns = (1 << A) | (1 << C);
          break;

        case 'R':
          ns = (1 << A) | (1 << G);
          break;

        case 'W':
          ns = (1 << A) | (1 << T);
          break;

        case 'S':
          ns = (1 << C) | (1 << G);
          break;

        case 'Y':
          ns = (1 << C) | (1 << T);
          break;

        case 'K':
          ns = (1 << G) | (1 << T);
          break;

        case 'B':
          ns = (1 << C) | (1 << G) | (1 << T);
          break;

        case 'D':
          ns = (1 << A) | (1 << G) | (1 << T);
          break;

        case 'H':
          ns = (1 << A) | (1 << C) | (1 << T);
          break;

        case 'V':
          ns = (1 << A) | (1 << C) | (1 << G);
          break;

        case 'N':
          ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
          break;

        case 'X':
          ns = (1 << A) | (1 << C) | (1 << G) | (1 << T);
          break;

        case '?':
          ns = (1 << A) | (1 << C) | (1 << G) | (1 << T) | (1 << O);
          break;

        case 'O':
          ns = 1 << O;
          break;

        case '-':
          ns = 1 << O;
          break;
      }
      ((dnapars_node*)t->nodep[i])->base[j] = ns;
      ((pars_node*)t->nodep[i])->numsteps[j] = 0;
    }
  }
}  /* makevalues */


// End.
