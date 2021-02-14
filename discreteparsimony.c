/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "discreteparsimony.h"

#define SAMPLES 1000

long nonodes, endsite, outgrno, nextree, which;
boolean interleaved, printdata, outgropt, treeprint, dotdiff;
steptr weight, category, alias, location, ally;
sequence inputSequences;

extern double *threshwt;
extern double **fsteps;
extern long minwhich;
extern boolean usertree, reusertree;
extern double nsteps[maxuser], minsteps;


tree* discretepars_tree_new(long nonodes, long spp)
{
  /* allocate and call initialization routines for a discrete
   * characters parsimony method tree */

  tree* t = generic_tree_new(nonodes, spp);
  parsimony_tree_init(t, nonodes, spp);
  discretepars_tree_init(t, nonodes, spp);
  return t;
} /* discretepars_tree_new */


void discretepars_tree_init(tree* t, long nonodes, long spp)
{
  /* initialize a tree for discrete-character parsimony methods */

  t->nuview = discretepars_tree_nuview;
  t->evaluate = discretepars_tree_evaluate;
  ((pars_tree*)t)->branchcollapsible = discretepars_tree_branchcollapsible;
} /* discretepars_tree_init */


node * discretepars_node_new(node_type type, long index) // RSGbugfix
{
  /* allocate and initialize a node for a discrete data parsimony method */

  node *n = Malloc(sizeof(discretepars_node));
  discretepars_node_init(n, type, index);
  return n;
} /* discretepars_node */


void discretepars_node_init(node* node, node_type type, long index)
{
  /* initialize a discrete characters parsimony method node */
  discretepars_node *n = (discretepars_node *)node;

  pars_node_init(node, type, index);
  node->free = discretepars_node_free;

  n->discbase = (discbaseptr)Malloc(endsite * sizeof(unsigned char));
  n->discnumnuc = (discnucarray *)Malloc(endsite * sizeof(discnucarray));
} /* discretepars_node_init */


void discretepars_node_free(node **n)
{
  /* free a discrete characters parsimony method node */
  pars_node_free(n);
} /* discretepars_node_free */


void inputdata(long chars)
{
  /* input the names and sequences for each species used by pars */
  long i, j, k, l;
  long basesread=0, basesnew=0, nsymbol=0, convsymboli=0;
  Char charstate;
  boolean allread, done, found;

  if (printdata)
    headings(chars, "Sequences", "---------");
  basesread = 0;
  allread = false;
  while (!(allread))
  {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp)
    {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i - 1);
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eoff(infile))
      {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile)))
        {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ')
            continue;
          if ((strchr("!\"#$%&'()*+,-./0123456789:;<=>?@\
                       ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`                \
                       abcdefghijklmnopqrstuvwxyz{|}~", charstate)) == NULL)
          {
            printf("\nERROR:  Bad symbol: %c at position %ld of species %ld.\n\n", charstate, j+1, i);
            exxit(-1);
          }
          j++;
          inputSequences[i - 1][j - 1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < chars)
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;

      scan_eoln(infile);

      if ((interleaved && j != basesnew) ||
          (!interleaved && j != chars))
      {
        printf("\nERROR:  Sequences out of alignment at position %ld.\n\n", j);
        exxit(-1);
      }
      i++;
    }
    if (interleaved)
    {
      basesread = basesnew;
      allread = (basesread == chars);
    }
    else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata)
  {
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++)
    {
      for (j = 1; j <= spp; j++)
      {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars)
          l = chars;
        for (k = (i - 1) * 60 + 1; k <= l; k++)
        {
          if (dotdiff && (j > 1 && inputSequences[j - 1][k - 1] == inputSequences[0][k - 1]))
            charstate = '.';
          else
            charstate = inputSequences[j - 1][k - 1];
          putc(charstate, outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 1; i <= chars; i++)
  {
    nsymbol = 0;
    for (j = 1; j <= spp; j++)
    {
      if ((nsymbol == 0) && (inputSequences[j - 1][i - 1] != '?'))
      {
        nsymbol = 1;
        convsymboli = 1;
        convtab[0][i-1] = inputSequences[j-1][i-1];
      }
      else if (inputSequences[j - 1][i - 1] != '?')
      {
        found = false;
        for (k = 1; k <= nsymbol; k++)
        {
          if (convtab[k - 1][i - 1] == inputSequences[j - 1][i - 1])
          {
            found = true;
            convsymboli = k;
          }
        }
        if (!found)
        {
          nsymbol++;
          convtab[nsymbol-1][i - 1] = inputSequences[j - 1][i - 1];
          convsymboli = nsymbol;
        }
      }
      if (nsymbol <= 8)
      {
        if (inputSequences[j - 1][i - 1] != '?')
          inputSequences[j - 1][i - 1] = (Char)('0' + (convsymboli - 1));
      }
      else
      {
        printf("\nERROR:  More than maximum of 8 symbols in column %ld.\n\n", i);
        exxit(-1);
      }
    }
  }
}  /* inputdata */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in pars */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;
  while (gap > 0)
  {
    for (i = gap + 1; i <= chars; i++)
    {
      j = i - gap;
      flip = true;
      while (j > 0 && flip)
      {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = true;
        k = 1;
        while (k <= spp && tied)
        {
          flip = (inputSequences[k - 1][jj - 1] > inputSequences[k - 1][jg - 1]);
          tied = (tied && inputSequences[k - 1][jj - 1] == inputSequences[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = weight[j - 1];
        weight[j - 1] = weight[j + gap - 1];
        weight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */


void sitecombine(long chars)
{
  /* combine sites that have identical patterns */
  /* used in pars */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < chars)
  {
    j = i + 1;
    tied = true;
    while (j <= chars && tied)
    {
      k = 1;
      while (k <= spp && tied)
      {
        tied = (tied &&
                inputSequences[k - 1][alias[i - 1] - 1] == inputSequences[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (tied)
      {
        weight[i - 1] += weight[j - 1];
        weight[j - 1] = 0;
        ally[alias[j - 1] - 1] = alias[i - 1];
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */


void sitescrunch(long chars)
{
  /* move so one representative of each pattern of
     sites comes first */
  /* used in pars */
  long i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done)
  {
    if (ally[alias[i - 1] - 1] != alias[i - 1])
    {
      if (j <= i)
        j = i + 1;
      if (j <= chars)
      {
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
        } while (!(found || j > chars));
        if (found)
        {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = weight[i - 1];
          weight[i - 1] = weight[j - 1];
          weight[j - 1] = itemp;
        }
        else
          done = true;
      }
      else
        done = true;
    }
    i++;
    done = (done || i >= chars);
  }
}  /* sitescrunch */


void makevalues(tree *t, boolean usertree)
{
  long i, j;
  unsigned char ns=0;

  for (j = 0; j < endsite; j++)
  {
    for (i = 0; i < spp; i++)
    {
      switch (inputSequences[i][alias[j] - 1])
      {
        case '0':
          ns = 1 << zero;
          break;

        case '1':
          ns = 1 << one;
          break;

        case '2':
          ns = 1 << two;
          break;

        case '3':
          ns = 1 << three;
          break;

        case '4':
          ns = 1 << four;
          break;

        case '5':
          ns = 1 << five;
          break;

        case '6':
          ns = 1 << six;
          break;

        case '7':
          ns = 1 << seven;
          break;

        case '?':
          ns = (1 << zero) | (1 << one) | (1 << two) | (1 << three) |
            (1 << four) | (1 << five) | (1 << six) | (1 << seven);
          break;
      }
      ((discretepars_node*)t->nodep[i])->discbase[j] = ns;
      ((pars_node*)t->nodep[i])->numsteps[j] = 0;
    }
  }
}  /* makevalues */


long discgetlargest(long *discnumnuc)
{
  /* find the largest in array numnuc */
  long i, largest;

  largest = 0;
  for (i = (long)zero; i <= (long)seven; i++)
    if (discnumnuc[i] > largest)
      largest = discnumnuc[i];
  return largest;
} /* discgetlargest */


void dischyptrav(tree* t, node *r_, discbaseptr hypset_, long b1, long b2, boolean bottom_)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  discbaseptr ancset;
  discnucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  ancset = Malloc(endsite * sizeof(unsigned char));
  tempnuc = (discnucarray *)Malloc(endsite * sizeof(discnucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    memset(((discretepars_node*)Vars.r)->discnumnuc, 0, endsite * sizeof(discnucarray));
  for (i = b1 - 1; i < b2; i++)
  {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip)
    {
      p = Vars.r->next;
      for (k = (long)zero; k <= (long)seven; k++)
        if (Vars.anc & (1 << k))
          ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]++;
      do {
        for (k = (long)zero; k <= (long)seven; k++)
          if (((discretepars_node*)p->back)->discbase[j - 1] & (1 << k))
            ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest =
        discgetlargest(((discretepars_node*)Vars.r)->discnumnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)zero; k <= (long)seven; k++)
      {
        if (((discretepars_node*)Vars.r)->discnumnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      ((discretepars_node*)Vars.r)->discbase[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = ((discretepars_node*)t->nodep[Vars.r->back->index - 1]) ->discbase[j - 1];
    Vars.nonzero = (Vars.nonzero || (((discretepars_node*)Vars.r)->discbase[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || ((discretepars_node*)Vars.r)->discbase[j - 1] != Vars.anc);
  }
  dischyprint(t, b1, b2, &Vars);
  Vars.bottom = false;
  if (!Vars.r->tip)
  {
    memcpy(tempnuc, ((discretepars_node*)Vars.r)->discnumnuc, endsite * sizeof(discnucarray));
    q = Vars.r->next;
    do {
      memcpy(((discretepars_node*)Vars.r)->discnumnuc, tempnuc, endsite * sizeof(discnucarray));
      for (i = b1 - 1; i < b2; i++)
      {
        j = location[ally[i] - 1];
        for (k = (long)zero; k <= (long)seven; k++)
          if (((discretepars_node*)q->back)->discbase[j - 1] & (1 << k))
            ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]--;
        largest = discgetlargest(((discretepars_node*)Vars.r)->
                                 discnumnuc[j - 1]);
        ancset[j - 1] = 0;
        for (k = (long)zero; k <= (long)seven; k++)
          if (((discretepars_node*)Vars.r)->discnumnuc[j - 1][k] == largest)
            ancset[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset[j - 1];
      }
      dischyptrav(t, q->back, ancset, b1, b2, Vars.bottom);
      q = q->next;
    } while (q != Vars.r);
  }
  free(ancset);
}  /* dischyptrav */


void disc_hypstates(tree* t, long chars)
{
  /* fill in and describe states at interior nodes */
  /* used in pars */
  long i, n;
  discbaseptr nothing;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                            ");
  if (dotdiff)
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (discbaseptr)Malloc(endsite * sizeof(unsigned char));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++)
  {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    dischyptrav(t, t->root, nothing, i * 40 - 39, n, true);
  }
  free(nothing);
}  /* hypstates */



void discinitbase(node *p, long sitei)
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
      memcpy(((discretepars_node*)q)->discnumnuc,
             ((discretepars_node*)p)->discnumnuc, endsite * sizeof(discnucarray));
      for (i = (long)zero; i <= (long)seven; i++)
      {
        if (((discretepars_node*)q->back)->discbase[sitei - 1] & (1 << i))
          ((discretepars_node*)q)->discnumnuc[sitei - 1][i]--;
      }
      if (p->back)
      {
        for (i = (long)zero; i <= (long)seven; i++)
        {
          if (((discretepars_node*)p->back)->discbase[sitei - 1] & (1 << i))
            ((discretepars_node*)q)->discnumnuc[sitei - 1][i]++;
        }
      }
      largest = discgetlargest(((discretepars_node*)q)->discnumnuc[sitei - 1]);
      ((discretepars_node*)q)->discbase[sitei - 1] = 0;
      for (i = (long)zero; i <= (long)seven; i++)
      {
        if (((discretepars_node*)q)->discnumnuc[sitei - 1][i] == largest)
          ((discretepars_node*)q)->discbase[sitei - 1] |= (1 << i);
      }
    }
    q = q->next;
  }
  q = p->next;
  while (q != p)
  {
    if (q->back != NULL)
      discinitbase(q->back, sitei);
    q = q->next;
  }
} /* initbase */


void inittreetrav(node *p, long sitei)
{
  /* traverse tree to clear boolean initialized and set up base */
  node *q;

  if (p == NULL)                                /* unless it's an empty tip */
    return;
  if (p->tip)
  {
    discinitmin((discretepars_node*)p, sitei, false);   /* initialize state */
    p->initialized = true;                       /* mark tip as initialized */
    return;
  }
  q = p->next;                                  /* for an interior fork ... */
  while (q != p)
  {
    if (q != NULL)
      inittreetrav(q->back, sitei);          /* ... go recursively outwards */
    q = q->next;
  }
  discinitmin((discretepars_node*)p, sitei, true); /* initializing state,,, */
  p->initialized = false;      /* ... marking them as needing to be updated */
  q = p->next;
  while (q != p)                 /* ,,, and continue around fork doing that */
  {
    discinitmin((discretepars_node*)q, sitei, true);
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
        disccompmin(p, q->back);
      q = q->next;
    }
  }
}  /* minpostorder */


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
  denom = 1;
  for (i = (long)zero; i <= (long)seven; i++)
  {
    for (j = (long)zero; j <= (long)seven; j++)
    {
      if (i == j)
        cost = 0;
      else
        cost = 1;
      if (((discretepars_node*)subtr1)->disccumlengths[i] != -1 &&
          (((discretepars_node*)subtr2)->disccumlengths[j] != -1))
      {
        if (((discretepars_node*)subtr1)->disccumlengths[i] + cost +
            ((discretepars_node*)subtr2)->disccumlengths[j] < minn)
        {
          minn = ((discretepars_node*)subtr1)->disccumlengths[i] + cost +
            ((discretepars_node*)subtr2)->disccumlengths[j];
          nom = 0;
          denom = 1;
        }
        if (((discretepars_node*)subtr1)->disccumlengths[i] + cost +
            ((discretepars_node*)subtr2)->disccumlengths[j] == minn)
        {
          nom += ((discretepars_node*)subtr1)->discnumreconst[i] *
            ((discretepars_node*)subtr2)->discnumreconst[j] * cost;
          denom += ((discretepars_node*)subtr1)->discnumreconst[i] *
            ((discretepars_node*)subtr2)->discnumreconst[j];
        }
      }
    }
  }
  *brlen = (double)nom/(double)denom;
} /* branchlength */


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
      q->v += (weight[sitei - 1]  * (*brlen));
      q->back->v += (weight[sitei - 1] * (*brlen));
      if (!q->back->tip)
        branchlentrav(q->back, root, sitei, chars, brlen, treenode);
    }
    q = q->next;
  } while (q != p);
}  /* branchlentrav */



void treeout(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in pars */
  node *q;
  long i, n;
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
      treeout(q->back, nextree, col, root);
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
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout */


void standev(long chars, long numtrees, long minwhich, double minsteps, double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH) on user trees */
   /* used in pars */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

  if (numtrees == 2)
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees)
    {
      fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1]);
      if (minwhich == which)
        fprintf(outfile, "  <------ best\n");
      else
      {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < endsite; i++)
        {
          if (weight[i] > 0)
          {
            wt = weight[i];
            sumw += wt;
            temp = (fsteps[which - 1][i] - fsteps[minwhich - 1][i]);
            sum += temp;
            sum2 += temp * temp / wt;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%9.1f %12.4f",
                (nsteps[which - 1] - minsteps), sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  }
  else
  {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n" , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    }
    else
    {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees * sizeof(double *));
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees * sizeof(double));
    sumw = 0.0;
    for (i = 0; i < endsite; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++)      /* compute covariances of trees */
    {
      sum = nsteps[i]/sumw;
      for (j = 0; j <=i; j++)
      {
        sum2 = nsteps[j]/sumw;
        temp = 0.0;
        for (k = 0; k < endsite; k++)
        {
          wt = weight[k];
          if (weight[k] > 0)
          {
            temp = temp + wt*(fsteps[i][k]/(wt)-sum) *(fsteps[j][k]/(wt)-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++)                /* Cholesky decomposition */
    {
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++)
      {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-23)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees * sizeof(double)); /* resampled sum steps */
    P = (double *)Malloc(numtrees * sizeof(double));      /* vector of P's  */
    r = (double *)Malloc(numtrees * sizeof(double));     /* normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0];               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i])
        sum2 = nsteps[i];
    for (i = 1; i <= SAMPLES; i++)      /* loop over resampled trees */
    {
      for (j = 0; j < numtrees; j++)    /* draw normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++)    /* compute vectors */
      {
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get min of vector */
        if (f[j] < sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (nsteps[j]-sum2 < f[j] - sum)
          P[j] += 1.0 / SAMPLES;
    }
    fprintf(outfile, "Tree    Steps   Diff Steps   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++)
    {
      fprintf(outfile, "%3ld%10.1f", i+1, nsteps[i]);
      if ((minwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        fprintf(outfile, "  %9.1f %10.3f", nsteps[i]-sum2, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
    fprintf(outfile, "\n");
    free(P);             /* free the variables we Malloc'ed */
    free(f);
    free(r);
    for (i = 0; i < numtrees; i++)
      free(covar[i]);
    free(covar);
  }
}  /* standev */


void disc_treelength(node *root, long chars, pointarray treenode)
{
  /*  calls branchlentrav at each site */
  long sitei;
  double trlen;

  initbranchlen(root);
  for (sitei = 1; sitei <= endsite; sitei++)
  {
    trlen = 0.0;
    discinitbase(root, sitei);
    inittreetrav(root, sitei);
    branchlentrav(root, root, sitei, chars, &trlen, treenode);
  }
} /* treelength */


void discinitmin(discretepars_node *p, long sitei, boolean internal)
{
  long i;

  if (internal)
  {
    for (i = (long)zero; i <= (long)seven; i++)
    {
      p->disccumlengths[i] = 0;
      p->discnumreconst[i] = 1;
    }
  }
  else
  {
    for (i = (long)zero; i <= (long)seven; i++)
    {
      if (p->discbase[sitei - 1] & (1 << i))
      {
        p->disccumlengths[i] = 0;
        p->discnumreconst[i] = 1;
      }
      else
      {
        p->disccumlengths[i] = -1;
        p->discnumreconst[i] = 0;
      }
    }
  }
} /* initdiscmin */


void dischyprint(tree* t, long b1, long b2, struct LOC_hyptrav *htrav)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k;
  boolean dot, found;

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
    htrav->tempset = ((discretepars_node*)htrav->r)->discbase[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc = ((discretepars_node*)t->nodep[htrav->r->back->index - 1]) ->discbase[j - 1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)
      putc('.', outfile);
    else
    {
      found = false;
      k = (long)zero;
      do {
        if (htrav->tempset == (1 << k))
        {
          putc(convtab[k][i - 1], outfile);
          found = true;
        }
        k++;
      } while (!found && k <= (long)seven);
      if (!found)
        putc('?', outfile);
    }
    if (i % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


void disccompmin(node *p, node *desc)
{
  /* computes minimum lengths up to p */
  long i, j, minn, cost, desclen, descrecon=0, maxx;

  maxx = 10 * spp;
  for (i = (long)zero; i <= (long)seven; i++)
  {
    minn = maxx;
    for (j = (long)zero; j <= (long)seven; j++)
    {
      if (i == j)
        cost = 0;
      else
        cost = 1;
      if (((discretepars_node*)desc)->disccumlengths[j] == -1)
      {
        desclen = maxx;
      }
      else
      {
        desclen = ((discretepars_node*)desc)->disccumlengths[j];
      }
      if (minn > cost + desclen)
      {
        minn = cost + desclen;
        descrecon = 0;
      }
      if (minn == cost + desclen)
      {
        descrecon += ((discretepars_node*)desc)->discnumreconst[j];
      }
    }
    ((discretepars_node*)p)->disccumlengths[i] += minn;
    ((discretepars_node*)p)->discnumreconst[i] *= descrecon;
  }
  p->initialized = true;
} /* compmin */


void discretepars_tree_nuview(tree* t, node*p)
{
  /* nuview for multistate discrete characters parsimony methods */
  node *q;
  discretepars_node *qback;
  long i, j, steps, largest;
  unsigned char base1, newbase;
  long numnuc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  boolean bif;
  long root = 0;

/* debug: already called?     generic_tree_nuview(t, p);  */
  bif = (count_sibs(p) == 2);

  for ( i = 0 ; i < endsite ; i++ ) /* do for each representative character */
  {
    newbase = 0xff;
    steps = 0;
    for ( q = p->next ; q != p ; q = q->next )
    {                                    /* go around the loop at this fork */
      qback = (discretepars_node*)q->back;
      if ( qback == NULL )
      {
        root = 1;
        continue;
      }
      /* if we root the tree we can safely ignore the root */
      base1 = qback->discbase[i];
      newbase &= base1;               /* intersection with base sets so far */
      steps += ((pars_node*)qback)->numsteps[i];
    }
    if ( newbase == 0 )
    {
      if ( !bif )
      {                                /* case where fork is multifurcating */
        memset(numnuc, 0, sizeof(numnuc));
        for (j = (long)zero; j <= (long)seven; j++)
        {
          for ( q = p->next ; q != p ; q = q->next )
          {
            qback = (discretepars_node*) q->back;
            if ( qback == NULL )
              continue;
            if ( qback->discbase[i] & (1 << j) )
              numnuc[j]++;
          }
        }
        largest = discgetlargest(numnuc);
        for (j = (long)zero; j <= (long)seven; j++) {
          if (numnuc[j] == largest )   /* make set of states that tie */
          {                             /* for most parsimonious here */
            newbase |= 1 << j;
          }
        }
        steps += (weight[i]) * (count_sibs(p) - largest - root);
      } /* above counts descendants that don't have most parsimonious */
      else
      {   /* optimized for bifurcation, code above still works though */
        newbase = ((discretepars_node*)p->next->back)->discbase[i] |
          ((discretepars_node*)p->next->next->back)->discbase[i];
        steps += weight[i];
      }
    }

    ((discretepars_node*)p)->discbase[i] = newbase;
    ((pars_node*)p)->numsteps[i] = steps;
  }
  p->initialized = true;
} /* discretepars_tree_nuview */


double discretepars_tree_evaluate(tree* t, node *n, boolean dummy)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, steps = 0;
  double term;
  double sum = 0.0;
  discretepars_node* p = (discretepars_node*)n;
  discretepars_node* q = (discretepars_node*)n->back;
  unsigned char base1, base2;

  generic_tree_evaluate(t, n, dummy);

  for (i = 0; i < endsite; i++)
  {
    if (p != NULL) {
      steps = ((pars_node*)p)->numsteps[i];
      base1 = p->discbase[i];
    }
    if (q != NULL) {
      steps += ((pars_node*)q)->numsteps[i];
      base2 = q->discbase[i];
    }
    if ((p != NULL) && (q != NULL)) {  /* if neither end of branch is empty */
      if ( ((base1 & base2) == 0) )
        steps += weight[i];       /* a step in this branch if sets disjunct */
    }
    if ( ((pars_tree*)t)->supplement)    /* add steps from polymorphic tips */
      steps += ((pars_tree*)t)->supplement(t, i);

    if (steps <= threshwt[i])                    /* for threshold parsimony */
      term = steps;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser && !reusertree) /* make array of steps */
      fsteps[which - 1][i] = term; /* ... for up to maximum number of trees */
  }
  if (usertree && which <= maxuser && !reusertree)
  {             /* find which tree has fewest steps and how many there were */
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
}  /* discretepars_tree_evaluate */


boolean discretepars_tree_branchcollapsible(tree* t, node* n)
{
  /* is this branch in a discrete characters parsimony tree
   * collapsible without altering the number of steps in the tree? */
  boolean collapsible = true;
  long i;
  node* q;

  if ( (n->tip == true) || (n->back->tip == true) )
    return false;                   /* not if either end of branch is a tip */
  if ( (n->index == t->root->index) || (n->back->index == t->root->index) )
    return false;        /* not if branch is connected to the rootmost fork */

  q = n->back;
  if ( q->initialized == false ) t->nuview(t, q);      /* get updated views */
  if ( n->initialized == false ) t->nuview(t, n);

  for ( i = 0 ; i < endsite ; i++ )
  {                        /* collapse only if all sites do not have a step */
    if ( (((discretepars_node*)q)->discbase[i] &
          ((discretepars_node*)n)->discbase[i] ) == 0)
      return false;
  }
  return collapsible;
} /* discretechars_tree_branchcollapsible */


// End.
