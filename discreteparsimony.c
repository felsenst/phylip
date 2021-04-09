/* Version 4.0.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* This file has functions for non-DNA discrete state parsimony methods */


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

  if (printdata)   /* debug: for discrete states maybe not call them "sequences"? */
    headings(chars, "Sequences", "---------");
  basesread = 0;
  allread = false;
  while (!(allread))
  {                  /* eat white space -- if the separator line has spaces on it */
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
          if (charstate == '\n' || charstate == '\t')  /* a newline or tab  */
            charstate = ' ';                     /* will be seen as a blank */
          if (charstate == ' ')                 /* if it's a blank, move on */
            continue;
          if ((strchr("!\"#$%&'()*+,-./0123456789:;<=>?@\
                       ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`                \
                       abcdefghijklmnopqrstuvwxyz{|}~", charstate)) == NULL)
          {
            printf(
               "\nERROR:  Bad symbol: %c at position %ld of species %ld.\n\n",
                charstate, j+1, i);
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
  checknames(spp);                      /* Check names array for duplicates */
  if (printdata)
  {                                        /* print table of data if needed */
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++)   /* for character blocks */
    {
      for (j = 1; j <= spp; j++)             /* print one species at a time */
      {
        for (k = 0; k < nmlngth; k++)             /* write out species name */
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars)
          l = chars;
        for (k = (i - 1) * 60 + 1; k <= l; k++)    /* print chracter states */
        {         /* doing dot-differencing if needed on subsequent species */
          if (dotdiff && (j > 1 && inputSequences[j - 1][k - 1]
                                    == inputSequences[0][k - 1]))
            charstate = '.';
          else
            charstate = inputSequences[j - 1][k - 1];
          putc(charstate, outfile);                  /* print the state out */
          if (k % 10 == 0 && k % 60 != 0)  /* add blank every 10 characters */
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 1; i <= chars; i++)               /* make for each character ... */
  {           /* a table indicating which symbols are replaced by 0, 1, ... */
    nsymbol = 0;
    for (j = 1; j <= spp; j++)
    {
      if ((nsymbol == 0) && (inputSequences[j - 1][i - 1] != '?'))
      {
        nsymbol = 1;
        convsymboli = 1;
        convtab[0][i-1] = inputSequences[j-1][i-1];
      }
      else if (inputSequences[j - 1][i - 1] != '?')    /* if a symbol there */
      {
        found = false;
        for (k = 1; k <= nsymbol; k++)
        {                          /* check whether it has been seen before */
          if (convtab[k - 1][i - 1] == inputSequences[j - 1][i - 1])
          {
            found = true;
            convsymboli = k;
          }
        }
        if (!found)              /* ... and if it has not been seen yet ... */
        {
          nsymbol++;
          convtab[nsymbol-1][i - 1] = inputSequences[j - 1][i - 1];
          convsymboli = nsymbol;
        }
      }
      if (nsymbol <= 8)          /* put that symbol in the conversion table */
      {
        if (inputSequences[j - 1][i - 1] != '?')
          inputSequences[j - 1][i - 1] = (Char)('0' + (convsymboli - 1));
      }
      else     /* unless are more than 8, in which case crash informatively */
      {
        printf("\nERROR:  More than maximum of 8 symbols in column %ld.\n\n", i);
        exxit(-1);
      }
    }
  }
}  /* inputdata */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order.  The Shell sort is
   * not optimally fast but is pretty good, being about O(n^(4/3)), and
   * is easy to program.  For much agonizing about speed and optimal
   * schemes of gap sizes, see Wikipedia.  The name "Shell" is not a
   * reference to sleight-of-hand but is the name of its inventor */
  /* used in pars */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;     /* start with a "gap" half the length of the array */
  while (gap > 0)
  {
    for (i = gap + 1; i <= chars; i++)    /* compare pairs that far apart */
    {
      j = i - gap;      /* i is the upper, j the lower member of the pair */
      flip = true;
      while (j > 0 && flip)
      {
        jj = alias[j - 1];          /* sort the pair in order of the states */
        jg = alias[j + gap - 1];                    /* ... in their aliases */
        tied = true;
        k = 1;
        while (k <= spp && tied)   /* proceed down column of the data table */
        {
          flip = (inputSequences[k - 1][jj - 1]
                   > inputSequences[k - 1][jg - 1]);
          tied = (tied && (inputSequences[k - 1][jj - 1]
                            == inputSequences[k - 1][jg - 1]));
          k++;                   /* ... until find a species where not tied */
        }
        if (!flip)              /* bail out of the loop if pair is in order */
          break;
        itemp = alias[j - 1];        /* exchange the order in array "alias" */
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = weight[j - 1];           /* ... and in the array of weights */
        weight[j - 1] = weight[j + gap - 1];
        weight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;       /* reduce gap size by half. sorted then it reaches zero */
  }
}  /* sitesort */


void sitecombine(long chars)
{
  /* combine sites that have identical patterns, starting with them sorted */
  /* used in pars */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < chars)                              /* for all characters ... */
  {
    j = i + 1;
    tied = true;
    while (j <= chars && tied)        /* look for a run of tied characters */
    {
      k = 1;
      while (k <= spp && tied)                       /* check whether tied */
      {
        tied = (tied && inputSequences[k - 1][alias[i - 1] - 1] == 
                          inputSequences[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (tied)        /* these two are to be combined, add up weights ... */
      {
        weight[i - 1] += weight[j - 1];
        weight[j - 1] = 0;
        ally[alias[j - 1] - 1] = alias[i - 1]; /* ... and bookkeep aliases */
      }
      j++;
    }
    i = j - 1;               /* go to next untied one and start from there */
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
          itemp = alias[i - 1];                         /* swap aliases ... */
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = weight[i - 1];                         /* ... and weights */
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
{  /* set up states at tip nodes by filling in bits in a byte
    * The "unknown" state has all seven bits filled in */
  long i, j;
  unsigned char ns=0;

  for (j = 0; j < endsite; j++)        /* for all representative characters */
  {
    for (i = 0; i < spp; i++)                               /* for all tips */
    {
      switch (inputSequences[i][alias[j] - 1])         /* look at the state */
      {
        case '0':                /* shift a 1 left that number of positions */
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
      ((discretepars_node*)t->nodep[i])->discbase[j] = ns; /* store the set */
      ((pars_node*)t->nodep[i])->numsteps[j] = 0.0;  /* no steps needed ... */
    }                                              /* ... at or above there */
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


void dischyptrav(tree* t, node *r_, discbaseptr hypset_, long b1, long b2,
                  boolean bottom_)
{
  /*  compute, print out states at one interior node for a range of sites */
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
  if (!Vars.r->tip)     /* if not a tip, step numbers conditional on states */
    memset(((discretepars_node*)Vars.r)->discnumnuc, 0,
             endsite * sizeof(discnucarray));
  for (i = b1 - 1; i < b2; i++)              /* for this range of sites ... */
  {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip)
    {
      p = Vars.r->next;
      for (k = (long)zero; k <= (long)seven; k++)  /* check possible states */
        if (Vars.anc & (1 << k))            /* if this state is not in  anc */ 
          ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]++;
      do {    /* go around fork circle checking for each state if not there */
        for (k = (long)zero; k <= (long)seven; k++)
          if (((discretepars_node*)p->back)->discbase[j - 1] & (1 << k))
            ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest =                         /* find the smallest count of steps */
        discgetlargest(((discretepars_node*)Vars.r)->discnumnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)zero; k <= (long)seven; k++)
      {                                      /* for each possible state ... */
        if (((discretepars_node*)Vars.r)->discnumnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);    /* ... add it into state set if tied */
      }
      ((discretepars_node*)Vars.r)->discbase[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)   /* debug: explain */
      Vars.anc =
     ((discretepars_node*)t->nodep[Vars.r->back->index - 1])->discbase[j - 1];
    Vars.nonzero = (Vars.nonzero ||
                     (((discretepars_node*)Vars.r)->discbase[j - 1]
                      & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || ((discretepars_node*)Vars.r)->discbase[j - 1]
                                   != Vars.anc);
  }
  dischyprint(t, b1, b2, &Vars); /* print the symbol for the possible state */
  Vars.bottom = false;
  if (!Vars.r->tip)
  {
    memcpy(tempnuc, ((discretepars_node*)Vars.r)->discnumnuc,
            endsite * sizeof(discnucarray));
    q = Vars.r->next;
    do {  /* debug: need to comment next part */
      memcpy(((discretepars_node*)Vars.r)->discnumnuc, tempnuc,
              endsite * sizeof(discnucarray));
      for (i = b1 - 1; i < b2; i++)       /* go through range of characters */
      {
        j = location[ally[i] - 1];      /* number of character representing */
        for (k = (long)zero; k <= (long)seven; k++)
          if (((discretepars_node*)q->back)->discbase[j - 1] & (1 << k))
            ((discretepars_node*)Vars.r)->discnumnuc[j - 1][k]--;
        largest = discgetlargest(((discretepars_node*)Vars.r)->
                                 discnumnuc[j - 1]);     /* find most steps */
        ancset[j - 1] = 0;
        for (k = (long)zero; k <= (long)seven; k++)  /* for possible states */
          if (((discretepars_node*)Vars.r)->discnumnuc[j - 1][k] == largest)
            ancset[j - 1] |= (1 << k); /* ... remove if tied for most steps */
        if (!Vars.bottom)
          Vars.anc = ancset[j - 1];
      }
      dischyptrav(t, q->back, ancset, b1, b2, Vars.bottom);   /* recurse on */
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
  if (dotdiff)                         /* if using dot-differencing, say so */
    fprintf(outfile, " ( . means same as in the node below it on tree)\n");
  nothing = (discbaseptr)Malloc(endsite * sizeof(unsigned char));
  for (i = 0; i < endsite; i++)            /* initialize for representative */
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++)    /* loop over char groups */
  {
    putc('\n', outfile);
    n = i * 40;                             /* last character of this group */
    if (n > chars)                    /* if  n  got past last character ... */
      n = chars;
    dischyptrav(t, t->root, nothing, i * 40 - 39, n, true);  /* reconstruct */
  }
  free(nothing);
}  /* hypstates */


void discinitbase(node *p, long sitei)
{
  /* traverse tree to initialize base at internal nodes */
  node *q;
  long i, largest;

  if (p->tip)                                    /* back out if it is a tip */
    return;
  q = p->next;
  while (q != p)                             /* loop over this fork's nodes */
  {
    if (q->back)
    {
      memcpy(((discretepars_node*)q)->discnumnuc,
         ((discretepars_node*)p)->discnumnuc, endsite * sizeof(discnucarray));
      for (i = (long)zero; i <= (long)seven; i++)
      {        /* if state  i  is at  q->back, reduce count of steps needed */
        if (((discretepars_node*)q->back)->discbase[sitei - 1] & (1 << i))
          ((discretepars_node*)q)->discnumnuc[sitei - 1][i]--;
      }
      if (p->back)
      {            /* if there is a node at the other end of the branch ... */
        for (i = (long)zero; i <= (long)seven; i++)
        {      /* if  p->back  has state  i  increase count of steps needed */
          if (((discretepars_node*)p->back)->discbase[sitei - 1] & (1 << i))
            ((discretepars_node*)q)->discnumnuc[sitei - 1][i]++;
        }
      }
      largest = discgetlargest(((discretepars_node*)q)->discnumnuc[sitei - 1]);
      ((discretepars_node*)q)->discbase[sitei - 1] = 0;
      for (i = (long)zero; i <= (long)seven; i++)
      {         /* set up state as those states that require fewest changes */
        if (((discretepars_node*)q)->discnumnuc[sitei - 1][i] == largest)
          ((discretepars_node*)q)->discbase[sitei - 1] |= (1 << i);
      }
    }
    q = q->next;
  }
  q = p->next;    /* go around circle again initializing bases for  q->back */
  while (q != p)
  {
    if (q->back != NULL)
      discinitbase(q->back, sitei);
    q = q->next;
  }
} /* initbase */


void inittreetrav(node *p, long sitei)
{
  /* traverse tree to clear boolean initialized and set up base.
   * unlike function  inittrav  this does not only set inward-looking
   * initialized booleans to false, but all of them at all forks */
  node *q;

  if (p == NULL)                                /* unless it's an empty tip */
    return;
  if (p->tip)
  {
    discinitmin((discretepars_node*)p, sitei, false);   /* initialize state */
 /* debug printf("initialize it in node %ld\n", p->index); */
    p->initialized = true;                       /* mark tip as initialized */
    return;
  } else {
    q = p->next;                                /* for an interior fork ... */
    while (q != p)
    {                            /* go around the fork circle, for each ... */
      if (q != NULL) {
        inittreetrav(q->back, sitei);        /* ... go recursively outwards */
        discinitmin((discretepars_node*)q, sitei, true);
/* debug  printf("initialize it in node %ld\n", q->index); */
        q->initialized = false;
      }
      q = q->next;
    }
    discinitmin((discretepars_node*)p, sitei, true);  /* initializing state */
/* printf("initialize it in node %ld\n", p->index); debug */
    p->initialized = false;    /* ... marking them as needing to be updated */
  }
} /* inittreetrav */


void disccompmin(node *p, node *desc)
{
  /* computes minimum lengths from  p  on beyond it, where we are going around
   * a fork circle and have got to the node in the circle whose back node
   * is  desc */
  long i, j, minn, cost, desclen, descrecon=0, maxx;

  maxx = 10 * spp;          /* a value bigger than number of steps could be */
  for (i = (long)zero; i <= (long)seven; i++) /* for all possible states... */
  {                                           /* ,,, at node  p */
    minn = maxx;
    for (j = (long)zero; j <= (long)seven; j++)        /* for state in desc */
    {
      if (i == j)               /* set cost of change to zero if same state */
        cost = 0;
      else                                       /* otherwise set it to one */
        cost = 1;
      if (((discretepars_node*)desc)->disccumlengths[j] == -1)
      {     /* if  state  i  would not be possible, make it too big a value */
        desclen = maxx;
      }
      else
      {                          /* number of steps from  desc  on outwards */
        desclen = ((discretepars_node*)desc)->disccumlengths[j];
      }
      if (minn > cost + desclen)         /* would it have too many changes? */
      {
        minn = cost + desclen;                 /* set  minn  to lower value */
        descrecon = 0;
      }
      if (minn == cost + desclen)         /* if  j  is a possible state ... */
      {        /* ... then increment the number of possible reconstructions */
        descrecon += ((discretepars_node*)desc->back)->discnumreconst[j];
      }
    }
    ((discretepars_node*)p)->disccumlengths[i] += minn;    /* add to length */
    ((discretepars_node*)p)->discnumreconst[i] *= descrecon;    /* multiply */
/*  debug printf("state %ld: minn, descrecon are %ld, %ld\n", i, minn, descrecon);*/
  }
  p->initialized = true;     /* mark stuff from here on out as already done */
} /* disccompmin */


void minpostorder(node *p, pointarray treenode)
{
  /* traverses an n-ary tree, computing minimum steps at each node */
  node *q;

  if (p == NULL)
    return;
  if (p->tip)
  {
    return;
  }
/* debug printf("minpostorder on node  %ld\n", p->index);  */
  if (!p->initialized) {
    q = p->next;    /* around fork ring, do minpostorder on backs as needed */
    while (q != p)
    {
      if ((q->back) != NULL) {
        if (!(q->back->initialized)) {
          minpostorder(q->back, treenode);
        }
        disccompmin(p, q->back);    /* part of conditional score for branch */
      }
      q = q->next;
    }
    p->initialized = true;  /* set node initialized once one has done those */
  }
}  /* minpostorder */
     

void branchlength(node *subtr1, node *subtr2, double *brlen, pointarray treenode)
{
  /* computes a branch length between two subtrees for a given site */
  long i, j, minn, cost, minsteps, numreconst, nom, denom;

  minpostorder(subtr1, treenode);   /* make sure have steps further out ... */
  minpostorder(subtr2, treenode);              /* ... at each end of branch */
  minn = 10 * spp;                 /* a value that is too big to be correct */
  nom = 0;                       /* the numerator that is being accumulated */
  denom = 0;                                     /* ... and the denominator */
  for (i = (long)zero; i <= (long)seven; i++) /* for all states at both ... */
  {
/* debug printf("for state %ld at %ld, disccumlengths = %ld\n", (long)i, subtr1->index, ((discretepars_node*)subtr1)->disccumlengths[i]); */
/* debug printf("for state %ld at %ld, disccumlengths = %ld\n", (long)i, subtr2->index, ((discretepars_node*)subtr2)->disccumlengths[i]); */
    for (j = (long)zero; j <= (long)seven; j++)       /* ... ends of branch */
    {
      if (i == j)                            /* count 1 for each difference */
        cost = 0;
      else
        cost = 1;
      if (((discretepars_node*)subtr1)->disccumlengths[i] != -1 &&
          (((discretepars_node*)subtr2)->disccumlengths[j] != -1))
      {                                  /* if both states are possible ... */
        minsteps = ((discretepars_node*)subtr1)->disccumlengths[i] + cost +
                    ((discretepars_node*)subtr2)->disccumlengths[j];
        if (minsteps < minn)
        {
          minn = minsteps;     /* update  minn  as smallest possible so far */
          nom = 0;               /* ... and restart accumulating tied steps */
          denom = 0;
        }
        if (minsteps == minn)
        {     /* add to numerator, denominator, weighted by reconstructions */
          numreconst = ((discretepars_node*)subtr1)->discnumreconst[i] *
                         ((discretepars_node*)subtr2)->discnumreconst[j];
          nom += numreconst * cost;
          denom += numreconst;
        }
      }
    }
  }
  if (denom == 0)
    *brlen = 0.0;
  else
    *brlen = (double)nom/(double)denom;    /* cast to doubles so get double */
/* debug printf("branch length for branch %ld:%ld is %6.2f\n", subtr1->index, subtr2->index, *brlen); */
} /* branchlength */


void branchlentrav(node *p, node *root, long sitei, long chars, double *brlen, pointarray treenode)
{
  /*  traverses the tree computing tree length at each branch */
  node *q;

/* debug if (p->back == NULL)
printf("computing branch length for branch %ld:(NULL) for site %ld\n", p->index, sitei); */
/* debug else */
/* debug printf("computing branch length for branch %ld:%ld for site %ld\n", p->index, p->back->index, sitei); */
  if (p->tip)         /* bail if  p  is a tip as already have branch length */
    return;
/* debug:  to avoid recursive infinite loop */
#if 0
  if (p->index == outgrno) /* if outgroup a tip go to nearest interior fork */
    p = p->back;
#endif
  q = p->next;
  do {
    if (q->back != NULL)        /* if this node does not connect to nothing */
    {
/* debug printf("computing branch length for branch %ld:%ld for site %ld\n", q->index, q->back->index, sitei); */
      branchlength(q, q->back, brlen, treenode);  /* get branch length here */
/* debug printf("branch length for site %ld branch %ld:%ld is: %f\n", sitei, q->index, q->back->index, *brlen); */
      q->v += (weight[sitei - 1]  * (*brlen));   /* set node branch lengths */
      q->back->v += (weight[sitei - 1] * (*brlen));
/* debug printf("branch length up to site %ld branch %ld:%ld is: %f\n", sitei, q->index, q->back->index, q->v); */
      if (!q->back->tip)                        /* traverse out of  q->back */
        branchlentrav(q->back, root, sitei, chars, brlen, treenode);
    }
    q = q->next;                              /* move on around fork circle */
  } while (q != p);
  if (p->back != NULL) {
    if (p == root) {
      branchlength(p, p->back, brlen, treenode); /* finally do for outgroup */
      p->v += (weight[sitei - 1]  * (*brlen));   /* set node branch lengths */
      p->back->v += (weight[sitei - 1] * (*brlen));
    }
/* debug printf("branch length up to site %ld branch %ld:%ld is: %f\n", sitei, p->index, p->back->index, q->v); */
  }
}  /* branchlentrav */


void treeout(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in pars */
  node *q;
  long i, n;
  Char c;

  if (p->tip)                                        /* if are at a tip ... */
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)         /* figure out length of tip name */
    {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)   /* write name, making blanks into underscores */
    {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;        /* this keeps track of which column on line we are at */
  }
  else
  {                                               /* if at an interior fork */
    putc('(', outtree);                       /* put a left parenthesis ... */
    (*col)++;
    q = p->next;
    while (q != p)
    {
      treeout(q->back, nextree, col, root);   /* ... then print subtree ... */
      q = q->next;
      if (q == p)                      /* ... if at last furc, no comma ... */
        break;
      putc(',', outtree);              /* ... if not, put a comma there ... */
      (*col)++;
      if (*col > 60)
      {                   /* ... go to a new line if got past column 60 ... */
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);           /* ... close with a right parenthesis ... */
    (*col)++;
  }
  if (p != root)                 /* if not done tree, bail out of recursion */
    return;                                    /* otherwise if are done ... */
  if (nextree > 2)    /* write tree weight after tree if more than one tree */
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

  if (numtrees == 2)                     /* for two species do the KHT test */
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees)                        /* loop for both trees */
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
        temp = sum / sumw;           /* following statement gets SD of mean */
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%9.1f %12.4f",
                (nsteps[which - 1] - minsteps), sd);
        if (sum > 1.95996 * sd)                /* if enough SD's above zero */
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  }
  else
  {                  /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile,
               "Shimodaira-Hasegawa test on first %d of %ld trees\n\n",
               MAXSHIMOTREES, numtrees);
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
            temp += wt*(fsteps[i][k]/(wt)-sum) * (fsteps[j][k]/(wt)-sum2);
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
    free(P);                             /* free the variables we Malloc'ed */
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
/* debug:  maybe not name this "root" since it can be different from the
   root pointer */
  long sitei;
  double trlen;

  initbranchlen(root);      /* debug:  initialize something-or-other ... */
  for (sitei = 0; sitei < endsite; sitei++)     /* for representative chars */
  {
    trlen = 0.0;
    discinitbase(root, sitei);    /* initialize the counts, reconstructions */
printf("initialize site %ld\n", sitei); /* debug */
    inittreetrav(root, sitei);                    /* traverse to initialize */
    inittreetrav(root->back, sitei);                   /* ... both ways out */
    branchlentrav(root, root, sitei, chars, &trlen, treenode);  /* go */
printf("numsteps[%ld]  = %10.6ld, %10.6f\n", sitei, ((pars_node*)root)->numsteps[sitei], trlen);  /* debug */
    ((pars_node*)root)->numsteps[sitei] = trlen;
printf("numsteps[%ld]  = %10.6ld, %10.6f\n", sitei, ((pars_node*)root)->numsteps[sitei], trlen);  /* debug */
  }
} /* disc_treelength */


void discinitmin(discretepars_node *p, long sitei, boolean internal)
{
  long i;

/* debug printf("doing discinitmin on node at %ld, for site %ld, internal = %d\n", ((node*)p)->index, sitei, internal);  */
  if (internal)                              /* for internal fork nodes ... */
  {
    for (i = (long)zero; i <= (long)seven; i++)   /* ... for each state ... */
    {
      p->disccumlengths[i] = 0;     /* initialize lengths beyond there zero */
      p->discnumreconst[i] = 1;     /* initialize number of reconstructions */
    }
  }
  else                                                /* for a tip node ... */
  {
    for (i = (long)zero; i <= (long)seven; i++)   /* ... for each state ... */
    {
      if (p->discbase[sitei - 1] & (1 << i))  /* if that state is there ... */
      {
        p->disccumlengths[i] = 0;   /* ... set length from there on to zero */
        p->discnumreconst[i] = 1;   /* ...  and there is one reconstruction */
      }
      else
      {
        p->disccumlengths[i] = -1;  /* ... if state is not there, signal it */
        p->discnumreconst[i] = 0;    /* ... and no possible reconstructions */
      }
    }
  }
} /* discinitmin */


void dischyprint(tree* t, long b1, long b2, struct LOC_hyptrav *htrav)
{
  /* print out states in sites b1 through b2 at node */
  long i, j, k;
  boolean dot, found;

  if (htrav->bottom)                            /* if at bottom of tree ... */
  {
    if (!outgropt)                       /* ... if there is no outgroup ... */
      fprintf(outfile, "       ");                          /* print blanks */
    else
      fprintf(outfile, "root   ");         /* ... otherwise the word "root" */
  }
  else                                                /* for most nodes ... */
    fprintf(outfile, "%4ld   ", htrav->r->back->index - spp); /* ... number */
  if (htrav->r->tip)                                  /* for tips, the name */
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[htrav->r->index - 1][i], outfile);
  }
  else                                          /* for non-tips, the number */
    fprintf(outfile, "%4ld      ", htrav->r->index - spp);
  if (htrav->bottom)       /* this prints whether or not can be a step here */
    fprintf(outfile, "          ");
  else if (htrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (htrav->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1; i <= b2; i++)               /* for a block of characters ... */
  {
    j = location[ally[i - 1] - 1];
    htrav->tempset = ((discretepars_node*)htrav->r)->discbase[j - 1];
    htrav->anc = htrav->hypset[j - 1];
    if (!htrav->bottom)
      htrav->anc =
       ((discretepars_node*)t->nodep[htrav->r->back->index-1])->discbase[j-1];
    dot = dotdiff && (htrav->tempset == htrav->anc && !htrav->bottom);
    if (dot)                 /* if dot-differencing, print a dot if need to */
      putc('.', outfile);
    else
    {
      found = false;
      k = (long)zero;                        /* look at all possible states */
      do {
        if (htrav->tempset == (1 << k))    /* if there's one possible state */
        {
          putc(convtab[k][i - 1], outfile);    /* write actual state symbol */
          found = true;
        }
        k++;
      } while (!found && k <= (long)seven);    /* see if need to keep going */
      if (!found)
        putc('?', outfile);        /* write ambiguity symbol if appropriate */
    }
    if (i % 10 == 0)                        /* a blank every ten characters */
      putc(' ', outfile);
  }
  putc('\n', outfile);             /* a newline character if at end of line */
}  /* hyprint */


void discretepars_tree_nuview(tree* t, node*p)
{
  /* nuview for multistate discrete characters parsimony methods */
  node *q;
  discretepars_node *qback;
  long i, j, steps, largest;
  unsigned char base1, newbase;
  long numnuc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  boolean bif;
  long atroot = true;

/* debug  generic_tree_nuview(t, p);   needed? */
  bif = (count_sibs(p) == 2);          /* boolean to indicate a bifurcation */

if (p->back != NULL)  /* debug     */
/* debug printf("update states at node %ld facing %ld\n", p->index, p->back->index);  */
/* debug  else  */
/* debug printf("update states at node %ld facing NULL\n", p->index);  */
  for ( i = 0 ; i < endsite ; i++ ) /* do for each representative character */
  {
    newbase = 0xff;
    steps = 0;
    for ( q = p->next ; q != p ; q = q->next )
    {                                    /* go around the loop at this fork */
      qback = (discretepars_node*)q->back;
      if ( qback == NULL )
      {                /* if we root the tree we can safely ignore the root */
        atroot = true;
        continue;
      }
      base1 = qback->discbase[i];
      newbase &= base1;               /* intersection with base sets so far */
      steps += ((pars_node*)qback)->numsteps[i];
/* debug printf("update site: %ld, steps = %ld\n", i, weight[i]*steps); */
    }
    if ( newbase == 0 )
    {
      if ( !bif ) 
      {                                /* case where fork is multifurcating */
        memset(numnuc, 0, sizeof(numnuc));
        for (j = (long)zero; j <= (long)seven; j++) /* for all possible ... */
        {                              /* ... states that might be here ... */
          for ( q = p->next ; q != p ; q = q->next )    /* go around circle */
          {
            qback = (discretepars_node*)q->back;
            if ( qback == NULL )             /* skip if nothing behind node */
              continue;
            if ( qback->discbase[i] & (1 << j) )       /* is bit  j  there? */
              numnuc[j]++;     /* increment count of how many states in set */
          }
        }
        largest = discgetlargest(numnuc);
        for (j = (long)zero; j <= (long)seven; j++) {
          if (numnuc[j] == largest )         /* make set of states that tie */
          {                                   /* for most parsimonious here */
            newbase |= 1 << j;
          }
        }
        steps += (weight[i]) * (count_sibs(p) - largest - root);
      }       /* above counts descendants that don't have most parsimonious */
      else
      {         /* optimized for bifurcation, code above still works though */
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
  /* determines the number of steps needed for a tree. This is
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
    if (p != NULL) {                     /* set up state set at one end ... */
      steps = ((pars_node*)p)->numsteps[i];
      base1 = p->discbase[i];
    }
    if (q != NULL) {                            /* ... and at the other end */
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
    else                              /* truncate it to the threshold value */
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser && !reusertree) /* make array of steps */
      fsteps[which - 1][i] = term; /* ... for up to maximum number of trees */
  }
  if (usertree && (which <= maxuser) && !reusertree)
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
      return false;    /* not collapsible since state sets do not share any */
  }              /* if we finish this loop then the branch has no steps ... */
  return collapsible;                    /* and we report it as collapsible */
} /* discretechars_tree_branchcollapsible */


/* End. */
