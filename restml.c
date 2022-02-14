/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "seq.h"
#include "ml.h"
#include "rest_common.h"

#define initialv        0.1     /* starting value of branch length          */
#define over            60   /* maximum width of a tree on screen */
#define NLR_TEMPS 5

typedef double* sitelike2;      /* Array of sitelength+1 likelihoods */
typedef sitelike2 *phenotype2;  /* Array of length ml_node->endsite */

typedef struct restml_tree _restml_tree;

typedef long (*get_trans_t)(_restml_tree *);
typedef void (*free_trans_t)(_restml_tree *, long);
typedef void (*free_all_trans_t)(void);

typedef struct restml_node {
  ml_node ml_node;      /* Base object, must be first */
  phenotype2 x;
  long branchnum;       /* Index of branch transition matrix entry */
} restml_node;

typedef struct restml_tree {
  struct ml_tree ml_tree;                     /* Base object, must be first */
  transptr trans;
  long *freetrans;
  long transindex;
  transmatrix *lr_temps;
  transmatrix travmatrix;
  get_trans_t get_trans;
  free_trans_t free_trans;
  free_all_trans_t free_all_trans;
} restml_tree;

#ifndef OLDC
/* function prototypes */
void   restml_inputnumbers(void);
void   getoptions(void);
void   allocrest(void);
void   setuppie(void);
void   doinit(void);
void   inputoptions(void);
void   restml_inputdata(void);
void   restml_makevalues(void);
void   getinput(void);
void   copymatrix(transmatrix, transmatrix);
void   maketrans(double, transmatrix, transmatrix, transmatrix);
void   branchtrans(long, double);
void   restml_tree_nuview(tree*t, node *p);
void   restml_tree_makenewv(tree* t, node *p);
void   restml_node_copy(node *c, node *d);
void   restml_buildnewtip(long , tree *);
void   restml_buildsimpletree(tree *);
void   restml_coordinates(node *, double, long *, double *, double *);
void   restml_printree(void);
double sigma(node *, double *);
void   describe(node *);
void   summarize(void);
void   restml_treeout(node *);
void   initialvtrav(restml_tree*, node *);
void   maketree(void);
void   adjust_lengths(tree *);
double adjusted_v(double v);
node * restml_node_new(node_type type, long index);
void   restml_node_init(node* n, node_type, long);
void   restml_node_print(node* n);
void   restml_node_allocx(ml_node* n, long endsite, long sitelength);
sitelike2 init_sitelike(long sitelength);
void   free_sitelike(sitelike2 sl);
void   copy_sitelike(sitelike2 dest, sitelike2 src, long sitelength);
void   reallocsites(void);
void   restml_initx(restml_node * n);
void   restml_init_forkring(restml_node * n);             // RSGbugfix
void   restml_tree_init(tree* t, long nonodes, long spp);
void   restml_node_freex(ml_node* n);

#if 0                                   // RSGbugfix: Never used.
void   restml_tree_re_move(tree* t, node **p, node **q);
#endif

void   restml_tree_smooth(tree* t, node* p);
double restml_tree_evaluate(tree* t, node *p, boolean saveit);
tree*  restml_tree_new(long nonodes, long spp);
void   alloctreetrans(tree *t, long sitelength);
transmatrix alloctrans(long sitelength);

#if 0                                   // RSGbugfix: Never used.
void   restml_tree_copy(restml_tree *a, restml_tree *b);
#endif

long   get_trans(restml_tree* t);
void   free_trans(restml_tree* t, long trans);
void   free_all_trans(restml_tree* t);
void   restml_tree_save_lr_nodes(tree* t, node* p, node* r);
void   restml_tree_restore_lr_nodes(tree* t, node* p, node* r);
void   restml_tree_save_traverses(tree* t, node* p, node* q);
void   restml_tree_restore_traverses(tree* t, node* p, node* q);
void   restml_tree_do_branchl_on_insert(restml_tree* t, node* p, node *q);
void   restml_tree_do_branchl_on_re_move(tree* t, node *p, node *q);
/* function prototypes */
#endif

long rcategs;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
long nonodes2, sites, enzymes, weightsum;
long sitelength;                /* Length of restriction site in bases default 6, overridden in menu. */
long datasets, ith, njumble, jumb = 0;
long inseed, inseed0;
boolean  global, jumble, lngths, weights, trout, trunc8, usertree, progress, reusertree, mulsets, justwts, firstset, improve, smoothit, inserting = false, polishing =false;
tree *curtree, *priortree, *bestree, *bestree2;
longer seed;
long *enterorder;
steptr aliasweight;
char *progname;
node *qwhere, *addwhere;

/* Local variables for maketree, propagated globally for C version: */
long       nextsp, numtrees, maxwhich, col, shimotrees;
double      maxlogl;
boolean     smoothed;
sitelike2    pie;
double      *l0gl;
double     **l0gf;
/* Char ch; */

/* variables added to keep treeread2() happy */
boolean goteof;
double trweight;
boolean haslengths;


long get_trans(restml_tree* t)
{
  long ret;
  ret = t->freetrans[t->transindex];
  t->transindex--;
  return ret;
}


void free_trans(restml_tree* t, long trans)
{
  /* FIXME This is a temporary workaround and probably slows things down a bit.
   * During rearrangements, this function is sometimes called more than once on
   * already freed nodes, causing the freetrans array to overrun other data. */
  long i;
  for ( i = 0 ; i < t->transindex; i++ )
  {
    if ( t->freetrans[i] == trans )
    {
#ifdef TRANS_DEBUG
      printf("\nERROR:  Trans %ld has already been freed!!\n", trans);
#endif
      return;
    }
  }
  /* end of temporary fix */

  t->transindex++;
  t->freetrans[t->transindex] = trans;
}


void free_all_trans(restml_tree* t)
{
  long i;

  for ( i = 0; i < nonodes2; i++ )
    t->freetrans[i] = i;
  t->transindex = nonodes2 - 1;
}


transmatrix alloctrans(long sitelength)
{
  transmatrix ret;
  long i;

  ret =  Malloc((sitelength + 1) * sizeof(double *));
  for (i = 0 ; i < sitelength + 1 ; i++)
    ret[i] = Malloc((sitelength + 1) * sizeof(double));
  return ret;
}


void alloctreetrans(tree *t, long sitelength)
{
  long i;

  ((restml_tree*)t)->trans = (transptr)Malloc(t->nonodes * sizeof(transmatrix));
  for (i = 0; i < t->nonodes; ++i)
    ((restml_tree*)t)->trans[i] = alloctrans(sitelength);
  ((restml_tree*)t)->travmatrix = alloctrans(sitelength);

  ((restml_tree*)t)->freetrans = Malloc(t->nonodes* sizeof(long));
  for ( i = 0; i < t->nonodes; i++ )
    ((restml_tree*)t)->freetrans[i] = i+1;
  ((restml_tree*)t)->transindex = t->nonodes - 1;
}  /* alloctreetrans */


void restml_initx(restml_node * n)
{ /* allocate views on a restml_node */
  long i, j;

  n->branchnum = -1;

  if(n->ml_node.endsite != endsite)
  {
#ifdef BUG_968
    printf("EWIFX -- realloc needed\n");
#endif
    n->ml_node.endsite = endsite;
    // BUG.968 -- free here
    n->x = NULL;
  }
  if(n->x == NULL)
  {
#ifdef BUG_968
    printf("BUG.968 -- realloc needed\n");
#endif
    n->x = (phenotype2)Malloc((1+endsite) * sizeof(sitelike2));
  }

  for ( i = 0 ; i <= endsite ; i++ )
  {
    if(n->x[i] == NULL)
      n->x[i] = init_sitelike(sitelength);
    if ( n->ml_node.node.tip == false )
      for ( j = 0 ; j <= sitelength; j++)
        n->x[i][j] = 1.0;
  }
} /* restml_initx */


void restml_init_forkring(restml_node * n) // RSGbugfix
{ /* initialize views at a restml_node */
#ifdef BUG_968
  printf("BUG.968 -- branch %ld x[0][0] %lf\n", n->branchnum, n->x[0][0]);
#endif
#if 0
  if(n->branchnum != -1)
  {
#endif
    restml_initx(n);
#if 0
  }
#endif
} /* restml_init_forkring */


tree* restml_tree_new(long nonodes, long spp)
{ /* initialize a new restml_tree */
  tree* t = Malloc(sizeof(restml_tree));
  restml_tree_init(t, nonodes, spp);
  return t;
} /* restml_tree_new */


void restml_tree_init(tree* t, long nonodes, long spp) // RSGbugfix
{
  long i, j;

  ml_tree_init(t, nonodes, spp);
  // don't add 1 to endsite here, it's taken care of in restml-specific routines
  allocx(nonodes, endsite, sitelength, (ml_node**)t->nodep);
  alloctreetrans(t, sitelength);
  ((restml_tree*)t)->lr_temps = Malloc(NLR_TEMPS * sizeof(transmatrix));
  for ( i = 0 ; i < NLR_TEMPS ; i++ )
  {
    ((restml_tree*)t)->lr_temps[i] = Malloc((sitelength + 1) * sizeof(double * ));
    for ( j = 0 ; j <= sitelength ; j++)
      ((restml_tree*)t)->lr_temps[i][j] = Malloc((sitelength + 1) * sizeof(double));
  }

  t->evaluate = restml_tree_evaluate;
  t->save_traverses = restml_tree_save_traverses;
  t->restore_traverses = restml_tree_restore_traverses;
  t->nuview = restml_tree_nuview;
  t->save_lr_nodes = restml_tree_save_lr_nodes;
  t->restore_lr_nodes = restml_tree_restore_lr_nodes;
  t->do_branchl_on_insert_f = (do_branchl_on_insert_t)restml_tree_do_branchl_on_insert;
  t->do_branchl_on_re_move_f = (do_branchl_on_re_move_t) restml_tree_do_branchl_on_re_move;
  ((ml_tree*)t)->makenewv = (makenewv_t)restml_tree_makenewv;

  ((restml_tree*)t)->get_trans = get_trans;
  ((restml_tree*)t)->free_trans = free_trans;
} /* restml_tree_init */


sitelike2 init_sitelike(long sitelength)
{ /* allocate site likelihood arrays at a site */
  return Malloc((sitelength+1) * sizeof(double));
} /* init_sitelike */


void free_sitelike(sitelike2 sl)
{ /* free the site likelihood arrays at a site */
  free(sl);
} /* free_sitelike */


void copy_sitelike(sitelike2 dest, sitelike2 src, long sitelength)
{ /* copy the site likelihood values */
  memcpy(dest, src, (sitelength+1) * sizeof(double));
} /* copy_sitelike */


node* restml_node_new(node_type type, long index) // RSGbugfix
{ /* make a new node */
  node* n;

  n = Malloc(sizeof(restml_node));
  restml_node_init(n, type, index);

  return n;
}  /* restml_node_new */


void restml_node_init(node* nn, node_type type, long index)
{ /* initialize the values and functions at a node */
  restml_node *n = (restml_node *)nn;
  n->branchnum = -1;
  ml_node_init(nn, type, index);
  n->ml_node.allocx = restml_node_allocx;
  n->ml_node.node.copy = restml_node_copy;
  n->ml_node.freex = restml_node_freex;
  n->ml_node.node.node_print_f = restml_node_print;
} /* restml_node_init */


void restml_node_print(node * n)
{ /* some sort of debugging code */
  ml_node_print(n);
  restml_node * rn = (restml_node*)n;
  printf(" restml(branchnum:%ld)", rn->branchnum);
} /* restml_node_print */


void restml_inputnumbers(void)
{ /* read and print out numbers of species and sites */
  if(fscanf(infile, "%ld%ld%ld", &spp, &sites, &enzymes) < 3)
  {
    printf("\nERROR reading input file.\n\n");
    exxit(-1);
  }
  nonodes2 = spp * 2 - 1;
}  /* restml_inputnumbers */


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  fprintf(outfile, "\nRestriction site Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n", VERSION);
  putchar('\n');
  sitelength = 6;
  trunc8 = true;
  global = false;
  improve = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  for (;;)
  {
    cleerhome();
    printf("\nRestriction site Maximum Likelihood");
    printf(" method, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n", (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree)
    {
      printf("  N          Use lengths from user trees?  %s\n", (lngths ? "Yes" : "No"));
    }
    printf("  A               Are all sites detected?  %s\n", (trunc8 ? "No" : "Yes"));
    printf("  W                       Sites weighted?  %s\n", (weights ? "Yes" : "No"));
    if (!usertree || reusertree)
    {
      printf("  S        Speedier but rougher analysis?  %s\n", (improve ? "No, not rough" : "Yes"));
      printf("  G                Global rearrangements?  %s\n", (global ? "Yes" : "No"));
      printf("  J   Randomize input order of sequences?  ");
      if (jumble)
        printf("Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("No. Use input order\n");
    }
    printf("  L                          Site length?%3ld\n", sitelength);
    printf("  O                        Outgroup root?  %s%3ld\n", (outgropt ? "Yes, at sequence number" : "No, use as outgroup species"), outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets, (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n", (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n", ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1    Print out the data at start of run  %s\n", (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n", (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n", (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n", (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("UNAWSGJLOTMI01234", ch) != NULL)) || (usertree && (strchr("UNAWSLOTMI01234", ch) != NULL)))
    {
      switch (ch)
      {
        case 'A':
          trunc8 = !trunc8;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'S':
          improve = !improve;
          break;

        case 'G':
          global = !global;
          break;

        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
          break;

        case 'L':
          loopcount2 = 0;
          do {
            printf("New Sitelength?\n");
            if(scanf("%ld%*[^\n]", &sitelength)) {} // Read number and scan to EOL.
            (void)getchar();
            if (sitelength < 1)
              printf("BAD RESTRICTION SITE LENGTH: %ld.\n", sitelength);
            countup(&loopcount2, 10);
          } while (sitelength < 1);
          break;

        case 'N':
          lngths = !lngths;
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          else outgrno = 1;
          break;

        case 'U':
          usertree = !usertree;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets)
          {
            printf("Multiple data sets or multiple weights?");
            loopcount2 = 0;
            do {
              printf(" (type D or W)\n");
              if(scanf("%c%*[^\n]", &ch2)) {} // Read char and scan to EOL.
              (void)getchar();
              if (ch2 == '\n')
                ch2 = ' ';
              uppercase(&ch2);
              countup(&loopcount2, 10);
            } while ((ch2 != 'W') && (ch2 != 'D'));
            justwts = (ch2 == 'W');
            if (justwts)
              justweights(&datasets);
            else
              initdatasets(&datasets);
            if (!jumble)
            {
              jumble = true;
              initjumble(&inseed, &inseed0, seed, &njumble);
            }
          }
          break;

        case 'I':
          interleaved = !interleaved;
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          printdata = !printdata;
          break;

        case '2':
          progress = !progress;
          break;

        case '3':
          treeprint = !treeprint;
          break;

        case '4':
          trout = !trout;
          break;
      }
    }
    else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
}  /* getoptions */


void reallocsites(void)
{
  long i;
  for (i = 0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (Char *)Malloc(sites * sizeof(Char));
  }
  free(weight);
  free(alias);
  free(aliasweight);
  weight = (steptr)Malloc((sites+1) * sizeof(long));
  alias = (steptr)Malloc((sites+1) * sizeof(long));
  aliasweight = (steptr)Malloc((sites+1) * sizeof(long));
} /* reallocsites */


void allocrest(void)
{
  long i;

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *)Malloc(sites * sizeof(Char));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
  weight = (steptr)Malloc((sites+1) * sizeof(long));
  alias = (steptr)Malloc((sites+1) * sizeof(long));
  aliasweight = (steptr)Malloc((sites+1) * sizeof(long));
}  /* allocrest */


void setuppie(void)
{
  /* set up equilibrium probabilities of being a given
     number of bases away from a restriction site */
  long i;
  double sum;

  pie = init_sitelike(sitelength);
  pie[0] = 1.0;
  sum = pie[0];
  for (i = 1; i <= sitelength; i++)
  {
    pie[i] = 3 * pie[i - 1] * (sitelength - i + 1) / i;
    sum += pie[i];
  }
  for (i = 0; i <= sitelength; i++)
    pie[i] /= sum;
}  /* setuppie */


void doinit(void)
{
  /* initializes variables */
  restml_inputnumbers();
  getoptions();
  if (!usertree)
    nonodes2--;
  if (printdata)
    fprintf(outfile, "%4ld Species, %4ld Sites, %4ld Enzymes\n",
            spp, sites, enzymes);
  setuppie();
  allocrest();
}  /* doinit */


void inputoptions(void)
{
  /* read the options information */
  long i, cursp, curst, curenz;

  if (!firstset && !justwts)            // see dnaml.c:inputoptions()
  {
    if (eoln(infile))
      scan_eoln(infile);
    if(fscanf(infile, "%ld%ld%ld", &cursp, &curst, &curenz) < 3)
    {
      printf("\nERROR reading input file.\n\n");
      exxit(-1);
    }
    if (cursp != spp)
    {
      printf("\nERROR:  INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld.\n", ith);
      exxit(-1);
    }
    if (curenz != enzymes)
    {
      printf("\nERROR:  INCONSISTENT NUMBER OF ENZYMES IN DATA SET %4ld.\n", ith);
      exxit(-1);
    }
    sites = curst;
    reallocsites();
  }

  for (i = 1; i <= sites; i++)
    weight[i] = 1;
  weightsum = sites;
  if (justwts || weights)
    inputweights2(1, sites+1, &weightsum, weight, &weights, "RESTML");

  fprintf(outfile, "\n  Recognition sequences all %ld bases long\n", sitelength);
  if (trunc8)
    fprintf(outfile, "\nSites absent from all species are assumed to have been omitted\n\n");
  if (weights && printdata)
    printweights(outfile, 1, sites, weight, "Sites");
}  /* inputoptions */


void restml_inputdata(void)
{
  /* read the species and sites data */
  long i, j, k, l, sitesread, sitesnew=0;
  Char ch;
  boolean allread, done;

  if (printdata)
    putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 39)
    j = 39;
  if (printdata)
  {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Sites\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "-----\n\n");
  }
  sitesread = 0;
  allread = false;
  while (!(allread))
  {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      ch = gettc(infile);
    } while (ch == ' ' || ch == '\t');
    ungetc(ch, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp )
    {
      if ((interleaved && sitesread == 0) || !interleaved)
        initname(i - 1);
      if (interleaved)
        j = sitesread;
      else
        j = 0;
      done = false;
      while (!done && !eoff(infile))
      {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eoff(infile)))
        {
          ch = gettc(infile);
          if (ch == '\n' || ch == '\t')
            ch = ' ';
          if (ch == ' ')
            continue;
          uppercase(&ch);
          if (ch != '1' && ch != '0' && ch != '+' && ch != '-' && ch != '?')
          {
            printf("\nERROR:  Bad symbol %c", ch);
            printf(" at position %ld of species %ld.\n", j+1, i);
            exxit(-1);
          }
          if (ch == '1')
            ch = '+';
          if (ch == '0')
            ch = '-';
          j++;
          inputSequences[i - 1][j - 1] = ch;
        }
        if (interleaved)
          continue;
        if (j < sites)
          scan_eoln(infile);
        else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        sitesnew = j;
      scan_eoln(infile);
      if ((interleaved && j != sitesnew ) || ((!interleaved) && j != sites))
      {
        printf("\nERROR:  SEQUENCES OUT OF ALIGNMENT.\n");
        exxit(-1);}
      i++;
    }
    if (interleaved)
    {
      sitesread = sitesnew;
      allread = (sitesread == sites);
    }
    else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata)
  {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++)
    {
      for (j = 0; j < spp; j++)
      {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > sites)
          l = sites;
        for (k = (i - 1) * 60 + 1; k <= l; k++)
        {
          putc(inputSequences[j][k - 1], outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* restml_inputdata */


void restml_makevalues(void)
{
  /* set up fractional likelihoods at tips */
  long i, j, k, l, m;
  boolean found;

  for (k = 1; k <= endsite; k++)
  {
    j = alias[k];
    for (i = 0; i < spp; i++)
    {
      for (l = 0; l <= sitelength; l++)
        ((restml_node*)curtree->nodep[i])->x[k][l] = 1.0;
      switch (inputSequences[i][j - 1])
      {
        case '+':
          for (m = 1; m <= sitelength; m++)
            ((restml_node*)curtree->nodep[i])->x[k][m] = 0.0;
          break;

        case '-':
          ((restml_node*)curtree->nodep[i])->x[k][0] = 0.0;
          break;

        case '?':
          /* blank case */
          break;
      }
    }
  }
  for (i = 0; i < spp; i++)
  {
    for (k = 1; k <= sitelength; k++)
      ((restml_node*)curtree->nodep[i])->x[0][k] = 1.0;
    ((restml_node*)curtree->nodep[i])->x[0][0] = 0.0;
  }

  if (trunc8)
    return;

  if (enzymes > 1)
  {
    /* Search for site with '-' for all species */
    found = false;
    for ( i = 1; i <= endsite; i++ )
    {
      found = true;
      for ( k = 0; k < spp; k++ )
      {
        if ( inputSequences[k][alias[i] - 1] != '-' )
        {
          found = false;
          break;
        }
      }
      if (found)
        break;
    }

    if (found)
    {
      weightsum += (enzymes - 1) * weight[i];
      weight[i] *= enzymes;
    }
  }
}  /* restml_makevalues */


void getinput(void)
{
  /* reads the input data */
  inputoptions();
  if(!justwts || firstset)
  {
    restml_inputdata();
  }
  rest_makeweights(sites, &endsite);
  inittrees(&curtree, &bestree, &priortree, &bestree2, nonodes2, spp);
  restml_makevalues();
}  /* getinput */


void copymatrix(transmatrix tomat, transmatrix frommat)
{
  /* copy a matrix the size of transition matrix */
  int i, j;

  for (i = 0; i <= sitelength; ++i)
  {
    for (j = 0; j <= sitelength; ++j)
      tomat[i][j] = frommat[i][j];
  }
} /* copymatrix */


void maketrans(double p, transmatrix tempmatrix, transmatrix tempslope, transmatrix tempcurve)
{
  /* make transition matrix: product matrix with change
     probability p. Return the results tempmatrix.
     If both tempslope and tempcurve are non-NULL, slope
     and curvature matrices are returned as well.  */
  long i, j, k, m1, m2;
  double sump, sums=0, sumc=0, pover3, pijk, term;
  sitelike2 binom1, binom2;
  boolean nr = false;

  assert( tempmatrix != NULL );

  if ( tempslope != NULL && tempcurve != NULL )
  {
    nr = true;
  }

  binom1 = init_sitelike(sitelength);
  binom2 = init_sitelike(sitelength);
  pover3 = p / 3;
  for (i = 0; i <= sitelength; i++)
  {
    if (p > 1.0 - epsilon)
      p = 1.0 - epsilon;
    if (p < epsilon)
      p = epsilon;
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= sitelength - i; k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; ++j)
    {
      sump = 0.0;
      if (nr)
      {
        sums = 0.0;
        sumc = 0.0;
      }
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++)
      {
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        if (nr)
        {
          term = (j-i+2*k)/p - (sitelength-j-k)/(1.0-p) - (i-k)/(3.0-p);
          sums += pijk * term;
          sumc += pijk * (term * term - (j-i+2*k)/(p*p) - (sitelength-j-k)/((1.0-p)*(1.0-p)) - (i-k)/((3.0-p)*(3.0-p)) );
        }
      }
      tempmatrix[i][j] = sump;
      if (nr)
      {
        tempslope[i][j] = sums;
        tempcurve[i][j] = sumc;
      }
    }
  }
  free_sitelike(binom1);
  free_sitelike(binom2);
}  /* maketrans */


void branchtrans(long i, double p)
{
  /* make branch transition matrix for branch i with probability of change p */
  transmatrix tm;

  tm = ((restml_tree*)curtree)->trans[i - 1];
  assert( tm != NULL );                                    // CrashPoint ASSERT
  maketrans(p, tm, NULL, NULL);
}  /* branchtrans */


double restml_tree_evaluate(tree *tr, node *p, boolean saveit)
{
  /* evaluates the likelihood, using info. at one branch */
  double sum, sum2, y, liketerm, like0, lnlike0=0, term;
  long i, j, k;
  node *q;
  sitelike2 x1, x2;
  static transmatrix tempmatrix = NULL;

  if ( tempmatrix == NULL )
    tempmatrix = alloctrans(sitelength);

  /* Update views */
  generic_tree_evaluate(tr, p, saveit);

  x1 = init_sitelike(sitelength);
  x2 = init_sitelike(sitelength);
  sum = 0.0;
  q = p->back;
  y = p->v;
  maketrans(y, tempmatrix, NULL, NULL);
  copy_sitelike(x1, ((restml_node*)p)->x[0], sitelength);
  copy_sitelike(x2, ((restml_node*)q)->x[0], sitelength);
  if (trunc8)
  {
    like0 = 0.0;
    for (j = 0; j <= sitelength; j++)
    {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
        like0 += liketerm * tempmatrix[j][k] * x2[k];
    }
    lnlike0 = log(enzymes * (1.0 - like0));
  }
  for (i = 1; i <= endsite; i++)
  {
    copy_sitelike(x1, ((restml_node*)p)->x[i], sitelength);
    copy_sitelike(x2, ((restml_node*)q)->x[i], sitelength);
    sum2 = 0.0;
    for (j = 0; j <= sitelength; j++)
    {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
      {
        sum2 += liketerm * tempmatrix[j][k] * x2[k];
      }
    }
    term = log(sum2);
    if (trunc8)
      term -= lnlike0;
    if (usertree && (which <= shimotrees))
      l0gf[which - 1][i - 1] = term;
    sum += weight[i] * term;
  }

  /* *** debug  put a variable "saveit" in evaluate as third argument as to whether to save the KHT suff, and start using it */
  if (usertree)
  {
    if(which <= shimotrees)
      l0gl[which - 1] = sum;
    if (which == 1)
    {
      maxwhich = 1;
      maxlogl = sum;
    }
    else if (sum > maxlogl)
    {
      maxwhich = which;
      maxlogl = sum;
    }
  }

  /* Log likelihoods should be less than 1.0 (IDR) */
  assert( sum <= 1.0 );
  tr->score = sum;
  free_sitelike(x1);
  free_sitelike(x2);

  return sum;
}  /* restml_tree_evaluate */


void restml_tree_nuview(tree *t, node *p)
{
  /* recompute fractional likelihoods for one part of tree */
  long i, j, k, lowlim;
  double sumq, sumr;
  node *q, *r;
  sitelike2 xq, xr, xp;
  transmatrix tempq, tempr;
  double *tq, *tr;

  generic_tree_nuview(t, p);

  xq = init_sitelike(sitelength);
  xr = init_sitelike(sitelength);
  xp = init_sitelike(sitelength);
  tempq = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  tempr = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0;i<=sitelength;++i)
  {
    tempq[i] = (double *)Malloc((sitelength+1) * sizeof (double));
    tempr[i] = (double *)Malloc((sitelength+1) * sizeof (double));
  }
  if (trunc8)
    lowlim = 0;
  else
    lowlim = 1;

  q = p->next->back;
  r = p->next->next->back;
  copymatrix(tempq, ((restml_tree*)curtree)->trans[((restml_node*)q)->branchnum - 1]);
  copymatrix(tempr, ((restml_tree*)curtree)->trans[((restml_node*)r)->branchnum - 1]);

  for (i = lowlim; i <= endsite; i++)
  {
    copy_sitelike(xq, ((restml_node*)q)->x[i], sitelength);
    copy_sitelike(xr, ((restml_node*)r)->x[i], sitelength);
    for (j = 0; j <= sitelength; j++)
    {
      sumq = 0.0;
      sumr = 0.0;
      tq = tempq[j];
      tr = tempr[j];
      for (k = 0; k <= sitelength; k++)
      {
        sumq += tq[k] * xq[k];
        sumr += tr[k] * xr[k];
      }
      xp[j] = sumq * sumr;
    }
    copy_sitelike(((restml_node*)p)->x[i], xp, sitelength);
  }
  for (i=0;i<=sitelength;++i)
  {
    free(tempq[i]);
    free(tempr[i]);
  }

  free(tempq);
  free(tempr);
  p->initialized = true;
  free_sitelike(xq);
  free_sitelike(xr);
  free_sitelike(xp);
}  /* restml_tree_nuview */


void restml_tree_makenewv(tree* t, node *p)
{
  /* Newton-Raphson algorithm improvement of a branch length */
  long i, j, k, lowlim, it, ite;
  double sum, sums, sumc, like, slope, curve, liketerm, liket, y, yold=0, yorig, like0=0, slope0=0, curve0=0, oldlike=0, temp;
  boolean done, firsttime, better;
  node *q;
  sitelike2 xx1, xx2;
  double *tm, *ts, *tc; /* pointers to rows in matrices */

  /* These will be allocated only once for efficiency */
  static transmatrix tempmatrix = NULL;
  static transmatrix tempslope = NULL;
  static transmatrix tempcurve = NULL;

  if ( tempcurve == NULL )
  {
    tempmatrix = alloctrans(sitelength);
    tempslope = alloctrans(sitelength);
    tempcurve = alloctrans(sitelength);
  }

  xx1 = init_sitelike(sitelength);
  xx2 = init_sitelike(sitelength);
  q = p->back;
  y = p->v;
  yorig = y;
  if (trunc8)
    lowlim = 0;
  else
    lowlim = 1;
  done = false;
  firsttime = true;
  it = 1;
  ite = 0;
  while ((it < iterations) && (ite < 20) && (!done))
  {
    like = 0.0;
    slope = 0.0;
    curve = 0.0;
    maketrans(y, tempmatrix, tempslope, tempcurve);

    for (i = lowlim; i <= endsite; i++)
    {
      copy_sitelike(xx1, ((restml_node*)p)->x[i], sitelength);
      copy_sitelike(xx2, ((restml_node*)q)->x[i], sitelength);
      sum = 0.0;
      sums = 0.0;
      sumc = 0.0;
      for (j = 0; j <= sitelength; j++)
      {
        liket = xx1[j] * pie[j];
        tm = tempmatrix[j];
        ts = tempslope[j];
        tc = tempcurve[j];
        for (k = 0; k <= sitelength; k++)
        {
          liketerm = liket * xx2[k];
          sum += tm[k] * liketerm;
          sums += ts[k] * liketerm;
          sumc += tc[k] * liketerm;
        }
      }
      if (i == 0)
      {
        like0 = sum;
        slope0 = sums;
        curve0 = sumc;
      }
      else
      {
        like += weight[i] * log(sum);
        slope += weight[i] * sums/sum;
        temp = sums/sum;
        curve += weight[i] * (sumc/sum-temp*temp);
      }
    }
    if (trunc8 && fabs(like0 - 1.0) > 1.0e-10)
    {
      like -= weightsum * log(enzymes * (1.0 - like0));
      slope += weightsum * slope0 /(1.0 - like0);
      curve += weightsum * (curve0 /(1.0 - like0) + slope0*slope0/((1.0 - like0)*(1.0 - like0)));
    }
    better = false;
    if (firsttime)
    {
      yold = y;
      oldlike = like;
      firsttime = false;
      better = true;
    }
    else
    {
      if (like > oldlike)
      {
        yold = y;
        oldlike = like;
        better = true;
        it++;
      }
    }
    if (better)
    {
      y = y + slope/fabs(curve);
      if (y < epsilon)
        y = 10.0 * epsilon;
      if (y > 0.75)
        y = 0.75;
    }
    else
    {
      if (fabs(y - yold) < epsilon)
        ite = 20;
      y = (y + yold) / 2.0;
    }
    ite++;
    done = fabs(y-yold) < epsilon;
  }
  smoothed = (fabs(yold-yorig) < epsilon) && (yorig > 1000.0*epsilon);
  p->v = yold;
  q->v = yold;
  branchtrans(((restml_node*)p)->branchnum, yold);
  t->score = oldlike;
  free_sitelike(xx1);
  free_sitelike(xx2);
}  /* restml_tree_makenewv */


void restml_tree_do_branchl_on_insert(restml_tree* t, node* forknode, node *where)
{
  long i;
  long j;

  if (where->v >= 0.75)
    where->v = 0.75;
  else
    where->v = 0.75 * (1 - sqrt(1 - 1.333333 * where->v));
  if ( where->v < epsilon )
    where->v = epsilon;

  where->back->v = where->v;
  forknode->next->next->v = where->v;
  forknode->next->next->back->v = forknode->next->next->v;

  for (i = 0; i <= endsite; i++)
  {
    for (j = 0; j <= sitelength; j++)
    {
      ((restml_node*)where->back)->x[i][j] = 1.0;
      ((restml_node*)where->back->next)->x[i][j] = 1.0;
      ((restml_node*)where->back->next->next)->x[i][j] = 1.0;
    }
  }

  // BUG.968 -- branchnum overwritten -- is it freed anywhere?
  ((restml_node*)where->back)->branchnum = ((restml_node*)where)->branchnum;
  ((restml_node*)forknode->next->next->back)->branchnum = ((restml_node*)forknode->next->next)->branchnum = t->get_trans(t);

  ((restml_node*)forknode)->branchnum = ((restml_node*)forknode->back)->branchnum = t->get_trans(t);
  branchtrans(((restml_node*)where)->branchnum, where->v);
  branchtrans(((restml_node*)forknode->next->next)->branchnum, forknode->next->next->v);
  branchtrans(((restml_node*)forknode)->branchnum, forknode->v);
}  /* restml_tree_do_branchl_on_insert */


void restml_tree_do_branchl_on_re_move(tree* t, node *p, node *q)
// does not call ml_tree_do_branchl_on_re_move
{
  free_trans((restml_tree*)t, ((restml_node*)q->back)->branchnum);
  free_trans(((restml_tree*)t), ((restml_node*)p)->branchnum);
  ((restml_node*)q->back)->branchnum = ((restml_node*)q)->branchnum;
  q->v =  0.75*(1 - (1 - 1.333333 * q->v) * (1 - 1.333333*p->v));

  if ( q->v > 1 - epsilon)
    q->v = 1 - epsilon;
  else if ( q->v < epsilon)
    q->v = epsilon;

  q->back->v = q->v;
  branchtrans(((restml_node*)q)->branchnum, q->v);

}  /* restml_tree_do_branchl_on_re_move */


void restml_node_copy(node *srcn, node *dstn)
{
  /* copy a node */
  restml_node *src = (restml_node *)srcn;
  restml_node *dst = (restml_node *)dstn;
  int i;
  long oldendsite = dst->ml_node.endsite;

  ml_node_copy(srcn, dstn);
  dst->branchnum = src->branchnum;

  /* If src->endsite differs from dst->endsite, free dst->x[] */
  if ( oldendsite != 0 && oldendsite != dst->ml_node.endsite)
  {
    for ( i = 0 ; i <= ((ml_node*)src)->endsite ; i++ )
      free(dst->x[i]);
    free(dst->x);
    oldendsite = 0;
  }
  /* If dst phenotype was just freed or never allocated, allocate one
   * the same size as src phenotype */
  if ( oldendsite == 0 )
  {
    dst->x = (phenotype2)Malloc(((((ml_node*)src)->endsite)+1) * sizeof(sitelike2));
    for ( i = 0 ; i <= ((ml_node*)src)->endsite ; i++ )
      dst->x[i] = init_sitelike(sitelength);
  }
  /* Now copy src phenotype to dst */
  for ( i = 0 ; i <= endsite ; i++ )
    copy_sitelike(dst->x[i], src->x[i], sitelength);
}  /* restml_node_copy */


#if 0                                   // RSGbugfix: Never used.
//
void restml_tree_copy(restml_tree *a, restml_tree *b)
{
  long i;

  /* Copy the entire tree */
  generic_tree_copy((tree*)a, (tree*)b);
  /* Copy the transition matrices */
  for ( i = 0 ; i<((tree*)a)->nonodes; i++ )
    copymatrix(b->trans[i], a->trans[i]);
}  /* restml_tree_copy */
//
#endif


void restml_buildsimpletree(tree *tr)
{
  /* set up and adjust branch lengths of a three-species tree */
  hookup(tr->nodep[enterorder[0] - 1], tr->nodep[enterorder[1] - 1]);
  tr->nodep[enterorder[0] - 1]->v = initialv;
  tr->nodep[enterorder[1] - 1]->v = initialv;
  branchtrans(enterorder[1], initialv);
  ((restml_node*)tr->nodep[enterorder[0] - 1])->branchnum = ((restml_node*)tr->nodep[enterorder[1] - 1])->branchnum = get_trans((restml_tree*)tr);
  tr->insert_(tr, tr->nodep[enterorder[2] - 1], tr->nodep[enterorder[1] - 1], false, false);
  tr->root = tr->nodep[enterorder[2]-1]->back;
}  /* restml_buildsimpletree */


void restml_tree_save_lr_nodes(tree *t, node *p, node *r)
{
  restml_tree *rt;

  unrooted_tree_save_lr_nodes((tree*)t, p, r);

  rt = (restml_tree *)t;
  copymatrix(rt->lr_temps[0], rt->trans[((restml_node *)r->back)->branchnum - 1]);
  copymatrix(rt->lr_temps[1], rt->trans[((restml_node*)r->back->next)->branchnum - 1]);
  copymatrix(rt->lr_temps[2], rt->trans[((restml_node*)r->back->next->next)->branchnum - 1]);
  copymatrix(rt->lr_temps[3], rt->trans[((restml_node*)p->next)->branchnum - 1]);
  copymatrix(rt->lr_temps[4], rt->trans [((restml_node*)p->next->next)->branchnum - 1]);
}


//void restml_tree_restore_lr_nodes(tree* t, node* p, node* r, node * restoredFork)

void restml_tree_restore_lr_nodes(tree* t, node* p, node* r)
{
  restml_tree *rt;

  unrooted_tree_restore_lr_nodes((tree*)t, p, r);

#if 0 // BUG.968
  restoredFork->branchnum = r->branchnum;
  restoredFork->next->branchnum = restoredFork->next->back->branchnum;
  restoredFork->next->next->branchnum = restoredFork->next->next->back->branchnum;
#endif

  // BUG.968 -- this relies on the re_move code in locrearrange_recurs
  // not over-writing the node info in r->back
  ((restml_node*)r)->branchnum = ((restml_node*)(r->back))->branchnum;
  ((restml_node*)r->back->next->back)->branchnum = ((restml_node*)(r->back->next))->branchnum;
  ((restml_node*)r->back->next->next->back)->branchnum = ((restml_node*)(r->back->next->next))->branchnum;

  rt = (restml_tree *)t;
  copymatrix(rt->trans[((restml_node*)r->back)->branchnum - 1], rt->lr_temps[0]);
  copymatrix(rt->trans[((restml_node*)r->back->next)->branchnum - 1], rt->lr_temps[1]);
  copymatrix(rt->trans[((restml_node*)r->back->next->next)->branchnum - 1], rt->lr_temps[2]);
  copymatrix(rt->trans[((restml_node*)p->next)->branchnum - 1], rt->lr_temps[3]);
  copymatrix(rt->trans[((restml_node*)p->next->next)->branchnum - 1], rt->lr_temps[4]);
}


void restml_coordinates(node *p, double lengthsum, long *tipy, double *tipmax, double *x)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

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
    (*x) = -0.75 * log(1.0 - 1.333333 * q->v);
    restml_coordinates(q->back, lengthsum + (*x), tipy, tipmax, x);
    q = q->next;
  } while ((p == curtree->root || p != q) && (p != curtree->root || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == curtree->root)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* restml_coordinates */


void restml_printree(void)
{
  /* prints out diagram of the tree */
  long tipy, i;
  double scale, tipmax, x;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  restml_coordinates(curtree->root, 0.0, &tipy, &tipmax, &x);
  scale = 1.0 / (tipmax + 1.000);
  for (i = 1; i <= tipy - down; i++)
    drawline2(i, scale, curtree);
  putc('\n', outfile);
}  /* restml_printree */


void restml_node_freex(ml_node* nn)
{
  restml_node *n = (restml_node *)nn;
  long i;

  for ( i = 0 ; i <= n->ml_node.endsite; i++ )
  {
    free(n->x[i]);
  }
  free(n->x);
}


void restml_node_allocx(ml_node* nn, long endsite, long sitelength)
{
  restml_node *n = (restml_node *)nn;
  long i, j;

  n->x = (phenotype2)Malloc((1+endsite) * sizeof(sitelike2));
  for ( i = 0 ; i <= endsite ; i++ )
  {
    n->x[i] = init_sitelike(sitelength);
    if ( n->ml_node.node.tip == false )
      for ( j = 0 ; j <= sitelength; j++)
        n->x[i][j] = 1.0;
  }
  n->ml_node.endsite = endsite;
}


void freex2(long nonodes, pointarray treenode)
{
  /* used in restml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++)
    free(((restml_node*)treenode[i])->x);
  for (i = spp; i < nonodes; i++)
  {
    p = treenode[i];
    for (j = 1; j <= 3; j++)
    {
      free(((restml_node*)p)->x);
      p = p->next;
    }
  }
}  /* freex2 */


double sigma(node *q, double *sumlr)
{
  /* get 1.95996 * approximate standard error of branch length */
  double sump, sumr, sums, sumc, p, pover3, pijk, Qjk, liketerm, f;
  double  slopef, curvef;
  long i, j, k, m1, m2;
  sitelike2 binom1, binom2;
  transmatrix Prob, slopeP, curveP;
  node *r;
  sitelike2 x1, x2;
  double term, TEMP;

  x1 = init_sitelike(sitelength);
  x2 = init_sitelike(sitelength);
  binom1 = init_sitelike(sitelength);
  binom2 = init_sitelike(sitelength);
  Prob   = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  slopeP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  curveP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0; i<=sitelength; ++i)
  {
    Prob[i]   = (double *)Malloc((sitelength+1) * sizeof(double));
    slopeP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
    curveP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
  }
  p = q->v;
  pover3 = p / 3;
  for (i = 0; i <= sitelength; i++)
  {
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= (sitelength - i); k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; j++)
    {
      sump = 0.0;
      sums = 0.0;
      sumc = 0.0;
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++)
      {
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        term = (j-i+2*k)/p - (sitelength-j-k)/(1.0-p) - (i-k)/(3.0-p);
        sums += pijk * term;
        sumc += pijk * (term * term
                        - (j-i+2*k)/(p*p)
                        - (sitelength-j-k)/((1.0-p)*(1.0-p))
                        - (i-k)/((3.0-p)*(3.0-p)) );
      }
      Prob[i][j] = sump;
      slopeP[i][j] = sums;
      curveP[i][j] = sumc;
    }
  }
  (*sumlr) = 0.0;
  sumc = 0.0;
  sums = 0.0;
  r = q->back;
  for (i = 1; i <= endsite; i++)
  {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    copy_sitelike(x1, ((restml_node*)q)->x[i], sitelength);
    copy_sitelike(x2, ((restml_node*)r)->x[i], sitelength);
    for (j = 0; j <= sitelength; j++)
    {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++)
      {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weight[i] * log(f / sumr);
    sums += weight[i] * slopef / f;
    TEMP = slopef / f;
    sumc += weight[i] * (curvef / f - TEMP * TEMP);
  }
  if (trunc8)
  {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    copy_sitelike(x1, ((restml_node*)q)->x[0], sitelength);
    copy_sitelike(x2, ((restml_node*)r)->x[0], sitelength);
    for (j = 0; j <= sitelength; j++)
    {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++)
      {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weightsum * log((1.0 - sumr) / (1.0 - f));
    sums += weightsum * slopef / (1.0 - f);
    TEMP = slopef / (1.0 - f);
    sumc += weightsum * (curvef / (1.0 - f) + TEMP * TEMP);
  }
  for (i=0;i<=sitelength;++i)
  {
    free(Prob[i]);
    free(slopeP[i]);
    free(curveP[i]);
  }
  free(Prob);
  free(slopeP);
  free(curveP);
  free_sitelike(x1);
  free_sitelike(x2);
  free_sitelike(binom1);
  free_sitelike(binom2);
  if (sumc < -1.0e-6)
    return ((-sums - sqrt(sums * sums - 3.841 * sumc)) / sumc);
  else
    return -1.0;
}  /* sigma */


void describe(node *p)
{
  /* print out information on one branch */
  double sumlr;
  long i;
  node *q;
  double s;

  q = p->back;
  fprintf(outfile, "%4ld      ", q->index - spp);
  fprintf(outfile, "    ");
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  if (q->v >= 0.75)
    fprintf(outfile, "     infinity");
  else
    fprintf(outfile, "%13.5f", -0.75 * log(1 - 1.333333 * q->v));
  if (p->iter)
  {
    s = sigma(q, &sumlr);
    if (s < 0.0)
      fprintf(outfile, "     (     zero,    infinity)");
    else
    {
      fprintf(outfile, "     (");
      if (q->v - s <= 0.0)
        fprintf(outfile, "     zero");
      else
        fprintf(outfile, "%9.5f", -0.75 * log(1 - 1.333333 * (q->v - s)));
      putc(',', outfile);
      if (q->v + s >= 0.75)
        fprintf(outfile, "    infinity");
      else
        fprintf(outfile, "%12.5f", -0.75 * log(1 - 1.333333 * (q->v + s)));
      putc(')', outfile);
    }
    if (sumlr > 1.9205)
      fprintf(outfile, " *");
    if (sumlr > 2.995)
      putc('*', outfile);
  }
  else
    fprintf(outfile, "            (not varied)");
  putc('\n', outfile);
  if (!p->tip)
  {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */


void summarize(void)
{
  /* print out information on branches of tree */

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n\n", curtree->score);
  fprintf(outfile, " \n");
  fprintf(outfile, " Between        And            Length");
  fprintf(outfile, "      Approx. Confidence Limits\n");
  fprintf(outfile, " -------        ---            ------");
  fprintf(outfile, "      ------- ---------- ------\n");
  describe(curtree->root->next->back);
  describe(curtree->root->next->next->back);
  describe(curtree->root->back);
  fprintf(outfile, "\n     *  = significantly positive, P < 0.05\n");
  fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n\n");
}  /* summarize */


void restml_treeout(node *p)
{
  /* write out file with representation of final tree */
  long i, n, w;
  Char c;
  double x;

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
    col += n;
  }
  else
  {
    putc('(', outtree);
    col++;
    restml_treeout(p->next->back);
    putc(',', outtree);
    col++;
    if (col > 45)
    {
      putc('\n', outtree);
      col = 0;
    }
    restml_treeout(p->next->next->back);
    if (p == curtree->root)
    {
      putc(',', outtree);
      col++;
      if (col > 45)
      {
        putc('\n', outtree);
        col = 0;
      }
      restml_treeout(p->back);
    }
    putc(')', outtree);
    col++;
  }
  if (p->v >= 0.75)
    x = -1.0;
  else
    x = -0.75 * log(1 - 1.333333 * p->v);
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree->root)
    fprintf(outtree, ";\n");
  else
  {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* restml_treeout */


/* restml_tree_save_traverses
 *
 * Save node q in t->temp, including the branch number and
 * a copy of the transmatrix
 */
void restml_tree_save_traverses(tree* t, node * p, node* q)
{
  long branchnum;

  restml_tree *rt = (restml_tree *)t;
  restml_node *rq = (restml_node *)q;

  generic_tree_save_traverses(t, p, q);
  branchnum = ((restml_node *)t->temp_q)->branchnum = rq->branchnum;
  copymatrix (rt->travmatrix, rt->trans[branchnum-1]);
}


void restml_tree_restore_traverses(tree *t, node * p, node *q)
{
  restml_tree *rt = (restml_tree *)t;
  restml_node *rq = (restml_node *)q;

  /* BUG.968 -- check that this should even happen at all */
  if(! q->back->tip)
  {
    rt->free_trans(rt, ((restml_node*)q->back->next)->branchnum);
    rt->free_trans(rt, ((restml_node*)q->back->next->next)->branchnum);
  }

  generic_tree_restore_traverses(t, p, q);
  rq->branchnum = ((restml_node*)t->temp_q)->branchnum;
  ((restml_node*)q->back)->branchnum = rq->branchnum;
  copymatrix (rt->trans[rq->branchnum-1], rt->travmatrix);
}


void initialvtrav(restml_tree* t, node *p)
{
  /* traverse tree to set initialized and v to initial values */
  node *q;

  if (((restml_node*)p)->branchnum == 0)
  {
    ((restml_node*)p)->branchnum = t->get_trans(t);
    ((restml_node*)p->back)->branchnum = ((restml_node*)p)->branchnum;
  }

  p->initialized = false;
  p->back->initialized = false;

  if ((!lngths) || p->iter)
  {
    branchtrans(((restml_node*)p)->branchnum, initialv);
    p->v = initialv;
    p->back->v = initialv;
  }
  else
  {
    branchtrans(((restml_node*)p)->branchnum, p->v);
  }

  if (!p->tip)
  {
    q = p->next;
    while ( q != p )
    {
      initialvtrav(t, q->back);
      q = q->next;
    }
  }
} /* initialvtrav */


double adjusted_v(double v)
{
  return 3.0/4.0 * (1.0-exp(-4.0/3.0 * v));
}


void adjust_lengths(tree* t)
{
  long i;
  for ( i = 0 ; i < spp ; i++)
  {
    t->nodep[i]->v = adjusted_v(t->nodep[i]->v);
  }
  for ( i = spp ; i < nonodes2 ; i++)
  {
    t->nodep[i]->v = adjusted_v(t->nodep[i]->v);
    t->nodep[i]->next->v = adjusted_v(t->nodep[i]->next->v);
    t->nodep[i]->next->next->v = adjusted_v(t->nodep[i]->next->next->v);
  }
}


void maketree(void)
{
  /* construct and rearrange tree */
  long i, j;
  double bestyet;

  if (usertree)
  {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *) Malloc(shimotrees * sizeof(double));
    l0gf = (double **) Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *)Malloc((1+endsite) * sizeof(double));
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;
    while (which <= numtrees)
    {
      preparetree(curtree);
      treeread2 (curtree, intree, &curtree->root, lngths, &trweight,
                  &goteof, &haslengths, &spp, false, nonodes2);
      fixtree(curtree);
      if ( outgropt )
        curtree->root = curtree->nodep[outgrno - 1]->back;
      ml_treevaluate(curtree, improve, reusertree, global, progress, priortree, bestree, (initialvtrav_t)initialvtrav );
      if ( reusertree && ( which == 1 || curtree->score > bestree2->score ))
      {
        curtree->copy(curtree, bestree2);
      }
      if ( reusertree && which == numtrees )
      {
        bestree2->copy(bestree2, curtree);
        curtree->root = curtree->nodep[0]->back;
        if ( outgropt )
          curtree->root = curtree->nodep[outgrno - 1]->back;
      }
      else if ( reusertree )
        continue;
      restml_printree();
      if (treeprint)
        summarize();
      if (trout)
      {
        col = 0;
        restml_treeout(curtree->root);
      }
      which++;
    }
    FClose(intree);
    if (numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite-1, maxlogl, l0gl, l0gf, aliasweight, seed);
  }
  else
  {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    if (progress)
    {
      printf("\nAdding species:\n");
      writename(0, 3, enterorder);
    }
    nextsp = 3;
    restml_buildsimpletree(curtree);
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    smoothit = improve;
    nextsp = 4;
    while (nextsp <= spp)
    {
      bestyet = - nextsp*sites*sitelength*log(4.0);
      if (smoothit)
        curtree->copy(curtree, priortree);
      curtree->addtraverse(curtree, curtree->nodep[enterorder[nextsp-1]-1],
                             curtree->root, further, &qwhere, &bestyet,
                             bestree, priortree, smoothit, NULL);
      if (smoothit)
        bestree->copy(bestree, curtree);
      else
      {
        curtree->insert_(curtree, curtree->nodep[enterorder[nextsp - 1] - 1], qwhere, false, false);
        curtree->smoothall(curtree, curtree->root);
        bestyet = curtree->score;
      }

      if (progress)
        writename(nextsp - 1, 1, enterorder);
      if (global && nextsp == spp)
      {
        if (progress)
        {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = spp ; j < nonodes2 ; j++)
            if ( (j - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
              putchar('-');
          putchar('!');
        }
        putchar('\n');
      }

      if (global && nextsp == spp)
        curtree->globrearrange(curtree, progress, smoothit);
      else
        curtree->locrearrange(curtree, curtree->nodep[enterorder[0] - 1], smoothit, priortree, bestree);

      nextsp++;
    }
    if (global && progress)
    {
      putchar('\n');
      fflush(stdout);
    }
    curtree->copy(curtree, bestree);
    if (njumble > 1)
    {
      if (jumb == 1)
        bestree->copy(bestree, bestree2);
      else
      {
        if (bestree2->score < bestree->score)
          bestree->copy(bestree, bestree2);
      }
    }
    if (jumb == njumble)
    {
      if (njumble > 1)
        bestree2->copy(bestree2, curtree);
      curtree->root = curtree->nodep[outgrno - 1]->back;
      restml_printree();
      summarize();
      if (trout)
      {
        col = 0;
        restml_treeout(curtree->root);
      }
    }
    destruct_tree(curtree);
  }

  if ( jumb < njumble )
    return;
  freex2(nonodes2, curtree->nodep);
  if (!usertree)
  {
    freex2(nonodes2, priortree->nodep);
    freex2(nonodes2, bestree->nodep);
    if (njumble > 1)
      freex2(nonodes2, bestree2->nodep);
  }
  else
  {
    free(l0gl);
    for (i=0;i<shimotrees;++i)
      free(l0gf[i]);
    free(l0gf);
  }
  if (jumb == njumble)
  {
    if (progress)
    {
      printf("\nOutput written to file \"%s\".\n\n", outfilename);
      if (trout)
        printf("Tree also written onto file \"%s\".\n", outtreename);
      putchar('\n');
    }
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* maximum likelihood phylogenies from restriction sites */
  initdata *funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Restml", "");        */
  argv[0] = "RestML";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = restml_node_new;
  funcs->tree_new = restml_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  for (ith = 1; ith <= datasets; ith++)
  {
    if (datasets > 1)
    {
      fprintf(outfile, "Data set # %ld:\n", ith);
      if (progress)
        printf("\nData set # %ld:\n", ith);
    }
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      maketree();
      fflush(outfile);
      fflush(outtree);
    }
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
  return 0;
}  /* maximum likelihood phylogenies from restriction sites */


// End.
