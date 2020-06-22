/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "dist.h"
#include "ml.h"

#define zsmoothings     10    /* number of zero-branch correction iterations */
#define epsilonf        0.000001   /* a very small but not too small number  */
#define delta           0.0001      /* a not quite so small number */

typedef struct fitch_tree{
  ml_tree ml_tree;
} fitch_tree;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   fitch_getinput(void);
void   secondtraverse(node *, double, long *, double *);
void   firsttraverse(node *, long *, double *);
double fitch_evaluate(tree *, node*, boolean);
void   nudists(node *, node *);
void   fitch_makenewv(tree* t, node *p);
void   makedists(node *);
void   makebigv(node *);
void   correctv(node *);
void   alter(node *, node *);
void   fitch_nuview(tree*, node *);
void   insert_(node *, node *, boolean);
void   fitch_setuptip(long, tree *);
void   fitch_buildnewtip(long, tree *, long);
void   fitch_buildsimpletree(tree *, long);
void   addtraverse(node *, node *, boolean, long *, boolean *);
void   rearrange(node *, long *, long *, boolean *);
void   describe(node *);
void   summarize(long);
void   nodeinit(node *);
void   initrav(node *);
void   treevaluate(void);
void   maketree(void);
void   globrearrange(long* numtrees, boolean* succeeded);
tree*  fitch_tree_new(long nonodes, long spp);
void   fitch_tree_init(tree* t, long nonodes, long spp);
void   fitchrun(void);
void   fitch(char * infilename, char * intreename, char * outfilename, char * outfileopt, char * outtreename,
             char * outtreeopt, char * Method, int BestTree, int UseLengths, double Power, int NegLengths,
             int LowerTMat, int UpperTMat, int Subreps, int GlobalRearr, int RandInput, int RandNum,
             int Njumble, int OutRoot, int OutNum, int MultData, int NumSets, int PrintData,
             int PrintInd, int PrintTree, int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH];
long outgrno, nums, col, datasets, ith, njumble, jumb = 0, nonodes = 0;
long inseed;
vector *x;
intvector *reps;
boolean minev, global, jumble, lngths, usertree, lower, upper, negallowed, outgropt, replicates, trout, printdata, progress, treeprint, mulsets, firstset, smoothit = true, smoothed = false, polishing;
double power;
double trweight; /* to make treeread happy */
boolean goteof, haslengths;  /* ditto ... */
boolean first; /* ditto ... */
node *addwhere;
longer seed, endsite, rcategs;
long *enterorder;
tree *curtree, *priortree, *bestree, *bestree2;
Char ch;
char *progname;


tree* fitch_tree_new(long nonodes, long spp)
{
  /* initialize a tree for Fitch */
  tree* t;

  t = Malloc(sizeof(fitch_tree));
  fitch_tree_init(t, nonodes, spp);
  return t;
} /* fitch_tree_new */


void fitch_tree_init(tree* t, long nonodes, long spp)
{
  /* set up functions for a tree for Fitch */

  fitch_tree *ft = (fitch_tree *)t;
  ml_tree_init(&(ft->ml_tree.tree), nonodes, spp);
  t->evaluate = fitch_evaluate;
  t->insert_ = ml_tree_insert_;
  t->re_move = ml_tree_re_move;
  t->nuview = fitch_nuview;
  ft->ml_tree.makenewv = fitch_makenewv;
} /* fitch_tree_init */


void getoptions(void)
{
  /* interactively set options */
  long inseed0=0, loopcount;
  Char ch;
  boolean done=false;

  putchar('\n');
  minev = false;
  global = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  lower = false;
  negallowed = false;
  outgrno = 1;
  outgropt = false;
  power = 2.0;
  replicates = false;
  trout = true;
  upper = false;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nFitch-Margoliash method version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  D      Method (F-M, Minimum Evolution)?  %s\n",
           (minev ? "Minimum Evolution" : "Fitch-Margoliash"));
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  N          Use lengths from user trees?  %s\n",
             (lngths ? "Yes" : "No"));
    }
    printf("  P                                Power?%9.5f\n", power);
    printf("  -      Negative branch lengths allowed?  %s\n",
           negallowed ? "Yes" : "No");
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  L         Lower-triangular data matrix?");
    if (lower)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  R         Upper-triangular data matrix?");
    if (upper)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  S                        Subreplicates?");
    if (replicates)
      printf("  Yes\n");
    else
      printf("  No\n");
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run");
    if (progress)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done)
    {
      if (((!usertree) && (strchr("DJOUNPG-LRSM01234", ch) != NULL))
          || (usertree && ((strchr("DOUNPG-LRSM01234", ch) != NULL))))
      {
        switch (ch)
        {
          case 'D':
            minev = !minev;
            if (minev && (!negallowed))
              negallowed = true;
            break;

          case '-':
            negallowed = !negallowed;
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
            lower = !lower;
            break;

          case 'N':
            lngths = !lngths;
            break;

          case 'O':
            outgropt = !outgropt;
            if (outgropt)
              initoutgroup(&outgrno, spp);
            break;

          case 'P':
            initpower(&power);
            break;

          case 'R':
            upper = !upper;
            break;

          case 'S':
            replicates = !replicates;
            break;

          case 'U':
            usertree = !usertree;
            break;

          case 'M':
            mulsets = !mulsets;
            if (mulsets)
              initdatasets(&datasets);
            jumble = true;
            if (jumble)
              initseed(&inseed, &inseed0, seed);
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
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
  if (lower && upper) {
    printf("ERROR:  Data matrix cannot be both uppeR and Lower triangular.\n");
    exxit(-1);
  }
}  /* getoptions */


void allocrest(void)
{
  /* allocate routing stuff */
  long i;

  x = (vector *)Malloc(spp * sizeof(vector));
  reps = (intvector *)Malloc(spp * sizeof(intvector));
  for (i=0;i<spp;++i)
  {
    x[i]=(vector)Malloc(nonodes * sizeof(double));
    reps[i]=(intvector)Malloc(spp * sizeof(long));
  }
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
} /* allocrest */


void doinit(void)
{
  /* initializes variables */

  inputnumbers2(&spp, &nonodes, 1);

  if (!javarun)
  {
    getoptions();
  }

  if (!usertree)
    nonodes--;

  curtree = functions.tree_new(nonodes, spp);
  if (!usertree) {
    bestree = functions.tree_new(nonodes, spp);
    priortree = functions.tree_new(nonodes, spp);
    if (njumble > 1)
      bestree2 = functions.tree_new(nonodes, spp);
  }

  allocrest();
}  /* doinit */


void inputoptions(void)
{
  /* print options information */
  if (!firstset)
    samenumsp2(ith);
  fprintf(outfile, "\nFitch-Margoliash method version %s\n\n", VERSION);
  if (minev)
    fprintf(outfile, "Minimum evolution method option\n\n");
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                   i  j        ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                             Obs\n\n");
  fprintf(outfile, "Negative branch lengths ");
  if (!negallowed)
    fprintf(outfile, "not ");
  fprintf(outfile, "allowed\n\n");
  if (global)
    fprintf(outfile, "global optimization\n\n");
}  /* inputoptions */


void fitch_getinput(void)
{
  /* reads the input data */
  inputoptions();
}  /* fitch_getinput */


void secondtraverse(node *q, double y, long *nx, double *sum)
{
  /* from each of those places go back to all others */
  /* nx comes from firsttraverse */
  /* sum comes from evaluate via firsttraverse */
  double z=0.0, TEMP=0.0;

  z = y + q->v;
  if (q->tip) {
    TEMP = ((dist_node*)q)->d[(*nx) - 1] - z;
    *sum += ((dist_node*)q)->w[(*nx) - 1] * (TEMP * TEMP);
  } else {
    secondtraverse(q->next->back, z, nx, sum);
    secondtraverse(q->next->next->back, z, nx, sum);
  }
}  /* secondtraverse */


void firsttraverse(node *p, long *nx, double *sum)
{
  /* go through tree calculating branch lengths */
  if (minev && (p != curtree->root))
    *sum += p->v;
  if (p->tip) {
    if (!minev) {
      *nx = p->index;
      secondtraverse(p->back, 0.0, nx, sum);
      }
  } else {
    firsttraverse(p->next->back, nx, sum);
    firsttraverse(p->next->next->back, nx, sum);
  }
}  /* firsttraverse */


double fitch_evaluate(tree *t, node* p, boolean dummy2)
{
  double sum=0.0;
  long nx=0;
  /* evaluate likelihood of a tree */
  generic_tree_evaluate(t, p, dummy2);
  firsttraverse(p->back, &nx, &sum);
  firsttraverse(p, &nx, &sum);
  if ((!minev) && replicates && (lower || upper))
    sum /= 2;
  t->score = -sum;
  return (-sum);
}  /* fitch_evaluate */


void nudists(node *x, node *y)
{
  /* compute distance between an interior node and tips */
  long nq=0, nr=0, nx=0, ny=0;
  double dil=0, djl=0, wil=0, wjl=0, vi=0, vj=0;
  node *qprime, *rprime;

  qprime = x->next;
  rprime = qprime->next->back;
  qprime = qprime->back;
  ny = y->index;
  dil = ((dist_node*)qprime)->d[ny - 1];
  djl = ((dist_node*)rprime)->d[ny - 1];
  wil = ((dist_node*)qprime)->w[ny - 1];
  wjl = ((dist_node*)rprime)->w[ny - 1];
  vi = qprime->v;
  vj = rprime->v;
  ((dist_node*)x)->w[ny - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    ((dist_node*)x)->d[ny - 1] = 0.0;
  else
    ((dist_node*)x)->d[ny - 1] = ((dil - vi) * wil + (djl - vj) * wjl) /
      (wil + wjl);
  nx = x->index;
  nq = qprime->index;
  nr = rprime->index;
  dil = ((dist_node*)y)->d[nq - 1];
  djl = ((dist_node*)y)->d[nr - 1];
  wil = ((dist_node*)y)->w[nq - 1];
  wjl = ((dist_node*)y)->w[nr - 1];
  ((dist_node*)y)->w[nx - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    ((dist_node*)y)->d[nx - 1] = 0.0;
  else
    ((dist_node*)y)->d[nx - 1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl);
}  /* nudists */


void makedists(node *p)
{
  /* compute distances among three neighbors of a node */
  long i=0, nr=0, ns=0;
  node *q, *r, *s;

  r = p->back;
  nr = r->index;
  for (i = 1; i <= 3; i++) {
    q = p->next;
    s = q->back;
    ns = s->index;
    if (((dist_node*)s)->w[nr - 1] + ((dist_node*)r)->w[ns - 1] <= 0.0)
      ((dist_node*)p)->dist = 0.0;
    else
      ((dist_node*)p)->dist =
        (((dist_node*)s)->w[nr - 1] * ((dist_node*)s)->d[nr - 1] +
         ((dist_node*)r)->w[ns - 1] * ((dist_node*)r)->d[ns - 1]) /
        (((dist_node*)s)->w[nr - 1] + ((dist_node*)r)->w[ns - 1]);
    p = q;
    r = s;
    nr = ns;
  }
}  /* makedists */


void makebigv(node *p)
{
  /* make new branch length */
  long i=0;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    if (p->iter) {
      p->v = (((dist_node*)p)->dist + ((dist_node*)r)->dist -
              ((dist_node*)q)->dist) / 2.0;
      p->back->v = p->v;
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */


void correctv(node *p)
{
  /* iterate branch lengths if some are to be zero */
  node *q, *r, *temp;
  long i=0, j=0, n=0, nq=0, nr=0, ntemp=0;
  double wq=0.0, wr=0.0;

  q = p->next;
  r = q->next;
  n = p->back->index;
  nq = q->back->index;
  nr = r->back->index;
  for (i = 1; i <= zsmoothings; i++) {
    for (j = 1; j <= 3; j++) {
      if (p->iter) {
        wr = ((dist_node*)r->back)->w[n - 1] +
          ((dist_node*)p->back)->w[nr - 1];
        wq = ((dist_node*)q->back)->w[n - 1] + ((dist_node*)p->back)->w[nq - 1];
        if (wr + wq <= 0.0 && !negallowed)
          p->v = 0.0;
        else
          p->v = ((((dist_node*)p)->dist - q->v) * wq +
                  (((dist_node*)r)->dist - r->v) * wr) / (wr + wq);
        if (p->v < 0 && !negallowed)
          p->v = 0.0;
        p->back->v = p->v;
      }
      temp = p;
      p = q;
      q = r;
      r = temp;
      ntemp = n;
      n = nq;
      nq = nr;
      nr = ntemp;
    }
  }
}  /* correctv */


void alter(node *x, node *y)
{
  /* traverse updating these views */
  nudists(x, y);
  if (!y->tip) {
    alter(x, y->next->back);
    alter(x, y->next->next->back);
  }
}  /* alter */


void fitch_nuview(tree* t, node *p)
{
  /* renew information about subtrees */
  node *q;

  alter(p, p->back);
  for (q = p->next ; q != p ; q = q->next )
    alter(q, q->back);
  p->initialized = true;
}  /* nuview */


void fitch_makenewv(tree* t, node *p)
{
  /* update branch lengths around a node */
  boolean iter;
  node* q;

  if (p->tip)
    return;
  makedists(p);

  iter = p->iter;
  for ( q = p->next ; p != q ; q = q->next )
    iter = iter || q->iter;

  if (iter) {
    if ( count_sibs(p) == 2)
      makebigv(p);
    correctv(p);
  }
  t->nuview(t, p);
}  /* update */


void fitch_setuptip(long m, tree *t)
{
  /* initialize branch lengths and views in a tip */
  long i=0;
  intvector n=(long *)Malloc(spp * sizeof(long));
  dist_node *WITH;

  WITH = (dist_node*)t->nodep[m - 1];
  memcpy(WITH->d, x[m - 1], (nonodes * sizeof(double)));
  memcpy(n, reps[m - 1], (spp * sizeof(long)));
  for (i = 0; i < spp; i++) {
    if (i + 1 != m && n[i] > 0) {
      if (WITH->d[i] < epsilonf)
        WITH->d[i] = epsilonf;
      WITH->w[i] = n[i] / exp(power * log(WITH->d[i]));
    } else {
      WITH->w[i] = 0.0;
      WITH->d[i] = 0.0;
    }
  }
  for (i = spp; i < nonodes; i++) {
    WITH->w[i] = 1.0;
    WITH->d[i] = 0.0;
  }
  WITH->node.index = m;
  if (WITH->node.iter) WITH->node.v = 0.0;
  free(n);
}  /* fitch_setuptip */


void fitch_buildnewtip(long m, tree *t, long nextsp)
{
  /* initialize and hook up a new tip */
  fitch_setuptip(m, t);
}  /* fitch_buildnewtip */


void fitch_buildsimpletree(tree *t, long nextsp)
{
  /* make and initialize a three-species tree */
  fitch_setuptip(enterorder[0], t);
  fitch_setuptip(enterorder[1], t);
  fitch_setuptip(enterorder[2], t);
  buildsimpletree(t, enterorder);
}  /* fitch_buildsimpletree */


void describe(node *p)
{
  /* print out information for one branch */
  long i=0;
  node *q;

  q = p->back;
  fprintf(outfile, "%4ld          ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.5f\n", q->v);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */


void summarize(long numtrees)
{
  /* print out branch lengths etc. */
  long i, j, totalnum;

  fprintf(outfile, "\nremember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
  if (!minev)
    fprintf(outfile, "Sum of squares = %11.5f\n\n", -curtree->score);
  else
    fprintf(outfile, "Sum of branch lengths = %11.5f\n\n", -curtree->score);
  if ((power == 2.0) && !minev) {
    totalnum = 0;
    for (i = 1; i <= nums; i++) {
      for (j = 1; j <= nums; j++) {
        if (i != j)
          totalnum += reps[i - 1][j - 1];
      }
    }
    fprintf(outfile, "Average percent standard deviation = ");
    fprintf(outfile, "%11.5f\n\n",
            100 * sqrt(-curtree->score / (totalnum - 2)));
  }
  fprintf(outfile, "Between        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  describe(curtree->root->next->back);
  describe(curtree->root->next->next->back);
  describe(curtree->root->back);
  fprintf(outfile, "\n\n");
  if (trout) {
    col = 0;
    treeout(curtree->root, &col, 0.43429445222, true, curtree->root);
  }
}  /* summarize */


void nodeinit(node *p)
{
  /* initialize a node */
  long i, j;

  for (i = 1; i <= 3; i++) {
    for (j = 0; j < nonodes; j++) {
      ((dist_node*)p)->w[j] = 1.0;
      ((dist_node*)p)->d[j] = 0.0;
    }
    p = p->next;
  }
  if ((!lngths) || p->iter)
    p->v = 1.0;
  if ((!lngths) || p->back->iter)
    p->back->v = 1.0;
}  /* nodeinit */


void initrav(node *p)
{
  /* traverse to initialize */
  if (p->tip)
    return;
  nodeinit(p);
  initrav(p->next->back);
  initrav(p->next->next->back);
}  /* initrav */


void treevaluate(void)
{
  /* evaluate user-defined tree, iterating branch lengths */
  long i;
  double oldlike;

  for (i = 1; i <= spp; i++)
    fitch_setuptip(i, curtree);
  unroot(curtree, nonodes);

  initrav(curtree->root);
  if (curtree->root->back != NULL) {
    initrav(curtree->root->back);
    curtree->evaluate(curtree, curtree->root, false);
    do {
      oldlike = curtree->score;
      curtree->smoothall(curtree, curtree->root);
      curtree->evaluate(curtree, curtree->root, false);
    } while (fabs(curtree->score - oldlike) > delta);
  }
  curtree->evaluate(curtree, curtree->root, false);
}  /* treevaluate */


void maketree(void)
{
  /* contruct the tree */
  long nextsp, numtrees=-1;
  boolean succeeded=false;
  long i, k, which;
  double bestyet;
  node *where, *p;
  boolean multf;

  if (usertree) {
    inputdata(replicates, printdata, lower, upper, x, reps);
    setuptree(curtree, nonodes);
    for (which = 1; which <= spp; which++)
      fitch_setuptip(which, curtree);
    if (eoln(infile))
      scan_eoln(infile);
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if (numtrees > MAXNUMTREES) {
      printf("\nERROR:  Number of input trees is read incorrectly from %s.\n", intreename);
      exxit(-1);
    }
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    first = true;
    which = 1;
    while (which <= numtrees) {
      treeread2 (intree, &curtree->root, curtree->nodep,
                 lngths, &trweight, &goteof, &haslengths, &spp, false, nonodes);
      nums = spp;
      treevaluate();
      printree(curtree->root, treeprint, false);
      if (treeprint)
        summarize(numtrees);
      destruct_tree(curtree);
      which++;
    }
    FClose(intree);
  } else {
    if (jumb == 1) {
      inputdata(replicates, printdata, lower, upper, x, reps);
      setuptree(curtree, nonodes);
      setuptree(priortree, nonodes);
      setuptree(bestree, nonodes);
      if (njumble > 1) setuptree(bestree2, nonodes);
    }
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    nextsp = 3;
    fitch_buildsimpletree(curtree, nextsp);
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    if (jumb == 1) numtrees = 1;
    nextsp = 4;
    if (progress) {
      sprintf(progbuf, "Adding species:\n");
      print_progress(progbuf);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }
    while (nextsp <= spp) {
      nums = nextsp;
      curtree->copy(curtree, priortree);
      fitch_buildnewtip(enterorder[nextsp - 1], curtree, nextsp);
      k = generic_tree_findemptyfork(curtree);
      p = curtree->get_fork(curtree, k);
      hookup(curtree->nodep[enterorder[nextsp-1]-1], p);
      p->v = initialv;
      p->back->v = initialv;
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
      curtree->root = curtree->nodep[enterorder[0] - 1]->back;
      curtree->addtraverse(curtree, p, curtree->root, true,
                            where, &bestyet, bestree, true);
      bestree->copy(bestree, curtree);
      if (progress) {
        writename(nextsp  - 1, 1, enterorder);
        phyFillScreenColor();
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        curtree->root = curtree->nodep[enterorder[0] - 1]->back;
        if (nextsp == spp  && global)
          curtree->globrearrange(curtree, progress, true);
        else{
          curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1], true,
                                priortree, bestree);
        }
        if (global && ((nextsp) == spp) && progress)
        {
          sprintf(progbuf, "\n   ");
          print_progress(progbuf);
        }
      }
      if (global && nextsp == spp) {
        putc('\n', outfile);
        if (progress)
        {
          sprintf(progbuf, "\n   ");
          print_progress(progbuf);
       }
      }
      curtree->copy(curtree, bestree);
      if (njumble > 1) {
        if (jumb == 1 && nextsp == spp)
          bestree->copy(bestree, bestree2);
        else if (nextsp == spp) {
          if (bestree2->score < bestree->score)
            bestree->copy(bestree, bestree2);
        }
      }
      nextsp++;
    }
    if (jumb == njumble) {
      if (njumble > 1) bestree2->copy(bestree2, curtree);
      curtree->root = curtree->nodep[outgrno - 1]->back;
      printree(curtree->root, treeprint, false);
      if (treeprint)
        summarize(numtrees);
    }
    destruct_tree(curtree);
  }
  if (jumb == njumble && progress) {
    sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout) {
      sprintf(progbuf, "Tree also written onto file \"%s\".\n", outtreename);
      print_progress(progbuf);
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
  }
}  /* maketree */


void fitchrun(void)
{
  /* do a single tree estimation */

  //printf("in fitchrun\n");
  fflush(stdout);

  // debug print JRMdebug
  /*
    printf("\nminev %i\n", minev);
    printf("global %i\n", global);
    printf("jumble %i\n", jumble);
    printf("njumble %li\n", njumble);
    printf("lngths %i\n", lngths);
    printf("lower %i\n", lower);
    printf("negallowed %i\n", negallowed);
    printf("outgrno %li\n", outgrno);
    printf("outgropt %i\n", outgropt);
    printf("power %f\n", power);
    printf("replicates %i\n", replicates);
    printf("trout %i\n", trout);
    printf("upper %i\n", upper);
    printf("usertree %i\n", usertree);
    printf("printdata %i\n", printdata);
    printf("progress %i\n", progress);
    printf("treeprint %i\n", treeprint);
    fflush(stdout);
  */

  int i;
  for (i=0;i<spp;++i) {
    enterorder[i]=0;}
  //printf("enterorder zeroed\n");//JRMdebug
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress) {
        sprintf(progbuf, "\nData set # %ld:\n\n", ith);
        print_progress(progbuf);
      }
    }

    fitch_getinput();
    //printf("fitch_getinput done, ith: %li\n", ith);  //JRMdebug
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      //printf("Calling maketree, jumb: %li\n", jumb);  //JRMdebug
      fflush(stdout);
      maketree();
    }

    firstset = false;
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
    fflush(outfile);
    fflush(outtree);
  }
} /* fitchrun */


void fitch(
  /* take parameters from Java interface */
  char * infilename,
  char * intreename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  char * Method,
  int BestTree,
  int UseLengths,
  double Power,
  int NegLengths,
  int LowerTMat,
  int UpperTMat,
  int Subreps,
  int GlobalRearr,
  int RandInput,
  int RandNum,
  int Njumble,
  int OutRoot,
  int OutNum,
  int MultData,
  int NumSets,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree
  )
{
  initdata *funcs;
  //printf("Hello from Fitch!\n");  // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Fitch";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dist_node_new;
  funcs->tree_new = fitch_tree_new;
  phylipinit(argc, argv, funcs, true);
  progname = argv[0];
  /*
  //minev = false;
  //global = false;
  //jumble = false;
  //njumble = 1;
  //lngths = false;
  //lower = false;
  //negallowed = false;
  //outgrno = 1;
  //outgropt = false;
  //power = 2.0;
  //replicates = false;
  //trout = true;
  //upper = false;
  //usertree = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;

  //char * infilename,
  //char * intreename,
  //char * outfilename,
  //char * outfileopt,
  //char * outtreename,
  //char * outtreeopt,
  //char * Method,
  //int BestTree,
  //int UseLengths,
  //double Power,
  //int NegLengths,
  //int LowerTMat,
  //int UpperTMat,
  //int Subreps,
  //int RandInput,
  //int RandNum,
  //int Njumble,
  //int OutRoot,
  //int OutNum,
  // int MultData,
  //int NumSets,
  //int PrintData,
  //int PrintInd,
  //int PrintTree,
  //int WriteTree
  */
  if (!strcmp(Method, "FM")) // Fitch-Margoliash
  {
    minev = false;
  }
  else //if (!strcmp(Method, "ME")) // Minimum Evolution
  {
    minev = true;
  }

  if (NegLengths != 0)
  {
    negallowed = true;
  }
  else
  {
    negallowed = false;
  }

  if (UseLengths != 0)
  {
    lngths = true;
  }
  else
  {
    lngths = false;
  }

  if (BestTree != 0)
  {
    usertree = false;
  }
  else
  {
    usertree = true;
  }

  if (LowerTMat != 0)
  {
    lower = true;
  }
  else
  {
    lower = false;
  }

  if (UpperTMat != 0)
  {
    upper = true;
  }
  else
  {
    upper = false;
  }

  if (Subreps != 0)
  {
    replicates = true;
  }
  else
  {
    replicates = false;
  }

  if (GlobalRearr != 0)
  {
    global = true;
  }
  else
  {
    global = false;
  }

  if (RandInput != 0)
  {
    jumble = true;
    inseed = RandNum;
    njumble = Njumble;
  }
  else
  {
    jumble = false;
    inseed = 1;
    njumble = 1;
  }

  if (OutRoot != 0)
  {
    outgropt = true;
    outgrno = OutNum;
  }
  else
  {
    outgropt = false;
    outgrno = 1;
  }

  if (MultData != 0)
  {
    mulsets = true;
    datasets = NumSets;
  }
  else
  {
    mulsets = false;
    datasets = 1;
  }

  if (WriteTree != 0)
  {
    trout = true;
  }
  else
  {
    trout = false;
  }

  if (PrintData != 0)
  {
    printdata = true;
  }
  else
  {
    printdata = false;
  }

  if (PrintInd != 0)
  {
    progress = true;
  }
  else
  {
    progress = false;
  }

  if (PrintTree != 0)
  {
    treeprint = true;
  }
  else
  {
    treeprint = false;
  }

  power = Power;

  // everything translated, start the run
  infile = fopen(infilename, "r");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  if (usertree)
  {
    intree = fopen(intreename, "r");
  }

  if (trout)
  {
    outtree = fopen(OuttreeName, outtreeopt);
    strcpy(outtreename, OuttreeName);
  }

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  firstset = true;

  doinit();
  fitchrun();

  if (trout)
    FClose(outtree);
  if (usertree)
    FClose(intree);
  FClose(outfile);
  FClose(infile);
  //printf("Done.\n\n");
} /* fitch */


int main(int argc, Char *argv[])
{
  initdata *funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Fitch", "");        */
  argv[0]="Fitch";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dist_node_new;
  funcs->tree_new = fitch_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;

  // do the work
  doinit();
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  fitchrun();

  if (trout)
    FClose(outtree);
  FClose(outfile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
} /* main */


// End.
