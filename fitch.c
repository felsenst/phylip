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
void   fitch_makenewv(tree*, node*);
void   makedists(node *);
void   makebigv(node *);
void   correctv(node *);
void   alter(node *, node *);
void   fitch_nuview(tree*, node *);
void   insert_(node *, node *, boolean);
void   fitch_setuptip(tree *, long);
void   fitch_buildnewtip(long, tree *, long);
void   initfitchnode(tree *, node **, long, long, long *, long *, initops,
                      pointarray, Char *, Char *, FILE *);
void   fitch_setupnewfork(tree *, long);
void   fitch_buildsimpletree(tree *, long);
void   rearrange(node *, long *, long *, boolean *);
node*  findroot(tree *, node *, boolean *);
void   describe(node *);
void   summarize(long);
void   nodeinit(node *);
void   initrav(node *);
void   treevaluate(void);
void   maketree(void);
void   globrearrange(long* numtrees, boolean* succeeded);
tree*  fitch_tree_new(long nonodes, long spp);
void   fitch_tree_init(tree*, long nonodes, long spp);
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
double trweight;                                  /* to make treeread happy */
boolean goteof, haslengths;                                    /* ditto ... */
boolean first;                                                 /* ditto ... */
node *addwhere, *there;
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
  t->do_newbl = true;
  fitch_tree_init(t, nonodes, spp);
  return t;
} /* fitch_tree_new */


void fitch_tree_init(tree* t, long nonodes, long spp)
{
  /* set up functions for a tree for Fitch */

  generic_tree_init((tree*)t, nonodes, spp);
  dist_tree_init((tree*)t, nonodes);   /* debug: need? */
  ((tree*)t)->evaluate = fitch_evaluate;
  ((tree*)t)->insert_ = ml_tree_insert_;
  ((tree*)t)->try_insert_ = ml_tree_try_insert_;
  ((tree*)t)->re_move = ml_tree_re_move;
  ((tree*)t)->nuview = fitch_nuview;
  ((tree*)t)->makenewv = fitch_makenewv;
  ((tree*)t)->smoothall = (tree_smoothall_t)ml_tree_smoothall;
  ((tree*)t)->do_newbl = true;
  ((tree*)t)->do_branchl_on_insert_f = ml_tree_do_branchl_on_insert;
  ((tree*)t)->do_branchl_on_re_move_f = ml_tree_do_branchl_on_re_move;
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
  fprintf(outfile, "                                    2\n");
  fprintf(outfile, "                  __ __  (Obs - Exp)\n");
  fprintf(outfile, "                  \\  \\   ------------\n");
  fprintf(outfile, "Sum of squares =  /_ /_        ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
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
  /* from each of those places go back to all others
   * nx   comes from firsttraverse
   * sum  comes from evaluate via firsttraverse */
  double z=0.0, TEMP=0.0;

  if (q) {
    z = y + q->v;
    if (q->tip) {
      TEMP = ((dist_node*)q)->d[(*nx) - 1] - z;
      *sum += ((dist_node*)q)->w[(*nx) - 1] * (TEMP * TEMP);
    } else {
      if (q->next->back != NULL)
        secondtraverse(q->next->back, z, nx, sum);
      if (q->next->next->back != NULL)
        secondtraverse(q->next->next->back, z, nx, sum);
    }
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
      if (p->back != NULL)
        secondtraverse(p->back, 0.0, nx, sum);
      }
  } else {
    if (p->next->back != NULL)
      firsttraverse(p->next->back, nx, sum);
    if (p->next->next->back != NULL)
      firsttraverse(p->next->next->back, nx, sum);
  }
}  /* firsttraverse */


double fitch_evaluate(tree *t, node* p, boolean dummy2)
{
  /* evaluate likelihood of a tree */
  double sum=0.0;
  long nx=0;

  generic_tree_evaluate(t, p, dummy2);
  if (p->back != NULL)
    firsttraverse(p->back, &nx, &sum);
  if (p != NULL)
  firsttraverse(p, &nx, &sum);
  if ((!minev) && replicates && (lower || upper))
    sum /= 2;
  t->score = -sum;
  return (-sum);
}  /* fitch_evaluate */


void nudists(node *x, node *y)
{
  /* compute distance between an interior node and tips 
   * y  is the node whose distance to the immediate descendants of  x
   * is being computed.  They are  qprime->back  and qprime-next->back.
   * This is the version for a bifurcation, so node has three neighbors */
  long nx=0, ny=0;    /* debug:  why this value? */
  double dil=0.0, djl=0.0, wil=0.0, wjl=0.0, vi=0.0, vj=0.0;
  node *qprime, *rprime, *qprimeback, *rprimeback;

  nx = x->index;
  ny = y->index;
  qprime = x->next;
  qprimeback = qprime->back;
  if (qprimeback != NULL) {
    vi = qprime->v;
    dil = ((dist_node*)qprimeback)->d[ny - 1];
    wil = ((dist_node*)qprimeback)->w[ny - 1];
  } else {
    dil = 0.0;
    wil = 0.0;
    qprime->v = 0.0;
  }
  rprime = qprime->next;
  rprimeback = rprime->back;
  if (rprimeback != NULL) {
    vj = rprime->v;
    djl = ((dist_node*)rprimeback)->d[ny - 1];
    wjl = ((dist_node*)rprimeback)->w[ny - 1];
  } else {
    djl = 0.0;
    wjl = 0.0;
    rprime->v = 0.0;
  }
  if (wil + wjl <= 0.0) {
    ((dist_node*)x)->d[ny - 1] = 0.0;
  }
  else {
    ((dist_node*)x)->d[ny - 1] = ((dil - vi) * wil + (djl - vj) * wjl) /
                                   (wil + wjl);
    ((dist_node*)x)->w[ny - 1] = wil + wjl;
  }
  ((dist_node*)y)->d[nx - 1] = ((dist_node*)x)->d[ny - 1];
  x->initialized = true;
}  /* nudists */


void makedists(node *p)
{
  /* compute distances among pairs of adjacent neighbors of an interior node.
   * (assumes a bifurcation so there are three */
  long npb=0, nqb=0, nrb=0;
  double d12, d23, d31;
  node *q, *r, *pb, *qb, *rb;

  if (p->back != NULL) {          /* node and number of the node behind  r  */
    pb = p->back;
    npb = pb->index;
  }
  q = p->next;
  if (q->back != NULL) {
    qb = q->back;
    nqb = qb->index;
  }
  r = q->next;
  if (r->back != NULL) {
    rb = r->back;
    nrb = rb->index;
  }
  if ((p->back != NULL) && (q->back != NULL)) {
    d12 = ((dist_node*)pb)->d[nqb - 1];
    }
  else {
    d12 = 0.0;
  }
  if ((q->back != NULL) && (r->back != NULL)) {
    d23 = ((dist_node*)qb)->d[nrb - 1];
    }
  else {
    d23 = 0.0;
  }
  if ((r->back != NULL) && (p->back != NULL)) {
    d31 = ((dist_node*)rb)->d[npb - 1];
  }
  else {
    d31 = 0.0;
  }
  ((dist_node*)p)->dist = d12;
  ((dist_node*)q)->dist = d23;
  ((dist_node*)r)->dist = d31;
}  /* makedists */


void makebigv(node *p)
{
  /* make new branch lengths around a bifurcating interior node
   * p->dist, q->dist, and r->dist are zero if near NULL root */
  long i=0;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    if (p->iter) {
      if (p->back != NULL) {
        p->v = (((dist_node*)p)->dist + ((dist_node*)r)->dist -
                ((dist_node*)q)->dist) / 2.0;
        p->back->v = p->v;
      }
      else {
        p->v = 0.0;
      }
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */


void correctv(node *p)
{
  /* iterate branch lengths if some are to be zero
   * Note this is only for a bifurcation with three neighbors of the
   * fork.  If any of these is null, don't do anything */
  node *q, *r, *temp;
  long i=0, j=0, n=0, nq=0, nr=0, ntemp=0;
  double wq=0.0, wr=0.0;

  q = p->next;
  r = q->next;
  if ((p->back != NULL) && (q->back != NULL)
        && (r->back != NULL)) {                /* skip all this if any NULL */
    n = p->back->index;
    nq = q->back->index;
    nr = r->back->index;
    for (i = 1; i <= zsmoothings; i++) {               /* do multiple times */
      for (j = 1; j <= 3; j++) {                   /* go around fork circle */
        if (p->iter) {
          wr = ((dist_node*)r->back)->w[n - 1] +
            ((dist_node*)p->back)->w[nr - 1];
          wq = ((dist_node*)q->back)->w[n - 1]
               + ((dist_node*)p->back)->w[nq - 1];
          if (((wr + wq) <= 0.0) && !negallowed)   /* if estimates megative */
            p->v = 0.0;
          else
            p->v = ((((dist_node*)p)->dist - q->v) * wq +
                    (((dist_node*)r)->dist - r->v) * wr) / (wr + wq);
          if (p->v < 0 && !negallowed)
            p->v = 0.0;
          p->back->v = p->v;
        }
        temp = p;                            /* move one step around circle */
        p = q;
        q = r;
        r = temp;
        ntemp = n;
        n = nq;
        nq = nr;
        nr = ntemp;
      }
    }
  }
}  /* correctv */


void alter(node *x, node *y)
{
  /* traverse updating these views */
  if (y != NULL) {
    nudists(x, y);  /* debug:  should this be after traversal? */
    if (!y->tip) {
      alter(x, y->next->back);
      alter(x, y->next->next->back);
    }
  }
}  /* alter */


void fitch_nuview(tree* t, node *p)
{
  double temp;
  /* renew information about subtrees */

  alter(p, p->back);
  if (p->back != NULL) {
    temp = ((dist_node*)(p))->d[p->back->index - 1];
    ((dist_node*)(p->back))->d[p->index - 1] = temp;
  }
}  /* fitch_nuview */


void fitch_makenewv(tree* t, node *p)
{
  /* update branch lengths around a node ,
   * if node is rootmost fork, be careful how you do this */
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
}  /* fitch_makenewv */


void fitch_setuptip(tree *t, long m)
{
  /* initialize branch lengths and views in a tip */
  long i=0;
  intvector n=(long *)Malloc(spp * sizeof(long));
  dist_node *which;

  which = (dist_node*)t->nodep[m - 1];
  memcpy(which->d, x[m - 1], (spp * sizeof(double)));  /* debug: too long? */
  memcpy(n, reps[m - 1], (spp * sizeof(long)));
  for (i = 0; i < spp; i++) {
    if (((i + 1) != m) && (n[i] > 0)) {
      if (which->d[i] < epsilonf)
        which->d[i] = epsilonf;
      which->w[i] = n[i] / exp(power * log(which->d[i]));
    } else {
      which->w[i] = 1.0;   /* debug: correct? */
      which->d[i] = 0.0;
    }
  }
  for (i = spp; i < nonodes; i++) {
    which->w[i] = 1.0;
    which->d[i] = 0.0;
  }
  which->node.index = m;
  if (which->node.iter) which->node.v = 0.0;
  free(n);
}  /* fitch_setuptip */


void fitch_buildnewtip(long m, tree *t, long nextsp)
{
  /* initialize and hook up a new tip */
  fitch_setuptip(t, m);
}  /* fitch_buildnewtip */


void initfitchnode(tree *treep, node **p, long len, long nodei, long *ntips,
                    long *parens, initops whichinit, pointarray nodep,
                    Char *str, Char *ch, FILE *intree)
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
/* debug:      ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs); */
      assert((*p)->index > 0);
      nodep[(*p)->index - 1] = (*p);
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
/* debug ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs);    */
      break;
    case tip:
      match_names_to_data (str, nodep, p, spp);
      break;
    case iter:
      (*p)->initialized = false;
      (*p)->v = initialv;
      (*p)->iter = true;
      if ((*p)->back != NULL)
      {
        (*p)->back->iter = true;
        (*p)->back->v = initialv;
        (*p)->back->initialized = false;
      }
      break;
    case length:
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      (*p)->v = valyew / divisor;
      (*p)->iter = false;
      if ((*p)->back != NULL)
      {
        (*p)->back->v = (*p)->v;
        (*p)->back->iter = false;
      }
      break;
    case hsnolength:
      haslengths = false;
      break;
    default:        /* cases hslength, treewt, unittrwt */
      break;        /* should never occur               */
  }
} /* initfitchnode */


void fitch_buildsimpletree(tree *t, long nextsp)
{
  /* make and initialize a three-species tree */
  fitch_setuptip(t, enterorder[0]);
  fitch_setuptip(t, enterorder[1]);
  fitch_setuptip(t, enterorder[2]);
  buildsimpletree(t, enterorder);
}  /* fitch_buildsimpletree */


void fitch_setupnewfork(tree *t, long m)
{ /* set up weights, distances for a new internal fork */
  long i;
  node* p;

  p = t->nodep[m-1];
  for (i = 0; i < nonodes; i++) {
    ((dist_node*)p)->w[i] = 1.0;
    ((dist_node*)p)->d[i] = 0.0;
  }
} /* setupnewtip */


void describe(node *p)
{
  /* print out information for one branch, number of fork or name of tip
   * at each end of the branch, and branch length */
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
  boolean start, found;
  node *p, *q;

  fprintf(outfile, "\nremember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
  if (!minev)                 /* print out quantity that is being minimized */
    fprintf(outfile, "Sum of squares = %11.5f\n\n", -curtree->score);
  else
    fprintf(outfile, "Sum of branch lengths = %11.5f\n\n", -curtree->score);
  if ((power == 2.0) && !minev) {  /* in original Fitch-Margoliash case ... */
    totalnum = 0;
    for (i = 1; i <= nums; i++) {                       /* compute APSD ... */
      for (j = 1; j <= nums; j++) {
        if (i != j)
          totalnum += reps[i - 1][j - 1];
      }
    }                                               /* ... and print it out */
    fprintf(outfile, "Average percent standard deviation = ");
    fprintf(outfile, "%11.5f\n\n",
            100 * sqrt(-curtree->score / (totalnum - 2)));
  }
  fprintf(outfile, "Between        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  q = curtree->root;
  q = findroot(curtree, q, &found);
  start = true;
  for (p = q; (start || (p != q)); p = p->next) {   /* around rootmost fork */
    start = false;
    if (p->back != NULL)    /* and if each node on circle has neighbors ... */
      describe(p->back);     /* recursively describe it and its descendants */
  }
  fprintf(outfile, "\n\n");
  if (trout) {             /* now write tree to output tree file if desired */
    col = 0;
    treeout(curtree->root, &col, 0.43429445222, true, curtree->root);
  }
}  /* summarize */


void nodeinit(node *p)
{
  /* initialize a node */
  long i, j;

  for (i = 1; i <= 3; i++) {        /* initialize its weights and distances */
    for (j = 0; j < nonodes; j++) {
      ((dist_node*)p)->w[j] = 1.0;
      ((dist_node*)p)->d[j] = 0.0;
    }
    p = p->next;
  }
  if ((!lngths) || p->iter)         /* and initial branch lengths if needed */
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
  double oldlike;

  unroot(curtree, nonodes);     /* debug: removes a root fork if bifurcating */

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
  boolean dummy_first=true, lastrearr, succeeded=false;
  long i, k, which, nextnode;
  double bestyet, *bestfound = NULL;
  node *p;

  if (usertree) {              /* the case where we are reading a user tree */
    nextnode = 0;
    inputdata(replicates, printdata, lower, upper, x, reps);
    dist_tree_init(curtree, nonodes);                  /* initialize a tree */
    for (which = 1; which <= spp; which++)               /* set up its tips */
      fitch_setuptip(curtree, which);
    if (eoln(infile))                 /* go to a new line if at end of line */
      scan_eoln(infile);
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree); /* get number of trees: count semicolons */
    if (numtrees > MAXNUMTREES) {
      printf("\nERROR:  Number of input trees is incorrect from %s.\n", intreename);
      exxit(-1);
    }
    if (treeprint) { /* print out that we are processing user-defined trees */
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    first = true;
    which = 1;
    while (which <= numtrees) {                           /* loop over them */
/* debug:      treeread2 (curtree, intree, &curtree->root, lngths, &trweight,
                  &goteof, &haslengths, &spp, false, nonodes);   debug */
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof,
                &dummy_first, &nextnode, &haslengths, initfitchnode,
                false, nonodes);
      nums = spp;
      treevaluate();  /* evaluate tree, if needed estimating branch lengths */
      printree(curtree->root, treeprint, false);            /* print it out */
      if (treeprint)             /* print table of nodes and branch lengths */
        summarize(numtrees);
      destruct_tree(curtree);                         /* dicombobulate tree */
      which++;
    }
    FClose(intree);
  } else {
    if (jumb == 1) {
      inputdata(replicates, printdata, lower, upper, x, reps);
    }
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    nextsp = 3;
    fitch_buildsimpletree(curtree, nextsp);
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    p = generic_newrootfork(curtree);
    fitch_setupnewfork(curtree, p->index);
    generic_insertroot(curtree, curtree->root, p); 
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
      fitch_setupnewfork(curtree, p->index);
      hookup(curtree->nodep[enterorder[nextsp-1]-1], p);
      p->v = initialv;
      p->back->v = initialv;
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
      curtree->root = curtree->nodep[enterorder[0] - 1]->back;
      curtree->addtraverse(curtree, p, curtree->root, further, there,
                             &bestyet, bestree, true, false, true, bestfound);
      bestree->copy(bestree, curtree);
      if (progress) {
        writename(nextsp  - 1, 1, enterorder);
        phyFillScreenColor();
      }
      succeeded = true;
      while (succeeded) {        /* as long as some rearrangement succeeded */
        succeeded = false;   /* debug: should move rootmost fork? */
        curtree->root = curtree->nodep[enterorder[0] - 1]->back;
        if ((nextsp == spp)  && global)    /* last time, SPR rearrangements */
          curtree->globrearrange(curtree, bestree, progress, true, bestfound);
        else {         /* earlier times, do nearest-neighbor rearrangements */ 
          curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1],
                    false, &bestyet, priortree, bestree, lastrearr, bestfound);
        }
        if (global && (nextsp == spp) && progress)
        {
          sprintf(progbuf, "\n   ");
          print_progress(progbuf);
        }
      }
      if (global && (nextsp == spp)) {
        {
          sprintf(progbuf, "\n   ");
          print_progress(progbuf);
        }
      }
      lastrearr = (nextsp == spp);
      if (lastrearr) {               /* if have finished adding all species */
        if (njumble > 1) {            /* if there is more than one jumbling */
          if (jumb == 1)
            bestree->copy(bestree, bestree2);  /* first tree is best so far */
          else {          /* if the tree is better than the best so far ... */
            if (bestree2->score < bestree->score)
              bestree->copy(bestree, bestree2); /* ... then put as best one */
          }
        }
      }
      nextsp++;
    }
    if (jumb == njumble) {       /* on last jumbling of species input order */
      if (njumble > 1) bestree2->copy(bestree2, curtree); /* bestree2 is it */
      putrootnearoutgroup(curtree, outgrno, true); /* root to near outgroup */
      printree(curtree->root, treeprint, false);    /* print the tree if OK */
      if (treeprint)
        summarize(numtrees);  /* print table of connections, branch lengths */
    }
    destruct_tree(curtree);                   /* recycle the tree structure */
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
  {
    FClose(outtree);
  }
  if (usertree)
  {
    FClose(intree);
  }
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
  functions.tree_new = fitch_tree_new;
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
  {
    FClose(outtree);
  }
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
