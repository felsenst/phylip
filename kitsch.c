/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <float.h>
#include "ml.h"
#include "phylip.h"
#include "dist.h"

#define epsilonk         0.000001   /* a very small but not too small number */

typedef struct kitsch_node {
  dist_node dist_node;
  double weight;
  boolean processed;
} kitsch_node;

typedef struct kitsch_tree {
  ml_tree ml_tree;
} kitsch_tree;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   doinit(void);
void   inputoptions(void);
void   getinput(void);
void   input_data(void);
void   scrunchtraverse(node *, node **, double *);
void   combine(node *, node *);
void   scrunch(node *);

void   secondtraverse(node *, node *, node *, node *, long, long, long, double *);
void   firstraverse(node *, node *, double *);
void   sumtraverse(node *, double *);
void   dtraverse(node *);
void   describe(void);
void   maketree(void);
node * kitsch_node_new(node_type type, long index);
void   kitsch_node_init(node* n, node_type type, long index);
void   kitsch_node_copy(node* src, node* dst);
tree * kitsch_tree_new(long nonodes, long spp);
void   kitsch_tree_init(tree* t, long nonodes, long spp);
double kitsch_tree_evaluate(tree* t, node *r, boolean dummy);
void   kitschrun(void);
void   kitsch(char * infilename, char * intreename, char * outfilename, char * outfileopt,
              char * outtreename, char * outtreeopt, char * Method, int BestTree, double Power,
              int NegLengths, int LowerTMat, int UpperTMat, int Subreps, int RandInput, int RandNum,
              int Njumble, int MultData, int NumSets, int PrintData, int PrintInd, int PrintTree,
              int WriteTree);
/* function prototypes */
#endif

/* to make ml.c happy */
boolean lngths, smoothit = true, smoothed = false, polishing  = false;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH];
long numtrees, col, datasets, ith, njumble, jumb = 0, nonodes = 0;

/*   numtrees is used by usertree option part of maketree */
long inseed;
tree *curtree, *bestree, *priortree;   /* pointers to all nodes in tree */
boolean minev, jumble, usertree, lower, upper, negallowed, replicates, trout, printdata, progress, treeprint, mulsets, firstset;
longer seed;
double power;
long *enterorder;

/* Local variables for maketree, propagated globally for C version: */
long examined;
double like, bestyet;
node *there;
boolean *names;
Char ch;
char *progname;
double trweight; /* to make treeread happy */
boolean goteof, haslengths, lengths;  /* ditto ... */
long rcategs;


node* kitsch_node_new(node_type type, long index)
{
  node* n;
  n = Malloc(sizeof(kitsch_node));
  kitsch_node_init(n, type, index);
  return n;
} /* kitsch_node_new */


void kitsch_node_init(node* n, node_type type, long index)
{
  dist_node_init(n, type, index);
  n->copy = kitsch_node_copy;
} /* kitsch_node_init */


void kitsch_node_copy(node* srcn, node* dstn)
{
  kitsch_node *src = (kitsch_node *)srcn;
  kitsch_node *dst = (kitsch_node *)dstn;
  dist_node_copy(srcn, dstn);
  dst->weight = src->weight;
  dst->processed = src->processed;
} /* kitsch_node_copy */


tree* kitsch_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(kitsch_tree));
  kitsch_tree_init(t, nonodes, spp);
  return t;
} /* kitsch_tree_new */


void kitsch_tree_init(tree* t, long nonodes, long spp)
{
  ml_tree_init(t, nonodes, spp);
  t->globrearrange = rooted_globrearrange;
  t->insert_ = rooted_tree_insert_;
  t->re_move = rooted_tree_re_move;
  t->locrearrange = rooted_locrearrange;
  t->save_lr_nodes = rooted_tree_save_lr_nodes;
  t->restore_lr_nodes = rooted_tree_restore_lr_nodes;
  t->evaluate = kitsch_tree_evaluate;
  t->smoothall = (tree_smoothall_t)no_op;
} /* kitsch_tree_init */


void getoptions(void)
{
  /* interactively set options */
  long inseed0, loopcount;
  Char ch;

  minev = false;
  jumble = false;
  njumble = 1;
  lower = false;
  negallowed = false;
  power = 2.0;
  replicates = false;
  upper = false;
  usertree = false;
  trout = true;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  for(;;) {
    cleerhome();
    printf("\nFitch-Margoliash method ");
    printf("with contemporary tips, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  D      Method (F-M, Minimum Evolution)?  %s\n",
           (minev ? "Minimum Evolution" : "Fitch-Margoliash"));
    printf("  U                 Search for best tree?  %s\n",
           usertree ? "No, use user trees in input file" : "Yes");
    printf("  P                                Power?%9.5f\n", power);
    printf("  -      Negative branch lengths allowed?  %s\n",
           (negallowed ? "Yes" : "No"));
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    if (!usertree) {
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
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("DJUP-LRSM12340", ch) != NULL))
        || (usertree && ((strchr("DUP-LRSM12340", ch) != NULL))))
    {
      switch (ch)
      {
        case 'D':
          minev = !minev;
          if (!negallowed)
            negallowed = true;
          break;

        case '-':
          negallowed = !negallowed;
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
    countup(&loopcount, 100);
  }
  if (upper && lower) {
    printf("ERROR:  Data matrix cannot be both uppeR and Lower triangular.\n");
    exxit(-1);
  }
}  /* getoptions */


void doinit(void)
{
  /* initializes variables */

  inputnumbers2(&spp, &nonodes, 1);
  if (!javarun)
  {
    getoptions();
  }
  curtree = functions.tree_new(nonodes, spp);
  if (!usertree && njumble > 1) {
    bestree = functions.tree_new(nonodes, spp);
    priortree = functions.tree_new(nonodes, spp);
  }
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
}  /* doinit */


void inputoptions(void)
{
  /* print options information */
  if (!firstset)
    samenumsp2(ith);
  fprintf(outfile, "\nFitch-Margoliash method ");
  fprintf(outfile, "with contemporary tips, version %s\n\n", VERSION);
  if (minev)
    fprintf(outfile, "Minimum evolution method option\n\n");
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                               ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
  fprintf(outfile, "negative branch lengths");
  if (!negallowed)
    fprintf(outfile, " not");
  fprintf(outfile, " allowed\n\n");
}  /* inputoptions */


void getinput(void)
{
  /* reads the input data */
  inputoptions();
}  /* getinput */


void input_data(void)
{
  /* read in distance matrix */
  long i, j, k, columns, n;
  boolean skipit, skipother;
  double x;
  columns = replicates ? 4 : 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  dist_tree_new(curtree, nonodes);
  if (!usertree && njumble > 1)
    dist_tree_new(bestree, nonodes);
  for (i = 0; i < spp; i++)
  {
    ((dist_node*)curtree->nodep[i])->d[i] = 0.0;
    ((dist_node*)curtree->nodep[i])->w[i] = 0.0;
    ((kitsch_node*)curtree->nodep[i])->weight = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 1; j <= spp; j++) {
      skipit = ((lower && j >= i + 1) || (upper && j <= i + 1));
      skipother = ((lower && i + 1 >= j) || (upper && i + 1 <= j));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        if(fscanf(infile, "%lf", &x) < 1)
        {
          printf("\nERROR reading input file.\n\n");
          exxit(-1);
        }
        ((dist_node*)curtree->nodep[i])->d[j - 1] = x;
        if (replicates)
        {
          if (eoln(infile))
            scan_eoln(infile);
          if(fscanf(infile, "%ld", &n) < 1)
          {
            printf("\nERROR reading input file.\n\n");
            exxit(-1);
          }
        }
        else
          n = 1;
        if (n > 0 && x < 0) {
          printf("NEGATIVE DISTANCE BETWEEN SPECIES%5ld AND %5ld\n",
                 i + 1, j);
          exxit(-1);
        }
        ((dist_node*)curtree->nodep[i])->w[j - 1] = n;
        if (skipother) {
          ((dist_node*)curtree->nodep[j - 1])->d[i] =
            ((dist_node*)curtree->nodep[i])->d[j - 1];
          ((dist_node*)curtree->nodep[j - 1])->w[i] =
            ((dist_node*)curtree->nodep[i])->w[j - 1];
        }
        if ((i == j) && (fabs(((dist_node*)curtree->nodep[i-1])->d[j-1]) > 0.000000001))
        {
          printf("\nERROR:  Diagonal element of row %ld of distance matrix ", i+2);
          printf("is not zero.\n");
          printf("        Is it a distance matrix?\n\n");
          exxit(-1);
        }
        if ((j < i) && (fabs(((dist_node*)curtree->nodep[i])->d[j-1]-
                             ((dist_node*)curtree->nodep[j-1])->d[i])
                        > 0.000000001)) {
          printf("ERROR:  Distance matrix is not symmetric:\n");
          printf("       (%ld,%ld) element and (%ld,%ld) element are unequal.\n", i+1, j+1, j+1, i+1);
          printf("        They are %10.6f and %10.6f, respectively.\n",
                 ((dist_node*)curtree->nodep[i])->d[j-1],
                 ((dist_node*)curtree->nodep[j])->d[i-1]);
          printf("        Is it a distance matrix?\n\n");
          exxit(-1);
        }
      }
    }
  }
  scan_eoln(infile);
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata)
  {
    for (i = 0; i < spp; i++)
    {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
      putc(' ', outfile);
      for (j = 1; j <= spp; j++)
      {
        fprintf(outfile, "%10.5f", ((dist_node*)curtree->nodep[i])->d[j - 1]);
        if (replicates)
          fprintf(outfile, " (%3ld)",
                  (long)((dist_node*)curtree->nodep[i])->w[j - 1]);
        if (j % columns == 0 && j < spp) {
          putc('\n', outfile);
          for (k = 1; k <= nmlngth + 1; k++)
            putc(' ', outfile);
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 0; i < spp; i++)
  {
    for (j = 0; j < spp; j++)
    {
      if (i + 1 != j + 1)
      {
        if (((dist_node*)curtree->nodep[i])->d[j] < epsilonk)
          ((dist_node*)curtree->nodep[i])->d[j] = epsilonk;
        ((dist_node*)curtree->nodep[i])->w[j] /=
          exp(power * log(((dist_node*)curtree->nodep[i])->d[j]));
      }
    }
  }
}  /* inputdata */


void scrunchtraverse(node *u, node **closest, double *tmax)
{
  /* traverse to find closest node to the current one */
  if (!((dist_node*)u)->sametime) {
    if (((dist_node*)u)->t > *tmax) {
      *closest = u;
      *tmax = ((dist_node*)u)->t;
    }
    return;
  }
  ((dist_node*)u)->t = ((dist_node*)curtree->nodep[u->back->index - 1])->t;
  if (!u->tip) {
    scrunchtraverse(u->next->back, closest, tmax);
    scrunchtraverse(u->next->next->back, closest, tmax);
  }
}  /* scrunchtraverse */


void combine(node *a, node *b)
{
  /* put node b into the set having the same time as a */
  if (((kitsch_node*)a)->weight + ((kitsch_node*)b)->weight <= 0.0)
    ((dist_node*)a)->t = 0.0;
  else
    ((dist_node*)a)->t = (((dist_node*)a)->t * ((kitsch_node*)a)->weight + ((dist_node*)b)->t *
                          ((kitsch_node*)b)->weight) / (((kitsch_node*)a)->weight +
                                                        ((kitsch_node*)b)->weight);
  ((kitsch_node*)a)->weight += ((kitsch_node*)b)->weight;
  ((dist_node*)b)->sametime = true;
}  /* combine */


void scrunch(node *s)
{
  /* see if nodes can be combined to prevent negative lengths */
  double tmax;
  node *closest;
  boolean found;

  closest = NULL;
  tmax = -1.0;
  do {
    if (!s->tip) {
      scrunchtraverse(s->next->back, &closest, &tmax);
      scrunchtraverse(s->next->next->back, &closest, &tmax);
    }
    found = (tmax > ((dist_node*)s)->t);
    if (found)
      combine(s, closest);
    tmax = -1.0;
  } while (found);
}  /* scrunch */


void secondtraverse(node *a, node *q, node *u, node *v, long i, long j, long k, double *sum)
{
  /* recalculate distances, add to sum */
  long l;
  double wil, wjl, wkl, wli, wlj, wlk, TEMP;

  if (!(((kitsch_node*)a)->processed || a->tip)) {
    secondtraverse(a->next->back, q, u, v, i, j, k, sum);
    secondtraverse(a->next->next->back, q, u, v, i, j, k, sum);
    return;
  }
  if (!(a != q && ((kitsch_node*)a)->processed))
    return;
  l = a->index;
  wil = ((dist_node*)u)->w[l - 1];
  wjl = ((dist_node*)v)->w[l - 1];
  wkl = wil + wjl;
  wli = ((dist_node*)a)->w[i - 1];
  wlj = ((dist_node*)a)->w[j - 1];
  wlk = wli + wlj;
  ((dist_node*)q)->w[l - 1] = wkl;
  ((dist_node*)a)->w[k - 1] = wlk;
  if (wkl <= 0.0)
    ((dist_node*)q)->d[l - 1] = 0.0;
  else
    ((dist_node*)q)->d[l - 1] = (wil * ((dist_node*)u)->d[l - 1] + wjl *
                                 ((dist_node*)v)->d[l - 1]) / wkl;
  if (wlk <= 0.0)
    ((dist_node*)a)->d[k - 1] = 0.0;
  else
    ((dist_node*)a)->d[k - 1] = (wli * ((dist_node*)a)->d[i - 1] + wlj *
                                 ((dist_node*)a)->d[j - 1]) / wlk;
  if (minev)
    return;
  if (wkl > 0.0) {
    TEMP = ((dist_node*)u)->d[l - 1] - ((dist_node*)v)->d[l - 1];
    (*sum) += wil * wjl / wkl * (TEMP * TEMP);
  }
  if (wlk > 0.0) {
    TEMP = ((dist_node*)a)->d[i - 1] - ((dist_node*)a)->d[j - 1];
    (*sum) += wli * wlj / wlk * (TEMP * TEMP);
  }
}  /* secondtraverse */


void firstraverse(node *q_, node *r, double *sum)
{  /* firsttraverse                              */
   /* go through tree calculating branch lengths */
  node *q;
  long i, j, k;
  node *u, *v;

  q = q_;
  if (q == NULL)
    return;
  ((dist_node*)q)->sametime = false;
  if (!q->tip)
  {
    firstraverse(q->next->back, r, sum);
    firstraverse(q->next->next->back, r, sum);
  }
  ((kitsch_node*)q)->processed = true;
  if (q->tip)
    return;
  u = q->next->back;
  v = q->next->next->back;
  i = u->index;
  j = v->index;
  k = q->index;
  if (((dist_node*)u)->w[j - 1] + ((dist_node*)v)->w[i - 1] <= 0.0)
    ((dist_node*)q)->t = 0.0;
  else
    ((dist_node*)q)->t = (((dist_node*)u)->w[j - 1] * ((dist_node*)u)->d[j - 1] +
                          ((dist_node*)v)->w[i - 1] * ((dist_node*)v)->d[i - 1]) /
      (2.0 * (((dist_node*)u)->w[j - 1] + ((dist_node*)v)->w[i - 1]));
  ((kitsch_node*)q)->weight = ((kitsch_node*)u)->weight +
    ((kitsch_node*)v)->weight + ((dist_node*)u)->w[j - 1] +
    ((dist_node*)v)->w[i - 1];
  if (!negallowed)
    scrunch(q);
  u->v = ((dist_node*)q)->t - ((dist_node*)u)->t;
  v->v = ((dist_node*)q)->t - ((dist_node*)v)->t;
  u->back->v = u->v;
  v->back->v = v->v;
  secondtraverse(r, q, u, v, i, j, k, sum);
}  /* firstraverse */


void sumtraverse(node *q, double *sum)
{
  /* traverse to finish computation of sum of squares */
  long i, j;
  node *u, *v;
  double TEMP, TEMP1;

  if (minev && (q != curtree->root))
    *sum += q->v;
  if (q->tip)
    return;
  sumtraverse(q->next->back, sum);
  sumtraverse(q->next->next->back, sum);
  if (!minev) {
    u = q->next->back;
    v = q->next->next->back;
    i = u->index;
    j = v->index;
    TEMP = ((dist_node*)u)->d[j - 1] - 2.0 * ((dist_node*)q)->t;
    TEMP1 = ((dist_node*)v)->d[i - 1] - 2.0 * ((dist_node*)q)->t;
    (*sum) += ((dist_node*)u)->w[j - 1] * (TEMP * TEMP) +
      ((dist_node*)v)->w[i - 1] * (TEMP1 * TEMP1);
  }
}  /* sumtraverse */


double kitsch_tree_evaluate(tree* t, node *r, boolean dummy)
{
  /* fill in times and evaluate sum of squares for tree
   * speed improvement could be made here.  we only evaluate the function at the
   * root and it means that we have to recalculate a bit of data during each evaluate*/
  double sum = 0.0;
  long i;

  (void)dummy;                          // RSGnote: Parameter never used.

  r = t->root;
  for (i = 0; i < (t->nonodes); i++)
    ((kitsch_node*)curtree->nodep[i])->processed = curtree->nodep[i]->tip;
  firstraverse(r, r, &sum);
  sumtraverse(r, &sum);
  examined++;
  if (replicates && (lower || upper))
    sum /= 2;
  like = -sum;
  t->score = like;
  return like;
}  /* evaluate */


void dtraverse(node *q)
{
  /* print table of lengths etc. */
  long i;

  if (!q->tip)
    dtraverse(q->next->back);
  if (q->back != NULL) {
    fprintf(outfile, "%4ld   ", q->back->index - spp);
    if (q->index <= spp) {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->index - 1][i], outfile);
    } else
      fprintf(outfile, "%4ld      ", q->index - spp);
    fprintf(outfile, "%13.5f", ((dist_node*)curtree->nodep[q->back->index - 1])->t - ((dist_node*)q)->t);
    q->v = ((dist_node*)curtree->nodep[q->back->index - 1])->t - ((dist_node*)q)->t;
    q->back->v = q->v;
    fprintf(outfile, "%16.5f\n", ((dist_node*)curtree->root)->t - ((dist_node*)q)->t);
  }
  if (!q->tip)
    dtraverse(q->next->next->back);
}  /* dtraverse */


void describe(void)
{
  /* prints table of lengths, times, sum of squares, etc. */
  long i, j;
  double totalnum;
  double TEMP;

  if (!minev)
    fprintf(outfile, "\nSum of squares = %10.3f\n\n", -like);
  else
    fprintf(outfile, "Sum of branch lengths = %10.3f\n\n", -like);
  if ((fabs(power - 2) < 0.01) && !minev) {
    totalnum = 0.0;
    for (i = 0; i < spp; i++)
    {
      for (j = 0; j < spp; j++)
      {
        if (i + 1 != j + 1 && ((dist_node*)curtree->nodep[i])->d[j] > 0.0)
        {
          TEMP = ((dist_node*)curtree->nodep[i])->d[j];
          totalnum += ((dist_node*)curtree->nodep[i])->w[j] * (TEMP * TEMP);
        }
      }
    }
    totalnum -= 2;
    if (replicates && (lower || upper))
      totalnum /= 2;
    fprintf(outfile, "Average percent standard deviation =");
    fprintf(outfile, "%10.5f\n\n", 100 * sqrt(-(like / totalnum)));
  }
  fprintf(outfile, "From     To            Length          Height\n");
  fprintf(outfile, "----     --            ------          ------\n\n");
  dtraverse(curtree->root);
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeoutr(curtree->root, &col, curtree);
  }
}  /* describe */


void maketree(void)
{
  /* constructs a binary tree from the pointers in curtree.nodep.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, which;
  double *bestfound, bestlike, bstlike2=0;
  boolean lastrearr;
  node *item;

  if (!usertree)
  {
    if (jumb == 1)
    {
      input_data();
      examined = 0;
    }
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree->root = curtree->nodep[enterorder[0] - 1];
    curtree->insert_(curtree, curtree->nodep[enterorder[1]-1], curtree->nodep[enterorder[0] - 1], false);
    if (progress)
    {
      sprintf(progbuf, "Adding species:\n");
      print_progress(progbuf);
      writename(0, 2, enterorder);
      phyFillScreenColor();
    }
    for (i = 3; i <= spp; i++)
    {
      bestyet = -DBL_MAX;
      item = curtree->nodep[enterorder[i - 1] - 1];
      curtree->addtraverse(curtree, item, curtree->root, false, there, &bestyet,
                             NULL, false, false, false, bestfound);
      curtree->insert_(curtree, item, there, true);
      like = bestyet;
      curtree->locrearrange(curtree, curtree->root, true, &bestyet,
                              bestree, priortree, lastrearr, bestfound);
      examined--;
      if (progress)
      {
        writename(i - 1, 1, enterorder);
        phyFillScreenColor();
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          sprintf(progbuf, "\nDoing global rearrangements\n");
          print_progress(progbuf);
          sprintf(progbuf, "  !");
          print_progress(progbuf);
          for (j = 1; j <= nonodes; j++)
          {
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
            {
                sprintf(progbuf, "-");
                print_progress(progbuf);
            }
          }
          sprintf(progbuf, "!\n");
          print_progress(progbuf);
          phyFillScreenColor();
        }
        bestlike = bestyet;
        curtree->globrearrange(curtree, bestree, progress, true, bestfound);
        if (njumble > 1) {
          if (jumb == 1 || (jumb > 1 && bestlike > bstlike2))
          {
            curtree->copy(curtree, bestree);
            bstlike2 = bestlike;
          }
        }
      }
    }
    if (njumble == jumb) {
      if (njumble > 1)
        bestree->copy(bestree, curtree);
      curtree->evaluate(curtree, curtree->root, false);
      printree(curtree->root, treeprint, true);
      if (treeprint)
        describe();
    }
  } else {
    input_data();
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if (treeprint)
      fprintf(outfile, "\n\nUser-defined trees:\n\n");
    names = (boolean *)Malloc(spp * sizeof(boolean));
    which = 1;
    while (which <= numtrees ) {
      treeread2 (intree, &curtree->root, curtree->nodep, lengths, &trweight,
                 &goteof, &haslengths, &spp, false, nonodes);
      if (curtree->root->back) {
        printf("Error:  Kitsch cannot read unrooted user trees\n");
        exxit(-1);
      }
      curtree->evaluate(curtree, curtree->root, false);
      printree(curtree->root, treeprint, true);
      if (treeprint)
        describe();
      which++;
    }
    FClose(intree);
    free(names);
  }
  if (jumb == njumble && progress) {
    sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree also written onto file \"%s\".\n", outtreename);
      print_progress(progbuf);
    }
    sprintf(progbuf,  "\n");
    print_progress(progbuf);
  }
  destruct_tree(curtree);

}  /* maketree */


void kitschrun(void)
{
  //printf("in kitschrun\n");
  fflush(stdout);

  // debug print JRMdebug
  /*
    printf("\nminev %i\n", minev);
    printf("jumble %i\n", jumble);
    printf("njumble %li\n", njumble);
    printf("lower %i\n", lower);
    printf("negallowed %i\n", negallowed);
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

  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %ld:\n", ith);
      if (progress) {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    getinput();
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    firstset = false;
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
    fflush(outfile);
    fflush(outtree);
  }
}


void kitsch(
  char * infilename,
  char * intreename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  char * Method,
  int BestTree,
  double Power,
  int NegLengths,
  int LowerTMat,
  int UpperTMat,
  int Subreps,
  int RandInput,
  int RandNum,
  int Njumble,
  int MultData,
  int NumSets,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree
  )
{
  initdata *funcs;
  //printf("Hello from Kitsch!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Kitsch";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = kitsch_node_new;
  funcs->tree_new = kitsch_tree_new;
  phylipinit(argc, argv, funcs, true);
  progname = argv[0];
  /*
  //minev = false;
  //jumble = false;
  //njumble = 1;
  //lower = false;
  //negallowed = false;
  //power = 2.0;
  //replicates = false;
  //upper = false;
  //usertree = false;
  //trout = true;
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
  //double Power,
  //int NegLengths,
  //int LowerTMat,
  //int UpperTMat,
  //int Subreps,
  //int RandInput,
  //int RandNum,
  //int Njumble,
  //int MultData,
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
  kitschrun();

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
  //printf("Done.\n\n");  // JRM debug
} /* kitsch */


int main(int argc, Char *argv[])
{  /* Fitch-Margoliash criterion with contemporary tips */
  initdata* funcs;

#ifdef MAC
  argc = 1;                /* macsetup("Kitsch", "");        */
  argv[0] = "Kitsch";
#endif

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = kitsch_node_new;
  funcs->tree_new = kitsch_tree_new;
  phylipinit(argc, argv, funcs, false);

  /* reads in spp, options, and the data, then calls maketree to construct the tree */
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  kitschrun();

  FClose(infile);
  FClose(outfile);
  FClose(outtree);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  phyRestoreConsoleAttributes();
  printf("Done.\n\n");
  return 0;
}  /* Fitch-Margoliash criterion with contemporary tips */


// End.
