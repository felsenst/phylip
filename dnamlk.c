/* Version 4.0. (c) Copyright 1986-2013 by the University of Washington
   and by Joseph Felsenstein.  Written by Joseph Felsenstein.  Permission is
   granted to copy and use this program provided no fee is charged for it
   and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "seq.h"
#include "ml.h"

#define over            60

typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1,
         xiy1xv;
  double *ww, *zz, *wwzz, *vvzz;
} valrec;

typedef struct dnamlk_tree{
  ml_tree ml_tree;
} dnamlk_tree;

typedef double contribarr[maxcategs];

valrec ***tbl;

#ifndef OLDC
/* function prototypes */
void    getoptions(void);
void    allocrest(void);
void    doinit(void);
void    inputoptions(void);
void    makeweights(void);
void    getinput(void);
void    inittable_for_usertree (FILE *);
void    inittable(void);
void    alloc_nvd(long, nuview_data *);
void    free_nvd(nuview_data *);
void    dnamlk_tree_nuview(tree* t, node *);
double  dnamlk_tree_evaluate(tree* t, node *p, boolean);
void    restoradd(node *, node *, node *, double);
void    initdnamlknode(tree *, node **, long, long, long *, long *, initops, pointarray, Char *, Char *, FILE *);
void    tymetrav(node *, double *);
void    dnamlk_coordinates(node *, long *);
void    dnamlk_drawline(long, double);
void    dnamlk_printree(void);
void    describe(node *);
void    reconstr(node *, long);
void    rectrav(node *, long, long);
void    summarize(void);
void    dnamlk_treeout(node *);
void    nodeinit(node *);
void    initrav(node *);
void    travinit(tree*, node *);
void    travsp(node *);
void    treevaluate(void);
void    maketree(void);
void    reallocsites(void);
void    save_tree_tyme(tree* save_tree, double tymes[]);
void    freetable(void);
tree *  dnamlk_tree_new(long nonodes, long spp);
void    dnamlk_tree_init(tree * t, long nonodes, long spp);
void    dnamlkrun(void);
void    dnamlk(char * infilename, char * intreename, char * wgtsfilename, char * catsfilename,
               char * outfilename, char * outfileopt, char * outtreename, char * outtreeopt,
               char * TreeUseMethod, int UseLengths, double TTratio, int useEmpBF, double BaseFreqA,
               double BaseFreqC, double BaseFreqG, double BaseFreqTU, int OneCat, int NumCats,
               double SiteRate1, double SiteRate2, double SiteRate3, double SiteRate4, double SiteRate5,
               double SiteRate6, double SiteRate7, double SiteRate8, double SiteRate9, char * RateVar,
               int AdjCor, double BlockLen, double CoeffVar, int NumRates, double HMMRate1, double HMMRate2,
               double HMMRate3, double HMMRate4, double HMMRate5, double HMMRate6, double HMMRate7, double HMMRate8,
               double HMMRate9, double HMMProb1, double HMMProb2, double HMMProb3, double HMMProb4, double HMMProb5,
               double HMMProb6, double HMMProb7, double HMMProb8, double HMMProb9, double InvarFract, int SitesWeight,
               int GlobalRe, int RandInput, int RandNum, int Njumble, int MultData, int MultDSet, int NumSeqs,
               int InputSeq, int PrintData, int PrintInd, int PrintTree, int WriteTree, int DotDiff, int RecHypo);
/* function prototypes */
#endif

long jumb = 0, nonodes = 0;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], catfilename[FNMLNGTH], weightfilename[FNMLNGTH];
double *rrate;
long sites, weightsum, categs, datasets, ith, njumble, numtrees, shimotrees;
/*  sites = number of sites in actual sequences numtrees = number of user-defined trees */
long inseed, inseed0, mx, mx0, mx1;
boolean freqsfrom, global, global2=0, jumble, trout, usertree, weights, rctgry, ctgry, ttr, auto_, progress, mulsets, firstset, hypstate, smoothit, polishing, justwts, gama, invar;
boolean inserting  = false;
boolean lengthsopt = false;             /* Use lengths in user tree option */
boolean lngths     = false;             /* Actually use lengths (depends on
                                           each input tree) */
tree *curtree, *bestree, *bestree2, *priortree;
node *qwhere;
double *tymes;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr, freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates, cv, alpha, lambda, lambda1, invarfrac;
long *enterorder;
steptr aliasweight;
double *rate;
double **term, **slopeterm, **curveterm;
longer seed;
double *probcat;
contribarr *contribution;
char *progname;
long rcategs;
long **mp;
char basechar[16]="acmgrsvtwyhkdbn";

/* Local variables for maketree, propagated globally for C version: */
long    k, maxwhich, col;
double  like, bestyet, maxlogl;
boolean lastsp, smoothed, succeeded;
double  *l0gl;
double  x[3], lnl[3];
double  expon1i[maxcategs], expon1v[maxcategs],
        expon2i[maxcategs], expon2v[maxcategs];
double  **l0gf;
Char ch, ch2;


tree* dnamlk_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(dnamlk_tree));
  dnamlk_tree_init(t, nonodes, spp);
  return t;
}


void dnamlk_tree_init(tree *t, long nonodes, long spp)
{
  ml_tree_init(t, nonodes, spp);                // calls allocx
  t->insert_  = mlk_tree_insert_;
  t->try_insert_ = ml_tree_try_insert_;
  t->re_move = mlk_tree_re_move;
  t->evaluate = dnamlk_tree_evaluate;
  t->globrearrange = rooted_globrearrange;
  t->locrearrange = rooted_locrearrange;
  ((ml_tree*)t)->makenewv = (makenewv_t)mlk_tree_makenewv;
  t->nuview = dnamlk_tree_nuview;
  t->save_lr_nodes = rooted_tree_save_lr_nodes;
  t->restore_lr_nodes = rooted_tree_restore_lr_nodes;
}


void getoptions(void)
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean done;
  boolean didchangecat, didchangercat;
  double probsum;

  putchar('\n');
  auto_ = false;
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  freqsfrom = true;
  gama = false;
  invar = false;
  global = false;
  hypstate = false;
  jumble = false;
  njumble = 1;
  lambda = 1.0;
  lambda1 = 0.0;
  lengthsopt = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nNucleic acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
      printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    if (usertree)
    {
      printf("  L           Use lengths from user tree?");
      if (lengthsopt)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  T        Transition/transversion ratio:");
    if (!ttr)
      printf("  2.0\n");
    else
      printf("  %8.4f\n", ttratio);
    printf("  F       Use empirical base frequencies?");
    if (freqsfrom)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  C   One category of substitution rates?");
    if (!ctgry)
      printf("  Yes\n");
    else
      printf("  %ld categories\n", categs);
    printf("  R           Rate variation among sites?");
    if (!rctgry)
      printf("  constant rate\n");
    else
    {
      if (gama)
        printf("  Gamma distributed rates\n");
      else
      {
        if (invar)
          printf("  Gamma+Invariant sites\n");
        else
          printf("  user-defined HMM of rates\n");
      }
      printf("  A   Rates at adjacent sites correlated?");
      if (!auto_)
        printf("  No, they are independent\n");
      else
        printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    if (!usertree)
    {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    if (!usertree)
    {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed = %8ld, %3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
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
    printf("  5   Reconstruct hypothetical sequences?  %s\n",
           (hypstate ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done)
    {
      uppercase(&ch);
      if (((!usertree) && (strchr("JUCRAFWGTMI012345", ch) != NULL))
          || (usertree && ((strchr("UCRAFWLTMI012345", ch) != NULL))))
      {
        switch (ch)
        {
          case 'C':
            ctgry = !ctgry;
            if (ctgry)
            {
              printf("\nSitewise user-assigned categories:\n\n");
              initcatn(&categs);
              if (rate)
              {
                free(rate);
              }
              rate    = (double *) Malloc(categs * sizeof(double));
              didchangecat = true;
              initcategs(categs, rate);
            }
            break;

          case 'R':
            if (!rctgry)
            {
              rctgry = true;
              gama = true;
            }
            else
            {
              if (gama)
              {
                gama = false;
                invar = true;
              }
              else
              {
                if (invar)
                  invar = false;
                else
                  rctgry = false;
              }
            }
            break;

          case 'A':
            auto_ = !auto_;
            if (auto_)
            {
              initlambda(&lambda);
              lambda1 = 1.0 - lambda;
            }
            break;

          case 'F':
            freqsfrom = !freqsfrom;
            if (!freqsfrom)
              initfreqs(&freqa, &freqc, &freqg, &freqt);
            break;

          case 'G':
            global = !global;
            break;

          case 'W':
            weights = !weights;
            break;

          case 'J':
            jumble = !jumble;
            if (jumble)
              initjumble(&inseed, &inseed0, seed, &njumble);
            else njumble = 1;
            break;

          case 'L':
            lengthsopt = !lengthsopt;
            break;

          case 'T':
            ttr = !ttr;
            if (ttr)
              initratio(&ttratio);
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
                phyFillScreenColor();
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

          case '5':
            hypstate = !hypstate;
            break;
        }
      }
      else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
  if (gama || invar)
  {
    loopcount = 0;
    do {
      printf("\nCoefficient of variation of substitution rate among sites (must be positive)\n");
      printf(" In gamma distribution parameters, this is 1/(square root of alpha)\n");
      phyFillScreenColor();
      if(scanf("%lf%*[^\n]", &cv)) {}   // Read number and scan to EOL.
      (void)getchar();
      countup(&loopcount, 10);
    } while (cv <= 0.0);
    alpha = 1.0 / (cv * cv);
  }
  if (!rctgry)
    auto_ = false;
  if (rctgry)
  {
    printf("\nRates in HMM");
    if (invar)
      printf(" (including one for invariant sites)");
    printf(":\n");
    initcatn(&rcategs);
    if (probcat)
    {
      free(probcat);
      free(rrate);
    }
    probcat = (double *) Malloc(rcategs * sizeof(double));
    rrate   = (double *) Malloc(rcategs * sizeof(double));
    didchangercat = true;
    if (gama)
      initgammacat(rcategs, alpha, rrate, probcat);
    else
    {
      if (invar)
      {
        loopcount = 0;
        do {
          printf("Fraction of invariant sites?\n");
          if(scanf("%lf%*[^\n]", &invarfrac)) {} // Read number and scan to EOL.
          (void)getchar();
          countup(&loopcount, 10);
        } while ((invarfrac <= 0.0) || (invarfrac >= 1.0));
        initgammacat(rcategs-1, alpha, rrate, probcat);
        for (i = 0; i < rcategs-1; i++)
          probcat[i] = probcat[i]*(1.0-invarfrac);
        probcat[rcategs-1] = invarfrac;
        rrate[rcategs-1] = 0.0;
      }
      else
      {
        initcategs(rcategs, rrate);
        initprobcat(rcategs, &probsum, probcat);
      }
    }
  }
  if (!didchangercat)
  {
    rrate      = Malloc(rcategs * sizeof(double));
    probcat    = Malloc(rcategs * sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat)
  {
    rate       = Malloc(categs * sizeof(double));
    rate[0]    = 1.0;
  }
}  /* getoptions */


void reallocsites(void)
{
  long i;

  for (i = 0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (char *)Malloc(sites * sizeof(char));
  }
  free(weight);
  free(category);
  free(alias);
  free(aliasweight);
  free(ally);
  free(location);

  weight      = (long *)Malloc(sites * sizeof(long));
  category    = (long *)Malloc(sites * sizeof(long));
  alias       = (long *)Malloc(sites * sizeof(long));
  aliasweight = (long *)Malloc(sites * sizeof(long));
  ally        = (long *)Malloc(sites * sizeof(long));
  location    = (long *)Malloc(sites * sizeof(long));
}


void allocrest(void)
{
  long i;

  inputSequences     = (Char **)Malloc(spp * sizeof(Char *));
  nayme  = (naym *)Malloc(spp * sizeof(naym));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (char *)Malloc(sites * sizeof(char));
  enterorder  = (long *)Malloc(spp * sizeof(long));
  weight      = (long *)Malloc(sites * sizeof(long));
  category    = (long *)Malloc(sites * sizeof(long));
  alias       = (long *)Malloc(sites * sizeof(long));
  aliasweight = (long *)Malloc(sites * sizeof(long));
  ally        = (long *)Malloc(sites * sizeof(long));
  location    = (long *)Malloc(sites * sizeof(long));
  tymes       = (double *)Malloc((nonodes - spp) * sizeof(double));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nNucleic acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n", VERSION);

  inputnumbers(&spp, &sites, &nonodes, 1);
  if (!javarun)
  {
    getoptions();
  }
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  allocrest();
}  /* doinit */


void inputoptions(void)
{
  long i;

  if (!firstset && !justwts)
  {
    samenumsp(&sites, ith);
    reallocsites();
  }

  for (i = 0; i < sites; i++)
    category[i] = 1;
  for (i = 0; i < sites; i++)
    weight[i] = 1;

  if (justwts || weights)
    inputweights(sites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if (ctgry && categs > 1)
  {
    inputcategs(0, sites, category, categs, "DnaMLK");
    if (printdata)
      printcategs(outfile, sites, category, "Site categories");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, weight, "Sites");
}  /* inputoptions */


void makeweights(void)
{
  /* make up bookkeeping vectors to avoid duplicate computations.
     for explanation of variables  alias, aliasweight, and ally see seq.c
     routines  sitesort2, sitecombine2, and sitescrunch2
     location[i-1] ends up as the position in the lexicographic order
     (lexicographic by site rate category and site pattern) of site i
     if site i is one whose category and site pattern represents the
     others that are tied with it.  Otherwise location[i-1] is 0
     endsite is the number (less 1) of such representative sites  */
  long i;

  for (i = 1; i <= sites; i++)
  {
    alias[i - 1] = i;
    ally[i - 1] = i;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }
  sitesort2(sites, aliasweight);
  sitecombine2(sites, aliasweight);
  sitescrunch2(sites, 1, 2, aliasweight);
  endsite = 0;
  for (i = 1; i <= sites; i++)
  {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
  {
    location[alias[i - 1] - 1] = i;
  }
  contribution = (contribarr *) Malloc(endsite * sizeof(contribarr));
}  /* makeweights */


void getinput(void)
{
  /* reads the input data */
  inputoptions();
  if (!freqsfrom)
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                 &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                 freqsfrom, true);
  if (!justwts || firstset)
    inputdata(sites);
  makeweights();
  inittrees(&curtree, &bestree, &priortree, &bestree2, nonodes, spp);
  makevalues2(rcategs, curtree->nodep, endsite, spp, inputSequences, alias);
  if (freqsfrom)
  {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight, curtree->nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy, &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange, freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


void inittable_for_usertree (FILE *intree)
{
  /* If there's a user tree, then the ww/zz/wwzz/vvzz elements need
     to be allocated appropriately. */
  long num_comma;
  long i, j;

  /* First, figure out the largest possible furcation, i.e. the number
     of commas plus one */
  countcomma(intree, &num_comma);
  num_comma++;

  for (i = 0; i < rcategs; i++)
  {
    for (j = 0; j < categs; j++)
    {
      /* Free the stuff allocated assuming bifurcations */
      free (tbl[i][j]->ww);
      free (tbl[i][j]->zz);
      free (tbl[i][j]->wwzz);
      free (tbl[i][j]->vvzz);

      /* Then allocate for worst-case multifurcations */
      tbl[i][j]->ww   = (double *) Malloc(num_comma * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc(num_comma * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc(num_comma * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc(num_comma * sizeof (double));
    }
  }
}  /* inittable_for_usertree */


void freetable(void)
{
  long i, j;

  for (i = 0; i < rcategs; i++)
  {
    for (j = 0; j < categs; j++)
    {
      free(tbl[i][j]->ww);
      free(tbl[i][j]->zz);
      free(tbl[i][j]->wwzz);
      free(tbl[i][j]->vvzz);
    }
  }
  for (i = 0; i < rcategs; i++)
  {
    for (j = 0; j < categs; j++)
      free(tbl[i][j]);
    free(tbl[i]);
  }
  free(tbl);
}


void inittable(void)
{
  /* Define a lookup table. Precompute values and print them out in tables */
  long i, j;
  double sumrates;

  tbl = (valrec ***) Malloc(rcategs * sizeof(valrec **));
  for (i = 0; i < rcategs; i++)
  {
    tbl[i] = (valrec **) Malloc(categs * sizeof(valrec *));
    for (j = 0; j < categs; j++)
      tbl[i][j] = (valrec *) Malloc(sizeof(valrec));
  }

  for (i = 0; i < rcategs; i++)
  {
    for (j = 0; j < categs; j++)
    {
      tbl[i][j]->rat = rrate[i]*rate[j];
      tbl[i][j]->ratxi = tbl[i][j]->rat * xi;
      tbl[i][j]->ratxv = tbl[i][j]->rat * xv;

      /* Allocate assuming bifurcations, will be changed later if
         neccesarry (i.e. there's a user tree) */
      tbl[i][j]->ww   = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc(2 * sizeof (double));
    }
  }
  sumrates = 0.0;
  for (i = 0; i < endsite; i++)
  {
    for (j = 0; j < rcategs; j++)
      sumrates += aliasweight[i] * probcat[j]
        * tbl[j][category[alias[i] - 1] - 1]->rat;
  }
  sumrates /= (double)sites;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++)
    {
      tbl[i][j]->rat /= sumrates;
      tbl[i][j]->ratxi /= sumrates;
      tbl[i][j]->ratxv /= sumrates;
    }

  if(jumb > 1)
    return;

  if (gama || invar)
  {
    fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
    fprintf(outfile, " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
  }
  if (rcategs > 1)
  {
    fprintf(outfile, "\nState in HMM    Rate of change    Probability\n\n");
    for (i = 0; i < rcategs; i++)
      if (probcat[i] < 0.0001)
        fprintf(outfile, "%9ld%16.3f%20.6f\n", i+1, rrate[i], probcat[i]);
      else if (probcat[i] < 0.001)
        fprintf(outfile, "%9ld%16.3f%19.5f\n", i+1, rrate[i], probcat[i]);
      else if (probcat[i] < 0.01)
        fprintf(outfile, "%9ld%16.3f%18.4f\n", i+1, rrate[i], probcat[i]);
      else
        fprintf(outfile, "%9ld%16.3f%17.3f\n", i+1, rrate[i], probcat[i]);
    putc('\n', outfile);
    if (auto_)
    {
      fprintf(outfile, "Expected length of a patch of sites having the same rate = %8.3f\n", 1/lambda);
      putc('\n', outfile);
    }
  }
  if (categs > 1)
  {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 0; i < categs; i++)
      fprintf(outfile, "%9ld%16.3f\n", i+1, rate[i]);
    fprintf(outfile, "\n\n");
  }
}  /* inittable */


void alloc_nvd(long num_sibs, nuview_data *local_nvd)
{
  /* Allocate blocks of memory appropriate for the number of siblings
     a given node has */
  local_nvd->yy     = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->wwzz   = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->vvzz   = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->vzsumr = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->vzsumy = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->sum    = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->sumr   = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->sumy   = (double *) Malloc(num_sibs * sizeof (double));
  local_nvd->xx     = (sitelike *) Malloc(num_sibs * sizeof (sitelike));
}  /* alloc_nvd */


void free_nvd(nuview_data *local_nvd)
{
  /* The natural complement to the alloc version */
  free (local_nvd->yy);
  free (local_nvd->wwzz);
  free (local_nvd->vvzz);
  free (local_nvd->vzsumr);
  free (local_nvd->vzsumy);
  free (local_nvd->sum);
  free (local_nvd->sumr);
  free (local_nvd->sumy);
  free (local_nvd->xx);
}  /* free_nvd */


void dnamlk_tree_nuview(tree *t, node *p)
/* current (modified dnaml) nuview */
{
  long i, j, k, l, num_sibs, sib_index;
  nuview_data *local_nvd;
  node *sib_ptr, *sib_back_ptr;
  sitelike p_xx;
  double lw;
  double correction;
  double maxx;

  num_sibs    = count_sibs (p);
  generic_tree_nuview(t, p);

  /* Allocate the structure and blocks therein for variables used in this function. */
  local_nvd = (nuview_data *) Malloc(sizeof (nuview_data));
  alloc_nvd (num_sibs, local_nvd);

  /* Loop 1: makes assignments to tbl based on some combination of
     what's already in tbl and the children's value of v */
  sib_ptr = p;
  for (sib_index=0; sib_index < num_sibs; sib_index++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (sib_back_ptr != NULL)
      lw = -fabs( get_tyme(p) - get_tyme(sib_back_ptr) );
    else
      lw = 0.0;

    for (i = 0; i < rcategs; i++)
      for (j = 0; j < categs; j++)
      {
        tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
        tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
        tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] * tbl[i][j]->zz[sib_index];
        tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
          tbl[i][j]->zz[sib_index];
      }
  }

  /* Loop 2: */
  for (i = 0; i < endsite; i++)
  {
    correction = 0;
    maxx = 0;
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      /* Loop 2.1 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        sib_ptr         = sib_ptr->next;
        sib_back_ptr    = sib_ptr->back;

        local_nvd->wwzz[sib_index] = tbl[j][k]->wwzz[sib_index];
        local_nvd->vvzz[sib_index] = tbl[j][k]->vvzz[sib_index];
        local_nvd->yy[sib_index]   = 1.0 - tbl[j][k]->zz[sib_index];
        if (sib_back_ptr != NULL)
        {
          memcpy(local_nvd->xx[sib_index],
                 ((dna_node*)sib_back_ptr)->x[i][j],
                 sizeof(sitelike));
          if ( j == 0)
            correction += ((ml_node*)sib_back_ptr)->underflows[i];
        }
        else
        {
          local_nvd->xx[sib_index][0] = 1.0;
          local_nvd->xx[sib_index][(long)C - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)G - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)T - (long)A] = 1.0;
        }
      }

      /* Loop 2.2 */
      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        local_nvd->sum[sib_index] =
          local_nvd->yy[sib_index] *
          (freqa * local_nvd->xx[sib_index][(long)A] +
           freqc * local_nvd->xx[sib_index][(long)C] +
           freqg * local_nvd->xx[sib_index][(long)G] +
           freqt * local_nvd->xx[sib_index][(long)T]);
        local_nvd->sumr[sib_index] =
          freqar * local_nvd->xx[sib_index][(long)A] +
          freqgr * local_nvd->xx[sib_index][(long)G];
        local_nvd->sumy[sib_index] =
          freqcy * local_nvd->xx[sib_index][(long)C] +
          freqty * local_nvd->xx[sib_index][(long)T];
        local_nvd->vzsumr[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumr[sib_index];
        local_nvd->vzsumy[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumy[sib_index];
      }

      /* Initialize to one, multiply incremental values for every
         sibling a node has */
      p_xx[(long)A] = 1 ;
      p_xx[(long)C] = 1 ;
      p_xx[(long)G] = 1 ;
      p_xx[(long)T] = 1 ;

      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        p_xx[(long)A] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)A] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)C] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)C] +
          local_nvd->vzsumy[sib_index];
        p_xx[(long)G] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)G] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)T] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)T] +
          local_nvd->vzsumy[sib_index];
      }

      for ( l = 0 ; l < ((long)T - (long)A + 1 ) ; l++ )
      {
        if (  p_xx[l] > maxx )
          maxx = p_xx[l];
      }

      /* And the final point of this whole function: */
      memcpy(((dna_node*)p)->x[i][j], p_xx, sizeof(sitelike));
    }
    ((ml_node*)p)->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE)
      fix_x(((dna_node*)p), i, maxx, rcategs);
    ((ml_node*)p)->underflows[i] += correction;

  }

  p->initialized = true;

  free_nvd (local_nvd);
  free (local_nvd);

}  /* dnamlk_tree_nuview */


double dnamlk_tree_evaluate(tree* t, node *p, boolean dummy)
{
  contribarr tterm;
  static contribarr like, nulike, clai;
  double sum, sum2, sumc=0, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3, sumterm, lterm;
  long i, j, k, lai;
  node *q, *r;
  sitelike x1, x2;
  sum = 0.0;

  (void)dummy;                                  // RSGnote: Parameter not used (probably correct here).

  if (p == t->root && (count_sibs(p) == 2))
  {
    r = p->next->back;
    q = p->next->next->back;
    y = get_tyme(r) + get_tyme(q) - 2 * get_tyme(p);
    if (!r->tip && !r->initialized) t->nuview (t, r);
    if (!q->tip && !q->initialized) t->nuview (t, q);
  }
  else if (p == t->root)
  {
    /* FIXME The following comment does not make sense because this is not
     * an internal node. */

    /* the next two lines copy tyme and x to p->next.  Normally they are
       not initialized for an internal node. */
    /* assumes bifurcation */
    /* FIXME This block is executed when root multifurcates! */

    set_tyme(p->next, get_tyme(p));
    t->nuview(t, p->next);
    r = p->next;
    q = p->next->back;
    y = fabs( get_tyme(p->next) - get_tyme(q) );
  }
  else
  {
    r = p;
    q = p->back;
    /* TODO Are tests for tip and initialized done in nuview() already? */
    if (!r->tip && !r->initialized) t->nuview (t, r);
    if (!q->tip && !q->initialized) t->nuview (t, q);
    y = fabs( get_tyme(r) - get_tyme(q) );
  }

  lz = -y;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++)
    {
      tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
      tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
      tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
      tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
    }
  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      if (y > 0.0)
      {
        y1 = 1.0 - tbl[j][k]->z1;
        z1zz = tbl[j][k]->z1zz;
        z1yy = tbl[j][k]->z1yy;
      }
      else
      {
        y1 = 0.0;
        z1zz = 1.0;
        z1yy = 0.0;
      }
      memcpy(x1, ((dna_node*)r)->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] + freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, ((dna_node*)q)->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] + freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) * (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) + (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) * (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] + freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] + freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] + freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, sizeof(contribarr));
    if (!auto_ && usertree && (which <= shimotrees))
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  if (auto_)
  {
    for (j = 0; j < rcategs; j++)
      like[j] = 1.0;
    for (i = 0; i < sites; i++)
    {
      sumc = 0.0;
      for (k = 0; k < rcategs; k++)
        sumc += probcat[k] * like[k];
      sumc *= lambda;
      if ((ally[i] > 0) && (location[ally[i]-1] > 0))
      {
        lai = location[ally[i] - 1];
        memcpy(clai, contribution[lai - 1], sizeof(contribarr));
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
      }
      else
      {
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc);
      }
      memcpy(like, nulike, sizeof(contribarr));
    }
    sum2 = 0.0;
    for (i = 0; i < rcategs; i++)
      sum2 += probcat[i] * like[i];
    sum += log(sum2);
  }

  t->score = sum;

  if (auto_ || !usertree)
    return sum;
  if(which <= shimotrees)
    l0gl[which - 1] = sum;
  if (which == 1)
  {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }

  if (sum > maxlogl)
  {
    maxwhich = which;
    maxlogl = sum;
  }
  return sum;
}  /* dnamlk_tree_evaluate */


void initdnamlknode(tree *treep, node **p, long len, long nodei, long *ntips, long *parens, initops whichinit, pointarray nodep, Char *str, Char *ch, FILE *intree)
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
      ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs);
      nodep[(*p)->index - 1] = (*p);
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs);
      (*p)->index = nodei;
      break;
    case tip:
      match_names_to_data (str, nodep, p, spp);
      break;
    case iter:
      (*p)->initialized = false;
      (*p)->v = initialv;
      (*p)->iter = true;
      if ((*p)->back != NULL)
        (*p)->back->iter = true;
      break;
    case length:
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      (*p)->v = valyew / divisor / fracchange;
      (*p)->iter = false;
      if ((*p)->back != NULL)
      {
        (*p)->back->v = (*p)->v;
        (*p)->back->iter = false;
      }
      break;
    case hslength:
      break;
    case hsnolength:
      if (usertree && lengthsopt && lngths)
      {
        printf("Warning: one or more lengths not defined in user tree number %ld.\n", which);
        printf("DNAMLK will attempt to optimize all branch lengths.\n\n");
        lngths = false;
      }
      break;
    case treewt:
      break;
    case unittrwt:
      break;
  }
} /* initdnamlknode */


void tymetrav(node *p, double *x)
{
  /* set up times of nodes */
  node *sib_ptr, *q;
  long i, num_sibs;
  double xmax;

  xmax = 0.0;
  if (!p->tip)
  {
    sib_ptr  = p;
    num_sibs = count_sibs(p);
    for (i=0; i < num_sibs; i++)
    {
      sib_ptr = sib_ptr->next;
      tymetrav(sib_ptr->back, x);
      if (xmax > (*x))
        xmax = (*x);
    }
  }
  else
    (*x)     = 0.0;
  set_tyme(p, xmax);
  if (!p->tip)
  {
    q = p;
    while (q->next != p)
    {
      q = q->next;
      set_tyme(q, get_tyme(p));
    }
  }
  (*x) = get_tyme(p) - p->v ;
}  /* tymetrav */


void dnamlk_coordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last, *pp1 =NULL, *pp2 =NULL;
  long num_sibs, p1, p2, i;

  if (p->tip)
  {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    dnamlk_coordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  num_sibs = count_sibs(p);
  p1 = (long)((num_sibs+1)/2.0);
  p2 = (long)((num_sibs+2)/2.0);
  i = 1;
  q = p->next;
  first  = q->back;
  do {
    if (i == p1) pp1 = q->back;
    if (i == p2) pp2 = q->back;
    last = q->back;
    q = q->next;
    i++;
  } while (q != p);
  p->xcoord = (long)(0.5 - over * get_tyme(p));
  p->ycoord = (pp1->ycoord + pp2->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* dnamlk_coordinates */


void dnamlk_drawline(long i, double scale)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = curtree->root;
  q = curtree->root;
  extra = false;
  if ((long)(p->ycoord) == i)
  {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
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
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * ((long)(p->xcoord) - (long)(q->xcoord)) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)(q->ycoord) == i && !done)
    {
      if (p->ycoord != q->ycoord)
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
      if ((long)(last->ycoord) > i && (long)(first->ycoord) < i &&
          i != (long)(p->ycoord))
      {
        putc('!', outfile);
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
    if (p != q)
      p = q;
  } while (!done);
  if ((long)(p->ycoord) == i && p->tip)
  {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* dnamlk_drawline */


void dnamlk_printree(void)
{
  /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;
  node *p;

  putc('\n', outfile);
  tipy = 1;
  dnamlk_coordinates(curtree->root, &tipy);
  p = curtree->root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)(get_tyme(p) - get_tyme(curtree->root) + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= tipy - down; i++)
    dnamlk_drawline(i, scale);
  putc('\n', outfile);
}  /* dnamlk_printree */


void describe(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;
  double v;

  if (p == curtree->root)
    fprintf(outfile, " root         ");
  else
    fprintf(outfile, "%4ld          ", p->back->index - spp);
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  if (p != curtree->root)
  {
    fprintf(outfile, "%11.5f", fracchange * (get_tyme(p) - get_tyme(curtree->root)));
    v = fracchange * ( get_tyme(p) - get_tyme(curtree->nodep[p->back->index - 1]));
    fprintf(outfile, "%13.5f", v);
  }
  putc('\n', outfile);
  if (!p->tip)
  {

    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++)
    {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      describe(sib_back_ptr);
    }
  }
}  /* describe */


void reconstr(node *p, long n)
{
  /* reconstruct and print out base at site n+1 at node p */
  long i, j, k, m, first, second, num_sibs;
  double f, sum, xx[4];
  node *q;

  j = location[ally[n]-1] - 1;
  for (i = 0; i < 4; i++)   /* do for each nucleotide */
  {
    if (p == curtree->root)    /* ... as  x  has not been computed there */
      f = 1.0;
    else
      f = ((dna_node*)p)->x[j][mx-1][i];
    num_sibs = count_sibs(p);
    q = p;
    for (k = 0; k < num_sibs; k++)
    {
      q = q->next;
      f *= ((dna_node*)q)->x[j][mx-1][i];
    }
    if (f > 0.0)
      f = exp(log(f)/(num_sibs-1.0));
    xx[i] = f;
  }
  xx[0] *= freqa;
  xx[1] *= freqc;
  xx[2] *= freqg;
  xx[3] *= freqt;
  sum = xx[0]+xx[1]+xx[2]+xx[3];
  for (i = 0; i < 4; i++)
    xx[i] /= sum;
  first = 0;
  for (i = 1; i < 4; i++)
    if (xx [i] > xx[first])
      first = i;
  if (first == 0)
    second = 1;
  else
    second = 0;
  for (i = 0; i < 4; i++)
    if ((i != first) && (xx[i] > xx[second]))
      second = i;
  m = 1 << first;
  if (xx[first] < 0.4999995)
    m = m + (1 << second);
  if (xx[first] > 0.95)
    putc(toupper(basechar[m - 1]), outfile);
  else
    putc(basechar[m - 1], outfile);
  if (rctgry && rcategs > 1)
    mx = mp[n][mx - 1];
  else
    mx = 1;
} /* reconstr */


void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  long num_sibs, i;
  node *sib_ptr;

  putc(' ', outfile);
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "  ");
  mx = mx0;
  for (i = m; i <= n; i++)
  {
    if ((i % 10 == 0) && (i != m))
      putc(' ', outfile);
    if (p->tip)
      putc(inputSequences[p->index-1][i], outfile);
    else
      reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip)
  {
    num_sibs = count_sibs(p);
    sib_ptr = p;
    for (i = 0; i < num_sibs; i++)
    {
      sib_ptr = sib_ptr->next;
      rectrav(sib_ptr->back, m, n);
    }
  }
  mx1 = mx;
}  /* rectrav */


void summarize(void)
{
  long i, j, mm;
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;

  mp = (long **)Malloc(sites * sizeof(long *));
  for (i = 0; i <= sites-1; ++i)
    mp[i] = (long *)Malloc(sizeof(long)*rcategs);
  fprintf(outfile, "\nLn Likelihood = %11.5f\n\n", curtree->score);
  fprintf(outfile, " Ancestor      Node      Node Height     Length\n");
  fprintf(outfile, " --------      ----      ---- ------     ------\n");
  describe(curtree->root);
  putc('\n', outfile);
  if (rctgry && rcategs > 1)
  {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        mp[i][j] = j + 1;
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
          {
            if (lambda * probcat[k - 1] * like[k - 1] > nulike[j])
            {
              nulike[j] = lambda * probcat[k - 1] * like[k - 1];
              mp[i][j] = k;
            }
          }
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= rcategs; i++)
    {
      if (probcat[i - 1] * like[i - 1] > mode)
      {
        mx = i;
        mode = probcat[i - 1] * like[i - 1];
      }
    }
    mx0 = mx;
    fprintf(outfile, "Combination of categories that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 1; i <= sites; i++)
    {
      fprintf(outfile, "%ld", mx);
      if (i % 10 == 0)
        putc(' ', outfile);
      if (i % 60 == 0 && i != sites)
      {
        putc('\n', outfile);
        for (j = 1; j <= nmlngth + 3; j++)
          putc(' ', outfile);
      }
      mx = mp[i - 1][mx - 1];
    }
    fprintf(outfile, "\n\n");
    marginal = (double **) Malloc(sites * sizeof(double *));
    for (i = 0; i < sites; i++)
      marginal[i] = (double *) Malloc(rcategs * sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
            nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] /= sum;
        marginal[i][j] = nulike[j];
      }
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = 0; i < sites; i++)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
            nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        marginal[i][j] *= like[j] * probcat[j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        sum += marginal[i][j];
      for (j = 0; j < rcategs; j++)
        marginal[i][j] /= sum;
    }
    fprintf(outfile, "Most probable category at each site if > 0.95 probability (\".\" otherwise)\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < sites; i++)
    {
      sum = 0.0;
      mm = 0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum)
        {
          sum = marginal[i][j];
          mm = j;
        }
      if (sum >= 0.95)
        fprintf(outfile, "%ld", mm+1);
      else
        putc('.', outfile);
      if ((i+1) % 60 == 0)
      {
        if (i != 0)
        {
          putc('\n', outfile);
          for (j = 1; j <= nmlngth + 3; j++)
            putc(' ', outfile);
        }
      }
      else if ((i+1) % 10 == 0)
        putc(' ', outfile);
    }
    putc('\n', outfile);
    for (i = 0; i < sites; i++)
      free(marginal[i]);
    free(marginal);
  }
  putc('\n', outfile);
  putc('\n', outfile);
  if (hypstate)
  {
    fprintf(outfile, "Probable sequences at interior nodes:\n\n");
    fprintf(outfile, "  node       ");
    for (i = 0; (i < 13) && (i < ((sites + (sites-1)/10 - 39) / 2)); i++)
      putc(' ', outfile);
    fprintf(outfile, "Reconstructed sequence (caps if > 0.95)\n\n");
    if (!rctgry || (rcategs == 1))
      mx0 = 1;
    for (i = 0; i < sites; i += 60)
    {
      k = i + 59;
      if (k >= sites)
        k = sites - 1;
      rectrav(curtree->root, i, k);
      putc('\n', outfile);
      mx0 = mx1;
    }
  }
  for (i = 0; i < sites; ++i)
    free(mp[i]);
  free(mp);
}  /* summarize */


void dnamlk_treeout(node *p)
{
  /* write out file with representation of final tree */
  node *sib_ptr;
  long i, n, w, num_sibs;
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
    sib_ptr = p;
    num_sibs = count_sibs(p);
    putc('(', outtree);
    col++;

    for (i=0; i < (num_sibs - 1); i++)
    {
      sib_ptr = sib_ptr->next;
      dnamlk_treeout(sib_ptr->back);
      putc(',', outtree);
      col++;
      if (col > 55)
      {
        putc('\n', outtree);
        col = 0;
      }
    }
    sib_ptr = sib_ptr->next;
    dnamlk_treeout(sib_ptr->back);
    putc(')', outtree);
    col++;
  }
  if (p == curtree->root)
  {
    fprintf(outtree, ";\n");
    return;
  }
  x = fracchange * ( get_tyme(p) - get_tyme(curtree->nodep[p->back->index - 1]));
  if (x > 0.0)
    w = (long)(0.4342944822 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.4342944822 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(outtree, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* dnamlk_treeout */


void nodeinit(node *p)
{
  /* set up times at one node */
  node *sib_ptr, *sib_back_ptr;
  long i, num_sibs;
  double lowertyme;

  sib_ptr = p;
  num_sibs = count_sibs(p);

  /* lowertyme = lowest of children's times */
  lowertyme = get_tyme(p->next->back);
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    if ( get_tyme(sib_back_ptr) < lowertyme)
      lowertyme = get_tyme(sib_back_ptr);
  }

  set_tyme(p, lowertyme - initialv);

  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    set_tyme(sib_ptr, get_tyme(p));
    sib_back_ptr->v = ( get_tyme(sib_back_ptr) - get_tyme(p) );
    sib_ptr->v = sib_back_ptr->v;
  }
}  /* nodeinit */


void initrav(node *p)
{

  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up times throughout tree */
  if (p->tip)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    initrav(sib_back_ptr);
  }

  nodeinit(p);
}  /* initrav */


void travinit(tree* t, node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up initial values */
  if (p == NULL)
    return;
  if (p->tip)
    return;
  if (p->initialized)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    travinit(t, sib_back_ptr);
  }

  t->nuview(t, p);
  p->initialized = true;
}  /* travinit */


void travsp(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to find tips */
  if (p == curtree->root)
    travinit(curtree, p);
  if (p->tip)
    travinit(curtree, p->back);
  else
  {
    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++)
    {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      travsp(sib_back_ptr);
    }
  }
}  /* travsp */


void treevaluate(void)
{
  /* evaluate likelihood of tree, after iterating branch lengths */
  long i, j,  num_sibs;
  node *sib_ptr;

  polishing = true;
  smoothit = true;
  if (lngths == 0 && usertree == 1)
  {
    for (i = 0; i < spp; i++)
      curtree->nodep[i]->initialized = false;
    for (i = spp; i < nonodes; i++)
    {
      sib_ptr = curtree->nodep[i];
      sib_ptr->initialized = false;
      num_sibs = count_sibs(sib_ptr);
      for (j=0 ; j < num_sibs; j++)
      {
        sib_ptr      = sib_ptr->next;
        sib_ptr->initialized = false;
      }
    }
    initrav(curtree->root);
    travsp(curtree->root);
  }

  if ( usertree && !lngths )
    for ( i = 0 ; i < smoothings ; i++ ) /* this converges slowly, there     */
      curtree->smoothall(curtree, curtree->root); /* already is a loop within smoothall, we need more! */
  else
    curtree->smoothall(curtree, curtree->root); /* we should already be close */

  curtree->evaluate(curtree, curtree->root, false);

}  /* treevaluate */


void maketree(void)
{
  /* constructs a binary tree from the pointers in curtree.nodep,
     adds each node at location which yields highest likelihood
     then rearranges the tree for greatest likelihood */

  long i, j, numtrees;
  double x;
  node *there, *item, *dummy, *q, *root=NULL;
  boolean dummy_haslengths, dummy_first, goteof, multf=false;
  long nextnode;

  inittable();

  if (!usertree)
  {
    /* enter species in order by default */
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder); /* shuffle enterorder[] */

    /* initialize tree */
    /* root is first to be popped from free_forkrings during first insert */
    curtree->root = curtree->nodep[enterorder[0]-1];
    /* FIXME is this necessary? */
    curtree->root = curtree->nodep[curtree->root->index - 1];

    /* FIXME back links should be zeroed already during remove */
    /* delete back links for tips */
    for (i = 0; i < spp; i++)
      curtree->nodep[i]->back = NULL;
    /* delete back links for all FORKNODES */
    for (i = spp; i < nonodes; i++)
    {
      q = curtree->nodep[i];
      q->back = NULL;
      while ((q = q->next) != curtree->nodep[i])
        q->back = NULL;
    }

    /* add first two tips */
    polishing = false;

    /* Inserting one tip at another creates a new fork between them, which here becomes the new root. */
    curtree->insert_(curtree, curtree->nodep[enterorder[1]-1], curtree->nodep[enterorder[0] - 1], false, false);

    if (progress)
    {
      sprintf(progbuf, "\nAdding species:\n");
      print_progress(progbuf);
      writename(0, 2, enterorder);
      phyFillScreenColor();
    }
    lastsp = false;
    smoothit = false;

    /* add remaining tips */
    for (i = 3; i <= spp; i++)
    {
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
      there = curtree->root;
      item = curtree->nodep[enterorder[i - 1] - 1];
      lastsp = (i == spp);
      curtree->addtraverse(curtree, item, curtree->root, true, &there, &bestyet, bestree, priortree, false, &multf);
      curtree->insert_(curtree, item, there, true, multf);
      curtree->smoothall(curtree, curtree->root);/* do we need this */
      curtree->locrearrange(curtree, curtree->root, false, priortree, bestree);
      if (progress)
      {
        writename(i - 1, 1, enterorder);
        phyFillScreenColor();
      }
    }
    if (global)
    {
      if (progress)
      {
        sprintf(progbuf, "Doing global rearrangements\n");
        print_progress(progbuf);
        sprintf(progbuf, "  !");
        print_progress(progbuf);
        for (j = 0; j < nonodes; j++)
        {
          if ( (j-spp) % (( nonodes / 72 ) + 1 ) == 0 )
          {
            sprintf(progbuf, "-");
            print_progress(progbuf);
          }
        }
        sprintf(progbuf, "!\n");
        print_progress(progbuf);
      }
      smoothit = true;
      curtree->globrearrange(curtree, progress, true);
      smoothit = false;
      curtree->copy(curtree, bestree);
    }
    curtree->copy(curtree, bestree);

    if (njumble > 1)                                    // RSGbugfix
    {
      // Release all FORKNODES by removing all but last tip.
      for (i = 0; i < spp - 1; i++ )
        curtree->re_move(curtree, curtree->nodep[i], &dummy, false);

      // Since root FORKNODE was removed, nullify ROOT (will be reset in upcoming COPY).
      curtree->root = NULL;

      if (jumb == 1 || bestree2->score < bestree->score)
        bestree->copy(bestree, bestree2);
    }

    if (jumb == njumble)
    {
      if (njumble > 1)
        bestree2->copy(bestree2, curtree);
      else
        bestree->copy(bestree, curtree);
      fprintf(outfile, "\n\n");
      treevaluate();
      curtree->score = curtree->evaluate(curtree, curtree->root, false);

      if (treeprint)
      {
        dnamlk_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        dnamlk_treeout(curtree->root);
      }
    }
  }
  else
  { /* usertree */
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    inittable_for_usertree (intree);
    numtrees = countsemic(intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *)Malloc(shimotrees * sizeof(double));
    l0gf = (double **)Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    fprintf(outfile, "\n\n");
    which = 1;
    while (which <= numtrees)
    {

      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      dummy_haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      lngths           = lengthsopt;

      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof, &dummy_first, &nextnode, &dummy_haslengths, initdnamlknode, false, nonodes);
      nonodes = nextnode;

      if (lngths)
        tymetrav(curtree->root, &x);

      if (goteof && (which <= numtrees))
      {
        /* if we hit the end of the file prematurely */
        printf ("\nERROR:  Trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit(-1);
      }
      treevaluate();
      if (treeprint)
      {
        dnamlk_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        dnamlk_treeout(curtree->root);
      }
      if(which < numtrees)
      {
        freex_notip(nonodes, curtree->nodep);
      }
      which++;
    }

    FClose(intree);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite, maxlogl, l0gl, l0gf,
               aliasweight, seed);
  }

  if (jumb == njumble)
  {
    if (progress)
    {
      sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
      print_progress(progbuf);
      if (trout)
      {
        sprintf(progbuf, "Tree also written onto file \"%s\".\n\n", outtreename);
        print_progress(progbuf);
      }
    }
    free(contribution);
    freex(nonodes, curtree->nodep);
    if (!usertree)
    {
      freex(nonodes, bestree->nodep);
      if (njumble > 1)
        freex(nonodes, bestree2->nodep);
    }
  }
  free(root);
  freetable();
} /* maketree */


/*?? Dnaml has a clean-up function for freeing memory, closing files, etc.
     Put one here too? */


void dnamlkrun(void)
{
  /*
  // debug printout // JRMdebug
  printf("\nctgry: %i\n", ctgry);
  printf("categs: %li\n", categs);
  printf("rctgry: %i\n", rctgry);
  printf("rcategs: %li\n", rcategs);
  printf("auto_: %i\n", auto_);
  printf("freqsfrom: %i\n", freqsfrom);
  printf("gama: %i\n", gama);
  printf("global: %i\n", global);
  printf("hypstate: %i\n", hypstate);
  printf("invar: %i\n", invar);
  printf("jumble: %i\n", jumble);
  printf("njumble: %li\n", njumble);
  printf("lngths: %i\n", lngths);
  printf("lambda: %f\n", lambda);
  printf("outgrno: %li\n", outgrno);
  printf("outgropt: %i\n", outgropt);
  printf("trout: %i\n", trout);
  printf("ttratio: %f\n", ttratio);
  printf("ttr: %i\n", ttr);
  printf("usertree: %i\n", usertree);
  printf("weights: %i\n", weights);
  printf("printdata: %i\n", printdata);
  printf("dotdiff: %i\n", dotdiff);
  printf("progress: %i\n", progress);
  printf("treeprint: %i\n", treeprint);
  printf("interleaved: %i\n", interleaved);
  printf("mulsets: %i\n", mulsets);
  printf("datasets: %li\n", datasets);
  */

  // do the work
  for (ith = 1; ith <= datasets; ith++)
  {
    ttratio = ttratio0;
    if (datasets > 1)
    {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    fflush(outfile);
    fflush(outtree);
  }
}


void dnamlk(
  char * infilename,
  char * intreename,
  char * wgtsfilename,
  char * catsfilename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  char * TreeUseMethod,
  int UseLengths,
  double TTratio,
  int useEmpBF,
  double BaseFreqA,
  double BaseFreqC,
  double BaseFreqG,
  double BaseFreqTU,
  int OneCat,
  int NumCats,
  double SiteRate1,
  double SiteRate2,
  double SiteRate3,
  double SiteRate4,
  double SiteRate5,
  double SiteRate6,
  double SiteRate7,
  double SiteRate8,
  double SiteRate9,
  char * RateVar,
  int AdjCor,
  double BlockLen,
  double CoeffVar,
  int NumRates,
  double HMMRate1,
  double HMMRate2,
  double HMMRate3,
  double HMMRate4,
  double HMMRate5,
  double HMMRate6,
  double HMMRate7,
  double HMMRate8,
  double HMMRate9,
  double HMMProb1,
  double HMMProb2,
  double HMMProb3,
  double HMMProb4,
  double HMMProb5,
  double HMMProb6,
  double HMMProb7,
  double HMMProb8,
  double HMMProb9,
  double InvarFract,
  int SitesWeight,
  int GlobalRe,
  int RandInput,
  int RandNum,
  int Njumble,
  int MultData,
  int MultDSet,
  int NumSeqs,
  int InputSeq,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree,
  int DotDiff,
  int RecHypo)
{
  initdata* funcs;

  //printf("Hello from DnaMLK!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Dnamlk";

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dna_node_new;
  funcs->tree_new = dnamlk_tree_new;
  progname = argv[0];

  phylipinit(argc, argv, funcs, true);

  //printf("init done\n"); // JRMdebug

  /* // JRMdebug
  // internal variables
  //ctgry = false;
  ##didchangecat = false;
  //rctgry = false;
  ##didchangercat = false;
  //categs = 1;
  //rcategs = 1;
  //auto_ = false;
  // freqsfrom = true;
  //gama = false;
  //global = false;
  //hypstate = false;
  //improve = false;
  //invar = false;
  //jumble = false;
  //njumble = 1;
  //lngths = false;
  //lambda = 1.0;
  //outgrno = 1;
  //outgropt = false;
  //trout = true;
  //ttratio = 2.0;
  ##ttr = false;
  //usertree = false;
  //weights = false;
  //printdata = false;
  //dotdiff = true;
  //progress = true;
  //treeprint = true;
  //interleaved = true;
  //mulsets = false;
  //datasets = 1;
  lambda1 = 0.0;
  lengthsopt = false;

  // from java
  //String infilename;
  //String intreename;
  //String wgtsfilename;
  //String catsfilename;
  //String outfilename;
  //String outfileopt;
  //String outtreename;
  //String outtreeopt;
  //String TreeUseMethod;
  //boolean UseLengths;
  //double TTratio;
  //boolean useEmpBF;
  //double BaseFreqA;
  //double BaseFreqC;
  //double BaseFreqG;
  //double BaseFreqTU;
  //boolean OneCat;
  //int NumCat;
  //double SiteRate1,
  //double SiteRate2,
  //double SiteRate3,
  //double SiteRate4,
  //double SiteRate5,
  //double SiteRate6,
  //double SiteRate7,
  //double SiteRate8,
  //double SiteRate9,
  //char * RateVar,
  //boolean AdjCor;
  //int BlockLen;
  //double CoeffVar,
  //int NumRates,
  //double HMMrate1,
  //double HMMrate2,
  //double HMMrate3,
  //double HMMrate4,
  //double HMMrate5,
  //double HMMrate6,
  //double HMMrate7,
  //double HMMrate8,
  //double HMMrate9,
  //double InvarFract,
  //boolean SitesWeight;
  //boolean Analysis;
  //boolean GlobalRe;
  //boolean RandInput;
  //int RandNum;
  //int Njumble;
  //boolean MultData;
  //boolean MultDSet;
  //int NumSeqs
  //boolean InputSeq;
  //boolean PrintData;
  //boolean PrintInd;
  //boolean PrintTree;
  //boolean WriteTree;
  //boolean DotDiff;
  //boolean RecHypo;
  */

  // transfer java data to local variables
  if (UseLengths != 0)
  {
    lngths = true;
  }
  else
  {
    lngths = false;
  }

  if (!strcmp(TreeUseMethod, "No"))
  {
    // No, use user trees in input file
    usertree   = true;
  }
  else
  {
    // Yes
    usertree   = false;
    lngths     = false; // shut it off no matter what the user entered
  }

  if (useEmpBF != 0)
  {
    freqsfrom = true;
  }
  else
  {
    freqsfrom = false;
  }

  if (OneCat != 0)
  {
    ctgry = false;
  }
  else
  {
    ctgry = true;
  }

  //printf("RateVar: %s\n", RateVar); // JRMdebug
  if (!strcmp(RateVar, "constant"))
  {
    // Constant rate
    rctgry = false;
    gama   = false;
    invar  = false;
    auto_  = false;
    //printf("Constant rate\n"); // JRMdebug
  }
  else if (!strcmp(RateVar, "gamma"))
  {
    // Gamma distance rates
    rctgry = true;
    gama   = true;
    invar  = false;
    //printf("Gamma distance\n"); // JRMdebug
  }
  else if (!strcmp(RateVar, "invar"))
  {
    // Gamma + invariant sites
    rctgry = true;
    gama   = false;
    invar  = true;
    //printf("Gamma +\n"); // JRMdebug
  }
  else
  {
    // User-defined HMM of rates
    rctgry = true;
    gama   = false;
    invar  = false;
    //printf("HMM\n"); // JRMdebug
  }

  if (AdjCor != 0)
  {
    auto_ = true;
    lambda = 1.0 / BlockLen;
  }
  else
  {
    auto_ = false;
    lambda = 1.0;
  }

  if (SitesWeight != 0)
  {
    weights = true;
  }
  else
  {
    weights = false;
  }

  if (GlobalRe != 0)
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

  if (MultData != 0)
  {
    mulsets = true;
    datasets = NumSeqs;
  }
  else
  {
    mulsets = false;
    datasets = 1;
  }

  if (MultDSet != 0)
  {
    justwts = false;
  }
  else
  {
    justwts = true;
  }

  if (InputSeq !=0)
  {
    interleaved = true;
  }
  else
  {
    interleaved = false;
  }

  if (WriteTree != 0)
  {
    trout = true;
  }
  else
  {
    trout = false;
  }

  if (DotDiff != 0)
  {
    dotdiff = true;
  }
  else
  {
    dotdiff = false;
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

  if (RecHypo != 0)
  {
    hypstate = true;
  }
  else
  {
    hypstate = false;
  }

  // transfer values
  ttratio   = TTratio;
  categs    = NumCats;
  rcategs   = NumRates;
  freqa     = BaseFreqA;
  freqc     = BaseFreqC;
  freqg     = BaseFreqG;
  freqt     = BaseFreqTU;
  cv        = CoeffVar;
  invarfrac = InvarFract;

  alpha = 1.0/(cv * cv);

  //printf("values transfered\n"); // JRMdebug
  //printf("maxcategs: %i\n", maxcategs); // JRMdebug

  // create arrays and initialize them
  // warning: if maxcategs ever changes, this has to change also
  rate = (double *) Malloc(maxcategs * sizeof(double));
  rate[0]   = SiteRate1;
  rate[1]   = SiteRate2;
  rate[2]   = SiteRate3;
  rate[3]   = SiteRate4;
  rate[4]   = SiteRate5;
  rate[5]   = SiteRate6;
  rate[6]   = SiteRate7;
  rate[7]   = SiteRate8;
  rate[8]   = SiteRate9;

  rrate = (double *) Malloc(maxcategs * sizeof(double));
  rrate[0]   = HMMRate1;
  rrate[1]   = HMMRate2;
  rrate[2]   = HMMRate3;
  rrate[3]   = HMMRate4;
  rrate[4]   = HMMRate5;
  rrate[5]   = HMMRate6;
  rrate[6]   = HMMRate7;
  rrate[7]   = HMMRate8;
  rrate[8]   = HMMRate9;

  probcat = (double *) Malloc(maxcategs * sizeof(double));
  probcat[0]   = HMMProb1;
  probcat[1]   = HMMProb2;
  probcat[2]   = HMMProb3;
  probcat[3]   = HMMProb4;
  probcat[4]   = HMMProb5;
  probcat[5]   = HMMProb6;
  probcat[6]   = HMMProb7;
  probcat[7]   = HMMProb8;
  probcat[8]   = HMMProb9;

  //printf("rates transfered\n"); // JRMdebug

  // conditional modification of arrays
  int i;
  if (gama)
  {
    initgammacat(rcategs, alpha, rrate, probcat);
  }
  else
  {
    if (invar)
    {
      initgammacat(rcategs-1, alpha, rrate, probcat);
      for (i = 0; i < rcategs-1; i++)
      {
        probcat[i] = probcat[i]*(1.0-invarfrac);
      }
      probcat[rcategs-1] = invarfrac;
      rrate[rcategs-1] = 0.0;
    }
  }

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

  if (ctgry)
  {
    //printf("catsfilename: %s\n", catsfilename); // JRMdebug
    catfile = fopen(catsfilename, "r");
    if (catfile == NULL)
      printf("catfile not found\n");
  }

  if (weights || justwts)
  {
    weightfile = fopen(wgtsfilename, "r");
  }

  if (trout)
  {
    outtree = fopen(OuttreeName, outtreeopt);
    strcpy(outtreename, OuttreeName);
  }

  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinit();
  ttratio0 = ttratio;

  //printf("calling dnamlrun\n");
  dnamlkrun();  // do the actual work

  FClose(infile);
  FClose(outfile);

  if (ctgry)
  {
    FClose(catfile);
  }

  if (weights || justwts)
  {
    FClose(weightfile);
  }

  if (usertree)
  {
    FClose(intree);
  }

  if (trout)
  {
    FClose(outtree);
  }

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood with molecular clock */
  initdata* funcs;

#ifdef MAC
  argc = 1;                /* macsetup("Dnamlk", "Dnamlk");        */
  argv[0] = "Dnamlk";
#endif

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dna_node_new;
  funcs->tree_new = dnamlk_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  mulsets = false;
  firstset = true;
  doinit();

  ttratio0    = ttratio;
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if (ctgry)
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);

  dnamlkrun();

  FClose(infile);
  FClose(outfile);
  FClose(outtree);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  printf("Done.\n\n");
  phyRestoreConsoleAttributes();

  return 0;
}  /* DNA Maximum Likelihood with molecular clock */


// End.
