/* Version 4.0.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Dan Fineman, Patrick Colacurcio, Elizabeth Walkup, and Jim McGill. */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PHYLIP_H
#include "phylip.h"
#define PHYLIP_H
#endif

#ifndef SEQ_H
#include "seq.h"
#endif

#ifndef BL_H
#include "bl.h"
#endif

#ifndef ML_H
#include "ml.h"
#endif

#ifndef MLDNA_H
#include "mldna.h"
#endif

struct tree *curtree, *bestree, *bestree2, *priortree;      /* global trees */

extern FILE *outfile, *infile, *intree, *outtree, *intree2, *workingplot;
extern FILE *weightfile, *catfile, *ancfile, *mixfile, *factfile;
extern FILE *progfile;
extern steptr weight, category, alias, location, ally;
extern long nodesize;

/* debug:  extern sequence inputSequences;  */

long outgrno, endsite;
long which;
boolean smoothed;

typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1, xiy1xv;
  double *ww, *zz, *wwzz, *vvzz;
} valrec;

typedef struct dnaml_tree {
  struct ml_tree ml_tree;
} dnaml_tree;

typedef long vall[maxcategs];
typedef double contribarr[maxcategs];

#ifndef OLDC
/* function prototypes */
void   dnaml_tree_new(struct tree**, long, long, long);
void   dnaml_tree_init(struct tree*, long, long);
void   dnaml_tree_setup(long, long);
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inittip(struct tree*, long);
void   inputoptions(void);
void   makeweights(void);
void   getinput(void);
void   inittable_for_usertree(FILE *);
void   inittable(void);
void   alloc_nvd (long, nuview_data *);
void   free_nvd (nuview_data *);
void   dnaml_tree_nuview(struct tree* , struct node *);
void   slopecurv(struct node *, double, double *, double *, double *);
void   dnaml_tree_makenewv(struct tree*, struct node *);
void   initdnamlnode(struct tree *, struct node *, long, long, long *,
                       long *, initops, pointarray, Char *, Char *, FILE *);
void   dnaml_coordinates(struct node *, double, long *, double *);
void   dnaml_printree(void);
void   sigma(struct node *, double *, double *, double *);
void   describe(struct node *);
void   reconstr(struct node *, long);
void   rectrav(struct node *, long, long);
void   summarize(void);
void   treevaluate(struct tree*);
void   maketree(void);
void   clean_up(void);
void   reallocsites(void);
void   dnaml_reroot(struct tree*);           // RSGbugfix: Name change.
void   dnaml_treeout(FILE *, struct tree*, struct node*);
double dnaml_tree_evaluate(struct tree*, struct node *, boolean);
void   freetable(void);
void   dnamlrun(void);
void   dnaml(char * infilename, char * intreename, char * wgtsfilename,
              char * catsfilename, char * outfilename, char * outfileopt,
              char * outtreename, char * outtreeopt, char * TreeUseMethod,
              int UseLengths, double TTratio, int useEmpBF, double BaseFreqA,
              double BaseFreqC, double BaseFreqG, double BaseFreqTU,
              int OneCat, int NumCats, double SiteRate1, double SiteRate2,
              double SiteRate3, double SiteRate4, double SiteRate5,
              double SiteRate6, double SiteRate7, double SiteRate8,
              double SiteRate9, char * RateVar, int AdjCor, double BlockLen,
              double CoeffVar, int NumRates, double HMMRate1, double HMMRate2,
              double HMMRate3, double HMMRate4, double HMMRate5,
              double HMMRate6, double HMMRate7, double HMMRate8,
              double HMMRate9, double HMMProb1, double HMMProb2,
              double HMMProb3, double HMMProb4, double HMMProb5,
              double HMMProb6, double HMMProb7, double HMMProb8,
              double HMMProb9, double InvarFract, int SitesWeight,
              int SpeedAn, int GlobalRe, int RandInput, int RandNum,
              int Njumble, int OutRoot, int OutNum, int MultData,
              int MultDSet, int NumSeqs, int InputSeq, int PrintData,
              int PrintInd, int PrintTree, int WriteTree, int DotDiff,
              int RecHypo);
/* function prototypes */
#endif

double fracchange;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff;
extern long rcategs;
boolean haslengths;

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH],
      outtreename[FNMLNGTH], catfilename[FNMLNGTH], weightfilename[FNMLNGTH];
double *rate, *rrate, *probcat;
long nonodes2, sites, weightsum, categs, datasets, ith, njumble, jumb = 0;
long parens;
boolean freqsfrom, global, jumble, weights, trout, usertree, 
         reusertree, ctgry, rctgry, auto_, hypstate, ttr, progress, mulsets,
         justwts, firstset, improve, thorough, smoothit, polishing, lngths,
         gama, invar;
struct tree *curtreee, *bestreee, *bestreee2, *priortreee;
struct tree **curtreep, **bestreep, **bestree2p, **priortreep;
struct node *qwhere;
initops whichinit;
double xi, xv, rho, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr, freqy,
        freqar, freqcy, freqgr, freqty, cv, alpha, lambda, invarfrac;
long *enterorder, inseed, inseed0;
steptr aliasweight;
contribarr *contribution, like, nulike, clai;
double **term, **slopeterm, **curveterm;
longer seed;
Char* progname;
char basechar[16]="acmgrsvtwyhkdbn";

/* Local variables for maketree, propagated globally for C version: */
long k, nextsp, numtrees, maxwhich, mx, mx0, mx1, shimotrees;
double dummy, maxlogl;
boolean succeeded;
double **l0gf;
double *l0gl;
valrec ***tbl;
Char ch, ch2;
long col;
vall *mp=NULL;


void dnaml_tree_new(struct tree** treep, long nonodes, long spp, 
		     long treesize)
{
  /* set up variables and then set up identities of functions */

  bl_tree_new(treep, nonodes, spp, sizeof(dnaml_tree));
  dnaml_tree_init(*treep, nonodes, spp);
} /* dnaml_tree_new */


void dnaml_tree_init(struct tree* t, long nonodes, long spp)
{
  /* set up functions for a dnaml_tree */

  t->evaluate = (tree_evaluate_t)dnaml_tree_evaluate;
  t->try_insert_ = (tree_try_insert_t)bl_tree_try_insert_;
  t->nuview = (tree_nuview_t)dnaml_tree_nuview;
  t->makenewv = (tree_makenewv_t)dnaml_tree_makenewv;
  t->get_fork = generic_tree_get_fork;
  t->insert_ = (tree_insert_t)bl_tree_insert_;
} /* dnaml_tree_init */


void dnaml_tree_setup(long nonodes, long spp)
{
  /* create and initialize the necessary trees */

  curtreep = (struct tree **)(&curtree);
  dnaml_tree_new(curtreep, nonodes, spp, sizeof(dnaml_tree));
  if (!usertree) {
    bestreep = (struct tree **)(&bestree);
    dnaml_tree_new(bestreep, nonodes, spp, sizeof(dnaml_tree));
    bestree2p = (struct tree **)(&bestree2);
    dnaml_tree_new(bestree2p, nonodes, spp, sizeof(dnaml_tree));
    priortreep = (struct tree **)(&priortree);
    dnaml_tree_new(priortreep, nonodes, spp, sizeof(dnaml_tree));
  }
} /* dnaml_tree_setup */


void inittip(struct tree* t, long m)
{ /* initialize and hook up a new tip;  m is the index of the tip */
/* debug:  not obvious that this is ever used */
  struct node *tmp;

  tmp = t->nodep[m - 1];
/* debug    memcpy(((mldna_node*)tmp)->x, x[m - 1], totalleles * sizeof(double));  */
  ((bl_node*)tmp)->deltav = 0.0;
}  /* inittip */


void getoptions(void)
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean didchangecat, didchangercat;
  double probsum;
  char* string;

  putchar('\n');
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  auto_ = false;
  freqsfrom = true;
  gama = false;
  global = false;
  hypstate = false;
  improve = true;
  invar = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  lambda = 1.0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  reusertree = false;
  weights = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  for (;;)
  {
    cleerhome();
    printf("Nucleic acid sequence Maximum Likelihood");
    printf(" method, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    if ( reusertree ) string = "Yes, rearrange on user tree";
    else if ( usertree ) string = "No, use user trees in input file";
    else string = "Yes";
    printf("  U                 Search for best tree?  %s\n", string);
    if (usertree && !reusertree)
    {
      printf(
  "  L          Use lengths from user trees?  %s\n", (lngths ? "Yes" : "No"));
    }
    printf(
   "  T        Transition/transversion ratio:%8.4f\n", (ttr ? ttratio : 2.0));
    printf("  F       Use empirical base frequencies?  %s\n", (freqsfrom ?
            "Yes" : "No"));
    printf("  C                One category of sites?");
    if (!ctgry || categs == 1)
      printf("  Yes\n");
    else
      printf("  %ld categories of sites\n", categs);
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
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    if (!usertree || reusertree)
    {
      printf("  S        Speedier but rougher analysis?  %s\n",
               (improve ? "No, not rough" : "Yes"));
      printf("  G                Global rearrangements?  %s\n",
               (global ? "Yes" : "No"));
    }
    if (!usertree)
    {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?  %s%3ld\n",
            (outgropt ? "Yes, at sequence number" :
              "No, use as outgroup species"), outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets, (justwts ? "sets of weights" :
              "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
            (interleaved ? "Yes" : "No, sequential"));
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
    printf("  5   Reconstruct hypothetical sequences?  %s\n",
            (hypstate ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("UTFCRAWSGJVOMI012345", ch) != NULL))
        || (usertree && ((strchr("ULTFCRAWSVOMI012345", ch) != NULL))))
    {
      switch (ch)
      {
        case 'F':
          freqsfrom = !freqsfrom;
          if (!freqsfrom)
          {
            initfreqs(&freqa, &freqc, &freqg, &freqt);
          }
          break;

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
            rate = (double *) Malloc(categs * sizeof(double));
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
            initlambda(&lambda);
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
          lngths = !lngths;
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          break;

        case 'T':
          ttr = !ttr;
          if (ttr)
          {
            initratio(&ttratio);
          }
          break;

        case 'U':
          if ((!usertree) && (!reusertree))
          {
            usertree = true;
            reusertree = false;
          }
          else if ((!reusertree) && usertree)
          {
            reusertree = true;
          }
          else
          {
            usertree = false;
            reusertree = false;
          }
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

        case '5':
          hypstate = !hypstate;
          break;
      }
    }
    else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
  if (gama || invar)
  {
    loopcount = 0;
    do {
      printf(
 "\nCoefficient of variation of substitution rate among sites (must be positive)\n");
      printf(
 " In gamma distribution parameters, this is 1/(square root of alpha)\n");
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
          countup (&loopcount, 10);
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
    rrate      = (double *) Malloc(rcategs * sizeof(double));
    probcat    = (double *) Malloc(rcategs * sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat)
  {
    rate       = (double *) Malloc(categs * sizeof(double));
    rate[0]    = 1.0;
  }
}  /* getoptions */


void reallocsites(void)
{
  long i;

  for (i=0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (Char *) Malloc(sites * sizeof(Char));
  }
  free(category);
  free(weight);
  free(alias);
  free(ally);
  free(location);
  free(aliasweight);

  category    = (long *) Malloc(sites * sizeof(long));
  weight      = (long *) Malloc(sites * sizeof(long));
  alias       = (long *) Malloc(sites * sizeof(long));
  ally        = (long *) Malloc(sites * sizeof(long));
  location    = (long *) Malloc(sites * sizeof(long));
  aliasweight = (long *) Malloc(sites * sizeof(long));

} /* reallocsites */


void allocrest(void)
{
  long i;

  inputSequences = (Char **) Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *) Malloc(sites * sizeof(Char));
  nayme       = (naym *) Malloc(spp * sizeof(naym));;
  enterorder  = (long *) Malloc(spp * sizeof(long));
  category    = (long *) Malloc(sites * sizeof(long));
  weight      = (long *) Malloc(sites * sizeof(long));
  alias       = (long *) Malloc(sites * sizeof(long));
  ally        = (long *) Malloc(sites * sizeof(long));
  location    = (long *) Malloc(sites * sizeof(long));
  aliasweight = (long *) Malloc(sites * sizeof(long));
}  /* allocrest */


void doinit(void)
{ /* initializes variables */
  fprintf(outfile, "\nNucleic acid sequence Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n", VERSION);

  inputnumbers(&spp, &sites, &nonodes2, 1);
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
  {
    inputweights(sites, weight, &weights);
  }
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if (ctgry && categs > 1)
  {
    inputcategs(0, sites, category, categs, "DnaML");
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
  sitesort2   (sites, aliasweight);
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
  term = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    term[i] = (double *) Malloc(rcategs * sizeof(double));
  slopeterm = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    slopeterm[i] = (double *) Malloc(rcategs * sizeof(double));
  curveterm = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    curveterm[i] = (double *) Malloc(rcategs * sizeof(double));
  mp = (vall *) Malloc(sites * sizeof(vall));
  contribution = (contribarr *) Malloc(endsite * sizeof(contribarr));
}  /* makeweights */


void getinput(void)
{
  /* reads the input data */
  inputoptions();
  if (!freqsfrom)
  {
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
           &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange, freqsfrom, true);
  }
  if (!justwts || firstset)
  {
    inputdata(sites);
  }
  makeweights();
  dnaml_tree_setup(nonodes2, spp);
  makevalues2(rcategs, ((struct tree*)curtree)->nodep, endsite, spp, 
		         inputSequences, alias);
  if (freqsfrom)
  {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight, 
		    ((struct tree*)curtree)->nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
           &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange, freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


void inittable_for_usertree(FILE *intree)
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
} /* freetable */


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
         necessary (i.e. there's a user tree) */
      tbl[i][j]->ww   = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc(2 * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc(2 * sizeof (double));
    }
  }
  if (!lngths)            /* restandardize rates */
  {
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
  }

  if(jumb > 1)
    return;

  if (rcategs > 1)
  {
    if (gama)
    {
      fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
      fprintf(outfile,
          " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
    }
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
      fprintf(outfile,
        "Expected length of a patch of sites having the same rate = %8.3f\n",
        1/lambda);
    putc('\n', outfile);
  }
  if (categs > 1)
  {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 0; i < categs; i++)
      fprintf(outfile, "%9ld%16.3f\n", i+1, rate[i]);
  }
  if ((rcategs  > 1) || (categs >> 1))
    fprintf(outfile, "\n\n");
}  /* inittable */


double dnaml_tree_evaluate(struct tree* t, struct node *p, boolean saveit)
{ /* dnaml version of evaluation of likelihood */
  contribarr tterm;
  double sum, sum2, sumc, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3,
          sumterm, lterm;
  long i, j, k, lai;
  struct node *q;
  sitelike x1, x2;

  if (p->back == NULL)                   /* ensure branch has non-null ends */
    p = p->next;
  generic_update(t, p);                                 /* views traversals */
  sum = 0.0;
  q = p->back;
  y = ((struct bl_node*)p)->v;
  lz = -y;
  for (i = 0; i < rcategs; i++)    /* get probabilities for different rates */
    for (j = 0; j < categs; j++)
    {
      tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
      tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
      tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
      tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
    }
  for (i = 0; i < endsite; i++)             /* do over all aliases of sites */
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)       /* ... and over all rate categories */
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
            /* the detailed putting-together of the probabilities of change */
      memcpy(x1, ((mldna_node*)p)->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
        freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, ((mldna_node*)q)->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
        freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
        (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
        (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
        (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
        freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
        freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
        freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i]
                         + ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, rcategs * sizeof(double));
    if (saveit && !auto_ && usertree && (which <= shimotrees) && !reusertree)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  for (j = 0; j < rcategs; j++)
    like[j] = 1.0;
  for (i = 0; i < sites; i++)              /* calculation for actual sites */
  {
    sumc = 0.0;
    for (k = 0; k < rcategs; k++)           /* ... and for rate categories */
      sumc += probcat[k] * like[k];
    sumc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0))
    {
      lai = location[ally[i] - 1];
      memcpy(clai, contribution[lai - 1], rcategs * sizeof(double));
      for (j = 0; j < rcategs; j++)    /* autocorrelated rates calculation */
        nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
    }
    else
    {
      for (j = 0; j < rcategs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc);
    }
    memcpy(like, nulike, rcategs * sizeof(double));
  }
  sum2 = 0.0;
  for (i = 0; i < rcategs; i++)
    sum2 += probcat[i] * like[i];
  sum += log(sum2);
  ((struct tree*)t)->score = sum;                /* putting into tree score */
  if (!saveit || auto_ || !usertree || reusertree)
    return sum;
  if (which <= shimotrees)             /* saving log likelihood in SH table */
    l0gl[which - 1] = sum;
  if (which == 1)
  {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl)                             /* upate which is besttree */
  {
    maxwhich = which;
    maxlogl = sum;
  }

  return sum;
}  /* dnaml_tree_evaluate */


void alloc_nvd (long num_sibs, nuview_data *local_nvd)
{
  /* Allocate blocks of memory appropriate for the number of siblings a node has */
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


void free_nvd (nuview_data *local_nvd)
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


void dnaml_tree_nuview(struct tree* t, struct node *p)
{ 
  /* recursive computation of conditional likelihoods for a node based
   * on those of its descendants */
  long i, j, k, l, num_sibs, sib_index;
  nuview_data *local_nvd;
  struct node *sib_ptr, *sib_back_ptr; /* for the ends of descendant branch */
  sitelike p_xx;    /* vector of site conditional likelihoods for that node */
  double lw;
  double correction;                  /* correction to cope with underflows */
  double maxx;                        /* maximum of conditional likelihoods */

  /* Allocate the structure and blocks for variables used in this function. */
  num_sibs = count_sibs((struct node*)p);
  local_nvd = (nuview_data *) Malloc(sizeof (nuview_data));
  alloc_nvd (num_sibs, local_nvd);

  /* Loop 1: makes assignments to tbl based on some combination of
     what's already in tbl and the children's values of v */
  sib_ptr = (struct node*)p;
  for (sib_index = 0; sib_index < num_sibs; sib_index++)
  {   /* for each descendant lineage tabulate some part of transition prob */
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (sib_back_ptr != NULL)
      lw = - ((struct bl_node*)sib_back_ptr)->v;
    else
      lw = 0.0;
    for (i = 0; i < rcategs; i++) { /* table of terms for transition probs */
      for (j = 0; j < categs; j++)
      {
        tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
        tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
        tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] *
                                      tbl[i][j]->zz[sib_index];
        tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
                                      tbl[i][j]->zz[sib_index];
      }
    }
  }
  for (i = 0; i < endsite; i++)                                 /* Loop 2: */
  {
    correction = 0;
    maxx = 0.0;
    k = category[alias[i]-1] - 1; /* get user-defined category for the site */
    for (j = 0; j < rcategs; j++)       /* Loop 2.1: for each rate category */
    {
      sib_ptr = (struct node*)p;
      for (sib_index = 0; sib_index < num_sibs; sib_index++)
      {
        sib_ptr         = sib_ptr->next;
        sib_back_ptr    = sib_ptr->back;

        if (sib_back_ptr != NULL) {        /* otherwise no table to fill in */
          if ( j == 0 )
            correction += ((struct ml_node*)sib_back_ptr)->underflows[i];

          local_nvd->wwzz[sib_index] = tbl[j][k]->wwzz[sib_index];
          local_nvd->vvzz[sib_index] = tbl[j][k]->vvzz[sib_index];
          local_nvd->yy[sib_index]   = 1.0 - tbl[j][k]->zz[sib_index];
          memcpy(local_nvd->xx[sib_index],
                        ((mldna_node*)sib_back_ptr)->x[i][j], sizeof(sitelike));
        }
      }

      sib_ptr = (struct node*)p;
      for (sib_index = 0; sib_index < num_sibs; sib_index++)
      {                    /* Loop 2.2:  pieces for each descendant lineage */
        sib_ptr = sib_ptr->next;
        if (sib_ptr != NULL) {
          sib_back_ptr = sib_ptr->back;
          if (sib_back_ptr != NULL) {
            local_nvd->sum[sib_index] =
              local_nvd->yy[sib_index] *
                (freqa * local_nvd->xx[sib_index][(long)A] +
                 freqc * local_nvd->xx[sib_index][(long)C] +
                 freqg * local_nvd->xx[sib_index][(long)G] +
                 freqt * local_nvd->xx[sib_index][(long)T]);
            local_nvd->sumr[sib_index] =
                                  freqar*local_nvd->xx[sib_index][(long)A]
                                   + freqgr*local_nvd->xx[sib_index][(long)G];
            local_nvd->sumy[sib_index] =
                                  freqcy*local_nvd->xx[sib_index][(long)C]
                                   + freqty*local_nvd->xx[sib_index][(long)T];
            local_nvd->vzsumr[sib_index] = local_nvd->vvzz[sib_index]
                                           * local_nvd->sumr[sib_index];
            local_nvd->vzsumy[sib_index] = local_nvd->vvzz[sib_index]
                                          * local_nvd->sumy[sib_index];
            }
         }
      }

        /* Initialize to one, multiply incremental values for every sibling */
      p_xx[(long)A] = 1.0 ;
      p_xx[(long)C] = 1.0 ;
      p_xx[(long)G] = 1.0 ;
      p_xx[(long)T] = 1.0 ;

      sib_ptr = (struct node*)p;
      for (sib_index = 0; sib_index < num_sibs; sib_index++)
      {    /* go over descendant lineages, multiplying pieces for each site */
        sib_ptr = sib_ptr->next;
        if (sib_ptr != NULL) {
          sib_back_ptr = sib_ptr->back;
          if (sib_back_ptr != NULL) {
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
         }
      }

      for ( l = 0 ; l < ((long)T - (long)A + 1); l++ )
      {
        if ( p_xx[l] > maxx )       /* find max of likelihood for each site */
          maxx = p_xx[l];
      }

      /* And the final point of this whole function: */
      memcpy(((mldna_node*)p)->x[i][j], p_xx, sizeof(sitelike));
    }

    ((ml_node*)p)->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE)
      fix_x((mldna_node*)p, i, maxx, rcategs);
    ((ml_node*)p)->underflows[i] += correction;
  }                                               /* end of loop over sites */

  /* debug */ if (p->back != 0) printf("nuview of %ld from %ld\n", p->index, p->back->index); else printf("nuvuiew of %ld from empty root\n", p->index);
  ((struct node*)p)->initialized = true;  /* mark node as has view updated */

  free_nvd (local_nvd);
  free (local_nvd);

}  /* dnaml_tree_nuview */


void slopecurv(struct node *p, double y, double *like, 
                 double *slope, double *curve)
{
  /* compute log likelihood, slope and curvature at node p */
  long i, j, k, lai;
  double sum, sumc, sumterm, lterm, sumcs, sumcc, sum2, slope2, curve2, temp;
  double lz, zz, z1, zzs, z1s, zzc, z1c, aa, bb, cc,
          prod1, prod2, prod12, prod3;
  contribarr thelike, nulike, nuslope, nucurve,
          theslope, thecurve, clai, cslai, cclai;
  struct node *q;
  sitelike x1, x2;

  q = ((struct node*)p)->back;
  sum = 0.0;
  lz = -y;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++)
    {
      tbl[i][j]->orig_zz = exp(tbl[i][j]->rat * lz);
      tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
    }
  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      if (y > 0.0)
      {
        zz = tbl[j][k]->orig_zz;
        z1 = tbl[j][k]->z1;
      }
      else
      {
        zz = 1.0;
        z1 = 1.0;
      }
      zzs = -tbl[j][k]->rat * zz ;
      z1s = -tbl[j][k]->ratxv * z1 ;
      temp = tbl[j][k]->rat;
      zzc = temp * temp * zz;
      temp = tbl[j][k]->ratxv;
      z1c = temp * temp * z1;
      memcpy(x1, ((mldna_node*)p)->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
        freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, ((mldna_node*)q)->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
        freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
        (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
        (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
        (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
        freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
        freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
        freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      aa = prod12 - prod3;
      bb = prod3 - prod1*prod2;
      cc = prod1 * prod2;
      term[i][j] = zz * aa + z1 * bb + cc;
      slopeterm[i][j] = zzs * aa + z1s * bb;
      curveterm[i][j] = zzc * aa + z1c * bb;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * term[i][j];
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i]
             + ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
    {
      term[i][j] = term[i][j] / sumterm;
      slopeterm[i][j] = slopeterm[i][j] / sumterm;
      curveterm[i][j] = curveterm[i][j] / sumterm;
    }
    sum += aliasweight[i] * lterm;
  }
  for (i = 0; i < rcategs; i++)
  {
    thelike[i] = 1.0;
    theslope[i] = 0.0;
    thecurve[i] = 0.0;
  }
  for (i = 0; i < sites; i++)
  {
    sumc = 0.0;
    sumcs = 0.0;
    sumcc = 0.0;
    for (k = 0; k < rcategs; k++)
    {
      sumc += probcat[k] * thelike[k];
      sumcs += probcat[k] * theslope[k];
      sumcc += probcat[k] * thecurve[k];
    }
    sumc *= lambda;
    sumcs *= lambda;
    sumcc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0))
    {
      lai = location[ally[i] - 1];
      memcpy(clai, term[lai - 1], rcategs * sizeof(double));
      memcpy(cslai, slopeterm[lai - 1], rcategs * sizeof(double));
      memcpy(cclai, curveterm[lai - 1], rcategs * sizeof(double));
      if (weight[i] > 1)
      {
        for (j = 0; j < rcategs; j++)
        {
          if (clai[j] > 0.0)
            clai[j] = exp(weight[i]*log(clai[j]));
          else clai[j] = 0.0;
          if (cslai[j] > 0.0)
            cslai[j] = exp(weight[i]*log(cslai[j]));
          else cslai[j] = 0.0;
          if (cclai[j] > 0.0)
            cclai[j] = exp(weight[i]*log(cclai[j]));
          else cclai[j] = 0.0;
        }
      }
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = ((1.0 - lambda) * thelike[j] + sumc) * clai[j];
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs) * clai[j]
          + ((1.0 - lambda) * thelike[j] + sumc) * cslai[j];
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc) * clai[j]
          + 2.0 * ((1.0 - lambda) * theslope[j] + sumcs) * cslai[j]
          + ((1.0 - lambda) * thelike[j] + sumc) * cclai[j];
      }
    }
    else
    {
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = ((1.0 - lambda) * thelike[j] + sumc);
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs);
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc);
      }
    }
    memcpy(thelike, nulike, rcategs * sizeof(double));
    memcpy(theslope, nuslope, rcategs * sizeof(double));
    memcpy(thecurve, nucurve, rcategs * sizeof(double));
  }
  sum2 = 0.0;
  slope2 = 0.0;
  curve2 = 0.0;
  for (i = 0; i < rcategs; i++)
  {
    sum2 += probcat[i] * thelike[i];
    slope2 += probcat[i] * theslope[i];
    curve2 += probcat[i] * thecurve[i];
  }
  sum += log(sum2);
  (*like) = sum;
  (*slope) = slope2 / sum2;

  /* compute curvature of log L to avoid overflow */
  (*curve) = curve2 / sum2 - *slope * *slope;
} /* slopecurv */


void dnaml_tree_makenewv(struct tree* t, struct node* p)
{
 /* Newton-Raphson algorithm / simple search improvement of a branch length */
  long it;
  double y=0.0, yold=0, like, slope, curve, oldlike, delta;
  boolean done, better, posslope;
  struct node *q;

/* debug  printf("smooth branch %ld:%ld \n", p->index, p->back->index);*/ 
/* debug  printf("((struct mldna_node*)p)->x[0][0][A] = %ld, %12.6f\n", p->index, ((struct mldna_node*)p)->x[0][0][0]);*/ 
/* debug printf("((struct mldna_node*)(p->back))->x[0][0][A] = %ld,%12.6f\n", p->back->index, ((struct mldna_node*)(p->back))->x[0][0][0]);*/ 
  if ((p->index == outgrno) || (p->back->index == outgrno)) {
    ((struct bl_node*)p)->v = epsilon;
    ((struct bl_node*)(p->back))->v = epsilon;
  } else {
    q = p->back;
    y = ((struct bl_node*)p)->v;
    yold = y;
    delta = y/2.0;        /* initial step size for non-Newton-Raphson steps */
    slopecurv (p, y, &like, &slope, &curve);
/* printf(" %ld:%ld v, like,  %10.6f %12.6f %12.6f %12.6f\n", p->index, q->index, y, like, slope, curve); debug */
    posslope = slope >= 0.0;
    oldlike = like;
    better = true;
    if (curve < 0.0) {
      delta = - slope / curve;
      y = y + delta;                            /* Newton-Raphson iteration */   
      }	  
    else {                                  /* when can't do Newton-Raphson */
      delta = yold/2.0;
      if (slope > 0.0)
        y = y + delta;
      else
        y = y - delta;
      }
    if (y <= epsilon) {
      y = epsilon;                          /* do not allow to go negative */
    }
    done = false;
    it = 0;
    while ((it < iterations) && (it < 20) && (!done))
    {
      slopecurv (p, y, &like, &slope, &curve);
/* printf(" %ld:%ld v, like,  %10.6f %12.6f %12.6f %12.6f\n", p->index, q->index, y, like, slope, curve); debug */
      better = like > oldlike;
      if (better) {                           /* if likelihood has improved */
        posslope = slope > 0.0;
        delta = y - yold;                               /* the step we made */
	oldlike = like;                            /* update likelihood ... */
	if (curve < 0.0) {
          delta = - slope/curve;                /* Newton-Raphson iteration */
          }
        else {
          if (((delta > 0.0) && posslope) || ((delta <= 0.0) && !posslope))
            delta = 2.0*delta;               /* step twice as far next time */
          else
            delta = - delta;            /* turn around and go the other way */
        }
        if (delta < -(y-epsilon))           /* if about to jump too far ... */
          delta = -(y-epsilon);             /* ... prepare to jump less far */
	yold = y;                                   /* update branch length */
        posslope = slope >= 0.0;
/* debug        printf("Better! next delta now %10.8f\n", delta); */
      } else {                                       /* if not better ... */
        delta = (y - yold)/2.0;                /* next time, a smaller step */
/* debug:  printf("Not better. y, yold now %10.8f, %10.8f, next delta now %10.8f\n", y, yold, delta); */
        }
      if (((delta < 0.0) && posslope) || ((delta >= 0.0) && !posslope))
        delta = - delta;
      y = yold + delta;                               /* take the next jump */
      if (fabs(delta) < epsilon)              /* if change is too small ... */
        it = 20;                        /* then don't do any more iterating */
      it++;
      done = fabs(delta) < 0.1*epsilon;  /* debug:  ? < epsilon ? */
    }
/* debug */ printf("dnaml_makenewv: %ld:%ld, now: %13.7f, was: %13.7f, like = %12.6f\n", p->index, q->index, y, ((struct bl_node*)p)->v, like);
    smoothed = (fabs(y-yold) < epsilon) && (yold > 10.0*epsilon);
    ((struct bl_node*)p)->v = yold;    /* the last one with better likelihood */
    ((struct bl_node*)(p->back))->v = yold;
    ((struct tree*)t)->score = oldlike;       /* score is the best likelihood */
  }
} /* dnaml_tree_makenewv */


void initdnamlnode(struct tree *treep, struct node *p, long len,
                    long nodei, long *ntips, long *parens, initops whichinit, 
                    pointarray nodep, Char *str, Char *ch, FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;
  struct bl_node *bln, *blnb;
  (void)len;                            // RSGnote: Parameter never used.
  (void)ntips;                          // RSGnote: Parameter never used.
  switch (whichinit)
  {
    case bottom:
      p = treep->get_forknode(treep, nodei);
      mldna_node_allocx(p, endsite, rcategs);
      assert(((struct node*)(p))->index > 0);
      nodep[((struct node*)(p))->index - 1] = (struct node*)(p);
      break;
    case nonbottom:
      p = treep->get_forknode(treep, nodei);
      mldna_node_allocx(p, endsite, rcategs);
      break;
    case tip:
      match_names_to_data (str, nodep, (struct node**)p, spp);
      break;
    case iter:
      bln = (struct bl_node*)p;
      ((struct node*)p)->initialized = false;
      bln->v = initialv;
      bln->iter = true;
      blnb = (struct bl_node*)(((struct node*)p)->back);
      if (blnb != NULL)
      {
        blnb->iter = true;
        blnb->v = initialv;
        ((struct node*)p)->initialized = false;
      }
      break;
    case length:
      blnb = (struct bl_node*)(((struct node*)p)->back);
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      ((struct bl_node*)p)->v = valyew / divisor / fracchange;
      ((struct bl_node*)p)->iter = false;
      if (((struct node*)p)->back != NULL)
      {
        ((struct bl_node*)((struct node*)p)->back)->v 
                              = ((struct bl_node*)p)->v;
        blnb->iter = false;
      }
      break;
    case hsnolength:
      haslengths = false;
      break;
    default:        /* cases hslength, treewt, unittrwt */
      break;        /* should never occur               */
  }
} /* initdnamlnode */


void dnaml_coordinates(struct node *p, double lengthsum, 
		        long *tipy, double *tipmax)
{
  /* establishes coordinates of nodes */
  struct node *q0, *q, *qprev, *first, *last;
  double xx;
  boolean atroot, dodo;

  if (p->tip)                             /* for tips, set ycoord, min, max */
  {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
  /* debug */ printf(" %ld: (x, y): (%10.4f, %10.4f), ymin, ymax   %ld, %ld\n", p->index, p->xcoord, p->ycoord, p->ymin, p->ymax );
    return;                            /* then bail if tip coordinates call */
  }
  atroot = (p == curtree->root);                       /* for interior node */
  q0 = p;                                        /* q0 starts at the node p */
  if (!atroot)
    q = q0->next;      /* unless at root node, starts at next one in circle */
  else
    q = q0;
  do {                /* go around ring, recursing into descendant subtrees */
    dodo = (atroot && (q->back != 0)) || (!atroot && (q != q0));
    if (dodo) {                     /* dodo is "do if not at end of circle" */
      xx = fracchange * ((struct bl_node*)q)->v;
      if (xx > 100.0)
        xx = 100.0;
      dnaml_coordinates(q->back,  lengthsum + xx, tipy, tipmax); /* recurse */
    }
    q = q->next;
  } while (dodo);
  if (atroot && (p->back != 0))  /* set first, last pointers to descendants */
    first = p->back;
  else
    first = p->next->back;
  q = p;                    /* find last immediate descendant and set "last"*/
  while (q->next != p) {      /* if we're all way around this interior node */
    qprev = q;
    q = q->next;
  }
  if (q->back == 0)                /* if we're at a node with an empty back */
    q = qprev;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);      /* how far our from root */
  if (p == curtree->root)     /* debug:  why this? */
    p->ycoord = p->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;              /* get leftmost descendant value of y */
  p->ymax = last->ymax;                        /* ... and rightmost one too */
  /* debug */ printf(" %ld: (x, y): (%10.4f, %10.4f), ymin, ymax   %ld, %ld,  first, last: %ld %ld\n", p->index, p->xcoord, p->ycoord, p->ymin, p->ymax, first->index, last->index);
}  /* dnaml_coordinates */


void dnaml_printree(void)
{
  /* prints out diagram of the tree2 */
  long tipy;
  double scale, tipmax;
  long i;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  if (curtree->root->tip)
    curtree->root = curtree->root->back;
  dnaml_coordinates(curtree->root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline2(i, scale, curtree);
  putc('\n', outfile);
}  /* dnaml_printree */


void sigma(struct node *p, double *sumlr, double *s1, double *s2)
{
  /* compute standard deviation */
  double tt, aa, like, slope, curv;
  struct bl_node *pp;

  pp = (struct bl_node*)p;
  slopecurv (p, pp->v, &like, &slope, &curv);
  tt = pp->v;
  pp->v = epsilon;
  ((struct bl_node*)(p->back))->v = epsilon;
  aa = curtree->evaluate(curtree, p, false);
  pp->v = tt;
  ((struct bl_node*)(p->back))->v = tt;
  (*sumlr) = curtree->evaluate(curtree, p, false) - aa;
  if (curv < -epsilon)
  {
    (*s1) = pp->v + (-slope - sqrt(slope * slope -  3.841 * curv)) / curv;
    (*s2) = pp->v + (-slope + sqrt(slope * slope -  3.841 * curv)) / curv;
  }
  else
  {
    (*s1) = -1.0;
    (*s2) = -1.0;
  }
}  /* sigma */


void describe(struct node *p)
{
  /* print out information for one branch */
  long i, num_sibs;
  struct node *q, *sib_ptr;
  struct bl_node *qq;
  double sumlr, sigma1, sigma2;

  if (p != NULL) {
    if (!p->tip && !p->initialized)
      generic_tree_nuview(curtree, p);
    if (!p->back->tip && !p->back->initialized)
      generic_tree_nuview(curtree, p->back);
    q = p->back;
    qq = (struct bl_node*)(p->back);

    assert(p->index > 0);                 // RSGdebug
    assert(q->index > 0);                 // RSGdebug

    if (q->tip)
    {
      fprintf(outfile, " ");
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->index-1][i], outfile);
      fprintf(outfile, "    ");
    }
    else
      fprintf(outfile, "  %4ld          ", q->index - spp);
    if (p->tip)
    {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[p->index-1][i], outfile);
    }
    else
      fprintf(outfile, "%4ld      ", p->index - spp);
    fprintf(outfile, "%15.5f", qq->v * fracchange);
    if (reusertree || !usertree || (usertree && !lngths) || 
		                        ((bl_node*)p)->iter )
    {
      sigma(q, &sumlr, &sigma1, &sigma2);
      if (sigma1 <= sigma2)
        fprintf(outfile, "     (     zero,    infinity)");
      else
      {
        fprintf(outfile, "     (");
        if (sigma2 <= 0.0)
          fprintf(outfile, "     zero");
        else
          fprintf(outfile, "%9.5f", sigma2 * fracchange);
        fprintf(outfile, ",%12.5f", sigma1 * fracchange);
        putc(')', outfile);
      }
      if (sumlr > 1.9205)
        fprintf(outfile, " *");
      if (sumlr > 2.995)
        putc('*', outfile);
    }
    putc('\n', outfile);
    if (!p->tip)
    {
      num_sibs = count_sibs(p);
      sib_ptr  = p;
      for (i=0; i < num_sibs; i++)
      {
        sib_ptr = sib_ptr->next;
        if (sib_ptr->back != NULL)
          describe(sib_ptr->back);
      }
    }
  }
}  /* describe */


void reconstr(struct node *p, long n)
{
  /* reconstruct and print out base at site n+1 at node p */
  long i, j, k, m, first, second, num_sibs;
  double f, sum, xx[4];
  struct node *q;

  j = location[ally[n]-1] - 1;
  for (i = 0; i < 4; i++)
  {
    f = ((mldna_node*)p)->x[j][mx-1][i];
    num_sibs = count_sibs(p);
    q = p;
    for (k = 0; k < num_sibs; k++)
    {
      q = q->next;
      f *= ((mldna_node*)q)->x[j][mx-1][i];
    }
    if (f > 0.0)   /* correct for overcounting of conditional likelihoods */
      f = exp(log(f)/num_sibs);
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


void rectrav(struct node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  long i;

  assert(p->index > 0);                 // RSGdebug

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
    rectrav(p->next->back, m, n);
    rectrav(p->next->next->back, m, n);
  }
  mx1 = mx;
}  /* rectrav */


void summarize(void)
{
  /* print out branch length information and node numbers */
  long i, j, num_sibs;
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;
  struct node   *sib_ptr;
  long mm = 0;   // RSGnote: Was formerly potentially used before being initialized.

  if (!treeprint)
    return;
  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree->score);
  fprintf(outfile, "\n Between        And             Length");
  if (!(!reusertree && usertree && lngths && haslengths))
    fprintf(outfile, "       Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, " -------        ---             ------");
  if (!(!reusertree && usertree && lngths && haslengths))
    fprintf(outfile, "       ------- ---------- ------");
  fprintf(outfile, "\n\n");
  for (i = spp; i < nonodes2; i++)
  {
    /* So this works with arbitrary multifurcations */
    if (curtree->nodep[i])
    {
      num_sibs = count_sibs (curtree->nodep[i]);
      sib_ptr  = curtree->nodep[i];
      for (j = 0; j < num_sibs; j++)
      {
        sib_ptr->initialized = false;
        sib_ptr              = sib_ptr->next;
      }
    }
  }

  describe(curtree->root->back);

  /* So this works with arbitrary multifurcations */
  num_sibs = count_sibs (curtree->root);
  sib_ptr  = curtree->root;
  for (i=0; i < num_sibs; i++)
  {
    sib_ptr = sib_ptr->next;
    describe(sib_ptr->back);
  }

  fprintf(outfile, "\n");
  if (!(!reusertree && usertree && lngths && haslengths))
  {
    fprintf(outfile, "     *  = significantly positive, P < 0.05\n");
    fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n");
  }
  dummy = curtree->evaluate(curtree, curtree->root, false);
  if (rctgry && rcategs > 1)
  {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
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
    fprintf(outfile,
 "Combination of categories that contributes the most to the likelihood:\n\n");
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
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
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
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
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
    fprintf(outfile, "Most probable category at each site\n\n");
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
      fprintf(outfile, "%ld", mm+1);
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
    fprintf(outfile,
"\n\nMost probable category at each site if > 0.95 probability (\".\" otherwise)\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < sites; i++)
    {
      sum = 0.0;
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
      rectrav(curtree->root->back, i, k);
      putc('\n', outfile);
    }
  }
}  /* summarize */


void dnaml_reroot(struct tree* t) 
{
  /* move root of tree */
  struct node *q;
  double newl;
  struct node *r = t->root;
  struct bl_node *rn, *rnn, *rnb, *rnnb;
  long numsibs = count_sibs(r);

  if ( numsibs > 2)
  {
    q = r;
    while ( q->next != r )
      q = q->next;
    q->next = r->next;
    t->release_forknode(t, r);
    t->nodep[spp] = q;
  }
  else
  {
    assert(r->back == NULL); // RSGnote: This assumes the FORKRING being manipulated
                             //has the ROOT FORKNODE pointing to NULL.


    rn = (struct bl_node*)(r->next);
    rnb = (struct bl_node*)(r->next->back);
    rnnb = (struct bl_node*)(r->next->next->back);
    rnn = (struct bl_node*)(r->next->next);

    newl = rnb->oldlen + rnn->oldlen;
    rnb->oldlen = newl;
    rnnb->oldlen = newl;

    newl = rn->v + rnn->v;
    rnb->v = newl;
    rnnb->v = newl;

    r->next->back->back = r->next->next->back;
    r->next->next->back->back = r->next->back;

   t->release_fork(t, r->index-1);
  }

  t->root = t->nodep[0]->back;
                // Reset ROOT; moved from line after DNAML_REROOT call.
} /* dnaml_reroot */


void dnaml_treeout(FILE  *outtree, struct tree* t, struct node* p) {
/* call bl_treeout to write tree out to tree output file */
  double bl_scale;

  bl_scale = fracchange;
  bl_treeout(outtree, t, p, bl_scale);
} /* dnaml_treeout */


void maketree(void)
{
  long i, k;
  boolean dummy_first, goteof;
  long nextnode, currentoutgrno;
  double bestyet;
  struct node *q;

  inittable();

  if (usertree)       /* case in which we are reading in user-defined trees */
  {
    if (!javarun)
    {   /* Open in binary: ftell() broken for UNIX line-endings under WIN32 */
      openfile(&intree, INTREE, "input tree file", "rb",progname, intreename);
    }
    //printf("calling inittable_for_usertree\n");
    inittable_for_usertree (intree);
    numtrees = countsemic(intree);
    //printf("numtrees: %li\n", numtrees);

    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2 && !reusertree)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *) Malloc(shimotrees * sizeof(double));
    l0gf = (double **) Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *) Malloc(endsite * sizeof(double));
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }

    /* This taken out of treeread, used to be [spp-1], but referring to [0]
     * produces output identical to what pre-modified dnaml produced. */
    whichinit = bottom;
    for (which = 1; which <= numtrees ; which++)    /* loop over user trees */
    {
      /* These initializations required each time through the loop since
       * multiple trees require re-initialization */
      haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
/* debug:        preparetree(curtree);   needed? */
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof,
                &dummy_first, &nextnode, &haslengths, whichinit,
                false, nonodes2);
/* debug:       fixtree(curtree);    needed? */
/* debug:      dnaml_reroot(curtree);   */                     // RSGbugfix: Name change.

      if (goteof && (which <= numtrees))
      {
        /* if we hit the end of the file prematurely */
        printf ("\nERROR:  Trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit(-1);
      }

      // RSGbugfix: Reset of ROOT moved to inside DNAML_REROOT.

      if ( outgropt )
        curtree->root = curtree->nodep[outgrno - 1]->back;

      improve = !lngths;    /* re-estimate lengths of branches in user tree */
      bl_treevaluate(curtree, improve, reusertree, global, progress,
                      priortree, bestree, (initialvtrav_t)bl_initialvtrav);
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

      dnaml_printree();
      summarize();

      if (trout) {
        dnaml_treeout(outtree, curtree, curtree->root);
      }
      if (which < numtrees)
      {
        freex_notip(nextnode, curtree->nodep);
      }
      else nonodes2 = nextnode;
    }
    FClose(intree);
    putc('\n', outfile);
    if (!auto_ && numtrees > 1 && weightsum > 1 && !reusertree )
    {
      standev2(numtrees, maxwhich, 0, endsite-1, maxlogl, l0gl, l0gf,
                aliasweight, seed);
    }
  }
  else  /* case in which make tree by sequential addition and rearrangement */
  {
    smoothit = improve;
    for (i = 1; i <= spp; i++)            /* If there's no input user tree, */
      enterorder[i - 1] = i;     /* will consider species in order, but ... */

    if (jumble)       /* ... if species to be in random order, permute them */
      randumize(seed, enterorder);

    if (progress)
    {
      sprintf(progbuf, "\nAdding species:\n");
      print_progress(progbuf);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }

    nextsp = 3;
    polishing = false;
    release_all_forks(curtree);                   /* make sure starts empty */
    buildsimpletree(curtree, enterorder);        /* make a fork with 3 tips */
    currentoutgrno = enterorder[0];
    if (enterorder[1] == outgrno)  /* set current outgroup to real outgroup */
      currentoutgrno = outgrno;
    if (enterorder[2] == outgrno)
      currentoutgrno = outgrno;
    curtree->root = curtree->nodep[currentoutgrno-1];
    generic_root_insert(curtree, curtree->nodep[currentoutgrno-1]->back);
    smoothit = improve;
    bl_initialvtrav(curtree, (struct bl_node*)(curtree->root));
    thorough = true;
    bl_update(curtree, curtree->root);
    bl_tree_smoothing(curtree, curtree->root);
    nextsp = 4;
    while (nextsp <= spp)
    {
      curtree->copy(curtree, priortree);                       /* save tree */
      k = generic_tree_findemptyfork(curtree);  /* connect next tip to fork */
      q = curtree->get_fork(curtree, k);
      curtree->nodep[k] = q;
      bl_hookup((struct node *)curtree->nodep[enterorder[nextsp-1]-1], 
                 (struct node *)q);
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
/* debug:  here make sure add rootmost fork in right place */
      curtree->addtraverse(curtree, (struct node*)q, curtree->root, further,
                             (struct node *)qwhere, &bestyet, bestree,
                             thorough, false, true, &bestyet);
      bestree->copy(bestree, curtree);
      if (currentoutgrno > enterorder[nextsp-1]) {
        currentoutgrno = enterorder[nextsp-1];
        curtree->root = curtree->nodep[currentoutgrno-1];
      }
      curtree->outgrno = currentoutgrno;
      if (progress)
      {
        writename(nextsp - 1, 1, enterorder);
        phyFillScreenColor();
      }
      if (global && nextsp == spp)
      {
        curtree->globrearrange(curtree, bestree, progress,
                                smoothit, &bestyet);
      }
      else
      {
        curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1],
                     smoothit, &bestyet, bestree, priortree, false, &bestyet);
      }
      nextsp++;
    }
    curtree->copy(curtree, bestree);
    if (global && progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
    if (njumble > 1)
    {
      if (jumb == 1)
        bestree->copy(bestree, bestree2);
      else
        if (bestree2->score < bestree->score)
          bestree->copy(bestree, bestree2);
    }
    if (jumb == njumble)
    {
      if (njumble > 1)
        bestree2->copy(bestree2, curtree);
      else
	bestree->copy(bestree, curtree);
      curtree->root = curtree->nodep[outgrno - 1]->back;
      for (i = 0; i < nonodes2; i++)
      {
        if (i < spp)
          curtree->nodep[i]->initialized = false;
        else
        {
          if (curtree->nodep[i] != NULL) {
            curtree->nodep[i]->initialized = false;
            curtree->nodep[i]->next->initialized = false;
            curtree->nodep[i]->next->next->initialized = false;
          }
        }
      }
      bl_treevaluate(curtree, improve, reusertree, global, progress,
                      priortree, bestree, (initialvtrav_t)bl_initialvtrav );
      bestree->copy(bestree, curtree);
      dnaml_printree();
      summarize();
      if (trout) {
        dnaml_treeout(outtree, curtree, curtree->root);
      }
    }
  }

  if (usertree)
  {
    free(l0gl);
    for (i=0; i < shimotrees; i++)
      free(l0gf[i]);
    free(l0gf);
  }
  freetable();
  if (jumb < njumble)
    return;
  free(contribution);
  free(mp);
  for (i=0; i < endsite; i++)
    free(term[i]);
  free(term);
  for (i=0; i < endsite; i++)
    free(slopeterm[i]);
  free(slopeterm);
  for (i=0; i < endsite; i++)
    free(curveterm[i]);
  free(curveterm);
  freex(nonodes2, curtree->nodep);
  if (!usertree)
  {
    freex(nonodes2, bestree->nodep);
    freex(nonodes2, priortree->nodep);
    if (njumble > 1)
      freex(nonodes2, bestree2->nodep);
  }
  if (progress)
  {
    sprintf(progbuf, "\n\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree also written onto file \"%s\".\n", outtreename);
      print_progress(progbuf);
    }
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
}  /* maketree */


void clean_up(void)
{
  /* Free and/or close stuff */
  long i;

  free (rrate);
  free (probcat);
  free (rate);
  /* Seems to require freeing every time... */
  for (i = 0; i < spp; i++)
  {
    free (inputSequences[i]);
  }
  free (inputSequences);
  free (nayme);
  free (enterorder);
  free (category);
  free (weight);
  free (alias);
  free (ally);
  free (location);
  free (aliasweight);

  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
}   /* clean_up */


void dnamlrun(void)
{
  /* debug printout  JRMdebug  */
  /*
    printf("\nctgry: %i\n", ctgry);
    printf("categs: %li\n", categs);
    printf("rctgry: %i\n", rctgry);
    printf("rcategs: %li\n", rcategs);
    printf("auto_: %i\n", auto_);
    printf("freqsfrom: %i\n", freqsfrom);
    printf("gama: %i\n", gama);
    printf("global: %i\n", global);
    printf("hypstate: %i\n", hypstate);
    printf("improve: %i\n", improve);
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
    printf("reusertree: %i\n", reusertree);
    printf("weights: %i\n", weights);
    printf("printdata: %i\n", printdata);
    printf("dotdiff: %i\n", dotdiff);
    printf("progress: %i\n", progress);
    printf("treeprint: %i\n", treeprint);
    printf("interleaved: %i\n", interleaved);
    printf("mulsets: %i\n", mulsets);
    printf("datasets: %li\n", datasets);
  */
  /* do the work  */
/* debug:  need?   if (!usertree) nonodes2--;  */
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n", ith);
    }
    ttratio = ttratio0;
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      maketree();
    }
    fflush(outfile);
    fflush(outtree);
  }
} /* dnamlrun */


void dnaml(
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
  int SpeedAn,
  int GlobalRe,
  int RandInput,
  int RandNum,
  int Njumble,
  int OutRoot,
  int OutNum,
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

  //printf("Hello from DnaML!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "DnaML";

  funcs.tree_new = (tree_new_t)dnaml_tree_new;
  funcs.tree_init = (tree_init_t)dnaml_tree_init;
  nodesize = (long)sizeof(mldna_node);
  funcs.node_new = (node_new_t)mldna_node_new;
  funcs.node_init = (node_init_t)mldna_node_init;
  progname = argv[0];

  phylipinit(argc, argv, &funcs, true);

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
  //reusertree = false;
  //weights = false;
  //printdata = false;
  //dotdiff = true;
  //progress = true;
  //treeprint = true;
  //interleaved = true;
  //mulsets = false;
  //datasets = 1;

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
  //boolean OutRoot;
  //int OutNum;
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

  //printf("TreeUseMethod: %s\n", TreeUseMethod); // JRMdebug
  if (!strcmp(TreeUseMethod, "No"))
  {
    // No, use user trees in input file
    reusertree = false;
    usertree   = true;
  }
  else if (!strcmp(TreeUseMethod, "rearrange"))
  {
    // Yes, rearrange on user tree
    reusertree = true;
    usertree   = true;
  }
  else
  {
    // Yes
    reusertree = false;
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

  if (SpeedAn != 0)
  {
    improve = true;
  }
  else
  {
    improve = false;
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

  /* printf("values transfered\n");    JRMdebug  */
  /* printf("maxcategs: %i\n", maxcategs);  JRMdebug */

  /* create arrays and initialize them  */
  /* warning: if maxcategs ever changes, this has to change also  */
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

  /* printf("rates transfered\n"); // JRMdebug  */

  /* conditional modification of arrays */
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

  /* everything translated, start the run */
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
      //printf("intreename: %s\n", intreename);
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
  doinit();                       /* read in the options and the data */
  ttratio0 = ttratio;

  dnamlrun();                                   /* do the actual work */

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
} /* dnaml */


int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood */
#ifdef MAC
  argc = 1;                /* macsetup("DnaML", "");        */
  argv[0] = "DnaML";
#endif
  funcs.tree_new = (tree_new_t)dnaml_tree_new;
  funcs.tree_init = (tree_init_t)dnaml_tree_init;
  nodesize = (long)sizeof(mldna_node);
  funcs.node_new = (node_new_t)mldna_node_new;
  funcs.node_init = (node_init_t)mldna_node_init;
  phylipinit(argc, argv, &funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  mulsets = false;
  datasets = 1;
  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinit();                                             /* present the menu */
  ttratio0 = ttratio;
  if (ctgry)
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r",
              argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  dnamlrun();     /* do the actual work */

  clean_up();
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* DNA Maximum Likelihood */

/* End. */
