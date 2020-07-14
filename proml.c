/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Lucas Mix, Akiko Fuseki, Sean Lamont,
   Andrew Keeffe, Dan Fineman, and Patrick Colacurcio.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "prom_common.h"
#include "seq.h"
#include "ml.h"

typedef long vall[maxcategs];
typedef double contribarr[maxcategs];

typedef struct proml_tree {
  ml_tree ml_tree;
} proml_tree;

#ifndef OLDC
/* function prototypes */
tree * proml_tree_new(long, long);
void   proml_tree_init(tree* t, long nonodes, long spp);
void   proml_tree_setup(long, long);
void   getoptions(void);
void   makeprotfreqs(void);
void   doinit(void);
void   inputoptions(void);
void   input_protdata(long);
void   more_allocation(void);
void   prot_makevalues(long, pointarray, long, long, sequence, steptr);
void   prot_inittable(void);
void   getinput(void);
void   prot_slopecurv(node *, double, double *, double *, double *);
void   proml_tree_makenewv(tree* t, node *p);
void   make_pmatrix(double **, double **, double **, long,
                     double, double, double *, double **);
boolean   rearrange(node *, node *, double *);
void   proml_coordinates(node *, double, long *, double *);
void   proml_printree(void);
void   sigma(node *, double *, double *, double *);
void   describe(node *);
void   prot_reconstr(node *, long);
void   rectrav(node *, long, long);
void   summarize(void);
void   initpromlnode(tree *, node **, long, long, long *, long *,
                      initops, pointarray, Char *, Char *, FILE *);
void   dnaml_treeout(node *);
void   maketree(void);
void   proml_reroot(tree*) ;            // RSGbugfix: Name change.
double proml_tree_evaluate(tree*, node *, boolean);
void   prot_freetable(void);
void   proml_tree_nuview(tree*, node*);
void   promlrun(void);
void   proml(char * infilename, char * intreename, char * wgtsfilename,
             char * catsfilename, char * outfilename, char * outfileopt,
             char * outtreename, char * outtreeopt, char * TreeUseMethod,
             char * ProbModel, int UseLengths, int OneCat, int NumCats,
             double SiteRate1, double SiteRate2, double SiteRate3,
             double SiteRate4, double SiteRate5, double SiteRate6,
             double SiteRate7, double SiteRate8, double SiteRate9,
             char * RateVar, int AdjCor, double BlockLen, double CoeffVar,
             int NumRates, double HMMRate1, double HMMRate2, double HMMRate3,
             double HMMRate4, double HMMRate5, double HMMRate6,
             double HMMRate7, double HMMRate8, double HMMRate9,
             double HMMProb1, double HMMProb2, double HMMProb3,
             double HMMProb4, double HMMProb5, double HMMProb6,
             double HMMProb7, double HMMProb8, double HMMProb9,
             double InvarFract, int SitesWeight, int SpeedAn,
             int GlobalRe, int RandInput, int RandNum, int Njumble,
             int OutRoot, int OutNum, int MultData, int MultDSet,
             int NumSeqs, int InputSeq, int PrintData, int PrintInd,
             int PrintTree, int WriteTree, int RecHypo);
/* function prototypes */
#endif

boolean haslengths;
Char infilename[100], outfilename[100], intreename[100], outtreename[100],
      catfilename[100], weightfilename[100];
long nonodes2, weightsum, datasets, ith, jumb = 0;
long inseed, inseed0, parens;
boolean global, jumble, weights, trout, usertree, reusertree, ctgry, rctgry,
         auto_, hypstate, progress, mulsets, justwts, firstset, thorough,
         improve, smoothit, polishing, lngths, gama, invar, usepam;
tree *curtree, *bestree, *bestree2, *priortree;
node *qwhere;
double cv, alpha, lambda, invarfrac;
contribarr *contribution, like, nulike, clai;
double **term, **slopeterm, **curveterm;
longer seed;
char *progname;
char aachar[26]="ARNDCQEGHILKMFPSTWYVBZX?*-";

/* Local variables for maketree, propagated globally for C version: */
long k, nextsp, numtrees, maxwhich, mx, mx0, mx1, shimotrees;
double dummy, maxlogl;
boolean smoothed;
double **l0gf;
double *l0gl;
double **tbl;
Char ch, ch2;
long col;
vall *mp;


tree* proml_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(proml_tree));
  proml_tree_init(t, nonodes, spp);
  return t;
} /* proml_tree_new */


void proml_tree_init(tree * t, long nonodes, long spp)
{
  ml_tree_init(t, nonodes, spp);
  t->evaluate = proml_tree_evaluate;
  t->nuview = proml_tree_nuview;
  ((ml_tree*)t)->makenewv = proml_tree_makenewv;
} /* proml_tree_init */


void getoptions(void)
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean didchangecat, didchangercat;
  double probsum;
  char* string;

  putchar('\n');

  categs = 1;
  njumble = 1;
  rcategs = 1;
  outgrno = 1;

  lambda = 1.0;

  loopcount = 0;

  auto_ = false;
  ctgry = false;
  didchangecat = false;
  didchangercat = false;
  gama = false;
  global = false;
  hypstate = false;
  improve = true;
  interleaved = true;
  invar = false;
  jumble = false;
  lngths = false;
  outgropt = false;
  printdata = false;
  progress = true;
  rctgry = false;
  reusertree = false;
  treeprint = true;
  trout = true;
  usejtt = true;
  usepam = false;
  usepmb = false;
  usertree = false;
  weights = false;

  for (;;)
  {
    cleerhome();
    printf("Amino acid sequence Maximum Likelihood");
    printf(" method, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    if ( reusertree ) string = "Yes, rearrange on user tree";
    else if ( usertree ) string = "No, use user trees in input file";
    else string = "Yes";
    printf("  U                 Search for best tree?  %s\n", string);
    if (usertree && !reusertree)
    {
      printf("  L          Use lengths from user trees?  %s\n",
              (lngths ? "Yes" : "No"));
    }
    printf("  P    JTT, PMB or PAM probability model?  %s\n",
            usejtt ? "Jones-Taylor-Thornton" :
            usepmb ? "Henikoff/Tillier PMB" : "Dayhoff PAM");
    printf("  C   One category of substitution rates?");
    if (!ctgry || categs == 1)
      printf("  Yes\n");
    else
      printf("  %ld categories\n", categs);
    printf("  R           Rate variation among sites?");
    if (!rctgry)
      printf("  constant rate of change\n");
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
    printf("  O                        Outgroup root?  %s%3ld\n", (outgropt ?
            "Yes, at sequence number" : "No, use as outgroup species"),
           outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n", (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n", (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n", (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n", (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n", (trout ? "Yes" : "No"));
    printf("  5   Reconstruct hypothetical sequences?  %s\n", (hypstate ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("UPCRAWSGJOMI012345", ch) != NULL))
        || (usertree && ((strchr("UPLCRAWSOMI012345", ch) != NULL))))
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

        case 'P':
          if (usejtt)
          {
            usejtt = false;
            usepmb = true;
          }
          else
          {
            if (usepmb)
            {
              usepmb = false;
              usepam = true;
            }
            else
            {
              usepam = false;
              usejtt = true;
            }
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

        case 'U':
          if ( !usertree && !reusertree)
          {
            usertree = true;
          }
          else
          {
            if (!reusertree && usertree )
            {
              reusertree = true;
            }
            else
            {
              usertree = false;
              reusertree = false;
            }
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
    countup(&loopcount, 100);
  }
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
  prom_init_protmats();
}  /* getoptions */


void makeprotfreqs(void)
{
  /* calculate amino acid frequencies based on eigmat */
  long i, mineig;

  mineig = 0;
  for (i = 0; i <= 19; i++)
    if (fabs(eigmat[i]) < fabs(eigmat[mineig]))
      mineig = i;
  memcpy(freqaa, probmat[mineig], 20 * sizeof(double));
  for (i = 0; i <= 19; i++)
    freqaa[i] = fabs(freqaa[i]);
} /* makeprotfreqs */


void doinit(void)
{ /* initializes variables */
  fprintf(outfile, "\nAmino acid sequence Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n", VERSION);

  inputnumbers(&spp, &sites, &nonodes2, 1);

  if (!javarun)
  {
    getoptions();
  }

  makeprotfreqs();

  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);

  prom_allocrest();
}  /* doinit */


void inputoptions(void)
{
  long i;

  if (!firstset)
  {
    samenumsp(&sites, ith);
    prom_reallocsites();
  }
  if (firstset)
  {
    for (i = 0; i < sites; i++)
      category[i] = 1;
    for (i = 0; i < sites; i++)
      weight[i] = 1;
  }
  if (justwts || weights)
    inputweights(sites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if ((ctgry && categs > 1) && (firstset || !justwts))
  {
    inputcategs(0, sites, category, categs, "ProML");
    if (printdata)
      printcategs(outfile, sites, category, "Site categories");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, weight, "Sites");
  fprintf(outfile, "%s model of amino acid change\n\n",
          (usejtt ? "Jones-Taylor-Thornton" :
           usepmb ? "Henikoff/Tillier PMB" : "Dayhoff PAM"));
}  /* inputoptions */


void input_protdata(long chars)
{
  /* input the names and sequences for each species */
  /* used by proml */
  long i, j, k, l, basesread, basesnew;
  Char charstate;
  boolean allread, done;

  if (printdata)
    headings(chars, "Sequences", "---------");
  basesread = 0;
  basesnew = 0;
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
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((strchr("ABCDEFGHIKLMNPQRSTVWXYZ*?-", charstate)) == NULL)
          {
            printf("\nERROR:  Bad amino acid: %c at position %ld of species %ld.\n", charstate, j+1, i);
            if (charstate == '.')
            {
              printf("        Periods (.) may not be used as gap characters.\n");
              printf("        The correct gap character is (-).\n");
            }
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
        printf("\nERROR:  SEQUENCES OUT OF ALIGNMENT AT POSITION %ld.\n", j);
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
  if (!printdata)
    return;
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
        if (j > 1 && inputSequences[j - 1][k - 1] == inputSequences[0][k - 1])
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
}  /* input_protdata */


void more_allocation(void)
{
  long i;

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
}


void prot_makevalues(long categs, pointarray treenode, long endsite, long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips   */
  /* a version of makevalues2 found in seq.c */
  /* used by proml                             */
  long i, j, k, l;
  long b;

  for (k = 0; k < endsite; k++)
  {
    j = alias[k];
    for (i = 0; i < spp; i++)
    {
      for (l = 0; l < categs; l++)
      {
        memset(((prot_node*)treenode[i])->x[k][l], 0, sizeof(double)*20);
        switch (y[i][j - 1])
        {

          case 'A':
            ((prot_node*)treenode[i])->x[k][l][0] = 1.0;
            break;

          case 'R':
            ((prot_node*)treenode[i])->x[k][l][(long)arginine   - (long)alanine] = 1.0;
            break;

          case 'N':
            ((prot_node*)treenode[i])->x[k][l][(long)asparagine - (long)alanine] = 1.0;
            break;

          case 'D':
            ((prot_node*)treenode[i])->x[k][l][(long)aspartic   - (long)alanine] = 1.0;
            break;

          case 'C':
            ((prot_node*)treenode[i])->x[k][l][(long)cysteine   - (long)alanine] = 1.0;
            break;

          case 'Q':
            ((prot_node*)treenode[i])->x[k][l][(long)glutamine  - (long)alanine] = 1.0;
            break;

          case 'E':
            ((prot_node*)treenode[i])->x[k][l][(long)glutamic   - (long)alanine] = 1.0;
            break;

          case 'G':
            ((prot_node*)treenode[i])->x[k][l][(long)glycine    - (long)alanine] = 1.0;
            break;

          case 'H':
            ((prot_node*)treenode[i])->x[k][l][(long)histidine  - (long)alanine] = 1.0;
            break;

          case 'I':
            ((prot_node*)treenode[i])->x[k][l][(long)isoleucine - (long)alanine] = 1.0;
            break;

          case 'L':
            ((prot_node*)treenode[i])->x[k][l][(long)leucine    - (long)alanine] = 1.0;
            break;

          case 'K':
            ((prot_node*)treenode[i])->x[k][l][(long)lysine     - (long)alanine] = 1.0;
            break;

          case 'M':
            ((prot_node*)treenode[i])->x[k][l][(long)methionine - (long)alanine] = 1.0;
            break;

          case 'F':
            ((prot_node*)treenode[i])->x[k][l][(long)phenylalanine - (long)alanine] = 1.0;
            break;

          case 'P':
            ((prot_node*)treenode[i])->x[k][l][(long)proline    - (long)alanine] = 1.0;
            break;

          case 'S':
            ((prot_node*)treenode[i])->x[k][l][(long)serine     - (long)alanine] = 1.0;
            break;

          case 'T':
            ((prot_node*)treenode[i])->x[k][l][(long)threonine  - (long)alanine] = 1.0;
            break;

          case 'W':
            ((prot_node*)treenode[i])->x[k][l][(long)tryptophan - (long)alanine] = 1.0;
            break;

          case 'Y':
            ((prot_node*)treenode[i])->x[k][l][(long)tyrosine   - (long)alanine] = 1.0;
            break;

          case 'V':
            ((prot_node*)treenode[i])->x[k][l][(long)valine     - (long)alanine] = 1.0;
            break;

          case 'B':
            ((prot_node*)treenode[i])->x[k][l][(long)asparagine - (long)alanine] = 1.0;
            ((prot_node*)treenode[i])->x[k][l][(long)aspartic   - (long)alanine] = 1.0;
            break;

          case 'Z':
            ((prot_node*)treenode[i])->x[k][l][(long)glutamine  - (long)alanine] = 1.0;
            ((prot_node*)treenode[i])->x[k][l][(long)glutamic   - (long)alanine] = 1.0;
            break;

          case 'X':                /* unknown aa                            */
            for (b = 0; b <= 19; b++)
              ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
            break;

          case '?':                /* unknown aa                            */
            for (b = 0; b <= 19; b++)
              ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
            break;

          case '*':                /* stop codon symbol                    */
            for (b = 0; b <= 19; b++)
              ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
            break;

          case '-':                /* deletion event-absent data or aa */
            for (b = 0; b <= 19; b++)
              ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
            break;
        }
      }
    }
  }
}  /* prot_makevalues */


void getinput(void)
{
  /* reads the input data */
  if (!justwts || firstset)
    inputoptions();
  if (!justwts || firstset)
    input_protdata(sites);
  prom_makeweights();
  inittrees(nonodes2, spp);
  prot_makevalues(rcategs, curtree->nodep, endsite, spp, inputSequences, alias);
}  /* getinput */


void prot_freetable(void)
{
  long i, j, k, l;
  for (j = 0; j < rcategs; j++)
  {
    for (k = 0; k < categs; k++)
    {
      for (l = 0; l < 20; l++)
        free(ddpmatrix[j][k][l]);
      free(ddpmatrix[j][k]);
    }
    free(ddpmatrix[j]);
  }
  free(ddpmatrix);

  for (j = 0; j < rcategs; j++)
  {
    for (k = 0; k < categs; k++)
    {
      for (l = 0; l < 20; l++)
        free(dpmatrix[j][k][l]);
      free(dpmatrix[j][k]);
    }
    free(dpmatrix[j]);
  }
  free(dpmatrix);

  for (j = 0; j < rcategs; j++)
    free(tbl[j]);
  free(tbl);

  for ( i = 0 ; i < max_num_sibs ; i++ )
    prom_free_pmatrix(i);
  free(pmatrices);
} /* prot_freetable */


void prot_inittable(void)
{
  /* Define a lookup table. Precompute values and print them out in tables */
  /* Allocate memory for the pmatrices, dpmatices and ddpmatrices          */
  long i, j, k, l;
  double sumrates;

  /* Allocate memory for pmatrices, the array of pointers to pmatrices     */

  pmatrices = (double *****) Malloc ( spp * sizeof(double ****));

  /* Allocate memory for the first 2 pmatrices, the matrix of conversion   */
  /* probabilities, but only once per run (aka not on the second jumble.   */

    prom_alloc_pmatrix(0);
    prom_alloc_pmatrix(1);

  /*  Allocate memory for one dpmatrix, the first derivative matrix        */

  dpmatrix = (double ****) Malloc(rcategs * sizeof(double ***));
  for (j = 0; j < rcategs; j++)
  {
    dpmatrix[j] = (double ***) Malloc(categs * sizeof(double **));
    for (k = 0; k < categs; k++)
    {
      dpmatrix[j][k] = (double **) Malloc(20 * sizeof(double *));
      for (l = 0; l < 20; l++)
        dpmatrix[j][k][l] = (double *) Malloc(20 * sizeof(double));
    }
  }

  /*  Allocate memory for one ddpmatrix, the second derivative matrix      */
  ddpmatrix = (double ****) Malloc(rcategs * sizeof(double ***));
  for (j = 0; j < rcategs; j++)
  {
    ddpmatrix[j] = (double ***) Malloc(categs * sizeof(double **));
    for (k = 0; k < categs; k++)
    {
      ddpmatrix[j][k] = (double **) Malloc(20 * sizeof(double *));
      for (l = 0; l < 20; l++)
        ddpmatrix[j][k][l] = (double *) Malloc(20 * sizeof(double));
    }
  }

  /* Allocate memory and assign values to tbl, the matrix of possible rates*/

  tbl = (double **) Malloc(rcategs * sizeof(double *));
  for (j = 0; j < rcategs; j++)
    tbl[j] = (double *) Malloc(categs * sizeof(double));

  for (j = 0; j < rcategs; j++)
    for (k = 0; k < categs; k++)
      tbl[j][k] = rrate[j]*rate[k];

  sumrates = 0.0;
  for (i = 0; i < endsite; i++)
  {
    for (j = 0; j < rcategs; j++)
      sumrates += aliasweight[i] * probcat[j]
        * tbl[j][category[alias[i] - 1] - 1];
  }
  sumrates /= (double)sites;
  for (j = 0; j < rcategs; j++)
    for (k = 0; k < categs; k++)
    {
      tbl[j][k] /= sumrates;
    }

  if(jumb > 1)
    return;

  if (gama)
  {
    fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
    fprintf(outfile,
    " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
  }
  if (rcategs > 1)
  {
    fprintf(outfile, "\nState in HMM   Rate of change    Probability\n\n");
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
      fprintf(outfile,
     "Expected length of a patch of sites having the same rate = %8.3f\n",
             1/lambda);
      putc('\n', outfile);
    }
  }
  if (categs > 1)
  {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (k = 0; k < categs; k++)
      fprintf(outfile, "%9ld%16.3f\n", k+1, rate[k]);
  }
  if ((rcategs  > 1) || (categs >> 1))
    fprintf(outfile, "\n\n");
}  /* prot_inittable */


void make_pmatrix(double **matrix, double **dmat, double **ddmat,
                  long derivative, double lz, double rat,
                  double *eigmat, double **probmat)
{
  /* Computes the R matrix such that matrix[m][l] is the joint probability */
  /* of m and l.                                                           */
  /* Computes a P matrix such that matrix[m][l] is the conditional         */
  /* probability of l given m.  This is accomplished by dividing all terms */
  /* in the R matrix by freqaa[m], the frequency of m.                     */

  long k, l, m;                 /* (l) original character state */
                                /* (m) final    character state */
                                /* (k) lambda counter           */
  double p0, p1, p2, q;
  double elambdat[20], delambdat[20], ddelambdat[20];
                                /* exponential term for matrix  */
                                /* and both derivative matrices */
  for (k = 0; k <= 19; k++)
  {
    elambdat[k] = exp(lz * rat * eigmat[k]);
    if(derivative != 0)
    {
        delambdat[k] = (elambdat[k] * rat * eigmat[k]);
        ddelambdat[k] = (delambdat[k] * rat * eigmat[k]);
    }
   }
  for (m = 0; m <= 19; m++)
  {
    for (l = 0; l <= 19; l++)
    {
      p0 = 0.0;
      p1 = 0.0;
      p2 = 0.0;
      for (k = 0; k <= 19; k++)
      {
        q = probmat[k][m] * probmat[k][l];
        p0 += (q * elambdat[k]);
        if(derivative !=0)
        {
          p1 += (q * delambdat[k]);
          p2 += (q * ddelambdat[k]);
        }
      }
      matrix[m][l] = p0 / freqaa[m];
      if(derivative != 0)
      {
        dmat[m][l] = p1 / freqaa[m];
        ddmat[m][l] = p2 / freqaa[m];
      }
    }
  }
}  /* make_pmatrix */


void proml_tree_nuview(tree* t, node *p)
{
  long i, j, k, l, m, num_sibs, sib_index;
  node *sib_ptr, *sib_back_ptr;
  psitelike prot_xx, x2;
  double lw, prod7;
  double **pmat;
  double maxx;
  double correction;

  /* Figure out how many siblings the current node has  */
  /* and be sure that pmatrices is large enough         */
  num_sibs = count_sibs(p);
  for (i = 0; i < num_sibs; i++)
    if (pmatrices[i] == NULL)
      prom_alloc_pmatrix(i);

  /* Make pmatrices for all possible combinations of category, rcateg      */
  /* and sib                                                               */
  sib_ptr = p;                                /* return to p */
  for (sib_index=0; sib_index < num_sibs; sib_index++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    lw = sib_back_ptr->v;

    for (j = 0; j < rcategs; j++)
      for (k = 0; k < categs; k++)
        make_pmatrix(pmatrices[sib_index][j][k], NULL, NULL, 0, lw, tbl[j][k], eigmat, probmat);
  }

  for (i = 0; i < endsite; i++)
  {
    maxx = 0;
    correction = 0;

    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {

      /* initialize to 1 all values of prot_xx */
      for (m = 0; m <= 19; m++)
        prot_xx[m] = 1;

      sib_ptr = p;                        /* return to p */
      /* loop through all sibs and calculate likelihoods for all possible*/
      /* amino acid combinations                                         */
      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        sib_ptr      = sib_ptr->next;
        sib_back_ptr = sib_ptr->back;

        if ( j == 0)
          correction += ((ml_node*)sib_back_ptr)->underflows[i];

        memcpy(x2, ((prot_node*)sib_back_ptr)->x[i][j], sizeof(psitelike));
        pmat = pmatrices[sib_index][j][k];
        for (m = 0; m <= 19; m++)
        {
          prod7 = 0;
          for (l = 0; l <= 19; l++)
            prod7 += (pmat[m][l] * x2[l]);
          prot_xx[m] *= prod7;
          if ( prot_xx[m] > maxx && sib_index == (num_sibs - 1))
            maxx = prot_xx[m];
        }
      }
      /* And the final point of this whole function: */
      memcpy(((prot_node*)p)->x[i][j], prot_xx, sizeof(psitelike));
    }
    ((ml_node*)p)->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE )
      fix_protx(((prot_node*)p), i, maxx, rcategs);
    ((ml_node*)p)->underflows[i] += correction;
  }

  p->initialized = true;
}  /* proml_tree_nuview */


void prot_slopecurv(node *p, double y, double *like, double *slope, double *curve)
{
  /* compute log likelihood, slope and curvature at node p */
  long i, j, k, l, m, lai;
  double sum, sumc, sumterm, lterm, sumcs, sumcc, sum2, slope2, curve2;
  double frexm = 0;                        /* frexm = freqaa[m]*x1[m] */
                                        /* frexml = frexm*x2[l]    */
  double prod4m, prod5m, prod6m;        /* elements of prod4-5 for */
                                        /* each m                   */
  double **pmat, **dpmat, **ddpmat;        /* local pointers to global*/
                                        /* matrices                   */
  double prod4, prod5, prod6;
  contribarr thelike, nulike, nuslope, nucurve,
    theslope, thecurve, clai, cslai, cclai;
  node *q;
  psitelike x1, x2;

  q = p->back;
  sum = 0.0;
  for (j = 0; j < rcategs; j++)
  {
    for (k = 0; k < categs; k++)
    {
      make_pmatrix(pmatrices[0][j][k], dpmatrix[j][k], ddpmatrix[j][k], 2, y, tbl[j][k], eigmat, probmat);
    }
  }
  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      memcpy(x1, ((prot_node*)p)->x[i][j], sizeof(psitelike));
      memcpy(x2, ((prot_node*)q)->x[i][j], sizeof(psitelike));
      pmat = pmatrices[0][j][k];
      dpmat = dpmatrix[j][k];
      ddpmat = ddpmatrix[j][k];
      prod4 = 0.0;
      prod5 = 0.0;
      prod6 = 0.0;
      for (m = 0; m <= 19; m++)
      {
        prod4m = 0.0;
        prod5m = 0.0;
        prod6m = 0.0;
        frexm = x1[m] * freqaa[m];
        for (l = 0; l <= 19; l++)
        {
          prod4m += x2[l] * pmat[m][l];
          prod5m += x2[l] * dpmat[m][l];
          prod6m += x2[l] * ddpmat[m][l];
        }
        prod4 += frexm * prod4m;
        prod5 += frexm * prod5m;
        prod6 += frexm * prod6m;
      }
      term[i][j] = prod4;
      slopeterm[i][j] = prod5;
      curveterm[i][j] = prod6;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * term[i][j];
    if (sumterm <= 0.0)
        sumterm = 0.00000001;        /* ?shouldn't get here?? */
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
    {
      term[i][j] = term[i][j] / sumterm;
      slopeterm[i][j] = slopeterm[i][j] / sumterm;
      curveterm[i][j] = curveterm[i][j] / sumterm;
    }
    sum += (aliasweight[i] * lterm);
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
        nulike[j]  = ((1.0 - lambda) * thelike[j]  + sumc) *  clai[j];
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs) * clai[j]
                   + ((1.0 - lambda) * thelike[j]  + sumc) *  cslai[j];
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc) * clai[j]
             + 2.0 * ((1.0 - lambda) * theslope[j] + sumcs) * cslai[j]
                   + ((1.0 - lambda) * thelike[j]  + sumc) *  cclai[j];
      }
    }
    else
    {
      for (j = 0; j < rcategs; j++)
      {
        nulike[j]  = ((1.0 - lambda) * thelike[j]  + sumc);
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

  /* Expressed in terms of *slope to prevent overflow */
  (*curve) = curve2 / sum2 - *slope * *slope;
} /* prot_slopecurv */


void proml_tree_makenewv(tree* t, node *p)
{
  /* Newton-Raphson algorithm improvement of a branch length */
  long it, ite;
  double y, yold=0, yorig, like, slope, curve, oldlike=0;
  boolean done, firsttime, better;
  node *q;

  q = p->back;
  y = p->v;
  yorig = y;
  done = false;
  firsttime = true;
  it = 1;
  ite = 0;
  while ((it < iterations) && (ite < 20) && (!done))
  {
    prot_slopecurv(p, y, &like, &slope, &curve);
    better = false;
    if (firsttime)             /* if no older value of y to compare with */
    {
      yold = y;
      oldlike = like;
      firsttime = false;
      better = true;
    }
    else
    {
      if (like > oldlike)      /* update the value of yold if it was better */
      {
        yold = y;
        oldlike = like;
        better = true;
        it++;
      }
    }
    if (better)
    {
      y = y + slope/fabs(curve);    /* Newton-Raphson, forced uphill-wards */
      if (y < epsilon)
        y = epsilon;
    }
    else
    {
      if (fabs(y - yold) < epsilon)
        ite = 20;
      y = (y + (7 * yold)) / 8;         /* retract 87% of way back */
    }
    ite++;
    done = fabs(y-yold) < 0.1*epsilon;
  }
  smoothed = (fabs(yold-yorig) < epsilon) && (yorig > 1000.0*epsilon);
  p->v = yold;                   /* the last one that had better likelihood */
  q->v = yold;
  ((tree*)t)->score = oldlike;
}  /* proml_tree_makenewv */


double proml_tree_evaluate(tree *t, node *p, boolean saveit)
{
  contribarr tterm;
  double sum, sum2, sumc, y, prod4, prodl, frexm, sumterm, lterm;
  double **pmat;
  long i, j, k, l, m, lai;
  node *q;
  psitelike x1, x2;

  generic_tree_evaluate(t, p, saveit);
  sum = 0.0;
  q = p->back;
  y = p->v;
  for (j = 0; j < rcategs; j++)
    for (k = 0; k < categs; k++)
      make_pmatrix(pmatrices[0][j][k], NULL, NULL, 0,
                    y, tbl[j][k], eigmat, probmat);
  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      memcpy(x1, ((prot_node*)p)->x[i][j], sizeof(psitelike));
      memcpy(x2, ((prot_node*)q)->x[i][j], sizeof(psitelike));
      prod4 = 0.0;
      pmat = pmatrices[0][j][k];
      for (m = 0; m <= 19; m++)
      {
        prodl = 0.0;
        for (l = 0; l <= 19; l++)
          prodl += (pmat[m][l] * x2[l]);
        frexm = x1[m] * freqaa[m];
        prod4 += (prodl * frexm);
      }
      tterm[j] = prod4;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    if (sumterm < 0.0)
      sumterm = 0.00000001;        /* ??? */
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, rcategs * sizeof(double));
    if (saveit && !auto_ && usertree && (which <= shimotrees) && !reusertree)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
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
      memcpy(clai, contribution[lai - 1], rcategs * sizeof(double));
      for (j = 0; j < rcategs; j++)
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
  curtree->score = sum;

  if (!saveit || auto_ || !usertree || reusertree)
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
}  /* proml_tree_evaluate */


void proml_coordinates(node *p, double lengthsum, long *tipy, double *tipmax)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

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
    proml_coordinates(q->back, lengthsum + xx, tipy, tipmax);
    q = q->next;
  } while ((p == curtree->root || p != q)
           && (p != curtree->root || p->next != q));
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
}  /* proml_coordinates */


void proml_printree(void)
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
  proml_coordinates(curtree->root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline2(i, scale, curtree);
  putc('\n', outfile);
}  /* proml_printree */


void sigma(node *p, double *sumlr, double *s1, double *s2)
{
  /* compute standard deviation */
  double tt, aa, like, slope, curv;

  prot_slopecurv(p, p->v, &like, &slope, &curv);
  tt = p->v;
  p->v = epsilon;
  p->back->v = epsilon;
  aa = curtree->evaluate(curtree, p, false);
  p->v = tt;
  p->back->v = tt;
  (*sumlr) = curtree->evaluate(curtree, p, false) - aa;
  if (curv < -epsilon)
  {
    (*s1) = p->v + (-slope - sqrt(slope * slope -  3.841 * curv)) / curv;
    (*s2) = p->v + (-slope + sqrt(slope * slope -  3.841 * curv)) / curv;
  }
  else
  {
    (*s1) = -1.0;
    (*s2) = -1.0;
  }
}  /* sigma */


void describe(node *p)
{
  /* print out information for one branch */
  long i, num_sibs;
  node *q, *sib_ptr;
  double sumlr, sigma1, sigma2;

  if (!p->tip && !p->initialized)
    curtree->nuview(curtree, p);
  if (!p->back->tip && !p->back->initialized)
    curtree->nuview(curtree, p->back);
  q = p->back;
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
  fprintf(outfile, "%15.5f", q->v);
  if (!usertree || (usertree && !lngths) || p->iter || reusertree)
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
        fprintf(outfile, "%9.5f", sigma2);
      fprintf(outfile, ",%12.5f", sigma1);
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
      describe(sib_ptr->back);
    }
  }
}  /* describe */


void prot_reconstr(node *p, long n)
{
  /* reconstruct and print out acid at site n+1 at node p */
  long i, j, k, first, num_sibs = 0;
  double f, sum, xx[20];
  node *q = NULL;

  if (p->tip)
    putc(inputSequences[p->index-1][n], outfile);
  else
  {
    num_sibs = count_sibs(p);
    if ((ally[n] == 0) || (location[ally[n]-1] == 0))
      putc('.', outfile);
    else
    {
      j = location[ally[n]-1] - 1;
      sum = 0;
      for (i = 0; i <= 19; i++)
      {
        f = ((prot_node*)p)->x[j][mx-1][i];
        if (!p->tip)
        {
          q = p;
          for (k = 0; k < num_sibs; k++)
          {
            q = q->next;
            f *= ((prot_node*)q)->x[j][mx-1][i];
          }
        }
        if (f > 0.0)
          f = exp(log(f)/(num_sibs-1.0));
        xx[i] = f * freqaa[i];
        sum += xx[i];
      }
      for (i = 0; i <= 19; i++)
        xx[i] /= sum;
      first = 0;
      for (i = 0; i <= 19; i++)
        if (xx[i] > xx[first])
          first = i;
      if (xx[first] > 0.95)
        putc(aachar[first], outfile);
      else
        putc(tolower(aachar[first]), outfile);
      if (rctgry && rcategs > 1)
        mx = mp[n][mx - 1];
      else
        mx = 1;
    }
  }
} /* prot_reconstr */


void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  node *q;
  long i;

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
    prot_reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip)
  {
    for (q = p->next; q != p; q = q->next)
    {
      rectrav(q->back, m, n);
    }
  }
  mx1 = mx;
}  /* rectrav */


void summarize(void)
{
  /* print out branch length information and node numbers */
  long i, j, num_sibs;
  long mm = 0;   // RSGnote: Formerly may have been referenced before being initialized.
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;
  node   *sib_ptr;

  if (!treeprint)
    return;
  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree->score);
  fprintf(outfile, "\n Between        And            Length");
  if (!(usertree && lngths && haslengths ) || reusertree)
    fprintf(outfile, "      Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, " -------        ---            ------");
  if (!(usertree && lngths && haslengths) || reusertree)
    fprintf(outfile, "      ------- ---------- ------");
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
  num_sibs = count_sibs(curtree->root);
  sib_ptr  = curtree->root;
  for (i=0; i < num_sibs; i++)
  {
    sib_ptr = sib_ptr->next;
    describe(sib_ptr->back);
  }

  fprintf(outfile, "\n");
  if (!(usertree && lngths && haslengths)|| reusertree)
  {
    fprintf(outfile, "     *  = significantly positive, P < 0.05\n");
    fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n");
  }
  curtree->evaluate(curtree, curtree->root, false);
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
    fprintf(outfile, "Most probable category at each site if > 0.95");
    fprintf(outfile, " probability (\".\" otherwise)\n\n");
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
      mx0 = mx1;
    }
  }
}  /* summarize */


void initpromlnode(tree *treep, node **p, long len, long nodei, long *ntips,
                    long *parens, initops whichinit, pointarray treenode,
                    Char *str, Char *ch, FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  (void)len;                            // RSGnote: Parameter never used.
  (void)ntips;                          // RSGnote: Parameter never used.
  (void)treenode;                       // RSGnote: Parameter never used.

  switch (whichinit)
  {
    case bottom:
      *p = treep->get_forknode(treep, nodei);
      (*p)->index = nodei;
      (*p)->tip = false;
      ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs);
      treep->nodep[(*p)->index - 1] = (*p);
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      ((ml_node*)*p)->allocx((ml_node*)*p, endsite, rcategs);
      (*p)->index = nodei;
      break;
    case tip:
      match_names_to_data(str, treep->nodep, p, spp);
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
    default:      /* cases hslength, treewt, unittrwt */
      break;      /* should never occur               */
  }
} /* initpromlnode */


void dnaml_treeout(node *p)
{
  /* write out file with representation of final tree2 */
  /* Only works for bifurcations! */
  long i, n, w;
  Char c;
  double x;
  node *q;
  boolean inloop;

  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)
    {
      if (nayme[p->index-1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)
    {
      c = nayme[p->index-1][i];
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

    inloop = false;
    q = p->next;
    do  {
      if (inloop)
      {
        putc(',', outtree);
        col++;
        if (col > 45)
        {
          putc('\n', outtree);
          col = 0;
        }
      }
      inloop = true;
      dnaml_treeout(q->back);
      q = q->next;
    } while ((p == curtree->root || p != q) && (p != curtree->root
             || p->next != q));
    putc(')', outtree);
    col++;
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
  if (p == curtree->root)
    fprintf(outtree, ";\n");
  else
  {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* dnaml_treeout */


void proml_reroot(tree* t)              // RSGbugfix: Name change.
{
  node *q;
  double newl;
  node *r = t->root;
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
    assert(r->back == NULL);       // RSGnote: This assumes the FORKRING being
                       //  manipulated has the ROOT FORKNODE pointing to NULL.

    newl = r->next->oldlen + r->next->next->oldlen;
    r->next->back->oldlen = newl;
    r->next->next->back->oldlen = newl;

    newl = r->next->v + r->next->next->v;
    r->next->back->v = newl;
    r->next->next->back->v = newl;

    r->next->back->back = r->next->next->back;
    r->next->next->back->back = r->next->back;

    t->release_fork(t, r);
  }

  t->root = t->nodep[0]->back; // Reset ROOT; moved from line just after call to PROML_REROOT.
} /* proml_reroot */


void maketree(void)
{
  long i;
  boolean dummy_first, goteof;
  pointarray dummy_treenode = NULL;
  long nextnode;
  double bestyet;
  node *p;

  prot_inittable();

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
      l0gf[i] = (double *) Malloc(endsite * sizeof(double));
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }

    /* This taken out of tree read, used to be [spp-1], but referring
       to [0] produces output identical to what the pre-modified dnaml
       produced. */
    for ( which = 1 ; which <= numtrees ; which++)
    {

      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      preparetree(curtree);
      treeread(curtree, intree, &curtree->root, dummy_treenode,
                &goteof, &dummy_first, &nextnode, &haslengths,
                initpromlnode, false, nonodes2);
      fixtree(curtree);
      proml_reroot(curtree);                     // RSGbugfix: Name change.

      if (goteof && (which <= numtrees))
      {
        /* if we hit the end of the file prematurely */
        printf ("\nERROR:  Trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit (-1);
      }

      // RSGbugfix: Reset of ROOT moved to inside PROML_REROOT.

      if ( outgropt )
        curtree->root = curtree->nodep[outgrno - 1]->back;
      ml_treevaluate(curtree, improve, reusertree, global, progress,
                      priortree, bestree, ml_initialvtrav);
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
      proml_printree();
      summarize();
      if (trout)
      {
        col = 0;
        dnaml_treeout(curtree->root);
      }
      if(which < numtrees)
      {
        prot_freex_notip(nextnode, curtree->nodep);
      }
      else
        nonodes2 = nextnode;
    }
    FClose(intree);
    putc('\n', outfile);
    if (!auto_ && numtrees > 1 && weightsum > 1 && !reusertree)
      standev2(numtrees, maxwhich, 0, endsite-1, maxlogl,
                l0gl, l0gf, aliasweight, seed);
  }
  else
  {
    /* If there's no input user tree, */
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    if (progress)
    {
      char *outstr = "\nAdding species:\n";
      print_progress(outstr);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }

    nextsp = 3;
    polishing = false;
    destruct_tree(curtree);
    buildsimpletree(curtree, enterorder);
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    smoothit = improve;
    thorough = true;
    nextsp = 4;

    while (nextsp <= spp)
    {
      curtree->copy(curtree, priortree);
      k = generic_tree_findemptyfork(curtree);
      p = curtree->get_fork(curtree, k);
      ml_hookup(curtree->nodep[enterorder[nextsp-1]-1], p);
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
      if (smoothit)
        curtree->copy(curtree, priortree);
      curtree->addtraverse(curtree, p, curtree->root, true, qwhere,
                            &bestyet, bestree, thorough);
      if (smoothit)
        bestree->copy(bestree, curtree);
      else
      {
        smoothit = true;
        curtree->insert_(curtree, curtree->nodep[enterorder[nextsp - 1] - 1], qwhere, false);
        smoothit = false;
        bestyet = curtree->score;
      }
      if (progress)
      {
        writename(nextsp - 1, 1, enterorder);
        phyFillScreenColor();
      }

      if (global && nextsp == spp)
      {
        curtree->globrearrange(curtree, progress, smoothit);
      }
      else
      {
        curtree->locrearrange(curtree, curtree->nodep[enterorder[0] - 1], smoothit, priortree, bestree) ;
      }

      if (!smoothit)
      {
        curtree->smoothall(curtree, curtree->root);
        bestyet = curtree->score;
      }
      nextsp++;
    }
    if (global && progress)
    {
      char *outstr = "\n";
      print_progress(outstr);
    }
    curtree->copy(curtree, bestree);
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
      curtree->root = curtree->nodep[outgrno - 1]->back;
      ml_treevaluate(curtree, improve, reusertree, global, progress, priortree, bestree, ml_initialvtrav );
      proml_printree();
      summarize();
      if (trout)
      {
        col = 0;
        dnaml_treeout(curtree->root);
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
  prot_freetable();
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
  prom_free_all_x(nonodes2, curtree->nodep);
  if (!usertree)
  {
    prom_free_all_x(nonodes2, bestree->nodep);
    prom_free_all_x(nonodes2, priortree->nodep);
    if (njumble > 1)
      prom_free_all_x(nonodes2, bestree2->nodep);
  }
  if (progress)
  {
    //progfile = fopen("progress.txt", "a");
    sprintf(progbuf, "\n\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree written onto file \"%s\".\n", outtreename);
      print_progress(progbuf);
    }
    sprintf(progbuf, "\n");
    print_progress(progbuf);
    fflush(progfile);
  }
}  /* maketree */


void promlrun(void)
{
  // debug printout // JRMdebug
  /*
    printf("\nctgry: %i\n", ctgry);
    printf("categs: %li\n", categs);
    printf("rctgry: %i\n", rctgry);
    printf("rcategs: %li\n", rcategs);
    printf("auto_: %i\n", auto_);
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
    printf("usertree: %i\n", usertree);
    printf("reusertree: %i\n", reusertree);
    printf("weights: %i\n", weights);
    printf("printdata: %i\n", printdata);
    printf("progress: %i\n", progress);
    printf("treeprint: %i\n", treeprint);
    printf("interleaved: %i\n", interleaved);
    printf("mulsets: %i\n", mulsets);
    printf("datasets: %li\n", datasets);
    printf("usejtt: %i\n", usejtt);
    printf("usepam: %i\n", usepam);
    printf("usepmb: %i\n", usepmb);
    fflush(stdout);
  */

  // do the work
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    getinput();
    more_allocation();

    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++) {
      max_num_sibs = 0;
      maketree();
    }
    fflush(outfile);
    fflush(outtree);
  }
}


void proml(
  char * infilename,
  char * intreename,
  char * wgtsfilename,
  char * catsfilename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  char * TreeUseMethod,
  char * ProbModel,
  int UseLengths,
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
  int RecHypo)
{
  initdata funcs;

  //printf("Hello from ProML!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Proml";

  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = prot_node_new;
  funcs.tree_new = proml_tree_new;

  progname = argv[0];

  phylipinit(argc, argv, &funcs, true);

  //printf("init done\n"); // JRMdebug

  /* // JRMdebug
  // internal variables
  //categs = 1;
  //njumble = 1;
  //rcategs = 1;
  //outgrno = 1;

  //lambda = 1.0;

  ##loopcount = 0;

  //auto_ = false;
  //ctgry = false;
  ##didchangecat = false;
  ##didchangercat = false;
  //gama = false;
  //global = false;
  //hypstate = false;
  //improve = false;
  //interleaved = true;
  //invar = false;
  //jumble = false;
  //lngths = false;
  //outgropt = false;
  //printdata = false;
  //progress = true;
  //rctgry = false;
  //reusertree = false;
  //treeprint = true;
  //trout = true;
  //usejtt = true;
  //usepam = false;
  //usepmb = false;
  //usertree = false;
  //weights = false;

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
  //String ProbModel;
  //boolean UseLengths;
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
  //boolean RecHypo;
  */

  // transfer java data to local variables
  if (!strcmp(ProbModel, "JTT"))
  {
    usejtt = true;
    usepmb = false;
    usepam = false;
  }
  else if (strcmp(ProbModel, "PMB"))
  {
    usejtt = false;
    usepmb = true;
    usepam = false;
  }
  else //"PAM"
  {
    usejtt = false;
    usepmb = false;
    usepam = true;
  }

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

  if (MultDSet != 0) //??
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
    progress =  true;
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
  categs    = NumCats;
  rcategs   = NumRates;
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
  //printf("rate transfered\n"); // JRMdebug

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
  //printf("rrate transfered\n"); // JRMdebug

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
  //printf("probcat transfered\n"); // JRMdebug

  //printf("rates transfered\n"); // JRMdebug
  // debug printout // JRMdebug
  /*
    printf("\nctgry: %i\n", ctgry);
    printf("categs: %li\n", categs);
    printf("rctgry: %i\n", rctgry);
    printf("rcategs: %li\n", rcategs);
    printf("auto_: %i\n", auto_);
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
    printf("usertree: %i\n", usertree);
    printf("reusertree: %i\n", reusertree);
    printf("weights: %i\n", weights);
    printf("printdata: %i\n", printdata);
    printf("progress: %i\n", progress);
    printf("treeprint: %i\n", treeprint);
    printf("interleaved: %i\n", interleaved);
    printf("mulsets: %i\n", mulsets);
    printf("datasets: %li\n", datasets);
    printf("usejtt: %i\n", usejtt);
    printf("usepam: %i\n", usepam);
    printf("usepmb: %i\n", usepmb);
    fflush(stdout);
  */

  // conditional modification of arrays
  int i;
  if (gama) {
    initgammacat(rcategs, alpha, rrate, probcat);
  }
  else {
    if (invar) {
      initgammacat(rcategs-1, alpha, rrate, probcat);
      for (i = 0; i < rcategs-1; i++)
      {
        probcat[i] = probcat[i]*(1.0-invarfrac);
      }
      probcat[rcategs-1] = invarfrac;
      rrate[rcategs-1] = 0.0;
    }
  }
  //printf("gamma initialized\n"); // JRMdebug
  //fflush(stdout); //JRMdebug

  prom_init_protmats();
  //printf("prom_init_protmat done\n"); // JRMdebug
  //fflush(stdout); //JRMdebug

  // everything translated, start the run
  infile = fopen(infilename, "r");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

  if (usertree)
  {
    intree = fopen(intreename, "r");
  }

  if (ctgry)
  {
    //printf("catsfilename: %s\n", catsfilename); // JRMdebug
    catfile = fopen(catsfilename, "r");
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

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  //printf("calling doinit\n");  //JRMdebug
  //fflush(stdout); //JRMdebug
  doinit();

  //printf("calling promlrun\n");  //JRMdebug
  //fflush(stdout); //JRMdebug
  promlrun();  // do the actual work

  FClose(infile);
  FClose(outfile);

  if (ctgry) {
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

  if (progress)
  {
    FClose(progfile);
  }

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* Protein Maximum Likelihood */
  initdata funcs;
#ifdef MAC
  argc = 1;             /* macsetup("ProML", "");        */
  argv[0] = "ProML";
#endif
  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = prot_node_new;
  funcs.tree_new = proml_tree_new;
  progname = argv[0];

  phylipinit(argc, argv, &funcs, false);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  mulsets = false;
  firstset = true;

  doinit();

  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if (ctgry)
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);

  promlrun();

  prom_clean_up();
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* Protein Maximum Likelihood */


// End.
