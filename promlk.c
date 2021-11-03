/* Version 4.0. (c) Copyright 1986-2013 by the University of Washington
  and by Joseph Felsenstein.  Written by Joseph Felsenstein and Lucas Mix.
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

typedef struct promlk_tree {
  ml_tree ml_tree;
} promlk_tree;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   makeprotfreqs(void);
void   doinit(void);
void   inputoptions(void);
void   input_protdata(long);
void   more_allocation(void);
void   prot_makevalues(long, pointarray, long, long, sequence, steptr);
void   getinput(void);
void   prot_inittable(void);
void   make_pmatrix(double **, double **, double **, long, double, double, double *, double **);
void   promlk_tree_nuview(tree*t, node *p);
void   restoradd(node *, node *, node *, double);
void   nodeinit(node *);
void   initrav(node *);
void   travinit(node *);
void   travsp(node *);
void   treevaluate(void);
void   promlk_coordinates(node *, long *);
void   promlk_drawline(long, double);
void   promlk_printree(void);
void   describe(node *);
void   prot_reconstr(node *, long);
void   rectrav(node *, long, long);
void   summarize(void);
void   promlk_treeout(node *);
void   initpromlnode(tree *, node **, long, long, long *, long *, initops, pointarray, Char *, Char *, FILE *);
void   tymetrav(node *, double *);
void   maketree(void);
void   prot_freetable(void);
tree * promlk_tree_new(long nonodes, long spp);
void   promlk_tree_init(tree *t, long nonodes, long spp);
double promlk_tree_evaluate(tree* t, node *p, boolean dummy);
void   promlk_tree_smooth(tree* t, node* p);
void   promlkrun(void);
void   promlk(char * infilename, char * intreename, char * wgtsfilename, char * catsfilename,
              char * outfilename, char * outfileopt, char * outtreename, char * outtreeopt,
              char * TreeUseMethod, char * ProbModel, int UseLengths, int OneCat, int NumCats,
              double SiteRate1, double SiteRate2, double SiteRate3, double SiteRate4, double SiteRate5,
              double SiteRate6, double SiteRate7, double SiteRate8, double SiteRate9, char * RateVar,
              int AdjCor, double BlockLen, double CoeffVar, int NumRates, double HMMRate1, double HMMRate2,
              double HMMRate3, double HMMRate4, double HMMRate5, double HMMRate6, double HMMRate7, double HMMRate8,
              double HMMRate9, double HMMProb1, double HMMProb2, double HMMProb3, double HMMProb4, double HMMProb5,
              double HMMProb6, double HMMProb7, double HMMProb8, double HMMProb9, double InvarFract,
              int SitesWeight, int GlobalRe, int RandInput, int RandNum, int Njumble,
              int MultData, int MultDSet, int NumSeqs, int InputSeq, int PrintData,
              int PrintInd, int PrintTree, int WriteTree, int RecHypo);
/* function prototypes */
#endif

long jumb = 0, nonodes = 0;
double **tbl;
Char infilename[100], outfilename[100], intreename[100], outtreename[100], catfilename[100], weightfilename[100];
long weightsum, datasets, ith, numtrees, shimotrees;

/*  sites = number of sites in actual sequences
  numtrees = number of user-defined trees */
long inseed, inseed0, mx, mx0, mx1;
boolean global, jumble, lngths, trout, usertree, weights, rctgry, ctgry, auto_, progress, mulsets, firstset, hypstate, smoothit, polishing, justwts, gama, invar, usepam;

boolean lengthsopt = false;             /* Use lengths in user tree option */
boolean lngths     = false;             /* Actually use lengths (depends on
                                           each input tree) */
node *qwhere;
double sumrates, cv, alpha, lambda, invarfrac;
longer seed;
contribarr *contribution;
char aachar[26]="ARNDCQEGHILKMFPSTWYVBZX?*-";
char *progname;
long nonodes2;

/* Local variables for maketree, propagated globally for C version: */
long    k, maxwhich, col;
double  like, bestyet, maxlogl;
boolean smoothed, succeeded;
double  *l0gl;
double  expon1i[maxcategs], expon1v[maxcategs], expon2i[maxcategs], expon2v[maxcategs];
node   *there;
double  **l0gf;
Char ch, ch2;
long **mp;


tree* promlk_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(promlk_tree));
  promlk_tree_init(t, nonodes, spp);
  return t;
}


void promlk_tree_init(tree * t, long nonodes, long spp)
{
  ml_tree_init(t, nonodes, spp);
  t->insert_  = mlk_tree_insert_;
  t->try_insert_ = ml_tree_try_insert_;
  t->re_move = mlk_tree_re_move;
  t->evaluate = promlk_tree_evaluate;
  t->globrearrange = rooted_globrearrange;
  t->locrearrange = rooted_locrearrange;
  ((ml_tree*)t)->makenewv = mlk_tree_makenewv;
  t->nuview = promlk_tree_nuview;
  t->save_lr_nodes = rooted_tree_save_lr_nodes;
  t->restore_lr_nodes = rooted_tree_restore_lr_nodes;
}


void getoptions(void)
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean didchangecat, didchangercat;
  double probsum;

  putchar('\n');

  categs = 1;
  njumble = 1;
  rcategs = 1;

  lambda = 1.0;

  loopcount = 0;

  auto_ = false;
  ctgry = false;
  didchangecat = false;
  didchangercat = false;
  gama = false;
  global = false;
  hypstate = false;
  interleaved = true;
  invar = false;
  jumble = false;
  lengthsopt = false;
  printdata = false;
  progress = true;
  rctgry = false;
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
    printf("\nAmino acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
      printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    printf("  P    JTT, PMB or PAM probability model?  %s\n",
           usejtt ? "Jones-Taylor-Thornton" :
           usepmb ? "Henikoff/Tillier PMB" : "Dayhoff PAM");
    if (usertree)
    {
      printf("  L           Use lengths from user tree?");
      if (lengthsopt)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
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
    if (((!usertree) && (strchr("UPCRJAFWGTMI012345", ch) != NULL))
        || (usertree && ((strchr("UPCRAFWLTMI012345", ch) != NULL))))
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
{
  /* initializes variables */
  fprintf(outfile, "\nAmino acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n", VERSION);

  inputnumbers(&spp, &sites, &nonodes, 1);

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
    inputcategs(0, sites, category, categs, "ProMLK");
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
            printf("\nERROR:  Bad amino acid: %c at position %ld of species %ld.\n", charstate, j, i);
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
  contribution = (contribarr *) Malloc(endsite * sizeof(contribarr));
}


void prot_makevalues(long categs, pointarray treenode, long endsite, long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips   */
  /* a version of makevalues2 found in seq.c */
  /* used by proml                           */
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

        case 'X':               /* unknown aa                       */
          for (b = 0; b <= 19; b++)
            ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
          break;

        case '?':               /* unknown aa                       */
          for (b = 0; b <= 19; b++)
            ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
          break;

        case '*':               /* stop codon symbol                */
          for (b = 0; b <= 19; b++)
            ((prot_node*)treenode[i])->x[k][l][b] = 1.0;
          break;

        case '-':               /* deletion event-absent data or aa */
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
#if 0                                   // RSGnote: Variable never used.
  long grcategs;
#endif

  /* reads the input data */
  if (!justwts || firstset)
    inputoptions();
  if (!justwts || firstset)
    input_protdata(sites);
  prom_makeweights();

#if 0                                   // RSGnote: Variable never used.
  grcategs = (categs > rcategs) ? categs : rcategs;
#endif

  inittrees(&curtree, &bestree, &priortree, &bestree2, nonodes, spp);
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
}


void prot_inittable(void)
{
  /* Define a lookup table. Precompute values and print them out in tables */
  /* Allocate memory for the pmatrices, dpmatices and ddpmatrices          */
  long i, j, k, l;
  double sumrates;

  /* Allocate memory for pmatrices, the array of pointers to pmatrices     */

  pmatrices = (double *****) Malloc (spp * sizeof(double ****));

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

  if (gama || invar)
  {
    fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
    fprintf(outfile,
    " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
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
    fprintf(outfile, "\n\n");
  }
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


void promlk_tree_nuview(tree* t, node *p)
{
  long b, i, j, k, l, m, num_sibs, sib_index;
  node *sib_ptr, *sib_back_ptr;
  psitelike prot_xx, x2;
  double lw, prod7;
  double **pmat;
  double maxx, correction;

  /* Figure out how many siblings the current node has  */
  /* and be sure that pmatrices is large enough         */
  num_sibs = count_sibs(p);
  for (i = 0; i < num_sibs; i++)
    if (pmatrices[i] == NULL)
      prom_alloc_pmatrix(i);

  /* Recursive calls, should be called for all children */
  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (!(sib_back_ptr == NULL))
      if (!sib_back_ptr->tip && !sib_back_ptr->initialized)
        t->nuview(t, sib_back_ptr);
  }

  /* Make pmatrices for all possible combinations of category, rcateg      */
  /* and sib                                                               */
  sib_ptr = p;                          /* return to p */
  for (sib_index=0; sib_index < num_sibs; sib_index++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (sib_back_ptr != NULL)
      lw =  fabs( get_tyme(p) - get_tyme(sib_back_ptr) );
    else
      lw = 0.0;

    for (j = 0; j < rcategs; j++)
      for (k = 0; k < categs; k++)
        make_pmatrix(pmatrices[sib_index][j][k], NULL, NULL, 0, lw,
                                        tbl[j][k], eigmat, probmat);
  }

  for (i = 0; i < endsite; i++)
  {
    correction = 0;
    maxx = 0;
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {

      /* initialize to 1 all values of prot_xx */
      for (m = 0; m <= 19; m++)
        prot_xx[m] = 1;

      sib_ptr = p;                      /* return to p */
      /* loop through all sibs and calculate likelihoods for all possible*/
      /* amino acid combinations                                         */
      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        sib_ptr      = sib_ptr->next;
        sib_back_ptr = sib_ptr->back;


        if (sib_back_ptr != NULL)
        {
          memcpy(x2, ((prot_node*)sib_back_ptr)->x[i][j], sizeof(psitelike));
          if ( j == 0 )
            correction += ((ml_node*)sib_back_ptr)->underflows[i];
        }
        else
          for (b = 0; b <= 19; b++)
            x2[b] = 1.0;
        pmat = pmatrices[sib_index][j][k];
        for (m = 0; m <= 19; m++)
        {
          prod7 = 0;
          for (l = 0; l <= 19; l++)
            prod7 += (pmat[m][l] * x2[l]);
          prot_xx[m] *= prod7;
          if ( prot_xx[m] > maxx && sib_index == (num_sibs - 1 ))
            maxx = prot_xx[m];
        }
      }
      /* And the final point of this whole function: */
      memcpy(((prot_node*)p)->x[i][j], prot_xx, sizeof(psitelike));
    }
    ((ml_node*)p)->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE )
      fix_protx((prot_node*)p, i, maxx, rcategs);
    ((ml_node*)p)->underflows[i] += correction;
  }

  p->initialized = true;
}  /* promlk_tree_nuview */


double promlk_tree_evaluate(tree* t, node *p, boolean dummy)
{
  contribarr tterm;
  static contribarr like, nulike, clai;
  double sum, sum2, sumc=0, y, prod4, prodl, frexm, sumterm, lterm;
  double **pmat;
  long i, j, k, l, m, lai;
  node *q, *r;
  psitelike x1, x2;

  (void)dummy;                          // RSGnote: Parameter not used.
  sum = 0.0;

  if (p == curtree->root && (count_sibs(p) == 2))
  {
    r = p->next->back;
    q = p->next->next->back;
    y = get_tyme(r) + get_tyme(q) - 2 * get_tyme(p);
    if (!r->tip && !r->initialized) t->nuview(t, r);
    if (!q->tip && !q->initialized) t->nuview(t, q);
  }
  else if (p == curtree->root)
  {
    /* the next two lines copy tyme and x to p->next.  Normally they are
       not initialized for an internal node. */
    /* assumes bifurcation */
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
    if (!r->tip && !r->initialized) t->nuview(t, r);
    if (!q->tip && !q->initialized) t->nuview(t, q);
    y = fabs( get_tyme(r) - get_tyme(q) );
  }

  for (j = 0; j < rcategs; j++)
    for (k = 0; k < categs; k++)
      make_pmatrix(pmatrices[0][j][k], NULL, NULL, 0, y, tbl[j][k], eigmat, probmat);
  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      memcpy(x1, ((prot_node*)r)->x[i][j], sizeof(psitelike));
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
        sumterm = 0.00000001;   /* ??? */
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, rcategs * sizeof(double));
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
}  /* promlk_tree_evaluate */


void restoradd(node *below, node *newtip, node *newfork, double prevtyme)
{
/* Restore "new" tip and fork to place "below".  Restore tymes. */
/* Assumes bifurcation. */
  hookup(newfork, below->back);
  hookup(newfork->next, below);
  hookup(newtip, newfork->next->next);
  curtree->nodep[newfork->index-1] = newfork;
  set_tyme(newfork, prevtyme);
/* assumes bifurcations */
  set_tyme(newfork->next, prevtyme);
  set_tyme(newfork->next->next, prevtyme);
} /* restoradd */


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
    if (get_tyme(sib_back_ptr) < lowertyme)
      lowertyme = get_tyme(sib_back_ptr);
  }

  set_tyme(p, lowertyme - 0.1);

  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++)
  {
    sib_ptr = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    set_tyme(sib_ptr, get_tyme(p) );
    sib_back_ptr->v = get_tyme(sib_back_ptr) - get_tyme(p);
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


void travinit(node *p)
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
    travinit(sib_back_ptr);
  }

  curtree->nuview(curtree, p);
  p->initialized = true;
}  /* travinit */


void travsp(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to find tips */
  if (p == curtree->root)
    travinit(p);
  if (p->tip)
    travinit(p->back);
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

  if (!lngths)
    initrav(curtree->root);
  travsp(curtree->root);
  if ( usertree && !lngths )
    for ( i = 0 ; i < smoothings ; i++ ) /* this converges slowly, there     */
      curtree->smoothall(curtree, curtree->root); /* already is a loop within smoothall, we need more! */
  else
    curtree->smoothall(curtree, curtree->root); /* we should already be close */

  curtree->evaluate(curtree, curtree->root, 0);
}  /* treevaluate */


void promlk_coordinates(node *p, long *tipy)
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
    promlk_coordinates(q->back, tipy);
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
  p->xcoord = (long)(0.5 - over * get_tyme(p) );
  p->ycoord = (pp1->ycoord + pp2->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* promlk_coordinates */


void promlk_drawline(long i, double scale)
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
}  /* promlk_drawline */


void promlk_printree(void)
{
 /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;
  node *p;

  putc('\n', outfile);
  tipy = 1;
  promlk_coordinates(curtree->root, &tipy);
  p = curtree->root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)( get_tyme(p) - get_tyme(curtree->root) + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= tipy - down; i++)
    promlk_drawline(i, scale);
  putc('\n', outfile);
}  /* promlk_printree */


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
    fprintf(outfile, "%11.5f", ( get_tyme(p) - get_tyme(curtree->root) ));
    v = ( get_tyme(p) - get_tyme(curtree->nodep[p->back->index -1]) );
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
    prot_reconstr(p, i);
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
  long i, j;
  long mm = 0;                          // RSGnote: Variable formerly may have been referenced before being initialized.
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
      {
        // RSGnote: Variable formerly may have been referenced before being initialized.
        fprintf(outfile, "%ld", mm+1);
      }
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
  for (i = 0; i <= sites-1; ++i)
    free(mp[i]);
  free(mp);
}  /* summarize */


void initpromlnode(tree *treep, node **p, long len, long nodei, long *ntips, long *parens, initops whichinit, pointarray treenode, Char *str, Char *ch, FILE *intree)
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
        (*p)->back->iter = true;
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
      if (usertree && lengthsopt && lngths)
      {
        printf("Warning: one or more lengths not defined in user tree number %ld.\n", which);
        printf("PROMLK will attempt to optimize all branch lengths.\n\n");
        lngths = false;
      }
      break;
    case unittrwt:
      treep->nodep[spp]->iter = false;
      break;
    default:      /* ignore cases hslength, treewt */
      break;      /* should never occur                 */
  }
} /* initpromlnode */


void promlk_treeout(node *p)
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
      promlk_treeout(sib_ptr->back);
      putc(',', outtree);
      col++;
      if (col > 55)
      {
        putc('\n', outtree);
        col = 0;
      }
    }
    sib_ptr = sib_ptr->next;
    promlk_treeout(sib_ptr->back);
    putc(')', outtree);
    col++;
  }
  if (p == curtree->root)
  {
    fprintf(outtree, ";\n");
    return;
  }
  x = get_tyme(p) - get_tyme(curtree->nodep[p->back->index - 1]);
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
}  /* promlk_treeout */


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
  (*x) = get_tyme(p) - p->v;
}  /* tymetrav */


void maketree(void)
{
  /* constructs a binary tree from the pointers in curtree->nodep,
     adds each node at location which yields highest likelihood
     then rearranges the tree for greatest likelihood */

  long i, j;
  long numtrees = 0;
  double x;
  node *item, *q, *root = NULL;
  boolean dummy_haslengths, dummy_first, goteof, multf;
  long nextnode;
#if 0                                   // RSGnote: Variable never used.
  long grcategs;
#endif
  pointarray dummy_treenode = NULL;

#if 0                                   // RSGnote: Variable never used.
  grcategs = (categs > rcategs) ? categs : rcategs;
#endif

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
    fprintf(outfile, "\n\n");
    for( which = 1 ; which <= numtrees ; which++)
    {

      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      dummy_haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      lngths           = lengthsopt;

      treeread(curtree, intree, &root, dummy_treenode, &goteof, &dummy_first, &nextnode, &dummy_haslengths, initpromlnode, false, nonodes);
      nonodes = nextnode;
      root = curtree->nodep[root->index - 1];
      curtree->root = root;

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
        promlk_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        promlk_treeout(curtree->root);
      }
      if(which < numtrees)
      {
        prot_freex_notip(nonodes, curtree->nodep);
      }
    }

    FClose(intree);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite, maxlogl, l0gl, l0gf, aliasweight, seed);
  }
  else
  {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree->root = curtree->nodep[spp];  /* debug: use enterorder?? */
    curtree->root->back = NULL;
    for (i = 0; i < spp; i++)
      curtree->nodep[i]->back = NULL;
    for (i = spp; i < nonodes; i++)
    {
      q = curtree->nodep[i];
      q->back = NULL;
      while ((q = q->next) != curtree->nodep[i])
        q->back = NULL;
    }
    polishing = false;
    curtree->insert_(curtree, curtree->nodep[enterorder[0]-1], curtree->nodep[enterorder[1]-1], false, false);

    if (progress)
    {
      sprintf(progbuf, "\nAdding species:\n");
      print_progress(progbuf);
      writename(0, 2, enterorder);
      phyFillScreenColor();
    }
    smoothit = false;
    for (i = 3; i <= spp; i++)
    {
      bestyet = UNDEFINED;
      there = curtree->root;
      item = curtree->nodep[enterorder[i - 1] - 1];
#if 0                                   // RSGnote: Variable never used.
      nufork = curtree->nodep[spp + i - 2];
#endif
      curtree->copy(curtree, priortree);
      like = curtree->evaluate(curtree, curtree->root, 0);
      curtree->addtraverse(curtree, item, curtree->root, further, &qwhere,
                            &bestyet, bestree, priortree, true, &multf);
      curtree->insert_(curtree, item, qwhere, false, multf);
      curtree->smoothall(curtree, curtree->root);
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
      }
      smoothit = true;
      curtree->globrearrange(curtree, progress, true);
      smoothit = false;
      curtree->copy(curtree, bestree);
    }
    curtree->copy(curtree, bestree);

    if (njumble > 1)                                    // RSGbugfix
    {
      // Free up FORKRINGS for use in next tree -- the order matters here as we want
      // the lowest-index forknode to be the first to be popped from the stack.
      for ( i = nonodes - 1 ; i >= spp ; i-- )
        curtree->release_fork(curtree, curtree->nodep[i]);

      // Since root FORKNODE was released, nullify ROOT (will be reset in upcoming COPY).
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
      curtree->evaluate(curtree, curtree->root, 0);

      if (treeprint)
      {
        promlk_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        promlk_treeout(curtree->root);
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
  prom_free_all_x(nonodes2, curtree->nodep);
  if (!usertree)
  {
    prom_free_all_x(nonodes2, bestree->nodep);
    if (njumble > 1)
      prom_free_all_x(nonodes2, bestree2->nodep);
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

  free(root);
} /* maketree */


void promlkrun(void)
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
    printf("invar: %i\n", invar);
    printf("jumble: %i\n", jumble);
    printf("njumble: %li\n", njumble);
    printf("lngths: %i\n", lngths);
    printf("lambda: %f\n", lambda);
    printf("trout: %i\n", trout);
    printf("usertree: %i\n", usertree);
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
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      max_num_sibs = 0;
      maketree();
    }
    fflush(outfile);
    fflush(outtree);
  }
}


void promlk(
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
  int RecHypo)
{
  initdata funcs;

  //printf("Hello from ProMLK!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Promlk";

  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = prot_node_new;
  funcs.tree_new = promlk_tree_new;

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
    usertree   = true;
  }
  else
  {
    // Yes
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
  //printf("calling doinit\n");  //JRMdebug
  //fflush(stdout); //JRMdebug
  doinit();

  //printf("calling promlrun\n");  //JRMdebug
  //fflush(stdout); //JRMdebug
  promlkrun();  // do the actual work

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

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* Protein Maximum Likelihood with molecular clock */
  initdata funcs;
#ifdef MAC
  argc = 1;             /* macsetup("Promlk", "");        */
  argv[0] = "Promlk";
#endif
  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = prot_node_new;
  funcs.tree_new = promlk_tree_new;
  phylipinit(argc, argv, &funcs, false);
  progname = argv[0];
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

  promlkrun();

  prom_clean_up();
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* Protein Maximum Likelihood with molecular clock */


// End.
