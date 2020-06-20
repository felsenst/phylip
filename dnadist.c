/* Phylip Version 4.0.
 * dnadist.c
 *
 * (c) Copyright 1993-2013 by the University of Washington.  Written by Joseph
 * Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.  Permission is
 * granted to copy and use this program provided no fee is charged for it and
 * provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "matrixd.h"
#include "pmatrix.h"
#include "seq.h"


#define ITERATIONSD     100   /* number of iterates of EM for each distance */

typedef struct valrec {
  double rat,
    ratxv,
    z1,
    y1,
    z1zz,
    z1yy,
    z1xv;
} valrec;

typedef enum dnadist_method {
  f84,
  kimura,
  jukes,
  logdet,
  similarity
} dnadist_method;

/* Default options */
dnadist_method  method  = f84;          /* method */
boolean ctgry           = false;        /* use more than one category */
long    categs          = 1;            /* number of categories */
double  rate[maxcategs] = { 1.0 };      /* default category rates */
boolean gama            = false;        /* use gamma distributed rates */
double  cvi             = 1.0;          /* coefficient of variation */
boolean invar           = false;        /* use gamma+invar rates */
double  invarfrac       = 0.0;          /* fraction invariant */
double  ttratio         = 2.0;          /* transition/transversion ratio */
boolean weights         = false;        /* use weights */
boolean freqsfrom       = true;         /* use freqs from data */
boolean mulsets         = false;        /* multiple data sets */
long    datasets        = 1;            /* number of sets */
boolean justwts         = false;        /* ...just weights */
boolean interleaved     = true;         /* data set interleaved */
boolean printdata       = false;        /* print data */
boolean dotdiff         = true;         /* use dot-difference format */
boolean progress        = true;         /* enable progress indicator */
boolean matrix_flags    = MAT_MACHINE;  /* Matrix output format */

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], catfilename[FNMLNGTH], weightfilename[FNMLNGTH];
long sites, weightsum, ith;

boolean firstset;
tree *curtree, *bestree, *priortree, *bestree2;       /* to make seq.h happy */
sitelike **x;    /* array of sitelike [0..spp-1][0..endsite-1] */
double xv, freqa, freqc, freqg, freqt, freqar, freqcy, freqgr, freqty, fracchange;
basefreq freq;
steptr oldweight;
double sumweightrat;                  /* these values were propagated  */
double *weightrat = NULL;             /* to global values from         */
valrec tbl[maxcategs];                /* function makedists.           */


#ifndef OLDC
/* function  prototypes */
void   dnadist_config(void);
void   allocsites(long nsites);
void   freesites(void);
void   reallocsites(long nsites);
void   dnadist_init(void);
void   inputcategories(void);
void   printcategories(void);
void   inputoptions(void);
void   dnadist_sitesort(void);
void   dnadist_sitecombine(void);
void   dnadist_sitescrunch(void);
void   makeweights(void);
void   dnadist_makevalues(void);
void   dnadist_empiricalfreqs(void);
void   getinput(void);
void   inittable(void);
double lndet(double (*a)[4]);
double makev(long, long);
void   makedists(double **d);
void   writedists(double **d);
void   dnadist_write_outfile(void);
void   dnadistrun(void);
void   dnadist(char *infilename, char *outfilename, char *outfileopt, char *catfilename, char *weightfilename,
               char *distanceoptions, char *gammadistributed, double coefvar, double pctinvar, int usettratio,
               double inttratio, int usesubrates, int innumcats, double inrate1, double inrate2, double inrate3,
               double inrate4, double inrate5, double inrate6, double inrate7, double inrate8, double inrate9,
               int useweights, int useEmpBF, double basefreqA, double basefreqC, double basefreqG, double basefreqTU,
               char *distmatrix, int usemultdataset, int useData, int numdatasets, int doseqinter, int doprintdata,
               int dodotdiff, int doprintind);
int    main(int argc, char **argv);
/* function  prototypes */
#endif


void dnadist_config(void)
{
  /* display menu and interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;
  boolean ttr = false, ok = false;
  char *str;

  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nNucleic acid sequence Distance Matrix program,");
    printf("PHYLIP version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    switch ( method )
    {
      case kimura:
        str = "Kimura 2-parameter";
        break;
      case jukes:
        str = "Jukes-Cantor";
        break;
      case logdet:
        str = "LogDet";
        break;
      case similarity:
        str = "Similarity table";
        break;
      case f84:
        str = "F84";
        break;
      default:
        str = "(unknown)";
        assert(0);
    }

    printf("  D  Distance (F84, Kimura, Jukes-Cantor, LogDet)?  %s\n", str);
    if ((method == kimura) ||
        (method == f84)    ||
        (method == jukes))
    {
      printf("  G          Gamma distributed rates across sites?  ");
      if (gama)
      {
        printf("Yes\n");
      }
      else
      {
        if (invar)
        {
          printf("Gamma+Invariant\n");
        }
        else
        {
          printf("No\n");
        }
      }
    }

    if ((method == kimura) ||
        (method == f84)) {
      printf("  T                 Transition/transversion ratio? ");
      if (!ttr)
        printf("  2.0\n");
      else
        printf("%7.4f\n", ttratio);
    }

    if ((method == kimura) ||
        (method == f84 )   ||
        (method == jukes)  ||
        gama               ||
        invar)
    {
      //if ( !(method & (logdet|similarity) || gama || invar) ) {
      printf("  C            One category of substitution rates?");
      if (!ctgry || categs == 1)
        printf("  Yes\n");
      else
        printf("  %ld categories\n", categs);
    }

    printf("  W                         Use weights for sites?");
    if (weights)
      printf("  Yes\n");
    else
      printf("  No\n");
    if ( method == f84 )
      printf("  F                Use empirical base frequencies?  %s\n",
             (freqsfrom ? "Yes" : "No"));
    printf("  L                       Form of distance matrix?  ");
    switch (matrix_flags)
    {
      case MAT_MACHINE:
        str = "Square";
        break;
      case MAT_LOWERTRI:
        str = "Lower-triangular";
        break;
      case MAT_HUMAN:
        str = "Human-readable";
        break;
      default:  /* shouldn't happen */
        assert( 0 );
        str = "(unknown)";
    }
    printf("%s\n", str);
    printf("  M                    Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I                   Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0            Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  1             Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2           Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    ok = true;
    switch (ch)
    {
      case 'D':
        switch (method)
        {
          case f84:
            method = kimura;
            freqsfrom = false;
            break;
          case kimura:
            method = jukes;
            freqsfrom = false;
            break;
          case jukes:
            method = logdet;
            freqsfrom = false;
            break;
          case logdet:
            method = similarity;
            freqsfrom = false;
            break;
          case similarity:
            method = f84;
            freqsfrom = true;
            break;
          default: /* shouldn't happen */
            assert(0);
            method = f84;
        }
        break;

      case 'G':
        if ( ctgry ||
             (method == logdet) ||
             (method == similarity)) {
          ok = false;
          break;
        }
        if (!(gama || invar))
          gama = true;
        else {
          if (gama) {
            gama = false;
            invar = true;
          } else {
            if (invar)
              invar = false;
          }
        }
        break;

      case 'C':
        if ((method == logdet) ||
            (method == similarity)) {
          ok = false;
          break;
        }
        ctgry = !ctgry;
        if (ctgry) {
          initcatn(&categs);
          initcategs(categs, rate);
        }
        break;

      case 'F':
        if ( method != f84 ) {
          ok = false;
          break;
        }
        freqsfrom = !freqsfrom;
        if (!freqsfrom)
          initfreqs(&freqa, &freqc, &freqg, &freqt);
        break;

      case 'W':
        weights = !weights;
        break;

      case 'L':
        /* square(machine-readable) -> lower-triangular -> square(human-readable) */
        switch ( matrix_flags )
        {
          case MAT_HUMAN:
            matrix_flags = MAT_MACHINE;
            break;
          case MAT_LOWERTRI:
            matrix_flags = MAT_HUMAN;
            break;
          case MAT_MACHINE:
            matrix_flags = MAT_LOWERTRI;
            break;
          default:
            assert(0);
            matrix_flags = MAT_HUMAN;
        }
        break;

      case 'T':
        ttr = !ttr;
        if (ttr)
          initratio(&ttratio);
        break;

      case 'M':
        if (mulsets) {
          mulsets = false;
          datasets = 1;
        }
        else {
          mulsets = true;
          printf("Multiple data sets or multiple weights?");
          loopcount2 = 0;
          do {
            printf(" (type D or W)\n");
            phyFillScreenColor();
            if(scanf("%c%*[^\n]", &ch2)) {} // Read char and scan to EOL.
            uppercase(&ch2);
            (void)getchar();
            countup(&loopcount2, 10);
          } while ((ch2 != 'W') && (ch2 != 'D'));
          justwts = (ch2 == 'W');
          if (justwts)
            justweights(&datasets);
          else
            initdatasets(&datasets);
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

      default:
        printf("Not a possible option!\n");
        printf("\nPress Enter or Return key to continue\n");
        (void)getchar();
    }
    if ( !ok ) {
      printf("That option is not allowed with these settings.\n");
      printf("\nPress Enter or Return key to continue\n");
      (void)getchar();
    }
    countup(&loopcount, 100);
  }

  /* Prevent similarity matrices from easily being used as input to other programs. */
  if (method == similarity) {
    if (matrix_flags == MAT_MACHINE)
      matrix_flags = MAT_HUMAN;
    else if (matrix_flags == MAT_LOWERTRI)
      matrix_flags = MAT_LOWER | MAT_HUMAN;
  }

  if (gama || invar) {
    loopcount = 0;
    do {
      printf("\n"
             "Coefficient of variation of substitution rate among sites (must be positive):\n"
             "  In gamma distribution parameters, this is 1/(square root of alpha)\n" );
      phyFillScreenColor();
      if(scanf("%lf%*[^\n]", &cvi)) {}  // Read number and scan to EOL.
      (void)getchar();
      countup(&loopcount, 10);
    } while (cvi <= 0.0);
    cvi = 1.0 / (cvi * cvi);
  }
  if (invar) {
    loopcount = 0;
    do {
      printf("Fraction of invariant sites?\n");
      if(scanf("%lf%*[^\n]", &invarfrac)) {} // Read number and scan to EOL.
      (void)getchar();
      countup (&loopcount, 10);
    } while ((invarfrac <= 0.0) || (invarfrac >= 1.0));
  }
}  /* dnadist_config */


void allocsites(long nsites)
{
  long i;

  assert(sites > 0);

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++) {
    inputSequences[i] = (Char *)Malloc(sites * sizeof(Char));
  }

  category  = (steptr)Malloc(nsites * sizeof(*category));
  oldweight = (steptr)Malloc(nsites * sizeof(*oldweight));
  weight    = (steptr)Malloc(nsites * sizeof(*weight));
  alias     = (steptr)Malloc(nsites * sizeof(*alias));
  ally      = (steptr)Malloc(nsites * sizeof(*ally));
  location  = (steptr)Malloc(nsites * sizeof(*location));
  weightrat = (double *)Malloc(nsites * sizeof(*weightrat));
}


void freesites(void)
{
  long i;

  if ( inputSequences != NULL ) {
    for (i = 0; i < spp; i++)
      free(inputSequences[i]);
    free(inputSequences);
    inputSequences = NULL;
  }

  free(category);
  category = NULL;
  free(oldweight);
  oldweight = NULL;
  free(weight);
  weight = NULL;
  free(alias);
  alias = NULL;
  free(ally);
  ally = NULL;
  free(location);
  location = NULL;
  free(weightrat);
  weightrat = NULL;
}


void reallocsites(long nsites)
{
  /* Reallocate arrays sized by number of sites, which may change between
   * datasets. */

  freesites();
  allocsites(nsites);
} /* reallocsites */


void dnadist_init(void)
{
  long nonodes_dummy;

  inputnumbers(&spp, &sites, &nonodes_dummy, 1);
}  /* dnadist_init */


void inputcategories(void)
{
  /* reads the categories for each site */
  long i;
  Char ch;

  for (i = 1; i < nmlngth; i++)
    gettc(infile);
  for (i = 0; i < sites; i++) {
    do {
      if (eoln(infile))
        scan_eoln(infile);
      ch = gettc(infile);
    } while (ch == ' ');
    category[i] = ch - '0';
  }
  scan_eoln(infile);
}  /* inputcategories */


void printcategories(void)
{ /* print out list of categories of sites */
  long i, j;

  fprintf(outfile, "Rate categories\n\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', outfile);
  for (i = 1; i <= sites; i++) {
    fprintf(outfile, "%ld", category[i - 1]);
    if (i % 60 == 0) {
      putc('\n', outfile);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', outfile);
    } else if (i % 10 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* printcategories */


void inputoptions(void)
{
  /* read options information */
  long i;

  if ( !(firstset || justwts) ) {
    samenumsp(&sites, ith);
    reallocsites(sites);
  }
  for (i = 0; i < sites; i++) {
    category[i] = 1;
    oldweight[i] = 1;
  }
  if (justwts || weights)
    inputweights(sites, oldweight, &weights);
  if (ctgry && categs > 1) {
    inputcategs(0, sites, category, categs, "DnaDist");
  }
  if (freqsfrom)
  {
    if ((method == jukes) ||
        (method == kimura) ||
        (method == logdet) )
    {
      printf(" WARNING: Cannot use empirical base frequencies\n"
             " with Jukes-Cantor, Kimura or logdet distances\n");
      exxit(-1);
    }
  }
  if (method == jukes)
    ttratio = 0.5000001;
}  /* inputoptions */


void dnadist_sitesort(void)
{
  /* Shell sort of sites lexicographically */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip;

  boolean tied = false;                 // RSGnote: Variable being initialized here only to silence compiler warning.

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip)
      {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        // RSGnote: BUG - "tied" referenced here before being initialized.
        flip = (oldweight[jj - 1] < oldweight[jg - 1] || (tied && category[jj - 1] > category[jg - 1]));
        tied = (category[jj - 1] == category[jg - 1]);
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
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* dnadist_sitesort */


void dnadist_sitecombine(void)
{
  /* combine sites that have identical patterns */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      tied = (category[alias[i - 1] - 1] == category[alias[j - 1] - 1]);
      k = 1;
      while (k <= spp && tied) {
        tied = (tied &&
                inputSequences[k - 1][alias[i - 1] - 1] == inputSequences[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (!tied)
        break;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* dnadist_sitecombine */


void dnadist_sitescrunch(void)
{
  /* move so one representative of each pattern of
     sites comes first */
  long i, j, itemp;
  boolean done, found, completed;

  done = false;
  i = 1;
  j = 2;
  while (!done) {
    if (ally[alias[i - 1] - 1] != alias[i - 1]) {
      if (j <= i)
        j = i + 1;
      if (j <= sites) {
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
          completed = (j > sites);
        } while (!(found || completed));
        if (found) {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
        } else
          done = true;
      } else
        done = true;
    }
    i++;
    done = (done || i >= sites);
  }
}  /* dnadist_sitescrunch */


void makeweights(void)
{
  /* make up weights vector to avoid duplicate computations */
  long i;
  double sumrates;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = i;
    weight[i - 1] = 0;
    location[i - 1] = 0;
  }
  dnadist_sitesort();
  dnadist_sitecombine();
  dnadist_sitescrunch();
  endsite = 0;
  for (i = 1; i <= sites; i++) {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += oldweight[i];
  sumrates = 0.0;
  for (i = 0; i < sites; i++)
    sumrates += oldweight[i] * rate[category[i] - 1];
  for (i = 0; i < categs; i++)
    rate[i] *= weightsum / sumrates;
  for (i = 0; i < sites; i++)
    if (oldweight[i] > 0)
      weight[location[ally[i] - 1] - 1] += oldweight[i];
}  /* makeweights */


void dnadist_makevalues(void)
{
  /* initialize fractional likelihoods */

  long sp, site;
  long sitealias;
  boolean ok = false;
  char basechar;

  x = (sitelike **)Malloc(spp * sizeof(sitelike *));
  for (sp = 0; sp < spp; sp++) {
    x[sp] = (sitelike *)Malloc(endsite * sizeof(sitelike));
    for (site = 0; site < endsite; site++) {
      x[sp][site][A] = 0.0;
      x[sp][site][C] = 0.0;
      x[sp][site][G] = 0.0;
      x[sp][site][T] = 0.0;

      sitealias = alias[site] - 1;
      basechar = inputSequences[sp][sitealias];
      ok = false;
      if ( strchr("AMRWDHVNX?O-", basechar) ) { /* A bases */
        x[sp][site][A] = 1.0;
        ok = true;
      }
      if ( strchr("CMYSBHVNX?O-", basechar) ) { /* C bases */
        x[sp][site][C] = 1.0;
        ok = true;
      }
      if ( strchr("GKRSBDVNX?O-", basechar) ) { /* G bases */
        x[sp][site][G] = 1.0;
        ok = true;
      }
      if ( strchr("TUKYWBDHNX?O-", basechar) ) { /* T bases */
        x[sp][site][T] = 1.0;
        ok = true;
      }
      if ( !ok ) {
        /* Shouldn't happen (bases are checked by inputdata()) */
        assert( 0 );
      }
    }
  }
}  /* dnadist_makevalues */


void dnadist_empiricalfreqs(void)
{
  /* Get empirical base frequencies from the data */
  long i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  freqa = 0.25;
  freqc = 0.25;
  freqg = 0.25;
  freqt = 0.25;
  for (k = 1; k <= 8; k++) {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < spp; i++) {
      for (j = 0; j < endsite; j++) {
        sum   =     0.0;
        sum  +=     freqa * x[i][j][A];
        sum  +=     freqc * x[i][j][C];
        sum  +=     freqg * x[i][j][G];
        sum  +=     freqt * x[i][j][T];

        w = weight[j];
        suma += w * freqa * x[i][j][A] / sum;
        sumc += w * freqc * x[i][j][C] / sum;
        sumg += w * freqg * x[i][j][G] / sum;
        sumt += w * freqt * x[i][j][T] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    freqa = suma / sum;
    freqc = sumc / sum;
    freqg = sumg / sum;
    freqt = sumt / sum;
  }
}  /* dnadist_empiricalfreqs */


void getinput(void)
{
  /* reads the input data */

  inputoptions();
  if (!justwts || firstset)
    read_sequences(sites);
  makeweights();
  dnadist_makevalues();
  if (method & (logdet|similarity))
    return;             /* freqs not used */
  if (freqsfrom) {
    /* f84 only */
    dnadist_empiricalfreqs();
    makebasefreq(&freq, freqa, freqc, freqg, freqt, ttratio);
  }
  else {
    /* Always use equal frequencies for Kimura and Jukes */
    if ( method & (kimura|jukes) ) {
      freqa = 0.25;
      freqc = 0.25;
      freqg = 0.25;
      freqt = 0.25;
    }
    makebasefreq(&freq, freqa, freqc, freqg, freqt, ttratio);

    if (freqa < 0.00000001) {
      freqa = 0.000001;
      freqc = 0.999999*freqc;
      freqg = 0.999999*freqg;
      freqt = 0.999999*freqt;
    }
    if (freqc < 0.00000001) {
      freqa = 0.999999*freqa;
      freqc = 0.000001;
      freqg = 0.999999*freqg;
      freqt = 0.999999*freqt;
    }
    if (freqg < 0.00000001) {
      freqa = 0.999999*freqa;
      freqc = 0.999999*freqc;
      freqg = 0.000001;
      freqt = 0.999999*freqt;
    }
    if (freqt < 0.00000001) {
      freqa = 0.999999*freqa;
      freqc = 0.999999*freqc;
      freqg = 0.999999*freqg;
      freqt = 0.000001;
    }
  }

  /* Unpack makebasefreq result */
  /* TODO use freq directly (make global?) */
  freqar     = freq.ar;
  freqcy     = freq.cy;
  freqgr     = freq.gr;
  freqty     = freq.ty;
  xv         = freq.xv;
  fracchange = freq.fracchange;
  ttratio    = freq.ttratio;
}  /* getinput */


void inittable(void)
{
  /* Define a lookup table. Precompute values and store in a table */
  long i;

  for (i = 0; i < categs; i++) {
    tbl[i].rat = rate[i];
    tbl[i].ratxv = rate[i] * xv;
  }
}  /* inittable */


double lndet(double (*a)[4])
{
  long i, j, k;
  double temp, ld;

  /*Gauss-Jordan reduction -- invert matrix a in place,
    overwriting previous contents of a.  On exit, matrix a
    contains the inverse, lndet contains the log of the determinant */
  ld = 1.0;
  for (i = 0; i < 4; i++) {
    ld *= a[i][i];
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < 4; j++)
      a[i][j] *= temp;
    for (j = 0; j < 4; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < 4; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
  if (ld <= 0.0)
    return(99.0);
  else
    return(log(ld));
}  /* lndet */


double makev(long m, long n)
{
  /* compute one distance */
  long i, j, k, l, it, num1, num2, idx;
  long numerator = 0, denominator = 0;
  double sum, sum1, sum2, sumyr, lz, aa, bb, cc, vv=0,
    p1, p2, p3, q1, q2, q3, tt, delta = 0, slope,
    xx1freqa, xx1freqc, xx1freqg, xx1freqt;
  double *prod = NULL;
  double *prod2 = NULL;
  double *prod3 = NULL;
  boolean quick = false;
  boolean jukesquick = false;
  boolean kimquick = false;
  boolean logdetquick = false;
  boolean similarityquick = false;
  boolean overlap;
  bases b;
  sitelike *xp = NULL, *xq = NULL;
  double *xx1 = NULL, *xx2 = NULL;
  double basetable[4][4];  /* for quick logdet */
  double basefreq1[4], basefreq2[4];

  xp = x[m-1];
  xq = x[n-1];

  /* check for overlap between sequences */
  overlap = false;
  for(i=0 ; i < sites ; i++)
  {
    if((strchr("NX?O-", inputSequences[m-1][i]) == NULL) &&
       (strchr("NX?O-", inputSequences[n-1][i]) == NULL))
    {
      overlap = true;
      break;
    }
  }
  if( !overlap )
  {
    printf("\nWARNING: No overlap between sequences "
           "%ld and %ld; -1.0 was written.\n", m, n );
    return -1.0;
  }

  /* Use quick method when only one category */
  if (!ctgry || categs == 1)
    quick = true;
  else
    quick = false;

  if ( method & (jukes|kimura|logdet|similarity) ) {
    numerator = 0;
    denominator = 0;
    for (i = 0; i < endsite; i++) {
      xx1 = (double *)(xp + i);
      xx2 = (double *)(xq + i);
      sum = 0.0;
      sum1 = 0.0;
      sum2 = 0.0;
      for (b = A; b <= T; b++) {
        sum1 += xx1[b];
        sum2 += xx2[b];
        sum += xx1[b] * xx2[b];
      }
      if (quick) {
        if ( sum1 != 1.0 && sum1 != 4.0 )
          quick = false;
        else if ( sum2 != 1.0 && sum2 != 4.0)
          quick = false;
      }
      if (sum1 == 1.0 && sum2 == 1.0) {
        numerator += (long)(weight[i] * sum);
        denominator += weight[i];
      }
    }
  }
  if (quick)
  {
    switch (method)
    {
      case f84:
        /* no quick f84 method */
        break;
      case kimura:
        kimquick = true;
        break;
      case jukes:
        jukesquick = true;
        break;
      case logdet:
        logdetquick = true;
        (void)logdetquick;              // RSGnote: Variable set but never referenced.
        break;
      case similarity:
        similarityquick = true;
        break;
      default:
        assert(0);      /* Shouldn't happen */
    }
  }
  if (method == similarity) {
    if (denominator < 1.0) {
      printf("\nWARNING: SPECIES %3ld AND %3ld HAVE NO BASES THAT", m, n);
      printf(" CAN BE COMPARED\n");
      printf("  -1.0 WAS WRITTEN\n");
      return -1.0;
    }
  }
  if (kimquick) {
    num1 = 0;
    num2 = 0;
    denominator = 0;
    for (i = 0; i < endsite; i++) {
      xx1 = (double *)(xp + i);
      xx2 = (double *)(xq + i);
      sum = 0.0;
      sum1 = 0.0;
      sum2 = 0.0;
      for (b = A; b <= T; b++) {
        sum1 += xx1[b];
        sum2 += xx2[b];
        sum += xx1[b] * xx2[b];
      }
      sumyr = (xx1[A] + xx1[G])
            * (xx2[A] + xx2[G]) +
              (xx1[C] + xx1[T]) *
              (xx2[C] + xx2[T]);
      if (sum1 == 1.0 && sum2 == 1.0) {
        num1 += (long)(weight[i] * sum);
        num2 += (long)(weight[i] * (sumyr - sum));
        denominator += weight[i];
      }
    }
    tt = ((1.0 - (double)num1 / denominator)-invarfrac)/(1.0-invarfrac);
    if (tt > 0.0) {
      delta = 0.1;
      tt = delta;
      it = 0;
      while (fabs(delta) > 0.0000002 && it < ITERATIONSD) {
        it++;
        if (!gama) {
          p1 = exp(-tt);
          p2 = exp(-xv * tt) - exp(-tt);
          p3 = 1.0 - exp(-xv * tt);
        } else {
          p1 = exp(-cvi * log(1 + tt / cvi));
          p2 = exp(-cvi * log(1 + xv * tt / cvi))
              - exp(-cvi * log(1 + tt / cvi));
          p3 = 1.0 - exp(-cvi * log(1 + xv * tt / cvi));
        }
        q1 = p1 + p2 / 2.0 + p3 / 4.0;
        q2 = p2 / 2.0 + p3 / 4.0;
        q3 = p3 / 2.0;
        q1 = q1 * (1.0-invarfrac) + invarfrac;
        q2 *= (1.0 - invarfrac);
        q3 *= (1.0 - invarfrac);
        if (!gama && !invar)
          slope = 0.5 * exp(-tt) * (num2 / q2 - num1 / q1) +
                  0.25 * xv * exp(-xv * tt) *
                 ((denominator - num1 - num2) * 2 / q3 - num2 / q2 - num1 / q1);
        else
          slope = 0.5 * (1 / (1 + tt / cvi)) * exp(-cvi * log(1 + tt / cvi)) *
                  (num2 / q2 - num1 / q1) + 0.25 * (xv / (1 + xv * tt / cvi)) *
                    exp(-cvi * log(1 + xv * tt / cvi)) *
                 ((denominator - num1 - num2) * 2 / q3 - num2 / q2 - num1 / q1);
        slope *= (1.0-invarfrac);
        if (slope < 0.0)
          delta = fabs(delta) / -2.0;
        else
          delta = fabs(delta);
        tt += delta;
      }
    }
    if ((delta >= 0.1) && (method != similarity)) {
      printf("\nWARNING: DIFFERENCE BETWEEN SPECIES %3ld AND %3ld", m, n);
      if (invar)
        printf(" TOO LARGE FOR INVARIABLE SITES\n");
      else
        printf(" TOO LARGE TO ESTIMATE DISTANCE\n");
      printf("  -1.0 WAS WRITTEN\n");
      return -1.0;
    }
    vv = fracchange * tt;
  }
  if (jukesquick || similarityquick) {
    if (jukesquick && (numerator * 4 <= denominator)) {
      printf("\nWARNING: INFINITE DISTANCE BETWEEN ");
      printf(" SPECIES %3ld AND %3ld\n", m, n);
      printf("  -1.0 WAS WRITTEN\n");
      return -1.0;
    }
    else if ( (jukesquick || similarityquick) && invar
              && (4 * (((double)numerator / denominator) - invarfrac)
                  <= (1.0 - invarfrac))) {
      printf("\nWARNING: DIFFERENCE BETWEEN SPECIES %3ld AND %3ld", m, n);
      printf(" TOO LARGE FOR INVARIABLE SITES\n");
      printf("  -1.0 WAS WRITTEN\n");
      return -1.0;
    }
    else {
      if (!gama && !invar)
        vv = -0.75 * log((4.0*((double)numerator / denominator) - 1.0) / 3.0);
      else if (!invar)
        vv = 0.75 * cvi * (exp(-(1/cvi)*
                               log((4.0 * ((double)numerator / denominator) - 1.0) / 3.0)) - 1.0);
      else
        vv = 0.75 * cvi * (exp(-(1/cvi)*
                               log((4.0 * ((double)numerator / denominator - invarfrac)/
                                    (1.0-invarfrac) - 1.0) / 3.0)) - 1.0);
    }
  }
  if (method == logdet) {
    if (!quick) {
      printf(" WARNING: Cannot calculate logdet distance\n");
      printf("  with partially ambiguous nucleotides.\n");
      printf("  -1.0 was written.\n");
      return -1.0;
    }
    else {
      /* compute logdet when no ambiguous nucleotides */
      for (i = 0; i < 4; i++) {
        basefreq1[i] = 0.0;
        basefreq2[i] = 0.0;
        for (j = 0; j < 4; j++)
          basetable[i][j] = 0.0;
      }
      for (i = 0; i < endsite; i++) {
        k = 0;
        while ( xp[i][k] == 0.0 )
          k++;
        basefreq1[k] += weight[i];
        l = 0;
        while ( xq[i][l] == 0.0 )
          l++;
        basefreq2[l] += weight[i];
        basetable[k][l] += weight[i];
      }
      vv = lndet(basetable);
      if (vv == 99.0) {
        printf("\nNegative or zero determinant for distance between species");
        printf(" %ld and %ld.\n", m, n);
        printf("  -1.0 WAS WRITTEN\n");
        return -1.0;
      }
      vv = -0.25*(vv - 0.5*(log(basefreq1[0])+log(basefreq1[1])
                            +log(basefreq1[2])+log(basefreq1[3])
                            +log(basefreq2[0])+log(basefreq2[1])
                            +log(basefreq2[2])+log(basefreq2[3])));
    }
  }
  else {
    /* For F84 and multi-category Jukes, similarity, and Kimura */
    prod = (double *)Malloc(sites * sizeof(double));
    prod2 = (double *)Malloc(sites * sizeof(double));
    prod3 = (double *)Malloc(sites * sizeof(double));
    for (i = 0; i < endsite; i++) {
      xx1 = (double *)(xp + i);
      xx2 = (double *)(xq + i);
      xx1freqa = xx1[A] * freqa;
      xx1freqc = xx1[C] * freqc;
      xx1freqg = xx1[G] * freqg;
      xx1freqt = xx1[T] * freqt;

      sum1 = xx1freqa + xx1freqc + xx1freqg + xx1freqt;
      sum2 = freqa * xx2[A] + freqc * xx2[C] +
             freqg * xx2[G] + freqt * xx2[T];
      prod[i] = sum1 * sum2;

      prod2[i] =
        (xx1freqa + xx1freqg) *
        (xx2[A] * freqar + xx2[G] * freqgr)
        +
        (xx1freqc + xx1freqt) *
        (xx2[C] * freqcy + xx2[T] * freqty);

      prod3[i] =
        xx1freqa * xx2[A] + xx1freqc * xx2[C] +
        xx1freqg * xx2[G] + xx1freqt * xx2[T];
    }
    tt = 0.1;
    delta = 0.1;
    it = 1;
    while (it < ITERATIONSD && fabs(delta) > 0.0000002) {
      slope = 0.0;
      if (tt > 0.0) {
        lz = -tt;
        for (i = 0; i < categs; i++) {
          if (!gama) {
            tbl[i].z1 = exp(tbl[i].ratxv * lz);
            tbl[i].z1zz = exp(tbl[i].rat * lz);
          }
          else {
            tbl[i].z1 = exp(-cvi*log(1.0-tbl[i].ratxv * lz/cvi));
            tbl[i].z1zz = exp(-cvi*log(1.0-tbl[i].rat * lz/cvi));
          }
          tbl[i].y1 = 1.0 - tbl[i].z1;
          tbl[i].z1yy = tbl[i].z1 - tbl[i].z1zz;
          tbl[i].z1xv = tbl[i].z1 * xv;
        }
        for (i = 0; i < endsite; i++) {
          idx = category[alias[i] - 1];
          cc = prod[i];
          bb = prod2[i];
          aa = prod3[i];
          if (!gama && !invar)
            slope += weightrat[i] * (tbl[idx - 1].z1zz * (bb - aa) +
                                     tbl[idx - 1].z1xv * (cc - bb)) /
                         (aa * tbl[idx - 1].z1zz + bb * tbl[idx - 1].z1yy +
                          cc * tbl[idx - 1].y1);
          else
            slope += (1.0-invarfrac) * weightrat[i] * (
                    ((tbl[idx-1].rat)/(1.0-tbl[idx-1].rat * lz/cvi))
                       * tbl[idx - 1].z1zz * (bb - aa) +
                    ((tbl[idx-1].ratxv)/(1.0-tbl[idx-1].ratxv * lz/cvi))
                       * tbl[idx - 1].z1 * (cc - bb)) /
                (aa * ((1.0-invarfrac)*tbl[idx - 1].z1zz + invarfrac)
                  + bb * (1.0-invarfrac)*tbl[idx - 1].z1yy
                  + cc * (1.0-invarfrac)*tbl[idx - 1].y1);
        }
      }
      if (slope < 0.0)
        delta = fabs(delta) / -2.0;
      else
        delta = fabs(delta);
      tt += delta;
      it++;
    }
    if ((delta >= 0.1) && (method != similarity)) {
      printf("\nWARNING: DIFFERENCE BETWEEN SPECIES %3ld AND %3ld", m, n);
      if (invar)
        printf(" TOO LARGE FOR INVARIABLE SITES\n");
      else
        printf(" TOO LARGE TO ESTIMATE DISTANCE\n");
      printf("  -1.0 WAS WRITTEN\n");
      return -1.0;
    }
    free(prod);
    free(prod2);
    free(prod3);
    if ( method == similarity ) {
      vv = (double)numerator / denominator;
    }
    else {
      vv = tt * fracchange;
    }
  }
  return fabs(vv);
}  /* makev */


void makedists(double **d)
{
  /* compute distance matrix */
  long i, j;
  double v;

  assert( d != NULL );
  inittable();
  for (i = 0; i < endsite; i++)
    weightrat[i] = weight[i] * rate[category[alias[i] - 1] - 1];
  if (progress) {
    sprintf(progbuf, "Distances calculated for species.\n");
    print_progress(progbuf);
    phyFillScreenColor();
  }

  /* set diagonals */
  for (i = 0; i < spp; i++) {
    if ( method == similarity )
      d[i][i] = 1.0;
    else
      d[i][i] = 0.0;
  }

  for (i = 1; i < spp; i++) {
    if (progress) {
      sprintf(progbuf, "    ");
      print_progress(progbuf);
      for (j = 0; j < nmlngth; j++)
      {
        sprintf(progbuf, "%c", nayme[i - 1][j]);
        print_progress(progbuf);
      }
      sprintf(progbuf, "   ");
      print_progress(progbuf);
    }
    for (j = i + 1; j <= spp; j++) {
      v = makev(i, j); /* compute value for (i, j) */
      d[i - 1][j - 1] = v;
      d[j - 1][i - 1] = v;
      if (progress) {
        sprintf(progbuf, ".");
        print_progress(progbuf);
      }
    }
    if (progress) {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
      phyFillScreenColor();
    }
  }
  if (progress) {
      sprintf(progbuf, "    ");
      print_progress(progbuf);
      for (j = 0; j < nmlngth; j++)
      {
        sprintf(progbuf, "%c", nayme[spp - 1][j]);
        print_progress(progbuf);
      }
      sprintf(progbuf, "\n");
      print_progress(progbuf);
  }
  for (i = 0; i < spp; i++) {
    free(x[i]);
  }
  free(x);
  x = NULL;
}  /* makedists */


void writedists(double **d)
{
  /* write out distances */
  char **names;

  assert( d != NULL );

  names = stringnames_new();
  output_matrix_d(outfile, d, spp, spp, names, names, matrix_flags);
  stringnames_delete(names);

  if (progress)
  {
    sprintf(progbuf, "\nDistances written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
  }
}  /* writedists */


void dnadist_write_outfile(void)
{
  long i;

  if (printdata) {
    fprintf(outfile, "\nNucleic acid sequence Distance Matrix program,");
    fprintf(outfile, " version %s\n\n", VERSION);

    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  }
  if ( method == similarity )
    fprintf(outfile, "  \n  Table of similarity between sequences\n");
  else if (printdata) {
    putc('\n', outfile);
    if (method == jukes)
      fprintf(outfile, "  Jukes-Cantor Distance\n");
    if (method == kimura)
      fprintf(outfile, "  Kimura 2-parameter Distance\n");
    if (method == f84)
      fprintf(outfile, "  F84 Distance\n");
    if (firstset)
    {
      if ((method == kimura) || (method == f84))
        fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n", ttratio);
    }
  }
  if (printdata) {
    if (ctgry && categs > 1) {
      printcategs(outfile, sites, category, "Site categories");
    }
    else if (printdata && (categs > 1)) {
      fprintf(outfile, "\nSite category   Rate of change\n\n");
      for (i = 1; i <= categs; i++)
        fprintf(outfile, "%12ld%13.3f\n", i, rate[i - 1]);
      putc('\n', outfile);
      printcategories();
    }
    if (weights)
      printweights(outfile, 0, sites, oldweight, "Sites");
    output_sequences(sites);
    if ((method == kimura) || (method == f84) || (method == jukes))
      print_basefreq(outfile, &freq, freqsfrom);
  }
}


void dnadistrun(void)
{
  // JRM debug printout
  /*
  printf("\nmethod: %i\n", method) ;
  printf("ctgry: %i\n", ctgry);
  printf("categs: %li\n", categs);
  for (ith = 1; ith <= maxcategs; ith++)
  {
    printf("rate[%li]: %f\n", ith, rate[ith]);
  }
  printf("gama: %i\n", gama);
  printf("cvi: %f\n", cvi);
  printf("invar: %i\n", invar);
  printf("invarfrac: %f\n", invarfrac);
  printf("ttratio: %f\n", ttratio);
  printf("weights: %i\n", weights);
  printf("freqsfrom: %i\n", freqsfrom);
  printf("mulsets: %i\n", mulsets);
  printf("datasets: %li\n", datasets);
  printf("justwts: %i\n", justwts);
  printf("interleaved: %i\n", interleaved);
  printf("printdata: %i\n", printdata);
  printf("dotdiff: %i\n", dotdiff);
  printf("progress: %i\n", progress);
  printf("matrix_flags: %i\n", matrix_flags);
  */

  double ttratio0;
  Matrix_double distance_matrix = NULL;
  long dummy;

  /* Read infile for number of species and sites */
  inputnumbers(&spp, &sites, &dummy, 1);

  /* Allocate global arrays */
  nayme = (naym *)Malloc(spp * sizeof(naym));
  allocsites(sites);
  distance_matrix = Matrix_double_new(spp, spp);

  /* Dataset loop */
  ttratio0 = ttratio;
  firstset = true;
  for (ith = 1; ith <= datasets; ith++) {
    ttratio = ttratio0;
    getinput();
    if (datasets > 1 && progress)
    {
      sprintf(progbuf, "Data set # %ld:\n\n", ith);
      print_progress(progbuf);
    }
    makedists(distance_matrix);
    dnadist_write_outfile();
    writedists(distance_matrix);
    firstset = false;
    fflush(outfile);
    fflush(outtree);
  }

  Matrix_double_delete(&distance_matrix, spp);
  freesites();
  free(nayme);
}


void dnadist(
  char *infilename,
  char *OutfileName,
  char *outfileopt,
  char *catfilename,
  char *weightfilename,
  char *distanceoptions,
  char *gammadistributed,
  double coefvar,
  double pctinvar,
  int usettratio,
  double inttratio,
  int usesubrates,
  int innumcats,
  double inrate1,
  double inrate2,
  double inrate3,
  double inrate4,
  double inrate5,
  double inrate6,
  double inrate7,
  double inrate8,
  double inrate9,
  int useweights,
  int useEmpBF,
  double basefreqA,
  double basefreqC,
  double basefreqG,
  double basefreqTU,
  char *distmatrix,
  int usemultdataset,
  int useData,
  int numdatasets,
  int doseqinter,
  int doprintdata,
  int dodotdiff,
  int doprintind)
{
  //printf("Hello from DnaDist!\n");  // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Dnadist";
  phylipinit(argc, argv, NULL, true); // Initialize PHYLIP

  //********
  ///dnadist_method  method  = f84;          /* method */
  ///boolean ctgry           = false;        /* use more than one category */
  ///long    categs          = 1;            /* number of categories */
  ///double  rate[maxcategs] = { 1.0 };      /* default category rates */
  ///boolean gama            = false;        /* use gamma distributed rates */
  ///double  cvi             = 1.0;          /* coefficient of variation */
  ///boolean invar           = false;        /* use gamma+invar rates */
  ///double  invarfrac       = 0.0;          /* fraction invariant */
  ///double  ttratio         = 2.0;          /* transition/transversion ratio */
  ///boolean weights         = false;        /* use weights */
  ///boolean freqsfrom       = true;         /* use freqs from data */
  ///boolean mulsets         = false;        /* multiple data sets */
  ///long    datasets        = 1;            /* number of sets */
  ///boolean justwts         = false;        /* ...just weights */
  ///boolean interleaved     = true;         /* data set interleaved */
  ///boolean printdata       = false;        /* print data */
  ///boolean dotdiff         = true;         /* use dot-difference format */
  ///boolean progress        = true;         /* enable progress indicator */
  ///boolean matrix_flags    = MAT_MACHINE;  /* Matrix output format */
  //*************
  // sort out booleans and enumerations

  if (!strcmp(distanceoptions, "F84"))
  {
    method = f84;
    if (useEmpBF != 0)
    {
      freqsfrom = true;
    }
    else
    {
      freqsfrom = false;
    }
  }
  else if (!strcmp(distanceoptions, "Jukes"))
  {
    method = jukes;
    freqsfrom = false;
  }
  else if (!strcmp(distanceoptions, "Kimura"))
  {
    method = kimura;
    freqsfrom = false;
  }
  else if (!strcmp(distanceoptions, "LogDet"))
  {
    method = logdet;
    freqsfrom = false;
  }
  else //(!strcmp(distanceoptions, "Similarity"))
  {
    method = similarity;
    freqsfrom = false;
  }

  if (!strcmp(gammadistributed, "Yes"))
  {
    gama = true;
    invar = false;
  }
  if (!strcmp(gammadistributed, "No"))
  {
    gama = false;
    invar = false;
  }
  if (!strcmp(gammadistributed, "Gamma"))
  {
    gama = false;
    invar = true;
  }

  if (usesubrates != 0)
  {
    ctgry = true;
  }
  else
  {
    ctgry = false;
  }

  if (useweights != 0)
  {
    weights = true;
  }
  else
  {
    weights = false;
  }

  if (usettratio != 0)
  {
    //printf("initratio called. ttratio: %f\n", ttratio); // JRMdebug
    initratio(&ttratio);
  }

  if (!strcmp(distmatrix, "Square"))
  {
    matrix_flags = MAT_MACHINE;
  }
  if (!strcmp(distmatrix, "Lower"))
  {
    matrix_flags = MAT_LOWERTRI;
  }
  if (!strcmp(distmatrix, "Human"))
  {
    matrix_flags = MAT_HUMAN;
  }

  if (useData != 0)
  {
    justwts = false;
  }
  else
  {
    justwts = true;
  }

  datasets = numdatasets;
  if (usemultdataset)
  {
    mulsets = true;
    if (justwts)
    {
      justweights(&datasets);
    }
    else
    {
      initdatasets(&datasets);
    }
  }
  else
  {
    mulsets = false;
  }

  if (doseqinter != 0)
  {
    interleaved = true;
  }
  else
  {
    interleaved = false;
  }

  if (doprintdata != 0)
  {
    printdata = true;
  }
  else
  {
    printdata = false;
  }

  if (dodotdiff != 0)
  {
    dotdiff = true;
  }
  else
  {
    dotdiff = false;
  }

  if (doprintind != 0)
  {
    progress = true;
  }
  else
  {
    progress = false;
  }

  // transfer values
  ttratio   = inttratio;
  categs    = innumcats;
  freqa     = basefreqA;
  freqc     = basefreqC;
  freqg     = basefreqG;
  freqt     = basefreqTU;
  cvi       = coefvar;
  invarfrac = pctinvar;

  // warning: if maxcategs ever changes, this has to change also
  rate[0]   = inrate1;
  rate[1]   = inrate2;
  rate[2]   = inrate3;
  rate[3]   = inrate4;
  rate[4]   = inrate5;
  rate[5]   = inrate6;
  rate[6]   = inrate7;
  rate[7]   = inrate8;
  rate[8]   = inrate9;

  // Prevent similarity matrices from easily being used as input to other programs
  if (method == similarity)
  {
    if (matrix_flags == MAT_MACHINE)
    {
      matrix_flags = MAT_HUMAN;
    }
    else if (matrix_flags == MAT_LOWERTRI)
    {
      matrix_flags = MAT_LOWER | MAT_HUMAN;
    }
  }

  // JRMdebug
  /*
  printf("infile: %s \n", infilename);
  printf("printMatrix: %i printcomp: %i\n", printMatrix, printcomp);
  printf("printTree: %i treeprint: %i\n", printTree, treeprint);
  printf("trout: %i\n", trout);
  printf("\nmsets: %li:\n", msets);
  */

  infile = fopen(infilename, "r");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

  if (ctgry)
  {
    catfile = fopen(catfilename, "r");
  }

  if (weights || justwts)
  {
    weightfile = fopen(weightfilename, "r");
  }

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  dnadistrun(); // do the actual work

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

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* DNA Distances by Maximum Likelihood */

#ifdef MAC
  argc = 1;                /* macsetup("Dnadist", "");        */
  argv[0] = "Dnadist";
#endif

  phylipinit(argc, argv, NULL, false); /* Initialize PHYLIP */

  /* Open input and output files */
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  dnadist_config(); /* Run configuration menu */

  /* Open optional files */
  if (ctgry)
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);

  dnadistrun();  // do the actual work

  if (catfile)
    FClose(catfile);
  if (weightfile)
    FClose(weightfile);
  if (infile)
    FClose(infile);
  if (outfile)
    FClose(outfile);

#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return EXIT_SUCCESS;
}  /* DNA Distances by Maximum Likelihood */


// End.
