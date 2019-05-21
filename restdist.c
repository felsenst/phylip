/* Version 4.0. (c) Copyright 1994-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "pmatrix.h"
#include "rest_common.h"
#include "seq.h"


#define initialv        0.1     /* starting value of branch length          */
#define iterationsr     20      /* how many Newton iterations per distance  */

#ifndef OLDC
/* function prototypes */
void restdist_inputnumbers(void);
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void restdist_inputdata(void);
void makev(long, long, double *);
void makedists(void);
void writedists(void);
void getinput(void);
void reallocsites(void);
void restdistrun(void);
void restdist(char * infilename, char * outfilename, char * outfileopt, int RestSites, int OrigNeiLi,
              int GammaDist, double CoeffVar, double TTRatio, double SiteLen, int SquareMat, int AnalyzeMult,
              int NumMult, int Sequential, int PrintData, int PrintInd);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH];
long sites, weightsum, datasets, ith;
boolean  restsites, neili, gama, weights, lower, progress, mulsets, firstset;
double ttratio, fracchange, cvi, sitelength, xi, xv;
double **d;
steptr aliasweight;
char *progname;
Char ch;


void restdist_inputnumbers(void)
{
  /* read and print out numbers of species and sites */
  if(fscanf(infile, "%ld%ld", &spp, &sites) < 2)
  {
    printf("\nERROR reading input file.\n\n");
    exxit(-1);
  }
}  /* restdist_inputnumbers */


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;

  putchar('\n');
  sitelength = 6.0;
  neili = false;
  gama = false;
  cvi = 0.0;
  weights = false;
  lower = false;
  printdata = false;
  progress = true;
  restsites = true;
  interleaved = true;
  ttratio = 2.0;
  loopcount = 0;

  for (;;)
  {
    cleerhome();
    printf("\nRestriction site or fragment distances, ");
    printf("version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  R           Restriction sites or fragments?  %s\n", (restsites ? "Sites" : "Fragments"));
    printf("  N        Original or modified Nei/Li model?  %s\n", (neili ? "Original" : "Modified"));
    if (!neili)
    {
      printf("  G  Gamma distribution of rates among sites?");
      if (!gama)
        printf("  No\n");
      else
        printf("  Yes\n");
      printf("  T            Transition/transversion ratio?  %f\n", ttratio);
    }
    printf("  S                              Site length? %4.1f\n", sitelength);
    printf("  L                  Form of distance matrix?  %s\n", (lower ? "Lower-triangular" : "Square"));
    printf("  M               Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I              Input sequences interleaved?  %s\n", (interleaved ? "Yes" : "No, sequential"));
    printf("  0       Terminal type (IBM PC, ANSI, none)?  %s\n", ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1       Print out the data at start of run?  %s\n", (printdata ? "Yes" : "No"));
    printf("  2     Print indications of progress of run?  %s\n", (progress ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("RDNGTSLMI012", ch) != NULL)
    {
      switch (ch)
      {
        case 'R':
          restsites = !restsites;
          break;

        case 'G':
          if (!neili)
            gama = !gama;
          break;

        case 'N':
          neili = !neili;
          break;

        case 'T':
          if (!neili)
            initratio(&ttratio);
          break;

        case 'S':
          loopcount2 = 0;
          do {
            printf("New Sitelength?\n");
            if(scanf("%lf%*[^\n]", &sitelength)) {} // Read number and scan to EOL.
            (void)getchar();
            if (sitelength < 1.0)
              printf("BAD RESTRICTION SITE LENGTH: %f\n", sitelength);
            countup(&loopcount2, 10);
          } while (sitelength < 1.0);
          break;

        case 'L':
          lower = !lower;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets)
            initdatasets(&datasets);
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
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }

  if (gama) {
    loopcount = 0;
    do {
      printf("\nCoefficient of variation of substitution rate among sites (must be positive)?\n");
      if(scanf("%lf%*[^\n]", &cvi)) {}  // Read number and scan to EOL.
      (void)getchar();
      countup(&loopcount, 100);
    } while (cvi <= 0.0);

    cvi = 1.0 / (cvi * cvi);
    printf("\n");
  }
  xi = (ttratio - 0.5)/(ttratio + 0.5);
  xv = 1.0 - xi;
  fracchange = xi*0.5 + xv*0.75;
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
  rest_makeweights(sites, &endsite);
}


void allocrest(void)
{
  long i;

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *)Malloc(sites * sizeof(Char));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  weight = (steptr)Malloc((sites+1) * sizeof(long));
  alias = (steptr)Malloc((sites+1) * sizeof(long));
  aliasweight = (steptr)Malloc((sites+1) * sizeof(long));
  d = (double **)Malloc(spp * sizeof(double *));
  for (i = 0; i < spp; i++)
    d[i] = (double*)Malloc(spp * sizeof(double));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nRestriction site or fragment distance matrix algorithm, version %s\n\n", VERSION);
  restdist_inputnumbers();
  if (!javarun)
  {
    getoptions();
  }
  if (printdata)
    fprintf(outfile, "\n %4ld Species, %4ld Sites\n", spp, sites);
  allocrest();
}  /* doinit */


void inputoptions(void)
{
  /* read the options information */
  Char ch;
  long i, extranum, cursp, curst;

  if (!firstset)
  {
    if (eoln(infile))
      scan_eoln(infile);
    if(fscanf(infile, "%ld%ld", &cursp, &curst) < 2)
    {
      printf("\nERROR reading input file.\n\n");
      exxit(-1);
    }
    if (cursp != spp)
    {
      printf("\nERROR:  INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld.\n", ith);
      exxit(-1);
    }
    sites = curst;
    reallocsites();
  }
  for (i = 1; i <= sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  if(fscanf(infile, "%*[ 0-9]") < 0)    // RSGnote: How many values must be read?
  {
    printf("\nERROR reading input file.\n\n");
    exxit(-1);
  }
  readoptions(&extranum, "W");
  for (i = 1; i <= extranum; i++)
  {
    matchoptions(&ch, "W");
    inputweights2(1, sites+1, &weightsum, weight, &weights, "RESTDIST");
  }
}  /* inputoptions */


void restdist_inputdata(void)
{
  /* read the species and sites data */
  long i, j, k, l, sitesread = 0, sitesnew = 0;
  Char ch;
  boolean allread, done;

  if (printdata)
    putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 39)
    j = 39;
  if (printdata) {
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
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      ch = gettc(infile);
    } while (ch == ' ' || ch == '\t');
    ungetc(ch, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp ) {
      if ((interleaved && sitesread == 0) || !interleaved)
        initname(i - 1);
      if (interleaved)
        j = sitesread;
      else
        j = 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eoff(infile))) {
          ch = gettc(infile);
          if (ch == '\n' || ch == '\t')
            ch = ' ';
          if (ch == ' ')
            continue;
          uppercase(&ch);
          if (ch != '1' && ch != '0' && ch != '+' && ch != '-' && ch != '?') {
            printf(" ERROR -- Bad symbol %c", ch);
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
        printf("ERROR:  SEQUENCES OUT OF ALIGNMENT.\n");
        exxit(-1);
      }
      i++;
    }
    if (interleaved) {
      sitesread = sitesnew;
      allread = (sitesread == sites);
    } else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata) {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
      for (j = 0; j < spp; j++) {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > sites)
          l = sites;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
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
}  /* restdist_inputdata */


void makev(long m, long n, double *v)
{
  /* compute one distance */
  long i, ii, it, numerator, denominator;
  double f, g=0, h, p1, p2, p3, q1, pp, tt, delta, vv;

  numerator = 0;
  denominator = 0;
  for (i = 0; i < endsite; i++) {
    ii = alias[i + 1];
    if ((inputSequences[m-1][ii-1] == '+') ||
        (inputSequences[n-1][ii-1] == '+')) {
      denominator += weight[i + 1];
      if ((inputSequences[m-1][ii-1] == '+') && (inputSequences[n-1][ii-1] == '+')) {
        numerator += weight[i + 1];
      }
    }
  }
  f = 2*numerator/(double)(denominator+numerator);
  if (restsites) {
    if (exp(-sitelength*1.38629436) > f) {
      printf("\nERROR:  Infinite distance between ");
      printf(" species %3ld and %3ld.\n", m, n);
      exxit(-1);
    }
  }
  if (!restsites) {
    if (!neili) {
      f = (sqrt(f*(f+8.0))-f)/2.0;
    }
    else {
      g = initialv;
      delta = g;
      it = 0;
      while (fabs(delta) > 0.00002 && it < iterationsr) {
        it++;
        h = g;
        g = exp(0.25*log(f * (3-2*g)));
        delta = g - h;
      }
    }
  }
  if ((!restsites) && neili)
    vv = - (2.0/sitelength) * log(g);
  else {
    if (neili && restsites) {
      pp = exp(log(f)/(2*sitelength));
      vv = -(3.0/2.0)*log((4.0/3.0)*pp - (1.0/3.0));
    } else {
      pp = exp(log(f)/sitelength);
      delta = initialv;
      tt = delta;
      it = 0;
      while (fabs(delta) > 0.000001 && it < iterationsr) {
        it++;
        if (gama) {
          p1 = exp(-cvi * log(1 + tt / cvi));
          p2 = exp(-cvi * log(1 + xv * tt / cvi))
            - exp(-cvi * log(1 + tt / cvi));
          p3 = 1.0 - exp(-cvi * log(1 + xv * tt / cvi));
        } else {
          p1 = exp(-tt);
          p2 = exp(-xv * tt) - exp(-tt);
          p3 = 1.0 - exp(-xv * tt);
        }
        q1 = p1 + p2 / 2.0 + p3 / 4.0;
        g = q1 - pp;
        if (g < 0.0)
          delta = fabs(delta) / -2.0;
        else
          delta = fabs(delta);
        tt += delta;
      }
      vv = fracchange * tt;
    }
  }
  *v = fabs(vv);
}  /* makev */


void makedists(void)
{
  /* compute distance matrix */
  long i, j;
  double v;

  if (progress)
  {
    sprintf(progbuf, "Distances calculated for species.\n");
    print_progress(progbuf);
  }
  for (i = 0; i < spp; i++)
    d[i][i] = 0.0;
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
      makev(i, j, &v);
      d[i - 1][j - 1] = v;
      d[j - 1][i - 1] = v;
      if (progress)
      {
        sprintf(progbuf, ".");
        print_progress(progbuf);
      }
    }
    if (progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
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
}  /* makedists */


void writedists(void)
{
  // RSGnote: Removed "custom" printing of number of species, since output_matrix_d() does so now.
  // if (!printdata)
  //   fprintf(outfile, "%5ld\n", spp);

  char **names;
  assert( d != NULL );
  names = stringnames_new();

  // Write out distances.
  // RSGnote: Substitution of output_matrix_d() for old custom code.
  output_matrix_d(outfile, d, spp, spp, names, names, lower ? MAT_LOWERTRI : MAT_MACHINE);

  stringnames_delete(names);

  if (progress)
  {
    sprintf(progbuf, "\nDistances written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
  }
}  /* writedists */


void getinput(void)
{
  /* reads the input data */
  inputoptions();
  restdist_inputdata();
  rest_makeweights(sites, &endsite);
}  /* getinput */


void restdistrun(void)
{
  //sitelength = 6.0;
  //neili = false;
  //gama = false;
  //cvi = 0.0;
  //weights = false;
  //lower = false;
  //printdata = false;
  //progress = true;
  //restsites = true;
  //interleaved = true;
  //ttratio = 2.0;
  // JRMdebug
  /*
    printf("\nsitelength %f\n", sitelength);
    printf("neili %i\n", neili);
    printf("gama %i\n", gama);
    printf("cvi %f\n", cvi);
    printf("weights %i\n", weights);
    printf("lower %i\n", lower);
    printf("printdata %i\n", printdata);
    printf("progress %i\n", progress);
    printf("restsites %i\n", restsites);
    printf("interleaved %i\n", interleaved);
    printf("ttratio %f\n", ttratio);
    fflush(stdout);
  */

  // do the actual work
  for (ith = 1; ith <= datasets; ith++) {
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1 && progress)
    {
        sprintf(progbuf, "\nData set # %ld:\n\n", ith);
        print_progress(progbuf);
    }
    makedists();
    writedists();
    fflush(outfile);
  }
}


void restdist(
  char * infilename,
  char * OutfileName,
  char * outfileopt,
  int RestSites,
  int OrigNeiLi,
  int GammaDist,
  double CoeffVar,
  double TTRatio,
  double SiteLen,
  int SquareMat,
  int AnalyzeMult,
  int NumMult,
  int Sequential,
  int PrintData,
  int PrintInd)
{
  //printf("Hello from RestDist!\n"); // JRMdebug
  //fflush(stdout);

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Restdist";
  phylipinit(argc, argv, NULL, true);
  progname = argv[0];

  /*
  //sitelength = 6.0;
  //neili = false;
  //gama = false;
  //cvi = 0.0;
  //weights = false;
  //lower = false;
  //printdata = false;
  //progress = true;
  //restsites = true;
  //interleaved = true;
  //ttratio = 2.0;
  //char * infile,
  //char * outfile,
  //char * outfileopt,
  //int RestSites,
  //int OrigNeiLi,
  //int GammaDist,
  //double CoeffVar,
  //double TTRatio,
  //double SiteLen,
  //int SquareMat,
  //int AnalyzeMult,
  //int NumMult,
  //int Sequential,
  //int PrintData,
  //int PrintInd)
  */
  weights = false; // never used
  sitelength = SiteLen;
  ttratio = TTRatio;

  if (RestSites != 0)
  {
    restsites = true;
  }
  else
  {
    restsites = false;
  }

  if (OrigNeiLi != 0)
  {
    neili = true;
  }
  else
  {
    neili = false;
  }

  if (GammaDist != 0)
  {
    gama = true;
    cvi = 1.0 / (CoeffVar * CoeffVar);
  }
  else
  {
    gama = false;
    cvi = 0.0;
  }

  if (SquareMat != 0)
  {
    lower = false;
  }
  else
  {
    lower = true;
  }

  if (Sequential != 0)
  {
    interleaved = false;
  }
  else
  {
    interleaved = true;
  }

  if (AnalyzeMult != 0)
  {
    mulsets = true;
    datasets = NumMult;
  }
  else
  {
    mulsets = false;
    datasets = 1;
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

  xi = (ttratio - 0.5)/(ttratio + 0.5);
  xv = 1.0 - xi;
  fracchange = xi*0.5 + xv*0.75;

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

  firstset = true;
  doinit();

  restdistrun();

  FClose(infile);
  FClose(outfile);

  //printf("\ndone\n"); // JRMdebug
  //fflush(stdout);
}


int main(int argc, Char *argv[])
{  /* distances from restriction sites or fragments */
#ifdef MAC
  argc = 1;                /* macsetup("Restdist", ""); */
  argv[0] = "Restdist";
#endif
  phylipinit(argc, argv, NULL, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input data file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();

  restdistrun();

  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Done.\n\n");
  return 0;
}  /* distances from restriction sites or fragments */


// End.
