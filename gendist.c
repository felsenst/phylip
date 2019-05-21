/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "pmatrix.h"

#define epsilong         0.02 /* a small number                            */

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   getalleles(void);
void   inputdata(void);
void   getinput(void);
void   makedists(void);
void   writedists(void);
void   freerest(void);
void   freex(void);
void gendistrun(void);
void gendist(char *infilename, char *outfilename, char *outfileopt, int AllAlleles,
             char *DistMeas, int SquareMat, int MulSets, int NumSets, int PrintInd);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH];
long loci, totalleles, df, datasets, ith, nonodes = 0;
long *alleles;
phenotype3 *x;
double **d;
boolean all, cavalli, lower, nei, reynolds, musat, mulsets, firstset, progress;


void getoptions(void)
{
  /* interactively set options */
  long loopcount;
  Char ch;

  all = false;
  cavalli = false;
  lower = false;
  nei = false;
  reynolds = false;
  musat = true;
  progress = true;
  loopcount = 0;

  for (;;)
  {
    cleerhome();
    printf("\nGenetic Distance Matrix program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  A   Input file contains all alleles at each locus?  %s\n", all ? "Yes" : "One omitted at each locus");
    printf("  N                        Use Nei genetic distance?  %s\n", nei ? "Yes" : "No");
    printf("  C                Use Cavalli-Sforza chord measure?  %s\n", cavalli ? "Yes" : "No");
    printf("  R                   Use Reynolds genetic distance?  %s\n", reynolds ? "Yes" : "No");
    printf("  U    Use delta-mu-squared microsatellite distance?  %s\n", musat? "Yes" : "No");
    printf("  L                         Form of distance matrix?  %s\n", lower ? "Lower-triangular" : "Square");
    printf("  M                      Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0              Terminal type (IBM PC, ANSI, none)?  %s\n", ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1            Print indications of progress of run?  %s\n", progress ? "Yes" : "No");
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("ACNMRL01", ch) != NULL)
    {
      switch (ch)
      {
        case 'A':
          all = !all;
          break;

        case 'C':
          cavalli = true;
          nei = false;
          reynolds = false;
          musat = false;
          break;

        case 'N':
          cavalli = false;
          nei = true;
          reynolds = false;
          musat = false;
          break;

        case 'R':
          reynolds = true;
          cavalli = false;
          nei = false;
          musat = false;
          break;

        case 'U':
          reynolds = false;
          cavalli = false;
          nei = false;
          musat = true;
          break;

        case 'L':
          lower = !lower;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets)
            initdatasets(&datasets);
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          progress = !progress;
          break;
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
  putchar('\n');
}  /* getoptions */


void allocrest(void)
{
  long i;

  x = (phenotype3 *)Malloc(spp * sizeof(phenotype3));
  d = (double **)Malloc(spp * sizeof(double *));
  for (i = 0; i < spp; i++)
    d[i] = (double *)Malloc(spp * sizeof(double));
  alleles = (long *)Malloc(loci * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
}  /* allocrest */


void freerest(void)
{
  long i;

  for (i = 0; i < spp; i++)
    free(d[i]);
  free(d);
  free(alleles);
  free(nayme);
}  /* freerest */


void freex(void)
{
  long i;

  for (i = 0; i < spp; i++)
    free(x[i]);
  free(x);
}


void doinit(void)
{
  /* initializes variables */

  inputnumbers(&spp, &loci, &nonodes, 1);
  if (!javarun)
  {
    getoptions();
  }
}  /* doinit */


void getalleles(void)
{
  long i;

  if (!firstset) {
    samenumsp(&loci, ith);
    free(alleles);
    alleles = (long *)Malloc(loci * sizeof(long));
  }
  totalleles = 0;
  scan_eoln(infile);
  for (i = 0; i < loci; i++) {
    if (eoln(infile))
      scan_eoln(infile);
    if(fscanf(infile, "%ld", &alleles[i]) < 1)
    {
      printf("\nERROR reading input file.\n\n");
      exxit(-1);
    }
    totalleles += alleles[i];
  }
  df = totalleles - loci;
}  /* getalleles */


void inputdata(void)
{
  /* read allele frequencies */
  long i, j, k, m, m1, n;
  double sum;

  for (i = 0; i < spp; i++)
    x[i] = (phenotype3)Malloc(totalleles * sizeof(double));
  for (i = 1; i <= spp; i++) {
    scan_eoln(infile);
    initname(i-1);
    m = 1;
    for (j = 1; j <= loci; j++) {
      sum = 0.0;
      if (all)
        n = alleles[j - 1];
      else
        n = alleles[j - 1] - 1;
      for (k = 1; k <= n; k++) {
        if (eoln(infile))
          scan_eoln(infile);
        if(fscanf(infile, "%lf", &x[i - 1][m - 1]) < 1)
        {
          printf("\nERROR reading input file.\n\n");
          exxit(-1);
        }
        sum += x[i - 1][m - 1];
        if (x[i - 1][m - 1] < 0.0) {
          printf("\nERROR:  Locus %ld in species %ld: an allele", j, i);
          printf(" frequency is negative.\n\n");
          exxit(-1);
        }
        m++;
      }
      if (all && fabs(sum - 1.0) > epsilong) {
        printf("\nERROR:  Locus %ld in species %ld: frequencies do not add up to 1.\n\n", j, i);
        for (m1 = 1; m1 <= n; m1 += 1) {
          if (m1 == 1)
            printf("%f", x[i-1][m-n+m1-2]);
          else {
            if ((m1 % 8) == 1)
              printf("\n");
            printf("+%f", x[i-1][m-n+m1-2]);
          }
        }
        printf(" = %f\n\n", sum);
        exxit(-1);
      }
      if (!all) {
        x[i - 1][m - 1] = 1.0 - sum;
        if (x[i-1][m-1] < -epsilong) {
          printf("\nERROR:  Locus %ld in species %ld: ", j, i);
          printf("frequencies add up to more than 1.\n\n");
          for (m1 = 1; m1 <= n; m1 += 1) {
            if (m1 == 1)
              printf("%f", x[i-1][m-n+m1-2]);
            else {
              if ((m1 % 8) == 1)
                printf("\n");
              printf("+%f", x[i-1][m-n+m1-2]);
            }
          }
          printf(" = %f\n\n", sum);
          exxit(-1);
        }
        m++;
      }
    }
  }
  checknames(spp);                      // Check NAYME array for duplicates.
}  /* inputdata */


void getinput(void)
{
  /* read the input data */
  getalleles();
  inputdata();
}  /* getinput */


void makedists(void)
{
  long i, j, k, m, locus;
  double s, s1, s2, s3, f;
  double temp, temp1, temp2;

  if (progress)
  {
    sprintf(progbuf, "Distances calculated for species.\n");
    print_progress(progbuf);
  }
  for (i = 0; i < spp; i++)
    d[i][i] = 0.0;
  for (i = 1; i <= spp; i++) {
    if (progress) {
      phyFillScreenColor();
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
    for (j = 0; j <= i - 2; j++) { /* ignore the diagonal */
      if (cavalli) {
        s = 0.0;
        for (k = 0; k < (totalleles); k++) {
          f = x[i - 1][k] * x[j][k];
          if (f > 0.0)
            s += sqrt(f);
        }
        d[i - 1][j] = 4 * (loci - s) / df;
      }
      if (nei) {
        s1 = 0.0;
        s2 = 0.0;
        s3 = 0.0;
        for (k = 0; k < (totalleles); k++) {
          s1 += x[i - 1][k] * x[j][k];
          temp = x[i - 1][k];
          s2 += temp * temp;
          temp = x[j][k];
          s3 += temp * temp;
        }
        if (s1 <= 1.0e-20) {
          d[i - 1][j] = -1.0;
          if (progress)
          {
              sprintf(progbuf, "\nWARNING: INFINITE DISTANCE BETWEEN SPECIES ");
              print_progress(progbuf);
              sprintf(progbuf, "%ld AND %ld; -1.0 WAS WRITTEN\n", i, j);
              print_progress(progbuf);
          }
          else
          {
              printf("\nWARNING: INFINITE DISTANCE BETWEEN SPECIES ");
              printf("%ld AND %ld; -1.0 WAS WRITTEN\n", i, j);
          }
        } else
          d[i - 1][j] = fabs(-log(s1 / sqrt(s2 * s3)));
      }
      if (reynolds) {
        s1 = 0.0;
        s2 = 0.0;
        for (k = 0; k < totalleles; k++) {
          temp = x[i - 1][k] - x[j][k];
          s1 += temp * temp;
          s2 += x[i - 1][k] * x[j][k];
        }
        d[i - 1][j] = s1 / (loci * 2 - 2 * s2);
      }
      if (musat) {    /* delta-mu-squared microsatellite distance */
        locus = 1;    /* counter for which locus we are in */
        m = 1;        /* counter for mobility of allele */
        s = 0.0;
        temp1 = 0.0;
        temp2 = 0.0;
        for (k = 0; k < totalleles; k++) {
          temp1 += m * x[i-1][k];   /* add up mean mobility of one species */
          temp2 += m * x[j][k];     /* ... and of other */
          if (m == alleles[locus-1]) {
            locus++;              /* next locus */
            m = 1;                /* reset mobility to 1 */
            s += (temp1-temp2)*(temp1-temp2);  /* add squared difference */
            temp1 = 0.0;
            temp2 = 0.0;
          }
          else m++;
        }
        d[i - 1][j] = s / loci;     /* ... and divide by the number of loci */
      }
      if (progress) {
        sprintf(progbuf, ".");
        print_progress(progbuf);
      }
      d[j][i - 1] = d[i - 1][j];
    }
    if (progress) {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
  }
  if (progress) {
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
}  /* makedists */


void writedists(void)
{
  // Write out distances.

  // RSGnote: Removed "custom" printing of number of species, since output_matrix_d() does so now.
  // fprintf(outfile, "%5ld\n", spp);

  char **names;
  assert( d != NULL );
  names = stringnames_new();

  // RSGnote: Substitution of output_matrix_d() for old custom code.
  output_matrix_d(outfile, d, spp, spp, names, names, lower ? MAT_LOWERTRI : MAT_MACHINE);

  stringnames_delete(names);

  if (progress)
  {
    sprintf(progbuf, "Distances written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
  }
}  /* writedists */


void gendistrun(void)
{
  // debug print JRMdebug
  /*
    printf("\nall %i\n", all);
    printf("cavalli %i\n", cavalli);
    printf("lower %i\n", lower);
    printf("nei %i\n", nei);
    printf("reynolds %i\n", reynolds);
    printf("musat %i\n", musat);
    printf("progress %i\n", progress);
    fflush(stdout);
  */

  // do the work
  for (ith = 1; ith <= (datasets); ith++) {
    allocrest();
    getinput();
    firstset = false;
    if ((datasets > 1) && progress) {
      sprintf(progbuf, "\nData set # %ld:\n\n", ith);
      print_progress(progbuf);
    }
    makedists();
    writedists();
    fflush(outfile);
    freerest();
    freex();
  }
}


void gendist(
  char *infilename,
  char *OutfileName,
  char *outfileopt,
  int AllAlleles,
  char *DistMeas,
  int SquareMat,
  int MulSets,
  int NumSets,
  int PrintInd)
{
  //printf("Hello from GenDist!\n"); // JRMdebug
  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Gendist";
  phylipinit(argc, argv, NULL, true);

  /*
  //  all = false;
  //  cavalli = false;
  //  lower = false;
  //  nei = false;
  //  reynolds = false;
  //  musat = true;
  //  progress = true;
  //    char *infile,
  //    char *outfile,
  //    char *outfileopt,
  //    char *DistMeas,
  //    int SquareMat,
  //    int AnalyzeMult,
  //    int NumMult,
  //    int PrintInd)

  */
  if (AllAlleles != 0)
  {
    all = true;
  }
  else
  {
    all = false;
  }

  cavalli = false;
  nei = false;
  reynolds = false;
  musat = false;
  if (!strcmp(DistMeas, "Nei")) // Nei genetic distance
  {
    nei = true;
  }
  else if (!strcmp(DistMeas, "CS")) // Cavalli-Sforza chord measure
  {
    cavalli = true;
  }
  else if (!strcmp(DistMeas, "Reynolds")) // Reynolds genetic distance
  {
    reynolds = true;
  }
  else //if (!strcmp(Method, "DeltaMu")) //Delta Mu Squared microsatellite distance
  {
    musat = true;
  }

  if (SquareMat != 0)
  {
    lower = false;
  }
  else
  {
    lower = true;
  }

  if (MulSets != 0)
  {
    mulsets = true;
    datasets = NumSets;
  }
  else
  {
    mulsets = false;
    datasets = 1;
  }

  if (PrintInd != 0)
  {
    progress = true;
  }
  else
  {
    progress = false;
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

  firstset = true;
  doinit();
  gendistrun();
  FClose(outfile);
  FClose(infile);
  //printf("Done.\n\n");

}


int main(int argc, Char *argv[])
{  /* main program */
#ifdef MAC
  argc = 1;                /* macsetup("Gendist", "");                */
  argv[0] = "Gendist";
#endif
  phylipinit(argc, argv, NULL, false);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  gendistrun();
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}


// End.
