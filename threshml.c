/* Version 4.0a. (c) Copyright 2004-2013 by the University of Washington.
   Written by Joseph Felsenstein.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/* TODO:  what to do if several adjacent interior nodes are 0 length apart?? */
/* TODO:  deallocate as well as allocate data arrays for each data set? */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cont.h"

typedef double** matrix;

#ifndef OLDC
/* function prototypes */
void getoptions (void);
void readdiscretedata (void);
void getthisname(void);
void readcontinuousdata (void);
void initthreshmlnode(tree *, node **, long, long, long *, long *, initops, pointarray, Char *, Char *, FILE *);
void readthetree (void);
void allocthings (void);
void readdata (void);
void doinit (void);
void initliabilities (node *);
void inittransforms (void);
void initcovariances (void);
void updatequadbranch (node *, boolean);
void updatequadratic (node *, boolean, boolean);
void updateliabilities (node *);
void iterateliabilities (long, boolean);
void cholesky (double **);
void zerocovariances(void);
void updatecovariances (node *);
void invert(matrix);
void correcttransform(boolean);
void givens(matrix, long, long, long, double, double, boolean);
void coeffs(double, double, double *, double *, double);
void tridiag(matrix, long, double);
void shiftqr(matrix, long, double);
void qreigen(matrix, long);
void reportmatrix (matrix);
void reportall(boolean);
double logdet(double**);
void inversesandlogdets(void);
void lrttest(void);
void freethings (void);
void domcmc (void);
/* end of function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], weightfilename[FNMLNGTH];
long inseed, inseed0;
longer seed;
boolean continuous, discrete, lengths, multrees, muldata, treeswithin, datawithin, cross, printdata, progress, firsttree, firstquad, nullchains;
long spp1, chars, chars1, chars2, nonodes, nbranches, proposed, accepted, burnin, cycles, steps, fsteps, sumsteps, numtrees, ndatas, nextnode;
long ith, ithwas, jth, jthwas, kth, sppwas, chars1was, chars2was;
tree curtree;
node *root;
long *zeros, *eigorder, *weightarray, weightsum;
char *progname;
double trweight;   /* added to make treeread happy */
boolean goteof, testing, lasttime, nexttolasttime;
boolean haslengths;   /* end of ones added to make treeread happy */
sequence yy;
initdata *funcs;
double stepsize, lwt, wt, c0, c1, c2, d0, d1, d2, d3, logdetMM0, logdetMM1;
double sumlrt, sumloglrt, sumdenominator;
double **x, *xx, **z, *zz, *eig;
double **AA, **BB, **BB0, **CC, **FF, **MM, **MM0, **MM1, **eigvecs;
double *sumnumerator;
char thisname[MAXNCH];


void getoptions (void)
{
  /* display menu and set the options */
  long loopcount;
  Char ch;
  boolean done;

  fprintf(outfile, "\nThreshold character Maximum Likelihood");
  fprintf(outfile, " method version %s\n\n", VERSION);
  putchar('\n');
  discrete = true;
  continuous = false;
  burnin = 1000;
  cycles = 20;
  steps = 100000;
  fsteps = 1000000;
  stepsize = 0.1;
  multrees = false;
  muldata = false;
  cross = false;
  treeswithin = false;
  datawithin = false;
  numtrees = 1;
  ndatas = 1;
  progress = true;
  lengths = false;
  printdata = false;
  testing = false;
  do {
    loopcount = 0;
    cleerhome();
    printf("\nThreshold character Maximum Likelihood");
    printf(" method version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  D                     Discrete characters?  ");
    if (discrete)
      printf("Yes\n");
    else
      printf("No\n");
    printf("  C                   Continuous characters?  ");
    if (continuous)
      printf("Yes\n");
    else
      printf("No\n");
    printf("  B                           Burn-in steps?  %ld\n", burnin);
    printf("  N                   How many chains to run  %ld\n", cycles);
    printf("  S            Length in steps of each chain  %ld\n", steps);
    if (testing)
      printf("  F   Length in steps of final testing chain  %ld\n", fsteps);
    printf("  P    Size of proposal in Metropolis update  %f\n", stepsize);
    printf("  T  LRT test of independence of characters?");
    if (testing)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  M     Multiple trees?  Multiple data sets?");
    if (muldata) {
      if (multrees) {
        if (cross) {
          if (datawithin)
            printf("  Trees x data sets\n");
          if (treeswithin)
            printf("  Data sets x trees\n");
        } else {
          if (datawithin)
            printf("  Multiple data sets per tree\n");
          if (treeswithin)
            printf("  Multiple trees per data set\n");
        }
      }
      else
        printf("  Multiple data sets, same tree\n");
    }
    else {
      if (multrees)
        printf("  One data set, multiple trees\n");
      else
        printf("  One data set, one tree\n");
    }
    printf("  0      Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1       Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2             Show progress of each chain?  %s\n",
           (progress ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done)
    {
      if (strchr("CDBNSFMTP012", ch) != NULL)
      {
        switch (ch)
        {
          case 'C':
            continuous = !continuous;
            break;

          case 'D':
            discrete = !discrete;
            break;

          case 'B':
            do {
              printf("Number of burn-in steps (must be nonnegative)?\n");
              if(scanf("%ld%*[^\n]", &burnin)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (burnin < 0);
            break;

          case 'N':
            do {
              printf("Number of sampling chains (must be positive)?\n");
              if(scanf("%ld%*[^\n]", &cycles)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (cycles < 0);
            break;

          case 'S':
            do {
              printf("Number of steps in final sampling chain (must be positive)?\n");
              if(scanf("%ld%*[^\n]", &steps)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (steps <= 0);
            break;

          case 'F':
            do {
              printf("Number of steps in final testing chain (must be positive)?\n");
              if(scanf("%ld%*[^\n]", &fsteps)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (fsteps <= 0);
            break;

          case 'M':
            if (multrees && muldata) {
              if (cross) {
                if (treeswithin) {
                  treeswithin = false;
                  datawithin = true;
                }
                else {
                  cross = false;
                  treeswithin = true;
                  datawithin = false;
                }
              } else {
                if (treeswithin) {
                  datawithin = true;
                  treeswithin = false;
                } else{
                  if (datawithin) {
                    muldata = false;
                    multrees = false;
                    datawithin = false;
                  }
                }
              }
            }
            else {
              if (multrees) {
                multrees = false;
                muldata = true;
              } else {
                if (muldata) {
                  multrees = true;
                  treeswithin = true;
                  cross = true;
                } else {
                  multrees = true;
                }
              }
            }
            break;

          case 'T':
            testing = !testing;
            break;

          case 'P':
            do {
              printf("Size of steps in Metropolis proposal (must be positive)?\n");
              printf("(bigger moves more effectively but rejects more often)\n");
              if(scanf("%lf%*[^\n]", &stepsize)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (stepsize <= 0.0);
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
    }
    countup(&loopcount, 100);
  } while (!done);
  if (multrees) {
    printf("Number of trees?\n");
    if(scanf("%ld%*[^\n]", &numtrees)) {} // Read number and scan to EOL.
    (void)getchar();
  }
  if (muldata) {
    printf("Number of data sets?\n");
    if(scanf("%ld%*[^\n]", &ndatas)) {} // Read number and scan to EOL.
    (void)getchar();
  }
  if (multrees && muldata && !cross) {
    if (datawithin && ((ndatas % numtrees) != 0)) {
      printf("ERROR:  In this case, the number of data sets\n");
      printf("        must be a multiple of the number of trees\n");
      printf("        but there are %ld data sets and %ld trees.\n",
             ndatas, numtrees);
      exxit(-1);
    }
    if (treeswithin && ((numtrees % ndatas) != 0)) {
      printf("ERROR:  In this case, the number of trees\n");
      printf("        must be a multiple of the number of data sets\n");
      printf("        but there are %ld data sets and %ld trees.\n",
             ndatas, numtrees);
      exxit(-1);
    }
  }
} /* getoptions */


void allocthings (void)
{
  /* allocate memory for liabilities, matrices etc. now that we know sizes */
  long i;

  if (!continuous) {  /* need them for the liabilities */
    x = (double **)Malloc((long)spp * sizeof(double *));
    for (i = 0; i < spp; i++)
      x[i] = (double *)Malloc((long)chars * sizeof(double));
  }
  z = (double **)Malloc(chars * sizeof(double *));    /* proposed liabilities */
  xx = (double *)Malloc(chars * sizeof(double));      /* the liabilities */
  zz = (double *)Malloc(chars * sizeof(double));    /* proposed liabilities */
  AA = (double **)Malloc(chars * sizeof(double *)); /* transform to dependence */
  for (i = 0; i < chars; i++)
    AA[i] = (double *)Malloc(chars * sizeof(double));
  BB = (double **)Malloc(chars * sizeof(double *)); /* temporary storage */
  for (i = 0; i < chars; i++)
    BB[i] = (double *)Malloc(chars * sizeof(double));
  BB0 = (double **)Malloc(chars * sizeof(double *)); /* temporary storage */
  for (i = 0; i < chars; i++)
    BB0[i] = (double *)Malloc(chars * sizeof(double));
  CC = (double **)Malloc(chars * sizeof(double *));     /* sampled covariances */
  for (i = 0; i < chars; i++)
    CC[i] = (double *)Malloc(chars * sizeof(double));
  FF = (double **)Malloc(chars * sizeof(double *));   /* diff of inverse mats */
  for (i = 0; i < chars; i++)
    FF[i] = (double *)Malloc(chars * sizeof(double));
  MM = (double **)Malloc(chars * sizeof(double *));     /* stored covariances */
  for (i = 0; i < chars; i++)
    MM[i] = (double *)Malloc(chars * sizeof(double));
  MM0 = (double **)Malloc(chars * sizeof(double *));/* constrained covariances */
  for (i = 0; i < chars; i++)
    MM0[i] = (double *)Malloc(chars * sizeof(double));
  MM1 = (double **)Malloc(chars * sizeof(double *));/* stored covariances */
  for (i = 0; i < chars; i++)
    MM1[i] = (double *)Malloc(chars * sizeof(double));
  eigvecs =  (double **)Malloc(chars * sizeof(double *));     /* eigenvectors */
  for (i = 0; i < chars; i++)
    eigvecs[i] = (double *)Malloc(chars * sizeof(double));
  eig = (double *)Malloc(chars * sizeof(double));        /* eigenvalues */
  eigorder = (long *)Malloc(chars * sizeof(long));      /* sort tags for them */
  weightarray = (long *)Malloc(chars * sizeof(long));        /* eigenvalues */
  alloctree(&curtree.nodep, 2*spp-1);
  setuptree(&curtree, 2*spp-1);
  /* allocview(&curtree, 2*spp-1, chars);  debug */
} /* allocthings */


void readcontinuousdata(void)
{
  /* read species data */
  long i, j, k, l;

  if (ith <= 1)
    nayme = (naym *)Malloc(spp * sizeof(naym));
  if (printdata) {
    fprintf(outfile, "Continuous characters\n\n");
    fprintf(outfile, "%4ld Species, %4ld Characters\n\n", spp, chars1);
    fprintf(outfile, "Name");
    fprintf(outfile, "                       Characters\n");
    fprintf(outfile, "----");
    fprintf(outfile, "                       ----------\n\n");
  }
  if (ith <= 1)
    x = (double **)Malloc((long)spp * sizeof(double *)); /* for characters */
  for (i = 0; i < spp; i++) {
    scan_eoln(infile);
    initname(i);                            /* read the name */
    k = i;
    if (ith <= 1)
      x[k] = (double *)Malloc((long)chars * sizeof(double)); /* room for */
    if (printdata) {
      for (j = 0; j < nmlngth; j++)
        fprintf(outfile, "%c", nayme[i][j]);
      fprintf(outfile, "  ");
    }
    for (j = 1; j <= chars1; j++)       /* ... the continuous characters */
    {
      if (eoln(infile))
        scan_eoln(infile);
      if(fscanf(infile, "%lf", &x[k][j - 1]) < 1) // read them in
      {
        printf("\nERROR reading input file.\n\n");
        exxit(-1);
      }
      if (printdata) {      /* print them and go to new line if out of room */
        fprintf(outfile, " %9.5f", x[k][j - 1]);
        if (j % 6 == 0) {
          putc('\n', outfile);
          for (l = 1; l <= nmlngth; l++)
            putc(' ', outfile);
        }
      }
    }
    if (printdata)
      putc('\n', outfile);
  }
  scan_eoln(infile);
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata)
    putc('\n', outfile);
}  /* readcontinuousdata */


void getthisname(void)
{
  /* read in species name and put it into array  thisname  */
  long j;

  for (j = 0; j < nmlngth; j++) {
    if (eoff(infile) || eoln(infile)) {
      printf("\nERROR:  End-of-Line or End-of-File in the middle of species name.\n\n");
      exxit(-1);
    }
    thisname[j] = gettc(infile);
    if ((thisname[j] == '(') || (thisname[j] == ')') || (thisname[j] == ':')
        || (thisname[j] == ',') || (thisname[j] == ';') || (thisname[j] == '[')
        || (thisname[j] == ']')) {
      printf("\nERROR:  Species name may not contain characters ( ) : ; , [ ] \n");
      printf("        In name of this species there is character %c.\n\n", thisname[j]);
      exxit(-1);
    }
  }
} /* gethisname */


void readdiscretedata (void)
{
  /* read in the 0/1 data */
  // RSGnote: Variable "basesread" never used; removed.
  long i, j, k, l;
  Char charstate;
  boolean allread, done, found;

  if (ith <= 1) {
    if (!continuous)
      nayme = (naym *)Malloc(spp * sizeof(naym));
    yy = (Char **)Malloc(spp * sizeof(Char *)); /* where store the characters */
    for (i = 0; i < spp; i++)
      yy[i] = (Char *)Malloc(chars * sizeof(Char));
  }
  if (printdata)
    headings(chars, "Characters", "-----------");
  allread = false;
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {                                   /* read, skipping blanks and tabs */
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {                            /* loop over species */
      if (!continuous) {  /* if there is no previous continuous data set */
        initname(i-1);   /* just take the name as is */
        k = i;
      }
      else {
        getthisname();   /* take the name for the (i+1)-th species */
        k = 0;
        found = false;
        while (!found) {    /* check to see if it is one of the species */
          found = true;
          for (l = 0; l < nmlngth; l++)
            found = found && (nayme[k][l] == thisname[l]);
          if (!found)
            k++;
          if (k >= spp)
            break;
        }
        if (!found) {
          printf("ERROR:  Name:  ");
          for (l = 0; l < nmlngth; l++)
            printf("%c", thisname[l]);
          printf(" not found in the discrete characters data set.\n\n");
          exxit(-1);
        }
      }
      j = 0;
      done = false;
      while (!done && !eoff(infile)) {
        while (j < chars2 && !(eoln(infile) || eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ')                   /* skip over blanks */
            continue;
          if ((strchr("01?", charstate)) == NULL) { /* look for bad characters */
            printf("\nERROR:  Bad symbol: %c at position %ld of species %ld.\n\n", charstate, j+1, i);
            exxit(-1);
          }
          j++;
          yy[i - 1][j - 1] = charstate;            /* store the character */
        }
        if (j < chars2)
          scan_eoln(infile);
        else if (j == chars2)
          done = true;
      }
      scan_eoln(infile);
      i++;
    }
    allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata) {                /* if requested, print out a table of data */
    for (i = 1; i <= ((chars2 - 1) / 60 + 1); i++) {  /* ranges of characters */
      for (j = 1; j <= spp; j++) {          /* loop over species */
        for (k = 0; k < nmlngth; k++)       /* print name */
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars2)
          l = chars2;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {     /* print characters */
          charstate = yy[j - 1][k - 1];
          putc(charstate, outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
} /* readdiscretedata */


void readdata (void)
{
  /* read in the data */
  chars = 0;
  chars1 = 0;
  chars2 = 0;
  if (continuous) {       /* if there's a data set of continuous characters */
    inputnumbers(&spp1, &chars1, &nonodes, 1);
    if (ith > 1) {
      if (sppwas != spp1) {
        printf("\nERROR:  Number of species in data set %ld is not the same\n", ith);
        printf("        as in all the previous data sets.\n");
      }
      if (chars1was != chars1) {
        printf("\nERROR:  Number of characters in continuous data set %ld is not the same\n", ith);
        printf("        as in all the previous data sets.\n");
      }
      if ((sppwas != spp) || (chars1was != chars1)) {
        printf("\n");
        exxit(-1);
      }
    }
    spp = spp1;
    chars = chars1;
    readcontinuousdata();
  }
  if (discrete) {          /* if there's a data set of discrete characters */
    inputnumbers(&spp, &chars2, &nonodes, 1);
    if (continuous && (spp1 != spp)) {
      printf("PROBLEM:  number of species in the continuous data set\n");
      printf(" does not match the number in the discrete data set.\n\n");
      exxit(-1);
    }
    if (ith > 1) {
      if (sppwas != spp) {
        printf("\nERROR:  Number of species in data set %ld is not the same\n", ith);
        printf("        as in all the previous data sets.\n");
      }
      if (chars2was != chars2) {
        printf("\nERROR:  Number of characters in discrete data set %ld is not the same\n", ith);
        printf("        as in all the previous data sets.\n");
      }
      if ((sppwas != spp) || (chars2was != chars2)) {
        printf("\n");
        exxit(-1);
      }
    }
    chars += chars2;
    readdiscretedata();
  }
  sppwas = spp;            /* record numbers of species, characters */
  chars1was = chars1;
  chars2was = chars2;
} /* readdata */


void doinit (void)
{
  /* read in the menu options and the data */
  getoptions();
  initseed(&inseed, &inseed0, seed);
} /* doinit */


void initthreshmlnode(tree * treep, node **p, long len,
                      long nodei, long *ntips, long *parens,
                      initops whichinit, pointarray treenode,
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
      treep->nodep[(*p)->index - 1] = (*p);
      ((cont_node_type*)(*p))->view = (phenotype3)Malloc((long)chars * sizeof(double));
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      (*p)->index = nodei;
      (*p)->tip = false;
      treep->nodep[(*p)->index - 1] = (*p);
      ((cont_node_type*)(*p))->view = (phenotype3)Malloc((long)chars * sizeof(double));
      break;
    case tip:
      match_names_to_data (str, treep->nodep, p, spp);
      ((cont_node_type*)(*p))->view = (phenotype3)Malloc((long)chars * sizeof(double));
      (*p)->deltav = 0.0;
      break;
    case length:
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      (*p)->v = valyew / divisor;
      if ((*p)->back != NULL)
      {
        (*p)->back->v = (*p)->v;
      }
      break;
    default:        /* cases of hslength, iter, hsnolength, treewt, unittrwt*/
      break;        /* not handled                                            */
  }
} /* initthreshmlnode */


void readthetree (void)
{
  /* read in the tree */
  nextnode = 0;
  haslengths = true;
  treeread(&curtree, intree, &root, curtree.nodep, &goteof, &firsttree,
           &nextnode, &haslengths, initthreshmlnode, false, nonodes);
} /* readthetree */


void inittransforms (void)
{
  /* set up the initial value of the transform */
  long i, j;

  for (i = 0; i < chars; i++) {     /* it is initially the identity matrix */
    for (j = 0; j < chars; j++)
      AA[i][j] = 0.0;
    AA[i][i] = 1.0;
  }
} /* inittransforms */


void initliabilities (node *p)
{
  /* set up the initial values of the liabilities */
  node *q, *qq;
  long i;
  double wt;

  if (p->tip) {        /* initially we assume cutoff is at 50% point, so ... */
    for (i = chars1; i < chars; i++) {                 /* for each character */
      if (yy[p->index - 1][i] == '0')
        ((cont_node_type*)p)->view[i] = -0.79788456;   /* mean of left half
                                                          of standard normal */
      else {
        if (yy[p->index - 1][i] == '1')
          ((cont_node_type*)p)->view[i] = 0.79788456;  /* mean of right half
                                                          of standard normal */
        else {
          if (yy[p->index - 1][i] == '?')
            ((cont_node_type*)p)->view[i] = 0.0;    /* mean of the whole
                                                       standard normal */
        }
      }
    }
    for (i = 0; i < chars1; i++)
      ((cont_node_type*)p)->view[i] = x[p->index - 1][i];    /* copy in the
                                                                character value */
  }
  else {
    for (i = 0; i < chars; i++)
      sumnumerator[i] = 0.0;
    sumdenominator = 0.0;
    q = p;
    qq = p->next;
    while (qq != q) {        /* for each branch descended from this node ... */
      initliabilities(qq->back);
      if (qq->v <= 0.0)
        wt = 1.0e9;
      else wt = 1/qq->v;
      for (i = 0; i < chars; i++)       /* get pieces for a weighted average */
        sumnumerator[i] += wt * ((cont_node_type*)qq->back)->view[i];
      sumdenominator += wt;
      qq = qq->next;
    }
    for (i = 0; i < chars; i++)        /* weighted average in each character */
      ((cont_node_type*)p)->view[i] = sumnumerator[i] / sumdenominator;
    q = p;
    qq = p->next;
    while (qq != q) {                /* for each nodelet on this circle ... */
      for (i = 0; i < chars; i++)         /* ... propagate the view */
        ((cont_node_type*)qq)->view[i] = ((cont_node_type*)q)->view[i];
      qq = qq->next;
    }
  }
} /* initliabilities */


void initcovariances (void)
{
  /* initialize covariances among liabilities */
  long i, j;

  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      CC[i][j] = 0.0;
  }
} /* initcovariances */


void updatequadbranch (node* p, boolean negative)
{
  long i, j;
  double sum, vv;
  /* subtract or add terms for one branch */

  vv = p->v;
  for (i = 0; i < chars; i++)    /* get change along each branch */
    zz[i] = ((cont_node_type*)p)->view[i] - ((cont_node_type*)(p->back))->view[i];
  sum = 0.0;
  for (i = 0; i < chars; i++)     /* quadratic form */
    for (j = 0; j < chars; j++)
      sum -= zz[i]*zz[j]*FF[i][j]/vv;
  if (negative)
    sumloglrt -= 0.5*sum;
  else
    sumloglrt += 0.5*sum;
} /* updatequadbranch */


void updatequadratic (node* p, boolean negative, boolean furstquad)
{
  /* subtract or add terms from the quadratic form near a node */
  node* q;

  if (p->tip) {   /* just do the ancestor branch for a tip */
    updatequadbranch (p, negative);
  } else {
    if (!root) {                  /* don't do the branch below the root */
      updatequadbranch (p, negative);
    }
    if (furstquad) {  /* for first time do not add descendant branches */
      q = p;
      do {                         /* loop among descendants */
        q = q->next;
        if (q != p)
          updatequadbranch (q, negative);
      } while (q != p);
    }
  }
} /* updatequadratic */


void updateliabilities (node* p)
{
  /* do a Gibbs update of interior node or a Metropolis update of a tip.
     Calls itself recursively to do this throughout the tree.  If this
     is the first tree of a chain (boolean "firstquad"), instead go through
     the tree (recursively) and sum up the quadratic terms of the
     likelihood ratio */
  double mean, var, sd, sum, sum2, wt, xold, xnew;
  long i, j;
  boolean rejected;
  node *q;

  if (nullchains && lasttime && !firstquad)
    updatequadratic (p, true, false);   /* subtract old terms from quadratic */
  if (p->tip) {                     /* weighted sampling of a tip value */
    /* draw new liabilities */
    if ((*p).v <= 0) {          /* if branch has zero length, just copy */
      for (i = 0; i < chars; i++)
        xx[i] = ((cont_node_type*)p->back)->view[i];
    }
    else {                           /* draw a new set of tip liabilities */
      rejected = false;
      proposed++;    /* counts only the Metropolis steps, not Gibbs steps */
      for (i = 0; i < chars1; i++) {
        xx[i] = ((cont_node_type*)p)->view[i];
        zz[i] = xx[i];   /* tip values left unchanged in continuous chars */
      }
      for (i = chars1; (i < chars) && (!rejected); i++) {
        xx[i] = ((cont_node_type*)p)->view[i];
        zz[i] = xx[i] + stepsize * normrand(seed); /* propose a liability */
        sum = 0.0;
        for (j = 0; j <= i; j++)          /* ... rejecting as soon as can */
          sum += AA[i][j] * zz[j];               /* compute new liability */
        rejected = ((yy[p->index - 1][i-chars1] == '0') && (sum > 0)) ||
          ((yy[p->index - 1][i-chars1] == '1') && (sum < 0));
      }
      sum = 0.0;
      var = (*p).v;
      if (!rejected) {
        for (i = chars1; i < chars; i++) {   /* or reject if density at new  */
          xold = xx[i];                        /* ... location is too low */
          xnew = zz[i];
          sum = sum + ((xnew-xold)  /* log of ratio of densities after, before */
                       *(xnew+xold-2.0*((cont_node_type*)p->back)->view[i]))
            /(2.0*var);
        }
        if (randum(seed) < exp(-sum)) {      /* else replace by new values */
          accepted++;
          for (i = chars1; i < chars; i++) { /* no need to change others */
            xx[i] = zz[i];
            ((cont_node_type*)p)->view[i] = xx[i];
          }
        }
        else rejected = true;
      }
    }
    if (nullchains && lasttime) {          /* add new terms to quadratic */
      updatequadratic (p, false, firstquad);
      if ((chars1 < chars) && !firstquad)
        sumsteps++;
    }
  }
  else {                          /* Gibbs sampler update of interior node */
    if ((p != root) && (p->v <= 0.0))      /* If too close to have changed */
      for (i = 0; i < chars; i++)
        ((cont_node_type*)p)->view[i] = ((cont_node_type*)(p->back))->view[i];
    else {                                             /* if not too close */
      sum2 = 0.0;                   /* this accumulates the sum of weights */
      q = p;
      if (p != root) {                                       /* at node  p */
        if ((*q).v <= 0.0)
          sum2 += 1.0e+9;
        else
          sum2 += 1.0/(*q).v;      /* sum of reciprocals of branch lengths */
      }
      q = q->next;
      while (q != p) {            /* at the others connected to its circle */
        if (q != root) {
          if ((*q->back).v <= 0.0)
            wt = 1.0e+9;
          else
            wt = 1.0/(*q->back).v;
          sum2 += wt;
        }
        q = q->next;
      }
      sd = sqrt(1.0/sum2);
      for (i = 0; i < chars; i++) {
        sum = 0.0;
        if (q != root) {      /* accumulates numerator of weighted average */
          if ((*q).v <= 0.0)
            sum = 1.0e+9*((cont_node_type*)(q->back))->view[i];
          else
            sum = ((cont_node_type*)(q->back))->view[i]/(*q->back).v;
        }
        q = q->next;
        while (q != p) {               /* ... for other nodes in its circle */
          if (q != root) {
            if ((*q->back).v <= 0.0)
              wt = 1.0e+9;
            else
              wt = 1.0/(*q->back).v;
            sum = sum + wt * ((cont_node_type*)(q->back))->view[i];
          }
          q = q->next;
        }
        mean = sum / sum2;             /* compute mean, standard deviation */
        ((cont_node_type*)p)->view[i] = xx[i] =
          mean + sd*normrand(seed);                   /* draw value */
        q = p;
        q = q->next;
        while (q != p) {               /* set all liabilities on that node */
          ((cont_node_type*)q)->view[i] = ((cont_node_type*)p)->view[i];
          q = q->next;
        }
      }
    }
    if (lasttime && testing) {
      updatequadratic (p, false, firstquad);   /* add new terms to quadratic */
      if (!firstquad) {
        sumlrt += exp(sumloglrt);
        sumsteps++;
      }
    }
    q = p;
    q = q->next;              /* if not, don't go out the one you came from */
    while (q != p) {                /* tree traversal to do this everywhere */
      updateliabilities(q->back);
      q = q->next;
    }
  }
} /* updateliabilities */


void zerocovariances (void)
{
  /* zero out the covariances */
  long i, j;

  for (i = 0; i < chars; i++)
    for (j = 0; j < chars; j++)
      CC[i][j] = 0.0;
} /* zerocovariances */


void updatecovariances (node* p)
{
  /* compute new values for the covariances */
  long i, j;
  node* q;

  if (p != root) {          /* compute covariances for branch below this */
    nbranches++;
    for (i = 0; i < chars; i++)
      for (j = 0; j < chars; j++)
        if ((*p).v > 0.0)
          CC[i][j] += (((cont_node_type*)p)->view[i]-((cont_node_type*)(p->back))->view[i])
            * (((cont_node_type*)p)->view[j]-((cont_node_type*)(p->back))->view[j]) / (*p).v;
  }
  if (!((*p).tip)) {
    q = p->next;
    while (q != p) {             /* tree traversal to do this on each branch */
      updatecovariances(q->back);
      q = q->next;
    }
  }
} /* updatecovariances */


void iterateliabilities (long timesteps, boolean estimatenewcovs)
{
  /* update the liabilities by Gibbs sampling in the interior nodes
     and rejection sampling at tips */
  long t, n, i, j, blanks;

  n = timesteps / 25;               /* how often we print out a dot */
  if (estimatenewcovs)
    initcovariances();
  proposed = 0;
  accepted = 0;
  if (progress) {
    blanks = 0;
    if (testing && !lasttime)
      blanks = (long)(log((double)fsteps/steps)/log(10.0)+0.01);
    for (i = 0; i < blanks; i++)
      printf(" ");
    fflush(stdout);
  }
  zerocovariances();
  if (lasttime && testing) {        /* get ready to add up LRT terms */
    sumlrt = 0.0;
    firstquad = true;
  }
  for (t = 0; t < timesteps; t++) { /* run a chain this long */
    if (progress) {
      if ((t%n) == 0) {
        printf("."); fflush(stdout);  /* print out a dot every so often */
      }
    }
    if (lasttime && nullchains)
      sumloglrt = 0.0;              /* get ready to get term for this tree */
    updateliabilities(root);        /* traverse tree to update liabilities */
    if (lasttime && nullchains && firstquad) {
      sumlrt += exp(sumloglrt);     /* add in term for this sampling */
      sumsteps++;                   /* increase count in number of samplings */
    }
    firstquad = false;              /* after first time make sure know that */
    if (estimatenewcovs) {
      nbranches = 0;
      updatecovariances(root);
    }
  }
  if (estimatenewcovs) {
    for (i = 0; i < chars; i++)     /* complete estimate of the covariances */
      for (j = 0; j < chars; j++)
        CC[i][j] /= (timesteps * nbranches);
  }
  if (progress) {
    if (discrete) /* report acceptance rate of Metropolis algortihm */
      printf("  %5.4f ", (double)accepted/proposed);
    else
      printf("  1.0000 ");   /* if only a Gibbs sampler, all are accepted */
  }
} /* iterateliabilities */


void cholesky (double **c)
{ /* Cholesky matrix square root in place,
     overwriting previous contents of c.  On exit, matrix c
     is the lower-triangular square root of its previous value */
  long i, j, k;
  double sum, temp;

  for (i = 0; i < chars; i++) { /* in-place Cholesky decomposition of C */
    sum = 0.0;
    for (j = 0; j <= i-1; j++)
      sum = sum + c[i][j] * c[i][j];
    if (c[i][i] <= sum)
      temp = 0.0;
    else
      temp = sqrt(c[i][i] - sum);
    c[i][i] = temp;
    for (j = i+1; j < chars; j++) {
      sum = 0.0;
      for (k = 0; k < i; k++)
        sum = sum + c[i][k] * c[j][k];
      if (fabs(temp) < 1.0E-12)
        c[j][i] = 0.0;
      else
        c[j][i] = (c[j][i] - sum)/temp;
      c[i][j] = 0.0;
    }
  } /* c is now the lower-triangular square root of its previous value */
} /* cholesky */


void invert(double **a)
{
  /* Gauss-Jordan reduction -- invert matrix a in place, overwriting previous
     contents of a.  On exit, matrix a contains the inverse. */
  long i, j, k;
  double temp;

  for (i = 0; i < chars; i++) {
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < chars; j++)
      a[i][j] *= temp;
    for (j = 0; j < chars; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < chars; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
}  /* invert */


void correcttransform(boolean nullchains)
{
  /* use the inferred covariances, square-rooting them, to correct the
   * transform from independence to the correlated liabilities */
  long i, j, k;
  double sum, norm;

  if (nullchains) {               /* hold any covariances zero as needed */
    for (i = 0; i < chars; i++)
      for (j = 0; j < chars; j++) {
        if (weightarray[i] != weightarray[j])
          CC[i][j] = 0.0;
      }
  }
  cholesky(CC);                 /* take matrix square root of covariances */
  for (i = 0; i < chars; i++)   /* multiply transform by CC to correct it */
    for (j = 0; j < chars; j++) {  /* note, the actual transform used on the */
      sum = 0;                     /*  original data is inverse of AA      */
      for (k = 0; k < chars; k++)
        sum += AA[i][k] * CC[k][j];
      BB[i][j] = sum;  /* don't update it right away as AA gets re-used */
    }
  norm = 0.0;                   /* get norm of change in transform */
  for (i = 0; i < chars; i++)
    for (j = 0; j < chars; j++) /* sum of squared changes of transform */
      norm += (AA[i][j]-BB[i][j])*(AA[i][j]-BB[i][j]);
  if (norm > 0.0)               /* norm is the root-sum-of-squares */
    norm = sqrt(norm);
  if (progress)
    printf("        %f\n", norm);  /* print out the norm so we can monitor */
  for (i = 0; i < chars; i++)   /* now update the transform AA */
    for (j = 0; j < chars; j++)
      AA[i][j] = BB[i][j];
  for (i = 0; i < chars; i++)   /* save a copy of overall covariances in MM */
    for (j = 0; j < chars; j++) {
      sum = 0.0;
      for (k = 0; k < chars; k++)
        sum += AA[i][k]*AA[j][k];
      MM[i][j] = sum;
    }
  /* get inverse of CC.  Multiply the tip and interior coordinates by it.  */
  invert(CC);                           /* do an in-place inverse of CC */
  for (i = 0; i < nonodes; i++) { /* transform coordinates to correct for CC */
    for (j = 0; j < chars; j++) {
      sum = 0.0;
      for (k = 0; k < chars; k++)
        sum += CC[j][k]*((cont_node_type*)curtree.nodep[i])->view[k];
      xx[j] = sum;
    }
    for (j = 0; j < chars; j++) {
      ((cont_node_type*)curtree.nodep[i])->view[j] = xx[j];
    }
  }
  if (testing && lasttime && !nullchains) {
    for (i = 0; i < chars; i++)   /* save a copy of covariances in MM1 */
      for (j = 0; j < chars; j++)
        MM1[i][j] = MM[i][j];
  }
  if (testing && nexttolasttime && nullchains) {
    for (i = 0; i < chars; i++)   /* save a copy of covariances in MM0 */
      for (j = 0; j < chars; j++)
        MM0[i][j] = MM[i][j];
  }
} /* correctransform */


void givens(double **a, long i, long j, long n, double ctheta,
            double stheta, boolean left)
{ /* Givens transform at i, j for 1..n with angle theta */
  long k;
  double d;

  for (k = 0; k < n; k++) {
    if (left) {
      d = ctheta * a[i - 1][k] + stheta * a[j - 1][k];
      a[j - 1][k] = ctheta * a[j - 1][k] - stheta * a[i - 1][k];
      a[i - 1][k] = d;
    } else {
      d = ctheta * a[k][i - 1] + stheta * a[k][j - 1];
      a[k][j - 1] = ctheta * a[k][j - 1] - stheta * a[k][i - 1];
      a[k][i - 1] = d;
    }
  }
}  /* givens */


void coeffs(double x, double y, double *c, double *s, double accuracy)
{ /* compute cosine and sine of theta */
  double root;

  root = sqrt(x * x + y * y);
  if (root < accuracy) {
    *c = 1.0;
    *s = 0.0;
  } else {
    *c = x / root;
    *s = y / root;
  }
}  /* coeffs */


void tridiag(double **a, long n, double accuracy)
{ /* Givens tridiagonalization */
  long i, j;
  double s, c;

  for (i = 2; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
      coeffs(a[i - 2][i - 1], a[i - 2][j - 1], &c, &s, accuracy);
      givens(a, i, j, n, c, s, true);
      givens(a, i, j, n, c, s, false);
      givens(eigvecs, i, j, n, c, s, true);
    }
  }
}  /* tridiag */


void shiftqr(double **a, long n, double accuracy)
{ /* QR eigenvalue-finder */
  long i, j;
  double approx, s, c, d, TEMP, TEMP1;

  for (i = n; i >= 2; i--) {
    do {
      TEMP = a[i - 2][i - 2] - a[i - 1][i - 1];
      TEMP1 = a[i - 1][i - 2];
      d = sqrt(TEMP * TEMP + TEMP1 * TEMP1);
      approx = a[i - 2][i - 2] + a[i - 1][i - 1];
      if (a[i - 1][i - 1] < a[i - 2][i - 2])
        approx = (approx - d) / 2.0;
      else
        approx = (approx + d) / 2.0;
      for (j = 0; j < i; j++)
        a[j][j] -= approx;
      for (j = 1; j < i; j++) {
        coeffs(a[j - 1][j - 1], a[j][j - 1], &c, &s, accuracy);
        givens(a, j, j + 1, i, c, s, true);
        givens(a, j, j + 1, i, c, s, false);
        givens(eigvecs, j, j + 1, n, c, s, true);
      }
      for (j = 0; j < i; j++)
        a[j][j] += approx;
    } while (fabs(a[i - 1][i - 2]) > accuracy);
  }
}  /* shiftqr */


void qreigen(double **prob, long n)
{ /* QR eigenvector/eigenvalue method for symmetric matrix. Leaves
     right eigenvectors as columns in eigvecs, and their eigenvalues in prob */
  double accuracy;
  long i, j;

  accuracy = 1.0e-6;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      eigvecs[i][j] = 0.0;
    eigvecs[i][i] = 1.0;
  }
  tridiag(prob, n, accuracy);
  shiftqr(prob, n, accuracy);
  for (i = 0; i < n; i++)
    eig[i] = prob[i][i];
}  /* qreigen */


void reportmatrix (matrix QQ) {
  /* print out the covariances and the axes that are independent */
  long i, j;

  fprintf(outfile, "       ");
  for (i = 0; i < chars; i++) {          /* print out column numbers */
    if ((i > 0) && (i%6 == 0))
      fprintf(outfile, "\n        ");
    fprintf(outfile, "     %2ld    ", i+1);
    if ((i+1) == chars)
      fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
  for (i = 0; i < chars; i++) {          /* print out the matrix */
    fprintf(outfile, " %2ld   ", i+1);
    for (j = 0; j < chars; j++) {
      fprintf(outfile, " %10.5f", QQ[i][j]);
      if ((j > 0) && (j%6 == 5) && ((j+1) != chars)) /* in six-column blocks */
        fprintf(outfile, "\n       ");
      if ((j+1) == chars)
        fprintf(outfile, "\n");
    }
  }
  fprintf(outfile, "\n\n");
} /* reportmatrix */


void reportall(boolean constrained)
{
  /* print out the transforms, covariance matrices and principal axes */
  long i, j;

  if (!constrained) {
    if (discrete && continuous) {
      fprintf(outfile, " In this output");
      if (chars1 == 1)
        fprintf(outfile, " the first character is continuous,\n");
      else
        fprintf(outfile, " the first %ld characters are continuous,\n", chars1);
      if (chars2 == 1)
        fprintf(outfile, " and the next character is discrete.\n\n");
      else
        fprintf(outfile, " and the next %ld characters are discrete.\n\n", chars2);
    }
  }
  if (constrained)
    fprintf(outfile, " Constrained covariance matrix of");
  else
    fprintf(outfile, " Covariance matrix of");
  if (continuous) {
    fprintf(outfile, " continuous character");
    if (chars1 > 1)
      fprintf(outfile, "s");
  }
  if (discrete) {
    if (continuous)
      fprintf(outfile, "\n and");
    if (chars2 > 1)
      fprintf(outfile, " liabilities of discrete characters");
    else
      fprintf(outfile, " liability of the discrete character");
  }
  fprintf(outfile, "\n");
  if (discrete && continuous)
    fprintf(outfile, " (the continuous characters are first)\n");
  fprintf(outfile, "\n");
  reportmatrix(MM);
  if (constrained) {
    fprintf(outfile, " In the constrained case,\n");
    fprintf(outfile, " transform from independent variables (columns)\n");
  }
  else
    fprintf(outfile, " Transform from independent variables (columns)\n");
  fprintf(outfile, " to");
  if (discrete)
    fprintf(outfile, " liabilities");
  if (continuous) {
    if (discrete)
      fprintf(outfile, " or");
    fprintf(outfile, " characters");
  }
  fprintf(outfile, " (rows)\n\n");
  reportmatrix(BB);
  qreigen (MM, chars);
  fprintf(outfile, "   For            Variance of\n");
  fprintf(outfile, "   variable       its change\n");
  fprintf(outfile, "   --------       -----------\n");
  for (i = 0; i < chars; i++)   /* set up sort tags for eigenvalue order */
    eigorder[i] = i;
  shellsort(eig, eigorder, chars);  /* sort eigenvalues in ascending order */
  for (i = 0; i < chars; i++)    /* print out the eigenvectors (variances) */
    fprintf(outfile, "  %6ld        %12.6f\n", i+1, eig[chars-i-1]);
  fprintf(outfile, "\n");
  for (i = 0; i < chars; i++)     /* copy eigenvectors into descending order */
    for (j = 0; j < chars; j++)
      BB[i][j] = eigvecs[eigorder[chars-i-1]][j];
  for (i = 0; i < chars; i++)     /* and copy back into array eigvecs */
    for (j = 0; j < chars; j++)
      eigvecs[i][j] = BB[i][j];
  if (constrained) {
    fprintf(outfile, " For the constrained estimates,\n");
    fprintf(outfile, " transform from liabilities or characters (columns)\n");
  }
  else
    fprintf(outfile, " Transform from liabilities or characters (columns)\n");
  fprintf(outfile, " to continuous variables or liabilities (rows)\n\n");
  reportmatrix (eigvecs);
} /* reportall */


double logdet(double **a)
{
  /* Gauss-Jordan log determinant calculation.
     in place, overwriting previous contents of a.  On exit,
     matrix a contains the inverse. Works only for positive definite A */
  long i, j, k;
  double temp, sum;

  sum = 0.0;
  for (i = 0; i < chars; i++) {
    if (fabs(a[i][i]) < 1.0E-37) {
      printf("ERROR:  Function logdet tried to invert singular matrix.\n");
      exxit(-1);
    }
    sum += log(a[i][i]);
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < chars; j++)
      a[i][j] *= temp;
    for (j = 0; j < chars; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < chars; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
  return(sum);
}  /* logdet */


void inversesandlogdets(void)
{
  /* compute for next-to-last chain, if testing, the log determinants of the
     estimated covariance matrix MM1 and the constrained one MM0, and their
     inverses */
  long i, j;

  logdetMM1 = logdet(MM1);   /* log det of covariances and invert them */
  logdetMM0 = logdet(MM0);   /* ditto for the constrained covariances */
  for (i = 0; i < chars; i++) /* difference between inverses of MM1, MM0 */
    for (j = 0; j < chars; j++)
      FF[i][j] = MM1[i][j]-MM0[i][j];
} /* inversesandlogdets */


void lrttest(void)
{
  /* carries out the Likelihood Ratio Test using the null hypothesis
     covariances and the alternative hypothesis covariances */
  double loglr;
  long i, j, df;

  loglr = -0.5*(spp-1)*(logdetMM1-logdetMM0) + log(sumlrt) - log((double)sumsteps);
  fprintf(outfile, "\n Likelihood Ratio Test:\n\n");
  fprintf(outfile,
          "   Testing whether these sets of characters evolve independently\n\n");
  fprintf(outfile, "    Set          consists of characters:\n");
  fprintf(outfile, "   -----         -----------------------\n");
  fprintf(outfile, "     1           ");
  for (i = 0; i < chars; i++)
    if (weightarray[i] == 0)
      fprintf(outfile, " %3ld", i+1);
  fprintf(outfile, "\n");
  fprintf(outfile, "     2           ");
  for (i = 0; i < chars; i++)
    if (weightarray[i] == 1)
      fprintf(outfile, " %3ld", i+1);
  fprintf(outfile, "\n\n   2 times log Likelihood Ratio = %lf", 2.0*loglr);
  df = 0;
  for (i = 0; i < chars-1; i++)
    for (j = i; j < chars; j++)
      if (weightarray[i] != weightarray[j])
        df += 1;
  fprintf(outfile, ",  df = %ld\n\n", df);
} /* LRT test of whether two sets of variables are independent */


void freethings (void)
{
  /* free memory set aside for character data, liabilities */
  long i;

  if (discrete) {
    for (i = 0; i < spp; i++)
      free(yy[i]);
    free(yy);
  }
  freeview(&curtree, 2*spp-1);    /* free the nodes on the tree */
  freetree(&curtree.nodep, 2*spp-1);
} /* freethings */


void domcmc (void)
{
  /* carry out the MCMC including testing runs if needed for one tree
     and one data set */
  long i, blanks = 0;                   // RSGnote: Formerly "blanks" may have been accessed before being initialized.
  boolean done;

  nullchains = false;                      /* estimating covariances */
  if (progress)
  {
    if (!nullchains)
      printf("\nMarkov chain Monte Carlo (MCMC) estimation of covariances:\n\n");
    else
      printf("\nMCMC estimation of constrained covariances:\n\n");
    blanks = 0;
    if (testing)
      blanks = (long)(log((double)fsteps/steps)/log(10.0));
  }
  do {
    if (progress)
    {
      printf("                                                 Acceptance    Norm of change\n");
      if (cycles < 10)
        printf("Chains (%1ld) ", cycles);
      else
        printf("Chains (%2ld)", cycles);
      printf("                                        rate        in transform\n");
      printf("------");
      for (i = 0; i < blanks; i++)
        printf(" ");
      printf("                                           ----------    --------------\n");
      fflush(stdout);
      printf("\nBurn-in: %ld updates\n", burnin);
    }
    inittransforms();                       /* set up initial liabilities */
    if (kth <= 1)
      sumnumerator = (double *)Malloc(chars * sizeof(double));
    initliabilities(root);
    for (i = 1; i <= cycles; i++)            /* for each chain to be run ... */
    {
      lasttime = (i == cycles);          /* is this the last one? */
      nexttolasttime = (i == (cycles-1));   /* ... the next-to-last one? */
      if (progress)
      {
        if (testing && lasttime)
          printf("Chain %ld: %ld updates ", i, fsteps);
        else
          printf("Chain %ld: %ld updates ", i, steps);
        if (i < 10)
          printf(" ");
      }
      sumsteps = 0;
      if (testing && lasttime)
        iterateliabilities(fsteps, true);     /* ... run the testing chain */
      else
        iterateliabilities(steps, true);      /* ... run the chain */
      correcttransform(nullchains); /* ... adjust transform to independence */
      if (testing && nexttolasttime && nullchains)
        inversesandlogdets();
    }
    if (!nullchains)
      done = true;
    else
      done = testing;       /* set whether doing a null-hypotheis run */
    reportall(nullchains);                           /* report results */
    nullchains = testing;
  } while (!done);
  if (testing)   /* after the last chain ... */
    lrttest(); /* infer the null hypothesis likelihood and do the LRT test */
} /* domcmc */


int main(int argc, Char *argv[])
{
  /* main program */
  long ncases, datasper, treesper;
  boolean datafirst, treesfirst;
  initdata *funcs;

#ifdef MAC
  argc = 1;
  argv[0] = "Threshml";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = cont_node_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree file", "rb", argv[0], intreename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  firsttree = true;
  doinit();                                 /* get settings from menu */
  if (testing) {  /* get indicators of independent characters if doing LRT */
    openfile(&weightfile, WEIGHTFILE, "weights file", "rb", argv[0],
             weightfilename);
    inputweights(chars, weightarray, &testing);
  }
  ncases = numtrees;    /* ncases will be how many tree-data pairs are done */
  if (ndatas > numtrees)
    ncases = ndatas;
  if (cross)
    ncases = ndatas * numtrees;
  datafirst = muldata && !treeswithin;    /* which one always gets advanced */
  treesfirst = multrees && !datawithin;
  treesper = ncases / ndatas;       /* how many trees per data set */
  datasper = ncases / numtrees;     /* how many data sets per tree */
  ith = 0;                              /* which data set was just done */
  jth = 0;                              /* which tree was just done */
  for (kth = 1; kth <= ncases; kth++) {
    ithwas = ith;
    jthwas = jth;
    if (datafirst || ((jth%treesper) == 0)) {
      ith++;
      readdata ();                      /* read a data set */
      if (kth == 1)
        allocthings();
    }
    if (treesfirst || ((ithwas%datasper) == 0)) {
      jth++;
      readthetree ();                   /* read a tree ... */
    }
    if (treesfirst && (ith > ithwas)) {
      if (ndatas > 1) {
        fprintf(outfile, "\nData set # %ld:\n\n", ith);
        printf("\n\nData set # %ld:\n", ith);
      }
    }
    if ((numtrees > 1) && (jth > jthwas)) {
      fprintf(outfile, "\nTree # %ld:\n\n", jth);
      printf("\n\nTree # %ld:\n", jth);
    }
    if (!treesfirst) {
      if (ndatas > 1) {
        fprintf(outfile, "\nData set # %ld:\n\n", ith);
        printf("\n\nData set # %ld:\n", ith);
      }
    }
    domcmc();           /* do an MCMC run for one tree and one data set */
    fflush(outfile);
    if ((cross && datawithin) && (ith == ndatas) && (kth < ncases)) {
      FClose(infile);
      openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
      ith = 0;
    }
    if ((cross && treeswithin) && (jth == numtrees) && (kth < ncases)) {
      FClose(intree);
      openfile(&intree, INTREE, "input tree file", "r", argv[0], intreename);
      jth = 0;
    }
  }
  printf("\nOutput written on output file \"%s\".\n\n", outfilename);
  FClose(outfile);
  FClose(infile);
  FClose(intree);
  if (testing)
    FClose(weightfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  freethings();
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
} /* threshml */


// End.
