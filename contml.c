/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cont.h"
#include "ml.h"

#define epsilon1        0.000001   /* small number */
#define epsilon2        0.02   /* not such a small number */

#define over            60

typedef struct contml_node {
  cont_node_type cont_node_var;
  double bigv;
  double dist;
} contml_node;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   getalleles(void);
void   inputdata(void);
void   transformgfs(void);
void   getinput(void);
void   contml_hookup(node *, node *);
void   sumlikely(node *, node *, double *);
double contml_tree_evaluate(tree *, node *, boolean);
double distance(node *, node *);
void   makedists(node *);
void   makebigv(contml_node *, boolean *);
void   correctv(node *);
void   littlev(node *);
void   contml_tree_makenewv(tree *, node *);
void   contml_tree_nuview(tree*, node *);
void   inittip(tree*,  long);
void   coordinates(node *, double, long *, double *);
void   drawline(long, double);
void   printree(void);
void   treeout(node *);
void   describe(node *, double, double);
void   summarize(void);
void   nodeinit(node *);
void   initrav(node *);
void   treevaluate(void);
void   maketree(void);
node*  contml_node_new(node_type, long);
void   contml_node_init(node*, node_type, long);
void   contml_node_reinit(node*);
void   contml_node_copy(node *, node *);
void   inittrees(void);
tree*  contml_tree_new(long, long);
void   contml_tree_init(tree*, long, long);
void   contml_tree_free(tree*);
void   contmlrun(void);
void   contml(char * infilename, char * intreename, char * OutfileName, char * outfileopt,
              char * OuttreeName, char * outtreeopt, int BestTree, int UseLengths, int GeneFreq,
              int AllAlleles, int OutRoot, int OutNum, int GlobalRearr, int RandInput,
              int RandNum, int Njumble, int MultData, int NumSets, int PrintData,
              int PrintInd, int PrintTree, int WriteTree);
/* function prototypes */
#endif

long rcategs;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH];
long nonodes2, loci, totalleles, df, outgrno, col, datasets, ith, njumble, jumb = 0;
long inseed, inseed0;
long *alleles, *locus, *weight;
phenotype3 *x;
boolean all, contchars, global, jumble, lngths, outgropt, trout, usertree, printdata, progress, treeprint, mulsets, firstset;
boolean smoothit = true;
longer seed;
long *enterorder;
tree *curtree, *priortree, *bestree, *bestree2;
long nextsp, numtrees, which, maxwhich, shimotrees;
/* From maketree, propagated to global */
boolean succeeded;
double maxlogl;
double l0gl[MAXSHIMOTREES];
double *pbar, *sqrtp, *l0gf[MAXSHIMOTREES];
Char ch;
char *progname;
double trweight;   /* added to make treeread happy */
boolean goteof;
boolean haslengths;   /* end of ones added to make treeread happy */
node *addwhere;

/* added to make ml.o happy, could use these later */
boolean smoothed = false;
boolean polishing = false;


tree * contml_tree_new(long nonodes, long spp)
{ /* create a contml_tree */
  tree* t = Malloc(sizeof(ml_tree));
  contml_tree_init(t, nonodes, spp);
  return t;
} /* contml_tree_new*/


void contml_tree_init(tree* t, long nonodes, long spp)
{ /* initialize a contml_tree */

  ml_tree_init(t, nonodes, spp);
  allocview(t, nonodes2, totalleles);
  t->evaluate = contml_tree_evaluate;
  t->nuview = contml_tree_nuview;
  t->makenewv = contml_tree_makenewv;
  t->free = contml_tree_free;
} /* contml_tree_init */

 
void contml_tree_free(tree* t)
{ /* free a contml_node */
  freeview(t, nonodes2);
  generic_tree_free(t);
} /* contml_tree_free */


void inittrees(void)
{ /* initialize the trees that will be needed */
  curtree = contml_tree_new(nonodes2, spp);
  if (usertree)
    return;
  bestree = contml_tree_new(nonodes2, spp);
  priortree = contml_tree_new(nonodes2, spp);
  if (njumble <= 1)
    return;
  bestree2 = contml_tree_new(nonodes2, spp);
} /* inittrees */


node * contml_node_new(node_type type, long index) // RSGbugfix
{ /* create a contml_node */
  node* n = Malloc(sizeof(contml_node));

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  contml_node_init(n, type, index);

  return n;

} /* contml_node_new */


void contml_node_init(node* n, node_type type, long index)
{ /* initialize a contml_node */
  contml_node *cn = (contml_node *)n;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  generic_node_init(&(cn->cont_node_var.node_var), type, index);

  cn->dist = 0;
  cn->cont_node_var.node_var.copy = contml_node_copy;
  cn->cont_node_var.node_var.reinit = contml_node_reinit;

} /* contml_node_init */


void contml_node_reinit(node* n)
{ /* re-init a contml_node */
  contml_node *cn = (contml_node *)n;
  generic_node_reinit(&(cn->cont_node_var.node_var));
  cn->dist = 0;
} /* contml_node_reinit */


void getoptions(void)
{ /* interactively set options */
  long loopcount;
  Char ch;
  boolean done;

  putchar('\n');
  global = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  outgrno = 1;
  outgropt = false;
  all = false;
  contchars = false;
  trout = true;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do
  {
    cleerhome();
    printf("\nContinuous character Maximum Likelihood");
    printf(" method version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  U                       Search for best tree?  %s\n", (usertree ? "No, use user trees in input" : "Yes"));
    if (usertree)
    {
      printf("  L                Use lengths from user trees?%s\n", (lngths ? "  Yes" : "  No"));
    }
    printf("  C  Gene frequencies or continuous characters?  %s\n", (contchars ? "Continuous characters" : "Gene frequencies"));
    if (!contchars)
      printf("  A   Input file has all alleles at each locus?  %s\n", (all ? "Yes" : "No, one allele missing at each"));
    printf("  O                              Outgroup root?  %s %ld\n", (outgropt ? "Yes, at species number" : "No, use as outgroup species"), outgrno);
    if (!usertree)
    {
      printf("  G                      Global rearrangements?  %s\n", (global ? "Yes" : "No"));
      printf("  J           Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed=%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M                 Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0         Terminal type (IBM PC, ANSI, none)?  %s\n", ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1          Print out the data at start of run  %s\n", (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n", (progress ? "Yes" : "No"));
    printf("  3                              Print out tree  %s\n", (treeprint ? "Yes" : "No"));
    printf("  4             Write out trees onto tree file?  %s\n", (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done)
    {
      if (((!usertree) && (strchr("JOUGACM12340", ch) != NULL)) || (usertree && ((strchr("LOUACM12340", ch) != NULL))))
      {
        switch (ch)
        {
          case 'A':
            if (!contchars)
              all = !all;
            break;

          case 'C':
            contchars = !contchars;
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
            usertree = !usertree;
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
      }
      else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /* getoptions */


void allocrest(void)
{  /* allocate arrays for number of alleles, the data coordinates, names etc */
  alleles = (long *)Malloc(loci * sizeof(long));
  if (contchars)
    locus = (long *)Malloc(loci * sizeof(long));
  x = (phenotype3 *)Malloc(spp * sizeof(phenotype3));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
}  /* allocrest */


void doinit(void)
{ /* initializes variables */
  inputnumbers(&spp, &loci, &nonodes2, 1);
  fprintf(outfile, "\nContinuous character Maximum Likelihood");
  fprintf(outfile, " method version %s\n\n", VERSION);

  if (!javarun)
  {
    getoptions();
  }

  if(!usertree)
    nonodes2--;
  if (printdata)
    fprintf(outfile, "\n%4ld Populations, %4ld Loci\n", spp, loci);
  allocrest();
}  /* doinit */


void getalleles(void)
{ /* set up number of alleles at loci */
  long i, j, m;

  if (!firstset)
    samenumsp(&loci, ith);
  if (contchars )
  {
    totalleles = loci;
    for (i = 1; i <= loci; i++)
    {
      locus[i - 1] = i;
      alleles[i - 1] = 1;
    }
    df = loci;
  }
  else
  {
    totalleles = 0;
    scan_eoln(infile);
    if (printdata)
    {
      fprintf(outfile, "\nNumbers of alleles at the loci:\n");
      fprintf(outfile, "------- -- ------- -- --- -----\n\n");
    }
    for (i = 1; i <= loci; i++)
    {
      if (eoln(infile))
        scan_eoln(infile);
      if (fscanf(infile, "%ld", &alleles[i - 1]) < 1)
      {
        printf("\nERROR:  Unable to read number of alleles at locus %ld.\n", i);
        exxit(-1);
      }
      if (alleles[i - 1] <= 0)
      {
        printf("\nERROR:  Bad number of alleles: %ld at locus %ld.\n", alleles[i-1], i);
        exxit(-1);
      }
      totalleles += alleles[i - 1];
      if (printdata)
        fprintf(outfile, "%4ld", alleles[i - 1]);
    }
    locus = (long *)Malloc(totalleles * sizeof(long));
    m = 0;
    for (i = 1; i <= loci; i++)
    {
      for (j = 0; j < alleles[i - 1]; j++)
        locus[m+j] = i;
      m += alleles[i - 1];
    }
    df = totalleles - loci;
  }
  for (i = 0; i < spp; i++)
    x[i] = (phenotype3)Malloc(totalleles * sizeof(double));
  pbar = (double *)Malloc(totalleles * sizeof(double));
  if (usertree)
    for (i = 0; i < MAXSHIMOTREES; i++)
      l0gf[i] = (double *)Malloc(totalleles * sizeof(double));
  if (printdata)
    putc('\n', outfile);
}  /* getalleles */


void inputdata(void)
{ /* read species data */
  long i, j, k, l, m, n, p;
  double sum;

  if (printdata)
  {
    fprintf(outfile, "\nName");
    if (contchars)
      fprintf(outfile, "                       Phenotypes\n");
    else
      fprintf(outfile, "                 Gene Frequencies\n");
    fprintf(outfile, "----");
    if (contchars)
      fprintf(outfile, "                       ----------\n");
    else
      fprintf(outfile, "                 ---- -----------\n");
    putc('\n', outfile);
    if (!contchars)
    {
      for (j = 1; j <= nmlngth - 8; j++)
        putc(' ', outfile);
      fprintf(outfile, "locus:");
      p = 1;
      for (j = 1; j <= loci; j++)
      {
        if (all)
          n = alleles[j - 1];
        else
          n = alleles[j - 1] - 1;
        for (k = 1; k <= n; k++)
        {
          fprintf(outfile, "%10ld", j);
          if (p % 6 == 0 && (all || p < df))
          {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth - 2; l++)
              putc(' ', outfile);
          }
          p++;
        }
      }
      fprintf(outfile, "\n\n");
    }
  }
  for (i = 0; i < spp; i++)
  {
    scan_eoln(infile);
    initname(i);
    if (printdata)
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
    m = 1;
    p = 1;
    for (j = 1; j <= loci; j++)
    {
      sum = 0.0;
      if (contchars)
        n = 1;
      else if (all)
        n = alleles[j - 1];
      else
        n = alleles[j - 1] - 1;
      for (k = 1; k <= n; k++)
      {
        if (eoln(infile))
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][m - 1]) < 1)
        {
          printf("\nERROR:  Unable to read allele frequency for species %ld, locus %ld.\n", i+1, j);
          exxit(-1);
        }
        sum += x[i][m - 1];
        if (!contchars && x[i][m - 1] < 0.0)
        {
          printf("\nERROR:  Locus %ld in species %ld: an allele", j, i+1);
          printf(" frequency is negative.\n");
          exxit(-1);
        }
        if (printdata)
        {
          fprintf(outfile, "%10.5f", x[i][m - 1]);
          if (p % 6 == 0 && (all || p < df))
          {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth; l++)
              putc(' ', outfile);
          }
        }
        p++;
        m++;
      }
      if (all && !contchars)
      {
        if (fabs(sum - 1.0) > epsilon2)
        {
          printf("\nERROR:  Locus %ld in species %ld: frequencies do not add up to 1.\n", j, i + 1);
          printf("\nFrequencies are:\n");
          for (l = 0; l <= m-3; l++)
            printf("%f+", x[i][l]);
          printf("%f = %f\n\n", x[i][m-2], sum);
          exxit(-1);
        }
        else
        {
          for (l = 0; l <= m-2; l++)
            x[i][l] /= sum;
        }
      }
      if (!all && !contchars)
      {
        x[i][m-1] = 1.0 - sum;
        if (x[i][m-1] < 0.0)
        {
          if (x[i][m-1] > -epsilon2)
          {
            for (l = 0; l <= m-2; l++)
              x[i][l] /= sum;
            x[i][m-1] = 0.0;
          }
          else
          {
            printf("\nERROR:  Locus %ld in species %ld: ", j, i + 1);
            printf("frequencies add up to more than 1.\n");
            printf("\nFrequencies are:\n");
            for (l = 0; l <= m-3; l++)
              printf("%f+", x[i][l]);
            printf("%f = %f\n\n", x[i][m-2], sum);
            exxit(-1);
          }
        }
        m++;
      }
    }
    if (printdata)
      putc('\n', outfile);
  }
  scan_eoln(infile);
  if (printdata)
    putc('\n', outfile);
  checknames(spp);                      // Check NAYME array for duplicates.
}  /* inputdata */


void transformgfs(void)
{ /* do stereographic projection transformation on gene frequencies to
     get variables that come closer to independent Brownian motions */
  long i, j, k, l, m, n, maxalleles;
  double f, sum;
  double *sumprod, *sqrtp, *pbar;
  phenotype3 *c;

  sumprod = (double *)Malloc(loci * sizeof(double));
  sqrtp = (double *)Malloc(totalleles * sizeof(double));
  pbar = (double *)Malloc(totalleles * sizeof(double));
  for (i = 0; i < totalleles; i++)      /* get mean gene frequencies */
  {
    pbar[i] = 0.0;
    for (j = 0; j < spp; j++)
      pbar[i] += x[j][i];
    pbar[i] /= spp;
    if (pbar[i] == 0.0)
      sqrtp[i] = 0.0;
    else
      sqrtp[i] = sqrt(pbar[i]);
  }
  for (i = 0; i < spp; i++)
  {
    for (j = 0; j < loci; j++)          /* for each locus, sum of root(p*x) */
      sumprod[j] = 0.0;
    for (j = 0; j < totalleles; j++)
      if ((pbar[j]*x[i][j]) >= 0.0)
        sumprod[locus[j]-1] += sqrtp[j]*sqrt(x[i][j]);
    for (j = 0; j < totalleles; j++)    /* the projection to tangent plane */
    {
      f = (1.0 + sumprod[locus[j]-1])/2.0;
      if (x[i][j] == 0.0)
        x[i][j] = (2.0/f - 1.0)*sqrtp[j];
      else
        x[i][j] = (1.0/f)*sqrt(x[i][j]) + (1.0/f - 1.0)*sqrtp[j];
    }
  }
  maxalleles = 0;
  for (i = 0; i < loci; i++)
    if (alleles[i] > maxalleles)
      maxalleles = alleles[i];
  c = (phenotype3 *)Malloc(maxalleles * sizeof(phenotype3));
  for (i = 0; i < maxalleles; i++)  /* enough room for any locus's contrasts */
    c[i] = (double *)Malloc(maxalleles * sizeof(double));
  m = 0;
  for (j = 0; j < loci; j++)            /* do this for each locus */
  {
    for (k = 0; k < alleles[j]-1; k++)  /* one fewer than # of alleles */
    {
      c[k][0] = 1.0;
      for (l = 0; l < k; l++)          /* for contrasts 1 to k make it ... */
      {
        sum = 0.0;
        for (n = 0; n <= l; n++)
          sum += c[k][n]*c[l][n];
        if (fabs(c[l][l+1]) > 0.000000001) /* ... orthogonal to those ones */
          c[k][l+1] = -sum / c[l][l+1];    /* set coeff to make orthogonal */
        else
          c[k][l+1] = 1.0;
      }
      sum = 0.0;
      for (l = 0; l <= k; l++) /* make it orthogonal to vector of sqrtp's */
        sum += c[k][l]*sqrtp[m+l];
      if (sqrtp[m+k+1] > 0.0000000001)
        c[k][k+1] = - sum / sqrtp[m+k+1];  /* ... setting last coeff */
      else
      {
        for (l = 0; l <= k; l++)
          c[k][l] = 0.0;
        c[k][k+1] = 1.0;
      }
      sum = 0.0;
      for (l = 0; l <= k+1; l++)
        sum += c[k][l]*c[k][l];
      sum = sqrt(sum);
      for (l = 0; l <= k+1; l++)
        if (sum > 0.0000000001)
          c[k][l] /= sum;
    }
    for (i = 0; i < spp; i++)     /* the orthonormal axes in the plane */
    {
      for (l = 0; l < alleles[j]-1; l++) /* compute the l-th one */
      {
        c[maxalleles-1][l] = 0.0;  /* temporarily store it ... */
        for (n = 0; n <= l+1; n++)
          c[maxalleles-1][l] += c[l][n]*x[i][m+n];
      }
      for (l = 0; l < alleles[j]-1; l++)
        x[i][m+l] = c[maxalleles-1][l];   /* replace the gene freqs by it */
    }
    m += alleles[j];
  }
  for (i = 0; i < maxalleles; i++)
    free(c[i]);
  free(c);
  free(sumprod);
  free(sqrtp);
  free(pbar);
} /* transformgfs */


void getinput(void)
{ /* reads the input data */
  getalleles();
  inittrees();
  inputdata();
  if (!contchars)
  {
    transformgfs();
  }
}  /* getinput */


void contml_hookup(node* p, node* q){
/* hook up two nodes, set branch length to initial value
   (one of the nodes may be in a fork circle) */

  hookup(p, q);
  p->v = initialv;
  q->v = initialv;
} /* contml_hookup */


void sumlikely(node *p, node *q, double *sum)
{ /* sum contribution to likelihood over forks in tree */
  long i, j, m;
  double term, sumsq, vee;
  double temp;

  if (!p->tip)
    sumlikely(p->next->back, p->next->next->back, sum);
  if (!q->tip)
    sumlikely(q->next->back, q->next->next->back, sum);
  if (p->back == q)
    vee = p->v;
  else
    vee = p->v + q->v;
  vee += p->deltav + q->deltav;
/* debug: this seems to give trouble so commenting it out ...
  if (vee <= 1.0e-10)
  {
    printf("\nERROR:  Check for two identical species and eliminate one from the data.\n");
    exxit(-1);
  }
debug:  */
  sumsq = 0.0;
  if (usertree && which <= MAXSHIMOTREES)
  {
    for (i = 0; i < loci; i++)
      l0gf[which - 1][i] += (1 - alleles[i]) * log(vee) / 2.0;
  }
  if (contchars)
  {
    m = 0;
    for (i = 0; i < loci; i++)
    {
      temp = ((cont_node_type*)p)->view[i] - ((cont_node_type*)q)->view[i];
      term = temp * temp;
      if (usertree && which <= MAXSHIMOTREES)
        l0gf[which - 1][i] -= term / (2.0 * vee);
      sumsq += term;
    }
  }
  else
  {
    m = 0;
    for (i = 0; i < loci; i++)
    {
      for (j = 1; j < alleles[i]; j++)
      {
        temp = ((cont_node_type*)p)->view[m+j-1] - ((cont_node_type*)q)->view[m+j-1];
        term = temp * temp;
        if (usertree && which <= MAXSHIMOTREES)
          l0gf[which - 1][i] -= term / (2.0 * vee);
        sumsq += term;
      }
      m += alleles[i];
    }
  }
  (*sum) += df * log(vee) / -2.0 - sumsq / (2.0 * vee);
}  /* sumlikely */


double contml_tree_evaluate(tree *t, node *p, boolean saveit)
{ /* evaluate likelihood of a tree */
  /* debug: (is saveit needed?) */
  long i;
  double sum;

  generic_tree_evaluate (t, p, saveit);    /* update views if needed */
  sum = 0.0;
  if (usertree && which <= MAXSHIMOTREES)
  {
    for (i = 0; i < loci; i++)
      l0gf[which - 1][i] = 0.0;
  }

  sumlikely(p, p->back, &sum);    /* this gets the likelihood, recursively */
  if (usertree && which <= MAXSHIMOTREES)
  {
    l0gl[which - 1] = sum;
    if (which == 1)
    {
      maxwhich = 1;
      maxlogl = sum;
    }
    else if (sum > maxlogl)
    {
      maxwhich = which;
      maxlogl = sum;
    }
  }
  t->score = sum;

  return sum;
}  /* contml_tree_evaluate */


double distance(node *p, node *q)
{ /* distance-squared between two nodes */
  long i, j, m;
  double sum, temp;

  sum = 0.0;
  if (!contchars)
  {
    m = 0;
    for (i = 0; i < loci; i++)
    {
      for (j = 0; j < alleles[i]-1; j++)
      {
        temp = ((cont_node_type*)p)->view[m+j] - ((cont_node_type*)q)->view[m+j];
        sum += temp * temp;
      }
      m += alleles[i];
    }
  }
  else
  {
    for (i = 0; i < totalleles; i++)
    {
      temp = ((cont_node_type*)p)->view[i] - ((cont_node_type*)q)->view[i];
      sum += temp * temp;
    }
  }
  return sum;
}  /* distance */


void makedists(node *p)
{ /* compute distances among three neighbors of a node */
  long i;
  node *q;

  if (p->tip)
    p = p->back;
  for (i = 1; i <= 3; i++)
  {
    q = p->next;
    ((contml_node*)p)->dist = distance(p->back, q->back);
    p = q;
  }
}  /* makedists */


void makebigv(contml_node *p, boolean *negatives)
{ /* make new branch lengths iteratively at a fork */
  long i;
  contml_node *temp, *q, *r;

  q = (contml_node*)((node*)p)->next;
  r = (contml_node*)((node*)q)->next;
  *negatives = false;
  for (i = 1; i <= 3; i++)
  {
    p->bigv = ((node*)p)->v + ((node*)p)->back->deltav;
    if (((node*)p)->iter)
    {
      p->bigv = (p->dist + r->dist - q->dist) / (df * 2);
      ((contml_node*)((node*)p)->back)->bigv = p->bigv;
      if (p->bigv < ((node*)p)->back->deltav)
        *negatives = true;
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */


void correctv(node *p)
{ /* iterate branch lengths if some are to be zero */
  node *q, *r, *temp;
  long i, j;
  double f1, f2, vtot;

  q = p->next;
  r = q->next;
  for (i = 1; i <= smoothings; i++)
  {
    for (j = 1; j <= 3; j++)
    {
      vtot = ((contml_node*)q)->bigv + ((contml_node*)r)->bigv;
      if (vtot > 0.0)
        f1 = ((contml_node*)q)->bigv / vtot;
      else
        f1 = 0.5;
      f2 = 1.0 - f1;
      ((contml_node*)p)->bigv = (f1 * ((contml_node*)r)->dist + f2 * ((contml_node*)p)->dist - f1 * f2 * ((contml_node*)q)->dist) / df;
      ((contml_node*)p)->bigv -= vtot * f1 * f2;
      if (((contml_node*)p)->bigv < p->back->deltav)
        ((contml_node*)p)->bigv = p->back->deltav;
      ((contml_node*)p->back)->bigv = ((contml_node*)p)->bigv;
      temp = p;
      p = q;
      q = r;
      r = temp;
    }
  }
}  /* correctv */


void littlev(node *p)
{ /* remove part of it that belongs to other branches 
   * This is a version that works only for bifurcating trees */
  long i;

  for (i = 1; i <= 3; i++)
  {
    if (p->iter)
      p->v = ((contml_node*)p)->bigv - p->back->deltav;
    if (p->back->iter)
      p->back->v = p->v;
    p = p->next;
  }
}  /* littlev */


void contml_tree_nuview(tree* t, node *p)
{ /* renew inward-looking view information in subtrees */
  long j, k, m;
  node *q, *r, *a, *b;
  double v1, v2, vtot, f1, f2;

  q = p->next;
  r = q->next;
  a = q->back;
  b = r->back;
  v1 = q->v + a->deltav;     /*  this is now   v1'   */
  v2 = r->v + b->deltav;     /*  this is now   v2'   */
  vtot = v1 + v2;            /*  this is now  v1' + v2'   */
  if (vtot > 0.0)
    f1 = v2 / vtot;
  else
    f1 = 0.5;
  f2 = 1.0 - f1;
  m = 0;
  for (j = 0; j < loci; j++)
  {
    for (k = 1; k <= alleles[j]; k++)
      ((cont_node_type*)p)->view[m+k-1] = f1*((cont_node_type*)a)->view[m+k-1]
                                     + f2 * ((cont_node_type*)b)->view[m+k-1];
    m += alleles[j];
  }
  p->deltav = v1 * f1;   /* so it is     v1' v2' / (v1' + v2')     */

}  /* contml_tree_nuview */


void contml_tree_makenewv(tree* t, node* p) {
/* Compute new branch length.  If after subtracting p->deltav it is negative,
 * then compute new branch lengths on the three branches connected to p
 * and do this iteratively, setting some to zero as needed. */
/* debug:    boolean negatives;    commented out maybe not needed  */

  p->v = distance(p, p->back);
  p->v = p->v - p->deltav - p->back->deltav;
  p->back->v = p->v;
  if (p->v < 0.0) {
    p->v = 0.0;    /* nearest legal value.  smoothing adjusts others */
    p->back->v = 0.0;
/* debug:    if (p->tip)
      p = p->back;
    makedists(p);     debug: probably need to do a loop around circle
    makebigv((contml_node*)p, &negatives);
    if (negatives)
      correctv(p);
    littlev(p);   debug */
  }
} /* contml_tree_makenewv */


void contml_node_copy(node *src, node *dst)
{ /* make a copy of a node */
  contml_node *c = (contml_node *)src;
  contml_node *d = (contml_node *)dst;

  cont_node_copy((node*)c, (node*)d);
  d->bigv = c->bigv;
  d->dist = c->dist;
}  /* contml_node_copy */


void inittip(tree* t, long m)
{ /* initialize and hook up a new tip;  m is the index of the tip */
  node *tmp;

  tmp = t->nodep[m - 1];
  memcpy(((cont_node_type*)tmp)->view, x[m - 1], totalleles * sizeof(double));
  tmp->deltav = 0.0;
}  /* inittip */


void coordinates(node *p, double lengthsum, long *tipy, double *tipmax)
{ /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip)
  {
    p->xcoord = lengthsum;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do
  {
    coordinates(q->back, lengthsum + q->v, tipy, tipmax);
    q = q->next;
  } while ((p == curtree->root || p != q) && (p != curtree->root || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = lengthsum;
  if (p == curtree->root)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(long i, double scale)
{ /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first = NULL, *last = NULL;
  boolean done;

  p = curtree->root;
  q = curtree->root;
  assert(p->index > 0);                 // RSGdebug
  extra = false;
  if (i == (long)p->ycoord && p == curtree->root)
  {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
  }
  else
    fprintf(outfile, "  ");
  do
  {
    if (!p->tip)
    {
      r = p->next;
      done = false;
      do
      {
        if (i >= (long)r->back->ymin && i <= (long)r->back->ymax)
        {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (p != curtree->root && r == p) || (p == curtree->root && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree->root)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done)
    {
      if ((long)p->ycoord != (long)q->ycoord)
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
      if ((long)last->ycoord > i && (long)first->ycoord < i && i != (long)p->ycoord)
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
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip)
  {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(void)
{ /* prints out diagram of the tree */
  long i;
  long tipy;
  double tipmax, scale;

  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree->root, 0.0, &tipy, &tipmax);
  scale = over / (tipmax + 0.0001);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale);
  putc('\n', outfile);
}  /* printree */


void treeout(node *p)
{ /* write out file with representation of final tree */
  long i, n, w;
  Char c;
  double x;

  assert(p->index > 0);                 // RSGdebug

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
    putc('(', outtree);
    col++;
    treeout(p->next->back);
    putc(',', outtree);
    col++;
    if (col > 55)
    {
      putc('\n', outtree);
      col = 0;
    }
    treeout(p->next->next->back);
    if (p == curtree->root)
    {
      putc(',', outtree);
      col++;
      if (col > 45)
      {
        putc('\n', outtree);
        col = 0;
      }
      treeout(p->back);
    }
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
    fprintf(outtree, ":%*.8f", (int)w + 7, x);
    col += w + 8;
  }
}  /* treeout */


void describe(node *p, double chilow, double chihigh)
{ /* print out information for one branch */
  long i;
  node *q;
  double bigv, delta;

  q = p->back;

  assert(p->index > 0);                 // RSGdebug
  assert(q->index > 0);                 // RSGdebug

  fprintf(outfile, "%3ld       ", q->index - spp);
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.8f", q->v);
  delta = p->deltav + p->back->deltav;
  bigv = p->v + delta;
  if (p->iter)
    fprintf(outfile, "   (%12.8f,%12.8f)", chilow * bigv - delta, chihigh * bigv - delta);
  fprintf(outfile, "\n");
  if (!p->tip)
  {
    describe(p->next->back, chilow, chihigh);
    describe(p->next->next->back, chilow, chihigh);
  }
}  /* describe */


void summarize(void)
{ /* print out branch lengths etc. */
  double chilow, chihigh;

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree->score);
  if (df == 1)
  {
    chilow = 0.000982;
    chihigh = 5.02389;
  }
  else if (df == 2)
  {
    chilow = 0.05064;
    chihigh = 7.3777;
  }
  else
  {
    chilow = 1.0 - 2.0 / (df * 9);
    chihigh = chilow;
    chilow -= 1.95996 * sqrt(2.0 / (df * 9));
    chihigh += 1.95996 * sqrt(2.0 / (df * 9));
    chilow *= chilow * chilow;
    chihigh *= chihigh * chihigh;
  }
  fprintf(outfile, "\nBetween     And             Length");
  if ((usertree && !lngths) || !usertree) 
    fprintf(outfile, "      Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, "-------     ---             ------");
  if ((usertree && !lngths) || !usertree) 
    fprintf(outfile, "      ------- ---------- ------");
  fprintf(outfile, "\n");
  describe(curtree->root->next->back, chilow, chihigh);
  describe(curtree->root->next->next->back, chilow, chihigh);
  describe(curtree->root->back, chilow, chihigh);
  fprintf(outfile, "\n\n");
  if (trout)
  {
    col = 0;
    treeout(curtree->root);
  }
}  /* summarize */


void nodeinit(node *p)
{ /* initialize a node -- for interior forks, use rough averages
   * debug: why do we need this?  Why not just use nodep? */
  node *q, *r;
  long i, j, m;

  if (p->tip)
    return;
  q = p->next->back;
  r = p->next->next->back;
  nodeinit(q);
  nodeinit(r);
  m = 0;
  for (i = 0; i < loci; i++)
  {
    for (j = 1; j < alleles[i]; j++)
      ((cont_node_type*)p)->view[m+j-1] = 0.5 * ((cont_node_type*)q)->view[m+j-1] + 0.5 * ((cont_node_type*)r)->view[m+j-1];
    m += alleles[i];
  }
  if ((!lngths) || p->iter)
    p->v = initialv;
  if ((!lngths) || p->back->iter)
    p->back->v = initialv;
}  /* nodeinit */


void initrav(node *p)
{ /* traverse to initialize */
  node *q;

  if (p->tip)
    nodeinit(p->back);
  else
  {
    q = p->next;
    while ( q != p )
    {
      initrav(q->back);
      q = q->next;
    }
    nodeinit(p);
  }
}  /* initrav */


void treevaluate(void)
{ /* evaluate user-defined tree, iterating branch lengths if needed */
  long i;
  double like;  /* to keep evaluate happy, not used */

  unroot(curtree, nonodes2);          /*  so root is at interior fork */
  inittravall (curtree, curtree->root);     /* set initializeds false */
  inittravall (curtree, curtree->root->back);
  curtree->donewbl = !lngths;
  if (!lngths) {        /* if no branch lengths, set them to initialv */
    ml_initialvtrav (curtree, curtree->root);
    ml_initialvtrav (curtree, curtree->root->back);
  }
  if ((!lngths) && curtree->donewbl) {
    for (i = 1; i <= smoothings * 4; i++)
      smooth(curtree, curtree->root);
  }
  else {
    inittravall(curtree, curtree->root);
    inittravall(curtree, curtree->root->back);
    ml_update(curtree, curtree->root);
  }
  like = curtree->evaluate(curtree, curtree->root, false);
}  /* treevaluate */


void maketree(void)
{ /* construct the tree */
  long i, k;
  node *p, *qwhere;
  double bestyet;

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
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      putc('\n', outfile);
    }
    for (which = 1; which <= spp; which++)
      inittip(curtree, which);
    which = 1;
    for (i = spp; i < nonodes2; i++) {
      p = curtree->get_fork(curtree, i);
      curtree->nodep[i] = p;
    }
    while (which <= numtrees)
    {
      treeread2 (intree, &curtree->root, curtree->nodep, lngths, &trweight,
                  &goteof, &haslengths, &spp, false, nonodes2);
      treevaluate();
      if (treeprint)
      {
        printree();
        summarize();
      }
      which++;
    }
    FClose(intree);
    if (numtrees > 1 && loci > 1 )
    {
      weight = (long *)Malloc(loci * sizeof(long));
      for (i = 0; i < loci; i++)
        weight[i] = 1;
      standev2(numtrees, maxwhich, 0, loci-1, maxlogl, l0gl, l0gf, seed);
      free(weight);
      fprintf(outfile, "\n\n");
    }
  }
  else                                  /* if we are constructing a tree */
  {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;

    if (jumble)
      randumize(seed, enterorder);

    nextsp = 3;

    // debug: RGS: destruct_tree() is ALWAYS called in all the other programs; doing same thing
    // here fixes SegFault bug due to something not getting initialized properly when using jumbling.
    destruct_tree(curtree);
    inittip(curtree, enterorder[1]);
    inittip(curtree, enterorder[2]);
    inittip(curtree, enterorder[3]);
    curtree->donewbl = true;
    buildsimpletree(curtree, enterorder);
    ml_initialvtrav (curtree, curtree->root);
    ml_initialvtrav (curtree, curtree->root->back);
    inittravall(curtree, curtree->root);
    inittravall(curtree, curtree->root->back);
    smooth(curtree, curtree->root);
    smooth(curtree, curtree->root->back);
    if (jumb == 1)
      numtrees = 1;
    nextsp = 4;

    if (progress)
    {
      sprintf(progbuf, "\nAdding species:\n");
      print_progress(progbuf);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }

    /* debug: make sure works properly when only 3 species */
    while (nextsp <= spp)    /* add rest of species to tree */
    {
      inittip(curtree, enterorder[nextsp - 1]);
      curtree->copy(curtree, priortree);
      bestree->score = UNDEFINED;
      bestyet = UNDEFINED;
      k = generic_tree_findemptyfork(curtree);
      p = curtree->get_fork(curtree, k);
      contml_hookup(curtree->nodep[enterorder[nextsp-1]-1],p);
      qwhere = NULL;
      curtree->addtraverse(curtree, p, curtree->root, true, qwhere,
                           &bestyet, bestree, true);
      bestree->copy(bestree, curtree);

      if (progress)
      {
        writename(nextsp - 1, 1, enterorder);
        phyFillScreenColor();
      }

      if ( global && nextsp == spp )
        curtree->globrearrange(curtree, progress, true);
      else
        curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1],
                              true, priortree, bestree);
      if (global && nextsp == spp)
        putc('\n', outfile);
      if (global && nextsp == spp && progress)
      {
        sprintf(progbuf, "\n");
        print_progress(progbuf);
      }
      if (njumble > 1)
      {
        if (jumb == 1 && nextsp == spp)
          bestree->copy(bestree, bestree2);
        else if (nextsp == spp)
        {
          if (bestree2->score < bestree->score)
            bestree->copy(bestree, bestree2);
        }
      }
      if (nextsp == spp && jumb == njumble)
      {
        if (njumble > 1) bestree2->copy(bestree2, curtree);
        curtree->root = curtree->nodep[outgrno - 1]->back;

        if (treeprint)
        {
          printree();
          summarize();
        }
      }
      nextsp++;
    }
  }

  if ( jumb < njumble)
    return;
  if (progress)
  {
    sprintf(progbuf, "\n\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree also written onto file \"%s\".\n\n", outtreename);
      print_progress(progbuf);
    }
  }
  freeview(curtree, nonodes2);
  if (!usertree)
  {
    freeview(bestree, nonodes2);
    freeview(priortree, nonodes2);
  }
  for (i = 0; i < spp; i++)
    free(x[i]);
  if (!contchars)
  {
    free(locus);
    free(pbar);
  }
}  /* maketree */


void contmlrun(void)
{
  /*
  // debug printout // JRMdebug
  printf("global: %i\n", global);
  printf("jumble: %i\n", jumble);
  printf("njumble: %li\n", njumble);
  printf("lngths: %i\n", lngths);
  printf("outgrno: %li\n", outgrno);
  printf("outgropt: %i\n", outgropt);
  printf("all: %i\n", all);
  printf("contchars: %i\n", contchars);
  printf("trout: %i\n", trout);
  printf("usertree: %i\n", usertree);
  printf("printdata: %i\n", printdata);
  printf("progress: %i\n", progress);
  printf("treeprint: %i\n", treeprint);
  printf("mulsets: %i\n", mulsets);
  printf("datasets: %li\n", datasets);
  */

  // do the work
  long i;
  for (ith = 1; ith <= datasets; ith++)
  {
    getinput();
    if (ith == 1)
    {
      firstset = false;
    }
    if (datasets > 1)
    {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      maketree();
    }
    fflush(outfile);
    fflush(outtree);
    if (usertree)
    {
      for (i = 0; i < MAXSHIMOTREES; i++)
      {
        free(l0gf[i]);
      }
    }
  }
} /* contmlrun */


void contml(
  char * infilename,
  char * intreename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  int BestTree,
  int UseLengths,
  int GeneFreq,
  int AllAlleles,
  int OutRoot,
  int OutNum,
  int GlobalRearr,
  int RandInput,
  int RandNum,
  int Njumble,
  int MultData,
  int NumSets,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree)
{
  initdata *funcs;
  //printf("Hello from ContML!\n"); // JRMdebug
  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Contml";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = contml_node_new;
  funcs->tree_new = contml_tree_new;
  phylipinit(argc, argv, funcs, true);
  progname = argv[0];

  /*
  //global = false;
  //jumble = false;
  //njumble = 1;
  //lngths = false;
  //outgrno = 1;
  //outgropt = false;
  //all = false;
  //contchars = false;
  //trout = true;
  //usertree = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  char * infile,
  char * intree,
  char * outfile,
  char * outfileopt,
  char * outtree,
  char * outtreeopt,
  //int BestTree,
  //int UseLengths,
  //int GeneFreq,
  //int AllAlleles,
  //int OutRoot,
  //int OutNum,
  //int GlobalRearr,
  //int RandInput,
  //int RandNum,
  //int Njumble,
  //int MultData,
  //int NumSets,
  //int PrintData,
  //int PrintInd,
  //int PrintTree,
  //int WriteTree)
  */

  if (AllAlleles != 0)
  {
    all = true;
  }
  else
  {
    all = false;
  }

  if (GeneFreq != 0)
  {
    contchars = false;
  }
  else
  {
    contchars = true;
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

  // everything translated, start the run
  infile = fopen(infilename, "r");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

  if (usertree)
  {
    intree = fopen(intreename, "r");
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

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  firstset = true;

  doinit();
  contmlrun();

  if (trout)
    FClose(outtree);
  if (usertree)
    FClose(intree);
  FClose(outfile);
  FClose(infile);
  //printf("Done.\n\n");
}


int main(int argc, Char *argv[])
{  /* main program */
  initdata *funcs;

#ifdef MAC
  argc = 1;                /* macsetup("Contml", "");                */
  argv[0] = "Contml";
#endif

  funcs = Malloc(sizeof(initdata));
  funcs->node_new = contml_node_new;
  funcs->tree_new = contml_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();

  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  contmlrun();

  FClose(outfile);
  FClose(outtree);
  FClose(infile);
  free(funcs);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();

  return 0;
}


// End.
