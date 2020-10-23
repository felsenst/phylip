/* Version 4.0. (c) Copyright 1993-2020
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "seq.h"
#include "dnaparsimony.h"

#define SAMPLES 1000

typedef boolean *boolptr;

typedef struct dnacomp_tree {
  dnapars_tree dnapars_tree;
} dnacomp_tree;

extern long nextree;                    /* parsimony.c */
extern long maxtrees;                   /* parsimony.c */

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   makeweights(void);
void   doinput(void);
void   mincomp(long );
double dnacomp_tree_evaluate(tree* t, node *, boolean);
void   localsavetree(void);
void   describe(void);
void   initboolnames(node *, boolean *);
void   maketree(void);
void   standev3(long, long, long, double, double *, long **, longer);
void   reallocchars(void);
tree*  dnacomp_tree_new(long nonodes, long spp);
void   dnacomp_tree_init(tree*, long nonodes, long spp);
void dnacomprun(void);
void dnacomp(char * infilename, char * intreename, char * OutfileName, char * outfileopt,
             char * weightsfilename, char * OuttreeName, char * outtreeopt, int searchbest,
             int TreeSave, int inputorder, int RandNum, int NumJumble, int OutRoot,
             int OutNum, int SitesWeighted, int AnalyzeMult, int MultDataSet, int NumMult,
             int InputSeq, int PrintData, int DotDiff, int PrintInd, int PrintTree,
             int PrintSteps, int PrintSeq, int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
node *p;
long chars, col, ith, njumble, jumb = 0, nonodes = 0, msets;
long inseed, inseed0;
boolean jumble, usertree, trout, weights, progress, stepbox, ancseq, firstset, mulsets, justwts;
steptr oldweight, necsteps;
tree *curtree, *bestree, *priortree;
long *enterorder;
Char basechar[32]="ACMGRSVTWYHKDBNO???????????????";
bestelm *bestrees;
boolean dummy;
longer seed;
Char ch;
Char * progname;

/* Local variables for maketree, propagated globally for C version: */
long maxwhich;
double like, maxsteps, bestyet, bestlike, bstlike2, bestfound;
boolean lastrearr, recompute;
double nsteps[maxuser];
long **fsteps;
node *there;
long *place;
boolptr in_tree;
baseptr nothing;
node *temp, *temp1;


tree* dnacomp_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(dnacomp_tree));
  dnacomp_tree_init(t, nonodes, spp);
  return t;
}


void dnacomp_tree_init(tree* t, long nonodes, long spp)
{
  dnapars_tree_init(t, nonodes, spp);
  t->evaluate = dnacomp_tree_evaluate;
}


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  putchar('\n');
  maxtrees = 10000;
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
  weights = false;
  justwts = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  interleaved = true;
  loopcount = 0;
  for (;;)
  {
    cleerhome();
    printf("\nDNA compatibility algorithm, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree)
    {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
      {
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      }
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets, (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4  Print steps & compatibility at sites  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts)
    {
      printf("WARNING:  W option and Multiple Weights options are both on.  ");
      printf("The W menu option is unnecessary and has no additional effect. \n");
    }
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WJOTUMI1234560", ch) != NULL))
        || (usertree && ((strchr("WOTUMI1234560", ch) != NULL))))
    {
      switch (ch)
      {
        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
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
              justweights(&msets);
            else
              initdatasets(&msets);
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

        case 'W':
          weights = !weights;
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
          stepbox = !stepbox;
          break;

        case '5':
          ancseq = !ancseq;
          break;

        case '6':
          trout = !trout;
          break;
      }
    }
    else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
}  /* getoptions */


void reallocchars(void)
{/* The amount of chars can change between runs
    this function reallocates all the variables
    whose size depends on the amount of chars */
  long i;

  for (i = 0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));
  }
  free(weight);
  free(oldweight);
  free(enterorder);
  free(necsteps);
  free(alias);
  free(ally);
  free(location);
  free(in_tree);

  weight = (steptr)Malloc(chars * sizeof(long));
  oldweight = (steptr)Malloc(chars * sizeof(long));
  enterorder = (long *)Malloc(spp * sizeof(long));
  necsteps = (steptr)Malloc(chars * sizeof(long));
  alias = (steptr)Malloc(chars * sizeof(long));
  ally = (steptr)Malloc(chars * sizeof(long));
  location = (steptr)Malloc(chars * sizeof(long));
  in_tree = (boolptr)Malloc(chars * sizeof(boolean));
}


void allocrest(void)
{
  long i;

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));
  bestrees = (bestelm *) Malloc(maxtrees * sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  weight = (steptr)Malloc(chars * sizeof(long));
  oldweight = (steptr)Malloc(chars * sizeof(long));
  enterorder = (long *)Malloc(spp * sizeof(long));
  necsteps = (steptr)Malloc(chars * sizeof(long));
  alias = (steptr)Malloc(chars * sizeof(long));
  ally = (steptr)Malloc(chars * sizeof(long));
  location = (steptr)Malloc(chars * sizeof(long));
  place = (long *)Malloc((2*spp-1) * sizeof(long));
  in_tree = (boolptr)Malloc(spp * sizeof(boolean));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nDNA compatibility algorithm, version %s\n\n", VERSION);

  inputnumbers(&spp, &chars, &nonodes, 1);
  if (!javarun)
  {
    getoptions();
  }
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, chars);
  allocrest();
}  /* doinit */


void makeweights(void)
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++)
  {
    alias[i - 1] = i;
    oldweight[i - 1] = weight[i - 1];
    ally[i - 1] = i;
  }
  sitesort(chars, weight);
  sitecombine(chars);
  sitescrunch(chars);
  endsite = 0;
  for (i = 1; i <= chars; i++)
  {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
}  /* makeweights */


void doinput(void)
{
  /* reads the input data */

  long i;

  if (justwts)
  {
    if (firstset)
      inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    inputweights(chars, weight, &weights);
    if (justwts)
    {
      fprintf(outfile, "\n\nWeights set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nWeights set # %ld:\n\n", ith);
        print_progress(progbuf);
      }
    }
    if (printdata)
      printweights(outfile, 0, chars, weight, "Sites");
  }
  else
  {
    if (!firstset)
    {
      samenumsp(&chars, ith);
      reallocchars();
    }
    inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (weights)
    {
      inputweights(chars, weight, &weights);
      if (printdata)
        printweights(outfile, 0, chars, weight, "Sites");
    }
  }
  makeweights();
  curtree = dnacomp_tree_new(nonodes, spp);
  bestree = dnacomp_tree_new(nonodes, spp);
  priortree = dnacomp_tree_new(nonodes, spp);
  dna_makevalues(curtree, usertree);
}  /* doinput */


void mincomp(long n)
{
  /* computes for each site the minimum number of steps
     necessary to accomodate those species already
     in the analysis, adding in species n */
  long i, j, k, l, m;
  bases b;
  long s;
  boolean allowable, deleted;

  in_tree[n - 1] = true;

  for (i = 0; i < endsite; i++)
    necsteps[i] = 3;

  for (m = 0; m <= 31; m++)
  {
    s = 0;
    l = -1;
    k = m;
    for (b = A; (long)b <= (long)O; b = (bases)((long)b + 1))
    {
      if ((k & 1) == 1)
      {
        s |= 1L << ((long)b);
        l++;
      }
      k /= 2;
    }
    for (j = 0; j < endsite; j++)
    {
      allowable = true;
      i = 1;
      while (allowable && (i <= spp))
      {
        if (in_tree[i - 1] && ((dnapars_node*)curtree->nodep[i - 1])->base[j]
            != 0)
        {
          if ((((dnapars_node*)curtree->nodep[i - 1])->base[j] & s) == 0)
            allowable = false;
        }
        i++;
      }
      if (allowable)
      {
        if (l < necsteps[j])
          necsteps[j] = l;
      }
    }
  }

  for (j = 0; j < endsite; j++)
  {
    deleted = false;
    for (i = 0; i < spp; i++)
    {
      if (in_tree[i] && ((dnapars_node*)curtree->nodep[i])->base[j] == 0)
        deleted = true;
    }
    if (deleted)
      necsteps[j]++;
  }

  for (i = 0; i < endsite; i++)
    necsteps[i] *= weight[i];

}  /* mincomp */


double dnacomp_tree_evaluate(tree* t, node *r, boolean dummy)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, term, steps;
  double sum;

  generic_tree_evaluate(t, r, dummy);

  sum = 0.0;
  for (i = 0; i < endsite; i++)
  {
    steps = ((pars_node*)r)->numsteps[i] + ((pars_node*)r->back)->numsteps[i];
    if ( (((dnapars_node*)r)->base[i] & ((dnapars_node*)r->back)->base[i] )
         == 0 )
      steps += weight[i];
    if (steps == necsteps[i])
      term = weight[i];
    else
      term = 0;
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser)
  {
    nsteps[which - 1] = sum;
    if (which == 1)
    {
      maxwhich = 1;
      maxsteps = sum;
    }
    else if (sum > maxsteps)
    {
      maxwhich = which;
      maxsteps = sum;
    }
  }
  like = sum;
  t->score = like;

  return like;
}  /* dnacomp_tree_evaluate */


void describe(void)
{
  /* prints ancestors, steps and table of numbers of steps in
     each site and table of compatibilities */
  long i, j, k;

  if (treeprint)
  {
    fprintf(outfile, "\ntotal number of compatible sites is ");
    fprintf(outfile, "%10.1f\n", like);
  }
  if (stepbox)
  {
    writesteps(curtree, chars, weights, oldweight);
    fprintf(outfile,
            "\n compatibility (Y or N) of each site with this tree:\n\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%ld", i);
    fprintf(outfile, "\n     *----------\n");
    for (i = 0; i <= (chars / 10); i++)
    {
      putc(' ', outfile);
      fprintf(outfile, "%3ld !", i * 10);
      for (j = 0; j <= 9; j++)
      {
        k = i * 10 + j;
        if (k > 0 && k <= chars)
        {
          if (((pars_node*)curtree->root)->numsteps
              [location[ally[k - 1] - 1] - 1] ==
              necsteps[location[ally[k - 1] - 1] - 1])
          {
            if (oldweight[k - 1] > 0)
              putc('Y', outfile);
            else
              putc('y', outfile);
          }
          else
          {
            if (oldweight[k - 1] > 0)
              putc('N', outfile);
            else
              putc('n', outfile);
          }
        }
        else
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
  }
  if (ancseq)
  {
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout)
  {
    col = 0;
    treeout(curtree->root, nextree, &col, curtree->root);
  }
}  /* describe */


void initboolnames(node *p, boolean *names)
{
  /* sets BOOLEANs that indicate tips */
  node *q;

  if (p->tip)
  {
    names[p->index - 1] = true;
    return;
  }
  q = p->next;
  while (q != p)
  {
    initboolnames(q->back, names);
    q = q->next;
  }
}  /* initboolnames */


void standev3(long chars, long numtrees, long maxwhich, double maxsteps, double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH) on user-defined trees */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

  if (numtrees == 2)
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Compatible  Difference  Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees)
    {
      fprintf(outfile, "%3ld  %11.1f", which, nsteps[which - 1]);
      if (maxwhich == which)
        fprintf(outfile, "  <------ best\n");
      else
      {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < chars; i++)
        {
          if (weight[i] > 0)
          {
            wt = weight[i];
            sumw += wt;
            temp = (fsteps[maxwhich - 1][i] - fsteps[which - 1][i]);
            sum += temp;
            sum2 += temp * temp / wt;
          }
        }
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum /sumw));
        fprintf(outfile, " %10.1f %11.4f",
                (maxsteps-nsteps[which - 1]), sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  }
  else
  {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n", MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    }
    else
    {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees * sizeof(double *));
    sumw = 0.0;
    for (i = 0; i < chars; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees * sizeof(double));
    for (i = 0; i < numtrees; i++)      /* compute covariances of trees */
    {
      sum = nsteps[i]/sumw;
      for (j = 0; j <=i; j++)
      {
        sum2 = nsteps[j]/sumw;
        temp = 0.0;
        for (k = 0; k < chars; k++)
        {
          if (weight[k] > 0)
          {
            wt = weight[k];
            temp = temp + (fsteps[i][k]- wt * sum)
              * (fsteps[j][k]- wt * sum2) / wt;
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++)      /* in-place Cholesky decomposition of trees x trees covariance matrix */
    {
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++)
      {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees * sizeof(double)); /* resamples sums */
    P = (double *)Malloc(numtrees * sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees * sizeof(double)); /* store normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0];             /* sum2 will be largest # of compat. sites */
    for (i = 1; i < numtrees; i++)
      if (sum2 < nsteps[i])
        sum2 = nsteps[i];
    for (i = 1; i <= SAMPLES; i++)      /* loop over resampled trees */
    {
      for (j = 0; j < numtrees; j++)    /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++)    /* compute vectors */
      {
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)    /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)    /* accumulate P's */
        if (sum2-nsteps[j] <= sum-f[j])
          P[j] += 1.0 / SAMPLES;
    }
    fprintf(outfile, "Tree   Compatible  Difference   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++)
    {
      fprintf(outfile, "%3ld  %10.1f", i+1, nsteps[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        fprintf(outfile, " %10.1f  %10.3f", sum2-nsteps[i], P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
    fprintf(outfile, "\n");
    free(P);             /* free the variables we Malloc'ed */
    free(f);
    free(r);
    for (i = 0; i < numtrees; i++)
      free(covar[i]);
    free(covar);
  }
}  /* standev2 */


void maketree(void)                     // RSGbugfix
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, k, numtrees, nextnode;
  boolean firsttree, goteof, haslengths, thorough = true;
  node *item, *dummy;
#if 0                                   // RSGbugfix
  pointarray nodep;
#endif
  boolean *names;
  boolean multif;

  if (!usertree)
  {
    recompute = true;
    for (i = 0; i < spp; i++)
      in_tree[i] = false;
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    buildsimpletree(curtree, enterorder);
    curtree->root = curtree->nodep[enterorder[0] - 1];
    if (progress)
    {
      sprintf(progbuf, "Adding species:\n");
      print_progress(progbuf);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }
    in_tree[0] = true;
    in_tree[1] = true;
    in_tree[2] = true;
    lastrearr = false;
    for (i = 4; i <= spp; i++)
    {
      mincomp(i);
      bestyet = -350.0 * spp * chars;
      item = curtree->nodep[enterorder[i - 1] - 1];
      there = curtree->root;
      k = generic_tree_findemptyfork(curtree);
      p = curtree->get_fork(curtree, k);
      hookup(item, p);
      curtree->addtraverse(curtree, item->back, curtree->root, true, there,
                        &bestyet, bestree, thorough, true, false, &bestfound);
      curtree->insert_(curtree, item->back, there, false);
      like = bestyet;
      curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1], true,
                       &bestyet,  bestree, priortree, (i == spp), &bestfound);
      if (progress)
      {
        writename(i - 1, 1, enterorder);
        phyFillScreenColor();
      }
      lastrearr = (i == spp);
      if (lastrearr)
      {
        if (progress)
        {
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
        if (jumb == 1)
        {
          bstlike2 = bestlike;
          nextree = 1;
        }
        curtree->globrearrange(curtree, bestree, progress,
                                 thorough, &bestfound);
      }
    }
    if (progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
    for (i = spp - 1; i >= 1; i--)
      curtree->re_move(curtree, curtree->nodep[i], &dummy, recompute);
    if (jumb == njumble)
    {
      if (treeprint)
      {
        putc('\n', outfile);
        if (nextree == 1)
          fprintf(outfile, "One most compatible tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree);
      }
      if (nextree > maxtrees + 1)
      {
        if (treeprint)
          fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      recompute = false;
      for (i = 0; i <= (nextree - 2); i++)
      {
        load_tree(curtree, i, bestrees);
        curtree->evaluate(curtree, curtree->root, 0);
        curtree->root = root_tree(curtree, curtree->root);
        curtree->nodep[curtree->root->index - 1] = curtree->root;
        printree(curtree);
        describe();
        reroot_tree(curtree, curtree->root); // RSGbugfix: Name change.
      }
    }
  }
  else
  {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    fsteps = (long **)Malloc(maxuser * sizeof(long *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (long *)Malloc(endsite * sizeof(long));
    names = (boolean *)Malloc(spp * sizeof(boolean));

#if 0                                   // RSGbugfix
    nodep = NULL;
#endif

    maxsteps = 0.0;
    which = 1;
    while (which <= numtrees)
    {
      firsttree = true;
      nextnode = 0;
      haslengths = true;
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof, &firsttree, &nextnode, &haslengths, initparsnode, false, nonodes);
      for (j = 0; j < spp; j++)
        names[j] = false;
      initboolnames(curtree->root, names);
      i = 0;
      for (j = 0; j < spp; j++)
      {
        in_tree[j] = names[j];
        if (in_tree[j])
          i = j;
      }
      mincomp(i+1);
      reroot_tree(curtree, curtree->root); // RSGbugfix: Name change.
      curtree->evaluate(curtree, curtree->root, false);
      curtree->root = root_tree(curtree, curtree->root);
      if (outgropt)
        reroot(curtree->nodep[outgrno - 1], curtree->root);
      printree(curtree);
      describe();
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 )
    {
      standev3(chars, numtrees, maxwhich, maxsteps, nsteps, fsteps, seed);
    }
    for (j = 1; j <= maxuser; j++)
      free(fsteps[j - 1]);
    free(fsteps);
    free(names);
  }
  if (jumb == njumble)
  {
    if (progress)
    {
      sprintf(progbuf, "\nOutput written to file \"%s\".\n", outfilename);
      print_progress(progbuf);
      if (trout)
      {
        sprintf(progbuf, "\nTrees also written onto file \"%s\".\n", outtreename);
        print_progress(progbuf);
      }
      sprintf(progbuf, "\nDone.\n\n");
      print_progress(progbuf);
    }
  }
}  /* maketree */


void dnacomprun(void)
{
  // JRMdebug
  /*
  printf("\njumble: %i\n", jumble);
  printf("njumble: %li\n", njumble);
  printf("outgrno %li\n", outgrno);
  printf("outgropt %i\n", outgropt);
  printf("trout %i\n", trout);
  printf("usertree %i\n", usertree);
  printf("weights %i\n", weights);
  printf("justwts %i\n", justwts);
  printf("printdata %i\n", printdata);
  printf("progress %i\n", progress);
  printf("treeprint %i\n", treeprint);
  printf("stepbox %i\n", stepbox);
  printf("ancseq %i\n", ancseq);
  printf("dotdiff %i\n", dotdiff);
  printf("interleaved %i\n", interleaved);
  fflush(stdout);
  */

  // do the work
  for (ith = 1; ith <= msets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (msets > 1 && !justwts) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "Data set # %ld:\n\n", ith);
        print_progress(progbuf);
      }
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    fflush(outfile);
    fflush(outtree);
  }
}


void dnacomp(
  char * infilename,
  char * intreename,
  char * OutfileName,
  char * outfileopt,
  char * weightsfilename,
  char * OuttreeName,
  char * outtreeopt,
  int searchbest,
  int TreeSave,
  int inputorder,
  int RandNum,
  int NumJumble,
  int OutRoot,
  int OutNum,
  int SitesWeighted,
  int AnalyzeMult,
  int MultDataSet,
  int NumMult,
  int InputSeq,
  int PrintData,
  int DotDiff,
  int PrintInd,
  int PrintTree,
  int PrintSteps,
  int PrintSeq,
  int WriteTree)
{
  initdata funcs;

  (void)searchbest;                     // RSGnote: Parameter never referenced.
  //printf("Hello from DnaComp!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Dnacomp";

  memset(&funcs, 0, sizeof(initdata));
  funcs.node_new = dnapars_node_new;
  funcs.tree_new = dnacomp_tree_new;

  phylipinit(argc, argv, &funcs, true);
  progname = argv[0];

  // JRMdebug
  /*
    maxtrees = 10000;
    jumble = false;
    njumble = 1;
    outgrno = 1;
    outgropt = false;
    trout = true;
    usertree = false;
    weights = false;
    justwts = false;
    printdata = false;
    dotdiff = true;
    progress = true;
    treeprint = true;
    stepbox = false;
    ancseq = false;
    interleaved = true;

    char * infile,
    char * intree,
    char * outfile,
    char * outfileopt,
    char * weightfile,
    char * outtree,
    char * outtreeopt,
    int searchbest,
    int TreeSave,
    int inputorder,
    int RandNum,
    int NumJumble,
    int OutRoot,
    int OutNum,
    int SitesWeighted,
    int AnalyzeMult,
    int multdataset,
    int NumSites,
    int inputseq,
    int PrintData,
    int DotDiff,
    int PrintInd,
    int PrintTree,
    int PrintSteps,
    int PrintSeq,
    int WriteTree)
  */

  maxtrees = TreeSave;

  if (inputorder != 0)
  {
    jumble = true;
  }
  else
  {
    jumble = false;
  }

  inseed =  RandNum;

  njumble = NumJumble;

  if (OutRoot != 0)
  {
    outgropt = true;
    outgrno = OutNum;
    initoutgroup(&outgrno, spp);
  }
  else
  {
    outgropt = false;
    outgrno = 1;
  }

  if (SitesWeighted != 0)
  {
    weights = true;
  }
  else
  {
    weights = false;
  }

  if (AnalyzeMult != 0)
  {
    mulsets = true;
    msets = NumMult;
  }
  else
  {
    mulsets = false;
    msets = 1;
  }

  if (MultDataSet != 0)
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

  if (PrintData != 0)
  {
    printdata = true;
  }
  else
  {
    printdata = false;
  }

  if (DotDiff != 0)
  {
    dotdiff = true;
  }
  else
  {
    dotdiff = false;
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
    trout = true;
  }
  else
  {
    trout = false;
  }

  if (PrintSteps != 0)
  {
    stepbox = true;
  }
  else
  {
    stepbox = false;
  }

  if (PrintSeq != 0)
  {
    ancseq = true;
  }
  else
  {
    ancseq = false;
  }

  if (WriteTree != 0)
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

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  firstset = true;
  doinit();

  if (weights || justwts)
  {
    weightfile = fopen(weightsfilename, "r");
  }

  if (trout)
  {
    outtree = fopen(OuttreeName, outtreeopt);
    strcpy(outtreename, OuttreeName);
  }

  if (jumble) {
    initjumble(&inseed, &inseed0, seed, &njumble);
  }

  dnacomprun();  // do the actual work

  FClose(infile);
  FClose(outfile);

  if (weights || justwts)
  {
    FClose(weightfile);
  }
  if (trout)
  {
    FClose(outtree);
  }
  if (usertree)
  {
    FClose(intree);
  }

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* DNA compatibility by uphill search */
  /* reads in spp, chars, and the data. Then calls maketree to
     construct the tree */
  initdata funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Dnacomp", "");        */
  argv[0]="Dnacomp";
#endif
  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = dnapars_node_new;
  funcs.tree_new = dnacomp_tree_new;
  progname = argv[0];
  phylipinit(argc, argv, &funcs, false);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  mulsets = false;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  msets = 1;
  firstset = true;
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  dnacomprun();

#if 0
  for (ith = 1; ith <= msets; ith++)
  {
    doinput();
    if (ith == 1)
      firstset = false;
    if (msets > 1 && !justwts)
    {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
        printf("Data set # %ld:\n\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    fflush(outfile);
    fflush(outtree);
  }
#endif

  FClose(infile);
  FClose(outfile);
  FClose(outtree);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  phyRestoreConsoleAttributes();
  exxit(0);
  return 0;
}  /* DNA compatibility by uphill search */


// End.
