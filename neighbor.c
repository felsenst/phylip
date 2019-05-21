/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Mary Kuhner, Jon Yamato, Joseph Felsenstein, Akiko Fuseki,
   Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <float.h>

#include "phylip.h"
#include "dist.h"

#ifndef OLDC
/* function prototypes */
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void getinput(void);
void describe(node *, double);
void summarize(void);
void nodelabel(boolean);
void jointree(void);
void maketree(void);
void freerest(void);
void neighborrun(void);
void neighbor(char * infilename, char * outfilename, char * outfileopt, char * outtreename,
              char * outtreeopt, char * Method, int LowerTMat, int UpperTMat, int Subreps,
              int RandInput, int RandNum, int OutRoot, int OutNum, int MultData, int NumSets,
              int PrintData, int PrintInd, int PrintTree, int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH];
long outgrno, col, datasets, ith, nonodes = 0;
long inseed;
vector *x;
intvector *reps;
boolean jumble, lower, upper, outgropt, replicates, trout, printdata, progress, treeprint, mulsets, njoin;
tree curtree;
longer seed;
long *enterorder;
Char progname[20];

/* Local variables for maketree, propagated globally for C version: */
node **cluster;


void getoptions(void)
{
  /* interactively set options */
  long inseed0 = 0, loopcount;
  Char ch;

  putchar('\n');
  jumble = false;
  lower = false;
  outgrno = 1;
  outgropt = false;
  replicates = false;
  trout = true;
  upper = false;
  printdata = false;
  progress = true;
  treeprint = true;
  njoin = true;
  loopcount = 0;
  for(;;) {
    cleerhome();
    printf("\nNeighbor-Joining/UPGMA method version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  N       Neighbor-joining or UPGMA tree?  %s\n",
           (njoin ? "Neighbor-joining" : "UPGMA"));
    if (njoin) {
      printf("  O                        Outgroup root?");
      if (outgropt)
        printf("  Yes, at species number%3ld\n", outgrno);
      else
        printf("  No, use as outgroup species%3ld\n", outgrno);
    }
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    printf("  J     Randomize input order of species?");
    if (jumble)
      printf("  Yes (random number seed =%8ld)\n", inseed0);
    else
      printf("  No. Use input order\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\n\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if  (ch == 'Y')
      break;
    if (strchr("NJOULRSM01234", ch) != NULL)
    {
      switch (ch)
      {
        case 'J':
          jumble = !jumble;
          if (jumble)
            initseed(&inseed, &inseed0, seed);
          break;

        case 'L':
          lower = !lower;
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          else
            outgrno = 1;
          break;

        case 'R':
          upper = !upper;
          break;

        case 'S':
          replicates = !replicates;
          break;

        case 'N':
          njoin = !njoin;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets)
            initdatasets(&datasets);
          jumble = true;
          if (jumble)
            initseed(&inseed, &inseed0, seed);
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
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
}  /* getoptions */


void allocrest(void)
{
  long i;

  x = (vector *)Malloc(spp * sizeof(vector));
  for (i = 0; i < spp; i++)
    x[i] = (vector)Malloc(spp * sizeof(double));
  reps = (intvector *)Malloc(spp * sizeof(intvector));
  for (i = 0; i < spp; i++)
    reps[i] = (intvector)Malloc(spp * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
  cluster = (node **)Malloc(spp * sizeof(node *));
}  /* allocrest */


void freerest(void)
{
  long i;

  for (i = 0; i < spp; i++)
    free(x[i]);
  free(x);
  for (i = 0; i < spp; i++)
    free(reps[i]);
  free(reps);
  free(nayme);
  free(enterorder);
  free(cluster);
}  /* freerest */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nNeighbor-Joining/UPGMA method version %s\n\n", VERSION);

  node *p;

  inputnumbers2(&spp, &nonodes, 2);
  nonodes += (njoin ? 0 : 1);
  if (!javarun)
  {
    getoptions();
  }
  alloctree(&curtree.nodep, nonodes+1);
  p = curtree.nodep[nonodes]->next;
  curtree.nodep[nonodes]->next = curtree.nodep[nonodes];
  free(p->next);
  free(p);
  allocrest();
}  /* doinit */


void inputoptions(void)
{
  /* read options information */

  if (ith != 1)
    samenumsp2(ith);
  putc('\n', outfile);
  if (njoin) {
    fprintf(outfile, " Neighbor-joining method\n");
    fprintf(outfile, "\n Negative branch lengths allowed\n\n");
  }
  else
    fprintf(outfile, " UPGMA method\n");
}  /* inputoptions */


void describe(node *p, double height)
{
  /* print out information for one branch */
  long i;
  node *q;

  q = p->back;
  if (njoin)
    fprintf(outfile, "%4ld          ", q->index - spp);
  else
    fprintf(outfile, "%4ld     ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
    putc(' ', outfile);
  } else {
    if (njoin)
      fprintf(outfile, "%4ld       ", p->index - spp);
    else {
      fprintf(outfile, "%4ld       ", p->index - spp);
    }
  }
  if (njoin)
    fprintf(outfile, "%12.5f\n", q->v);
  else
    fprintf(outfile, "%10.5f      %10.5f\n", q->v, q->v+height);
  if (!p->tip) {
    describe(p->next->back, height+q->v);
    describe(p->next->next->back, height+q->v);
  }
}  /* describe */


void summarize(void)
{
  /* print out branch lengths etc. */
  putc('\n', outfile);
  if (njoin) {
    fprintf(outfile, "remember:");
    if (outgropt)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  if (njoin) {
    fprintf(outfile, "\nBetween        And            Length\n");
    fprintf(outfile, "-------        ---            ------\n");
  } else {
    fprintf(outfile, "From     To            Length          Height\n");
    fprintf(outfile, "----     --            ------          ------\n");
  }
  describe(curtree.root->next->back, 0.0);
  describe(curtree.root->next->next->back, 0.0);
  if (njoin)
    describe(curtree.root->back, 0.0);
  fprintf(outfile, "\n\n");
}  /* summarize */


void nodelabel(boolean isnode)
{
  if (isnode)
  {
    sprintf(progbuf, "node   ");
    print_progress(progbuf);
  }
  else
  {
    sprintf(progbuf, "species");
    print_progress(progbuf);
  }
}  /* nodelabel */


void jointree(void)
{
  /* calculate the tree */
  long nc, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, iter;
  double fotu2, total, tmin, dio, djo, bi, bj, bk, dmin=0, da;
  long el[3];
  vector av;
  intvector oc;

  double *R;   /* added in revisions by Y. Ina */
  R = (double *)Malloc(spp * sizeof(double));

  for (i = 0; i <= spp - 2; i++) {
    for (j = i + 1; j < spp; j++) {
      da = (x[i][j] + x[j][i]) / 2.0;
      x[i][j] = da;
      x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = spp - 2.0;
  nextnode = spp + 1;
  av = (vector)Malloc(spp * sizeof(double));
  oc = (intvector)Malloc(spp * sizeof(long));
  for (i = 0; i < spp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  if (njoin)
    iter = spp - 3;
  else
    iter = spp - 1;
  for (nc = 1; nc <= iter; nc++) {
    for (j = 2; j <= spp; j++) {
      for (i = 0; i <= j - 2; i++)
        x[j - 1][i] = x[i][j - 1];
    }
    tmin = DBL_MAX;
    /* Compute sij and minimize */
    if (njoin) {     /* many revisions by Y. Ina from here ... */
      for (i = 0; i < spp; i++)
        R[i] = 0.0;
      for (ja = 2; ja <= spp; ja++) {
        jj = enterorder[ja - 1];
        if (cluster[jj - 1] != NULL) {
          for (ia = 0; ia <= ja - 2; ia++) {
            ii = enterorder[ia];
            if (cluster[ii - 1] != NULL) {
              R[ii - 1] += x[ii - 1][jj - 1];
              R[jj - 1] += x[ii - 1][jj - 1];
            }
          }
        }
      }
    } /* ... to here */
    for (ja = 2; ja <= spp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            if (njoin) {
              total = fotu2 * x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1];
               /* this statement part of revisions by Y. Ina */
            } else
              total = x[ii - 1][jj - 1];
            if (total < tmin) {
              tmin = total;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    /* compute lengths and print */
    if (njoin) {
      dio = 0.0;
      djo = 0.0;
      for (i = 0; i < spp; i++) {
        dio += x[i][mini - 1];
        djo += x[i][minj - 1];
      }
      dmin = x[mini - 1][minj - 1];
      dio = (dio - dmin) / fotu2;
      djo = (djo - dmin) / fotu2;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini - 1];
      bj -= av[minj - 1];
    } else {
      bi = x[mini - 1][minj - 1] / 2.0 - av[mini - 1];
      bj = x[mini - 1][minj - 1] / 2.0 - av[minj - 1];
      av[mini - 1] += bi;
    }
    if (progress) {
      sprintf(progbuf, "Cycle %3ld: ", iter - nc + 1);
      print_progress(progbuf);
      if (njoin)
        nodelabel((boolean)(av[mini - 1] > 0.0));
      else
        nodelabel((boolean)(oc[mini - 1] > 1.0));
      sprintf(progbuf, " %ld (%10.5f) joins ", mini, bi);
      print_progress(progbuf);
      if (njoin)
        nodelabel((boolean)(av[minj - 1] > 0.0));
      else
        nodelabel((boolean)(oc[minj - 1] > 1.0));
      sprintf(progbuf, " %ld (%10.5f)\n", minj, bj);
      print_progress(progbuf);
      phyFillScreenColor();
    }
    hookup(curtree.nodep[nextnode - 1]->next, cluster[mini - 1]);
    hookup(curtree.nodep[nextnode - 1]->next->next, cluster[minj - 1]);
    cluster[mini - 1]->v = bi;
    cluster[minj - 1]->v = bj;
    cluster[mini - 1]->back->v = bi;
    cluster[minj - 1]->back->v = bj;
    cluster[mini - 1] = curtree.nodep[nextnode - 1];
    cluster[minj - 1] = NULL;
    nextnode++;
    if (njoin)
      av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
    for (j = 0; j < spp; j++) {
      if (cluster[j] != NULL) {
        if (njoin) {
          da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
          if (mini - j - 1 < 0)
            x[mini - 1][j] = da;
          if (mini - j - 1 > 0)
            x[j][mini - 1] = da;
        } else {
          da = x[mini - 1][j] * oc[mini - 1] + x[minj - 1][j] * oc[minj - 1];
          da /= oc[mini - 1] + oc[minj - 1];
          x[mini - 1][j] = da;
          x[j][mini - 1] = da;
        }
      }
    }
    for (j = 0; j < spp; j++) {
      x[minj - 1][j] = 0.0;
      x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= spp; i++) {
    if (cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }
  if (!njoin) {
    curtree.root = cluster[el[0] - 1];
    curtree.root->back = NULL;
    free(av);
    free(oc);
    return;
  }
  bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = x[el[0] - 1][el[1] - 1] - bi;
  bk = x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];
  if (progress) {
    sprintf(progbuf, "last cycle:\n");
    print_progress(progbuf);
    sprintf(progbuf, " ");
    print_progress(progbuf);
    nodelabel((boolean)(av[el[0] - 1] > 0.0));
    sprintf(progbuf, " %ld  (%10.5f) joins ", el[0], bi);
    print_progress(progbuf);
    nodelabel((boolean)(av[el[1] - 1] > 0.0));
    sprintf(progbuf, " %ld  (%10.5f) joins ", el[1], bj);
    print_progress(progbuf);
    nodelabel((boolean)(av[el[2] - 1] > 0.0));
    sprintf(progbuf, " %ld  (%10.5f)\n", el[2], bk);
    print_progress(progbuf);
    phyFillScreenColor();
  }
  hookup(curtree.nodep[nextnode - 1], cluster[el[0] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next, cluster[el[1] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next->next, cluster[el[2] - 1]);
  cluster[el[0] - 1]->v = bi;
  cluster[el[1] - 1]->v = bj;
  cluster[el[2] - 1]->v = bk;
  cluster[el[0] - 1]->back->v = bi;
  cluster[el[1] - 1]->back->v = bj;
  cluster[el[2] - 1]->back->v = bk;
  curtree.root = cluster[el[0] - 1]->back;
  free(av);
  free(oc);
  free(R);
}  /* jointree */


void maketree(void)
{
  /* construct the tree */
  long i ;

  inputdata(replicates, printdata, lower, upper, x, reps);
  if (njoin && (spp < 3)) {
    printf("\nERROR:  Neighbor-Joining runs must have at least 3 species.\n\n");
    exxit(-1);
  }
  if (progress)
  {
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
  if (ith == 1)
    setuptree(&curtree, nonodes + 1);
  for (i = 1; i <= spp; i++)
    enterorder[i - 1] = i;
  if (jumble)
    randumize(seed, enterorder);
  for (i = 0; i < spp; i++)
    cluster[i] = curtree.nodep[i];
  jointree();
  if (njoin)
    curtree.root = curtree.nodep[outgrno - 1]->back;
  printree(curtree.root, treeprint, !njoin);
  if (treeprint)
    summarize();
  if (trout) {
    col = 0;
    if (njoin)
      treeout(curtree.root, &col, 0.43429448222, njoin, curtree.root);
    else
      curtree.root = curtree.root,
      treeoutr(curtree.root, &col, &curtree);
  }
  if (progress) {
    sprintf(progbuf, "\nOutput written on file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree written on file \"%s\".\n\n", outtreename);
      print_progress(progbuf);
    }
  }
}  /* maketree */


void neighborrun(void)
{
  // JRMdebug
  /*
  printf("\njumble %i\n", jumble);
  printf("lower  %i\n", lower);
  printf("outgrno  %li\n", outgrno);
  printf("outgropt  %i\n", outgropt);
  printf("mulsets %i\n", mulsets);
  printf("datasets %li\n", datasets);
  printf("replicates  %i\n", replicates);
  printf("trout  %i\n", trout);
  printf("upper  %i\n", upper);
  printf("printdata  %i\n", printdata);
  printf("progress  %i\n", progress);
  printf("treeprint  %i\n", treeprint);
  printf("njoin  %i\n", njoin);
  */

  ith = 1;
  while (ith <= datasets) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    inputoptions();
    maketree();
    fflush(outfile);
    fflush(outtree);
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
    ith++;
  }
}


void neighbor(
  char * infilename,
  char * OutfileName,
  char * outfileopt,
  char * OuttreeName,
  char * outtreeopt,
  char * Method,
  int LowerTMat,
  int UpperTMat,
  int Subreps,
  int RandInput,
  int RandNum,
  int OutRoot,
  int OutNum,
  int MultData,
  int NumSets,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree)
{
  initdata *funcs;

  //printf("Hello from Neighbor!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Neighbor";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dist_node_new;
  phylipinit(argc, argv, funcs, true);

  /*
  //jumble = false;
  //lower = false;
  //outgrno = 1;
  //outgropt = false;
  //replicates = false;
  //trout = true;
  //upper = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  //njoin = true;
  String infile,
  String outfile,
  String outfileopt,
  String outtree,
  String outtreeopt,
  String Method,
  //boolean LowerTMat,
  //boolean UpperTMat,
  //boolean Subreps,
  //boolean RandInput,
  //int RandNum,
  //boolean OutRoot,
  //int OutNum,
  //boolean MultData,
  //int NumSets,
  boolean PrintData,
  boolean PrintInd,
  boolean PrintTree,
  boolean WriteTree);
  */

  if (!strcmp(Method, "Neighbor"))
  {
    njoin = true;
  }
  else //"UPGMA"
  {
    njoin = false;
  }

  if (LowerTMat != 0)
  {
    lower = true;
  }
  else
  {
    lower = false;
  }

  if (UpperTMat != 0)
  {
    upper = true;
  }
  else
  {
    upper = false;
  }

  if (Subreps != 0)
  {
    replicates = true;
  }
  else
  {
    replicates = false;
  }

  if (RandInput != 0)
  {
    jumble = true;
    inseed = RandNum;
  }
  else
  {
    jumble = false;
    inseed = 1;
  }

  if (OutRoot != 0)
  {
    outgrno = OutNum;
    outgropt = true;
  }
  else
  {
    outgrno = 1;
    outgropt = false;
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

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  if (trout)
  {
    outtree = fopen(OuttreeName, outtreeopt);
    strcpy(outtreename, OuttreeName);
  }

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinit();
  neighborrun();

  FClose(infile);
  FClose(outfile);
  if (trout)
  {
    FClose(outtree);
  }
  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* main program */
  initdata *funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Neighbor", "");                */
  argv[0] = "Neighbor";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = dist_node_new;
  phylipinit(argc, argv, funcs, false);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  doinit();
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  neighborrun();

  FClose(infile);
  FClose(outfile);
  FClose(outtree);
  freetree(&curtree.nodep, nonodes+1);
  freerest();
  /* FIXME -- put in frees for funcs/initdata */
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}


// End.
