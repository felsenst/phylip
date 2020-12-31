/* Version 4.0 (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "discreteparsimony.h"

extern long maxtrees;   /* from parsimony.c */
extern long nextree;    /* from parsimony.c */

#ifndef OLDC
/* function prototypes */
void   pars_tree_setup(long, long);
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   makeweights(void);
void   doinput(void);
void   evaluate(node *);
void   tryadd(node *, node *, node *);         /* debug: used at all? */
void   addpreorder(node *, node *, node *);         /* debug: used at all? */
void   trydescendants(node *, node *, node *, node *, boolean);         /* debug: used at all? */
void   trylocal(node *, node *);         /* debug: used at all? */
void   trylocal2(node *, node *, node *);         /* debug: used at all? */
void   tryrearr(node *p, boolean *);         /* debug: used at all? */
void   repreorder(node *p, boolean *);         /* debug: used at all? */
void   rearrange(node **);         /* debug: used at all? */
void   describe(void);
void   pars_coordinates(node *, double, long *, double *);
void   pars_printree(void);
void   maketree(void);
void   freerest(void);
void   reallocchars(void);
void parsrun(void);
void pars(char *infilename, char *intreename, char *weightsfilename, char *outfilename,
          char *outfileopt, char *outtreename, char *outtreeopt, char *TreeSearchMethod,
          char *SearchOpts, int TreeSave, int InputOrder, int RandNum, int NumJumble,
          int OutRoot, int OutNum, int ThreshPars, double ThreshVal, int SitesWeighted,
          int AnalyzeMult, int MultDataset, int NumMult, int InputSeq, int PrintData,
          int DotDiff, int PrintInd, int PrintTree, int PrintSteps, int PrintSeq, int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
long chars, col, msets, ith, njumble, jumb = 0, nonodes = 0;
/*   chars = number of sites in actual sequences */
long inseed, inseed0;
double threshold, bestfound;
boolean jumble, usertree, reusertree, thresh, weights, thorough, rearrfirst,
    trout, progress, stepbox, ancseq, mulsets, justwts, firstset, mulf, multf;
steptr oldweight;
longer seed;
long *enterorder;
tree *curtree, *bestree, *priortree; /* use bestelm in final rearrangements */
char *progname;

/* Local variables for Pascal maketree, propagated globally for C version: */
long minwhich;
double minsteps, bestyet, bestlike, bstlike2, rebestyet;
boolean lastrearr, recompute;
double nsteps[maxuser];
long **fsteps;
node *there;
long *place;
extern bestelm *bestrees, **rebestrees;
double *threshwt;
discbaseptr nothing;
boolean *names;


void pars_tree_setup(long nonodes, long spp)
{
  /* call allocation and initialization of three new tree(s) */

  curtree = discretepars_tree_new(nonodes, spp);
  bestree = discretepars_tree_new(nonodes, spp);
  priortree = discretepars_tree_new(nonodes, spp);
} /* pars_tree_setup */


void getoptions(void)
{
  /* interactively set options */
  long inseed0, loopcount, loopcount2;
  Char ch, ch2;
  char* string;

  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  thorough = true;
  rearrfirst = false;
  maxtrees = 10000;
  trout = true;
  usertree = false;
  reusertree = false;
  weights = false;
  mulsets = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  dotdiff = true;
  interleaved = true;
  loopcount = 0;
  for (;;)
  {
    cleerhome();
    printf("\nDiscrete character parsimony algorithm, version %s\n\n", VERSION);
    printf("Setting for this run:\n");
    if ( reusertree ) string = "Yes, rearrange on user tree";
    else if ( usertree ) string = "No, use user trees in input file";
    else string = "Yes";
    printf("  U                 Search for best tree?  %s\n", string);
    if (!usertree)
    {
      printf("  S                        Search option?  ");
      if (thorough)
        printf("More thorough search\n");
      else if (rearrfirst)
        printf("Rearrange on one best tree\n");
      else
        printf("Less thorough\n");
      printf("  V              Number of trees to save?  %ld\n", maxtrees);
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number %ld\n", outgrno);
    else
      printf("  No, use as outgroup species %ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I            Input species interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"  : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           progress ? "Yes" : "No");
    printf("  3                        Print out tree  %s\n",
           treeprint ? "Yes" : "No");
    printf("  4          Print out steps in each site  %s\n",
           stepbox ? "Yes" : "No");
    printf("  5  Print character at all nodes of tree  %s\n",
           ancseq ? "Yes" : "No");
    if (ancseq || printdata)
      printf("  .  Use dot-differencing to display them  %s\n",
             dotdiff ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
           trout ? "Yes" : "No");
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WSVJOTUMI12345.60", ch) != NULL))
        || (usertree && (strchr("WSVOTUMI12345.60", ch) != NULL)))
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

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
          break;

        case 'W':
          weights = !weights;
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

        case 'U':
          if ( !usertree && !reusertree)
            usertree = true;
          else if (!reusertree && usertree)
            reusertree = true;
          else
          {
            usertree = false;
            reusertree = false;
          }

          break;

        case 'S':
          thorough = !thorough;
          if (!thorough)
          {
            printf("Rearrange on just one best tree?");
            loopcount2 = 0;
            do {
              printf(" (type Y or N)\n");
              phyFillScreenColor();
              if(scanf("%c%*[^\n]", &ch2)) {} // Read char and scan to EOL.
              (void)getchar();
              if (ch2 == '\n')
                ch2 = ' ';
              uppercase(&ch2);
              countup(&loopcount2, 10);
            } while ((ch2 != 'Y') && (ch2 != 'N'));
            rearrfirst = (ch2 == 'Y');
          }
          break;

        case 'V':
          loopcount2 = 0;
          do {
            printf("type the number of trees to save\n");
            phyFillScreenColor();
            if(scanf("%ld%*[^\n]", &maxtrees)) {} // Read number and scan to EOL.
            if (maxtrees > MAXNUMTREES)
              maxtrees = MAXNUMTREES;
            (void)getchar();
            countup(&loopcount2, 10);
          } while (maxtrees < 1);
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
          stepbox = !stepbox;
          break;

        case '5':
          ancseq = !ancseq;
          break;

        case '.':
          dotdiff = !dotdiff;
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
{
  long i;

  for (i = 0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));
  }
  for (i = 0; i < spp; i++)
  {
    free(convtab[i]);
    convtab[i] = (Char *)Malloc(chars * sizeof(Char));
  }

  free(weight);
  free(oldweight);
  free(alias);
  free(ally);
  free(location);

  weight = (long *)Malloc(chars * sizeof(long));
  oldweight = (long *)Malloc(chars * sizeof(long));
  alias = (long *)Malloc(chars * sizeof(long));
  ally = (long *)Malloc(chars * sizeof(long));
  location = (long *)Malloc(chars * sizeof(long));
}


void allocrest(void)
{
  long i;

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));
  convtab = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    convtab[i] = (Char *)Malloc(chars * sizeof(Char));
  if ( !reusertree )
    bestrees = allocbestree();
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
  place = (long *)Malloc(nonodes * sizeof(long));
  weight = (long *)Malloc(chars * sizeof(long));
  oldweight = (long *)Malloc(chars * sizeof(long));
  alias = (long *)Malloc(chars * sizeof(long));
  ally = (long *)Malloc(chars * sizeof(long));
  location = (long *)Malloc(chars * sizeof(long));
}  /* alocrest */


void doinit(void)
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  fprintf(outfile, "\nDiscrete character parsimony algorithm, version %s\n\n",
            VERSION);
  if (!javarun)
  {
    getoptions();
  }
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n\n", spp, chars);
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
  if (!thresh)
    threshold = spp;
  threshwt = (double *)Malloc(endsite * sizeof(double));
  for (i = 0; i < endsite; i++)
  {
    threshwt[i] = (threshold * weight[i]);
  }
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
  pars_tree_setup(nonodes, spp);                  /* set up the three trees */
  makevalues(curtree, usertree);   /* put information on characters at tips */
}  /* doinput */


void describe(void)
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  long indent;

  if (treeprint)
  {
    fprintf(outfile, "\nrequires a total of %10.3f\n", -(curtree->score));
    fprintf(outfile, "\n  between      and       length\n");
    fprintf(outfile, "  -------      ---       ------\n");
    printbranchlengths(curtree->root);
  }
  if (stepbox)
    writesteps(curtree, chars, weights, oldweight);
  if (ancseq)
  {
    disc_hypstates(curtree, chars);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout)
  {
    col = 0;
    indent = 0;
    treeout3(curtree->root, nextree, &col, indent, curtree->root);
  }
}  /* describe */


void pars_coordinates(node *p, double lengthsum, long *tipy,
                      double *tipmax)
{
  /* establishes coordinates of nodes where tree is drawn left to right
   * so y coordinate is across tips and x coordinate is along branches */
  node *q, *first, *last;
  double xx;

  if (p == NULL)
    return;
  if (p->tip)
  {                       /* if it is a tip, assign ycoord next amount down */
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))     /* keep track of how far to right tip is */
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do { 
    xx = q->v;
    if (xx > 100.0)     /* make sure tree doesn't stick out too far on line */
      xx = 100.0;
    pars_coordinates(q->back, lengthsum + xx, tipy, tipmax); /* recurse out */
    q = q->next;
  } while (p != q);
  first = p->next->back;             /* find immediate first descendant ,,, */
  q = p;
  while (q->next != p) {
    q = q->next;
    if (q->back != NULL)
      last = q->back;                                   /* ... and last one */
  }
  p->xcoord = (long)(over * lengthsum + 0.5);
  if ( (count_sibs(p) > 2) || ((p == curtree->root) 
                               && (p->next->next->back != NULL)) ) {
      p->ycoord = p->next->next->back->ycoord;
  }
  else   /* y coordinate is halfway between that of first, last descendants */
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;                   /* y coordinates of leftmost ... */
  p->ymax = last->ymax;   /* ... and rightmost tips descended from this one */
}  /* pars_coordinates */


void pars_printree(void)
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
  pars_coordinates(curtree->root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline3(i, scale, curtree->root);
  putc('\n', outfile);
}  /* pars_printree */


void maketree(void)
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, nextnode, outCount, missedCount;
  boolean firsttree, goteof, haslengths;

  // RSGnote: This was formerly uninitialized and potentially referenced before being set
  // below, depending on which branch of the first IF below was taken.
  long numtrees = 0;

  if (!usertree)
  {                            /* if sequentially adding, rearranging trees */
    lastrearr = false;
    hsbut(curtree, bestree, priortree, false, jumble, jumb, seed,
           progress, &bestfound);      /* call adding and local rearranging */

    if (progress)
    {      /* announce doing SPR rearrangements, print head of progress bar */
      sprintf(progbuf, "\nDoing global rearrangements");
      print_progress(progbuf);
      if (rearrfirst)
        sprintf(progbuf, " on the first of the trees tied for best\n");
      else
        sprintf(progbuf, " on all trees tied for best\n");
      print_progress(progbuf);
      sprintf(progbuf, "  !");
      print_progress(progbuf);
      for (j = 0; j < nonodes; j++)
      {
        if (j % ((nonodes / 72) + 1) == 0)
        {
          sprintf(progbuf, "-");
          print_progress(progbuf);
        }
      }
      sprintf(progbuf, "!\n");
      print_progress(progbuf);
    }

    phyFillScreenColor();
    grandrearr(curtree, bestree, progress, rearrfirst, &bestfound);
        /* call function doing collapsing of best trees, SPR rearrangements */
    if (progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
      phyFillScreenColor();
    }
    recompute = false;
    if (jumb == njumble)
    {
      missedCount = 0;
      outCount = nextree;
      collapsebestrees(curtree, bestrees, place,
                        chars, progress, &outCount);
      missedCount = nextree - 1 - maxtrees;
      if (treeprint)
      {
        putc('\n', outfile);
        if (outCount == 0)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
        {
          if (missedCount > 0)
          {
            fprintf(outfile, "as many as %ld trees may have been found\n",
                               missedCount + outCount + 1);
            fprintf(outfile, "here are the first %4ld of them\n",
                               outCount + 1 );
          }
          else
          {
            fprintf(outfile, "%6ld trees in all found\n", outCount);
          }
        }
      }
      if (treeprint)
        putc('\n', outfile);
      for (i = 0; i < outCount ; i++)           /* print out the best trees */
      {
        load_tree(curtree, i, bestrees);
        curtree->root = root_tree(curtree, curtree->root);
        initializetrav(curtree, curtree->root);
        initializetrav(curtree, curtree->root->back);
        curtree->evaluate(curtree, curtree->root, false);
        curtree->root = root_tree(curtree, curtree->root);
        curtree->nodep[curtree->root->index - 1] = curtree->root;
        disc_treelength(curtree->root, chars, curtree->nodep);
        pars_printree();
        describe();
      }
    }
  }
  else
  {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if (numtrees > MAXNUMTREES)
    {
      printf("\nERROR:  Number of input trees is read incorrectly from %s.\n\n", intreename);
      exxit(-1);
    }
    if ( reusertree )
      rebestrees = allocbestrees();
    if (numtrees > 2 && !reusertree)
      initseed(&inseed, &inseed0, seed);
    if (treeprint)
    {
      if (!reusertree)
        fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    if (!reusertree)       /* if not rearranging the trees that are read in */
    {
      fsteps = (long **)Malloc(maxuser * sizeof(long *));
      for (j = 1; j <= maxuser; j++)
        fsteps[j - 1] = (long *)Malloc(endsite * sizeof(long));
    }

#if 0                                   // RSGbugfix
    nodep = NULL;
#endif

    which = 1;
    if ( reusertree)
      rebestyet = UNDEFINED;
    while (which <= numtrees)
    {
      firsttree = true;
      nextnode = 0;
      haslengths = true;
/* debug:      preparetree(curtree);     need it?  */
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof, &firsttree, &nextnode, &haslengths, initparsnode, false, nonodes);
      fixtree(curtree);
      reroot_tree(curtree);
      curtree->evaluate(curtree, curtree->root, false);

      if ( reusertree )                      /* if rearrange the user trees */
      {                          /* but don't print them out till done that */
        rebestyet = curtree->score;
        bestrees = rebestrees[0];
        grandrearr(curtree, bestree, progress, rearrfirst, &bestfound);
        which++;
      }
      else
      {                                         /* print out this user tree */
        curtree->root = root_tree(curtree, curtree->root);
        if (treeprint)
          fprintf(outfile, "\n\n");
        if (outgropt)
          reroot(curtree->nodep[outgrno - 1], curtree->root);
        disc_treelength(curtree->root, chars, curtree->nodep);
        pars_printree();
        describe();
        which++;
      }
    }
    if ( reusertree )
    {        /* if user trees were rearranged, now, finally, print them out */
      for (i = 0; i <= (nextree - 2); i++)
      {
        load_tree(curtree, i, bestrees);
        curtree->evaluate(curtree, curtree->root, 0);
        curtree->root = root_tree(curtree, curtree->root);
        disc_treelength(curtree->root, chars, curtree->nodep);
        pars_printree();
        describe();
      }
    }

    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1  && !reusertree)     /* do KHT or SH Test */
      standev(chars, numtrees, minwhich, minsteps, nsteps, fsteps, seed);
    if ( !reusertree)
    {
      for (j = 1; j <= maxuser; j++)
        free(fsteps[j - 1]);
      free(fsteps);
    }
  }

  if (jumb == njumble)        /* if at end of last jumble or if no jumbling */
  {
    if (progress)    /* announce that the trees were written to output file */
    {
      sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
      print_progress(progbuf);
      if (trout)          /* ... and if needed, to the output tree file too */
      {
        sprintf(progbuf, "Tree");
        print_progress(progbuf);
        // RSGnote: If IF branch of first conditional is taken, "numtrees" would have been uninitialized here.
        if ((usertree && numtrees > 1) || (!usertree && nextree != 1))
        {
          sprintf(progbuf, "s");
          print_progress(progbuf);
        }
        sprintf(progbuf, " also written onto file \"%s\".\n\n", outtreename);
        print_progress(progbuf);
      }
    }
  }
}  /* maketree */


void freerest(void)
{
/* debug    free(threshwt);    for some reasons blows up so commented out, as in Dnapars */
}  /* freerest*/


void parsrun(void)
{
  /* run the inference of trees, for all data sets and all addition orders
   * of species */

  /*
    printf("jumble: %i\n", jumble);
    printf("njumble: %li\n", njumble);
    printf("outgrno: %li\n", outgrno);
    printf("outgropt: %i\n", outgropt);
    printf("thresh: %i\n", thresh);
    printf("thorough: %i\n", thorough);
    printf("rearrfirst: %i\n", rearrfirst);
    printf("maxtrees: %li\n", maxtrees);
    printf("trout: %i\n", trout);
    printf("usertree: %i\n", usertree);
    printf("reusertree: %i\n", reusertree);
    printf("weights: %i\n", weights);
    printf("mulsets: %i\n", mulsets);
    printf("printdata: %i\n", printdata);
    printf("progress: %i\n", progress);
    printf("treeprint: %i\n", treeprint);
    printf("stepbox: %i\n", stepbox);
    printf("ancseq: %i\n", ancseq);
    printf("dotdiff: %i\n", dotdiff);
    printf("interleaved: %i\n", interleaved);
    printf("justwts: %i\n", justwts);
  */

  for (ith = 1; ith <= msets; ith++) {                /* for each data set */
    if (msets > 1 && !justwts) {
      fprintf(outfile, "\nData set # %ld:\n\n", ith);
      if (progress)
      {                        /* print notification of progress on screen */
        sprintf(progbuf, "\nData set # %ld:\n\n", ith);
        print_progress(progbuf);
      }
    }
    doinput();                /* get input and set up tips of tree with it */
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++) /* for jumbling addition order */
      maketree(); /* reconstruct tree for one order of addition of species */
    fflush(outfile);
    fflush(outtree);
    freerest();
  }
} /* parsrun */


void pars(
  char *infilename,
  char *intreename,
  char *weightsfilename,
  char *OutfileName,
  char *outfileopt,
  char *OuttreeName,
  char *outtreeopt,
  char *TreeSearchMethod,
  char *SearchOpts,
  int TreeSave,
  int InputOrder,
  int RandNum,
  int NumJumble,
  int OutRoot,
  int OutNum,
  int ThreshPars,
  double ThreshVal,
  int SitesWeighted,
  int AnalyzeMult,
  int MultDataset,
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
  /* function that uses data sent from Java interface */

  initdata *funcs;
  //printf("Hello from Pars!\n"); // JRMdebug
  //fflush(stdout);

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Pars";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = discretepars_node_new;
  funcs->tree_new = discretepars_tree_new;
  phylipinit(argc, argv, funcs, true);
  progname = argv[0];
  /*
  //jumble = false;
  //njumble = 1;
  //outgrno = 1;
  //outgropt = false;
  //thresh = false;
  //thorough = true;
  //rearrfirst = false;
  //maxtrees = 10000;
  //trout = true;
  //usertree = false;
  //reusertree = false;
  //weights = false;
  //mulsets = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  //stepbox = false;
  //ancseq = false;
  //dotdiff = true;
  //interleaved = true;
  char *infile,
  char *intree,
  char *weightfile,
  char *outfile,
  char *outfileopt,
  char *outtree,
  char *outtreeopt,
  //char *TreeSearchMethod,
  //char *SearchOpts,
  //int TreeSave,
  //int InputOrder,
  //int RandNum,
  //int NumJumble,
  //int OutRoot,
  //int OutNum,
  //int ThreshPars,
  //double ThreshVal,
  //int SitesWeighted,
  //int AnalyzeMult,
  //int MultDataset,
  //int NumMult,
  //int InputSeq,
  //int PrintData,
  //int DotDiff,
  //int PrintInd,
  //int PrintTree,
  //int PrintSteps,
  //int PrintSeq,
  //int WriteTree)
  */

  if (!strcmp(TreeSearchMethod, "Best"))
  {
    usertree = false;
    reusertree = false;
  }
  else if (!strcmp(TreeSearchMethod, "User"))
  {
    usertree = true;
    reusertree = false;
  }
  else // Rearrange
  {
    usertree = true;
    reusertree = true;
  }

  if (!strcmp(SearchOpts, "More"))
  {
    thorough   = true;
    rearrfirst = false;
  }
  else if (!strcmp(SearchOpts, "Less"))
  {
    thorough   = false;
    rearrfirst = false;
  }
  else // Rearrange
  {
    thorough   = false;
    rearrfirst = true;
  }

  maxtrees = TreeSave;

  if (InputOrder != 0)
  {
    jumble = true;
    njumble = NumJumble;
  }
  else
  {
    jumble = false;
    njumble = 1;
  }

  inseed =  RandNum;

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

  if (ThreshPars != 0)
  {
    thresh = true;
    threshold = ThreshVal;
  }
  else
  {
    thresh = false;
    threshold = spp;
  }

  if (SitesWeighted != 0)
  {
    weights = true;
  }
  else
  {
    weights = false;
  }

  if (AnalyzeMult !=  0)
  {
    mulsets = true;
    msets = NumMult;
  }
  else
  {
    mulsets = false;
    msets = 1;
  }

  if (MultDataset != 0)
  {
    justwts = false;
  }
  else
  {
    justwts = true;
  }

  if (InputSeq != 0)
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
  //printf("Data in\n"); // JRMdebug
  //fflush(stdout);

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
    intree = fopen(intreename, "rb");
  }

  firstset = true;
  //printf("calling doinit\n"); // JRMdebug
  //fflush(stdout);
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
  //printf("files opened\n"); // JRMdebug
  //fflush(stdout);

  if (jumble) {
    initjumble(&inseed, &inseed0, seed, &njumble);
  }

  //printf("calling parsrun\n"); // JRMdebug
  //fflush(stdout);

  parsrun();  // do the actual work

  FClose(infile);
  if (usertree)
  {
    FClose(intree);
  }
  if (weights || justwts)
  {
    FClose(weightfile);
  }
  FClose(outfile);
  if (trout)
  {
    FClose(outtree);
  }
  //printf("\ndone\n"); // JRMdebug
} /* pars */


int main(int argc, Char *argv[])
{  /* Discrete character parsimony by uphill search */

  /* reads in spp, chars, and data. Calls maketree to construct the tree */
  initdata *funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Pars", "");                */
  argv[0]="Pars";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = discretepars_node_new;
  funcs->tree_new = discretepars_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  msets = 1;
  firstset = true;
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  parsrun();

  FClose(infile);
  FClose(outfile);
  if (weights || justwts) {
    FClose(weightfile);
  }
  if (trout) {
    FClose(outtree);
  }
  if (usertree) {
    FClose(intree);
  }
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  if (progress)
    printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* Discrete character parsimony by uphill search */


// End.
