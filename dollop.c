/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "disc.h"
#include "dollo.h"

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   doinput(void);
void   dollop_count(node *, steptr, steptr);
void   preorder(node *, steptr, steptr, long, boolean, long, bitptr, pointptr);
void   evaluate(node *);
void   savetree(void);
void   dollop_addtree(long *);

void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, node **, boolean *);
void   repreorder(node *, node **, boolean *);
void   rearrange(node **);
void   describe(void);
void   initdollopnode(tree *, node **, long, long, long *,
                      long *, initops, pointarray, Char *, Char *, FILE *);
void   maketree(void);
void   reallocchars(void);
void dolloprun(void);
void dollop(char *infilename, char *intreename, char *weightsfilename, char *ancfilename, char *outfilename,
            char *outfileopt, char *outtreename, char *outtreeopt, int BestTree, int UseDollo, int TreeMax,
            int InputOrder, int RandNum, int NumJumble, int ThreshDollop, double ThreshVal, int UseAncStates,
            int SitesWeighted, int AnalyzeMult, int MultDataset, int NumMult, int PrintData, int PrintInd,
            int PrintTree, int PrintSteps, int PrintSeq, int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH], ancfilename[FNMLNGTH];
tree * curtree;
long col, msets, ith, j, l, njumble, jumb = 0, nonodes = 0;
long maxtrees;
long inseed, inseed0;
boolean jumble, usertree, weights, thresh, ancvar, questions, dollo, trout,  progress, treeprint, stepbox, ancseq, mulsets, firstset, justwts;
boolean *ancone, *anczero, *ancone0, *anczero0;
double threshold;
double *threshwt;
longer seed;
long *enterorder;
double **fsteps;
steptr numsteps;
bestelm *bestrees;
Char *guess;
gbit *garbage;
char *progname;

/* Variables for treeread */
boolean goteof, firsttree, haslengths, phirst;
#if 0                                   // RSGbugfix: Global variable never used.
pointarray nodep;
#endif

/* Local variables for maketree, propagated globally for C version: */
long minwhich;
double like, bestyet, bestlike, bstlike2, minsteps;
boolean lastrearr;
double nsteps[maxuser];
node *there;
long fullset;
bitptr zeroanc, oneanc;
long *place;
Char ch;
boolean *names;
steptr numsone, numszero;
bitptr steps;


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  putchar('\n');
  maxtrees = 10000;
  ancvar = false;
  dollo = true;
  jumble = false;
  njumble = 1;
  thresh = false;
  threshold = 2 * spp;
  trout = true;
  usertree = false;
  goteof = false;
  weights = false;
  justwts = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  loopcount = 0;
  for (;;)
  {
    cleerhome();
    printf("\nDollo and polymorphism parsimony algorithm, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    printf("  P                     Parsimony method?  %s\n",
           dollo ? "Dollo" : "Polymorphism");
    if (!usertree)
    {
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
           ancvar ? "Yes" : "No");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           printdata ? "Yes" : "No");
    printf("  2  Print indications of progress of run  %s\n",
           progress ? "Yes" : "No");
    printf("  3                        Print out tree  %s\n",
           treeprint ? "Yes" : "No");
    printf("  4     Print out steps in each character  %s\n",
           stepbox ? "Yes" : "No");
    printf("  5     Print states at all nodes of tree  %s\n",
           ancseq ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
           trout ? "Yes" : "No");
    if(weights && justwts)
    {
      printf("WARNING:  W option and Multiple Weights options are both on.  ");
      printf("The W menu option is unnecessary and has no additional effect. \n");
    }
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WAPJTUM1234560", ch) != NULL))
        || (usertree && ((strchr("WAPTUM1234560", ch) != NULL))))
    {
      switch (ch)
      {
        case 'A':
          ancvar = !ancvar;
          break;

        case 'P':
          dollo = !dollo;
          break;

        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
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
{
  long i;

  free(extras);
  free(weight);
  free(threshwt);
  free(numsteps);
  free(ancone);
  free(anczero);
  free(ancone0);
  free(anczero0);
  free(numsone);
  free(numszero);
  free(guess);

  if (usertree)
  {
    for (i = 1; i <= maxuser; i++)
    {
      free(fsteps);
      fsteps[i - 1] = (double *)Malloc(chars * sizeof(double));
    }
  }

  extras = (steptr)Malloc(chars * sizeof(long));
  weight = (steptr)Malloc(chars * sizeof(long));
  threshwt = (double *)Malloc(chars * sizeof(double));
  numsteps = (steptr)Malloc(chars * sizeof(long));
  ancone = (boolean *)Malloc(chars * sizeof(boolean));
  anczero = (boolean *)Malloc(chars * sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars * sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars * sizeof(boolean));
  numsone = (steptr)Malloc(chars * sizeof(long));
  numszero = (steptr)Malloc(chars * sizeof(long));
  guess = (Char *)Malloc(chars * sizeof(Char));
}


void allocrest(void)
{
  long i;

  extras = (steptr)Malloc(chars * sizeof(long));
  weight = (steptr)Malloc(chars * sizeof(long));
  threshwt = (double *)Malloc(chars * sizeof(double));
  if (usertree)
  {
    fsteps = (double **)Malloc(maxuser * sizeof(double *));
    for (i = 1; i <= maxuser; i++)
      fsteps[i - 1] = (double *)Malloc(chars * sizeof(double));
  }
  bestrees = (bestelm *) Malloc(maxtrees * sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes * sizeof(long));
  numsteps = (steptr)Malloc(chars * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
  place = (long *)Malloc(nonodes * sizeof(long));
  ancone = (boolean *)Malloc(chars * sizeof(boolean));
  anczero = (boolean *)Malloc(chars * sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars * sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars * sizeof(boolean));
  numsone = (steptr)Malloc(chars * sizeof(long));
  numszero = (steptr)Malloc(chars * sizeof(long));
  guess = (Char *)Malloc(chars * sizeof(Char));
  zeroanc = (bitptr)Malloc(words * sizeof(long));
  oneanc = (bitptr)Malloc(words * sizeof(long));
  steps = (bitptr)Malloc(words * sizeof(long));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */
  inputnumbers(&spp, &chars, &nonodes, 1);
  words = chars / bits + 1;
  fprintf(outfile, "\nDollo and polymorphism parsimony algorithm,");
  fprintf(outfile, " version %s\n\n", VERSION);
  if(javarun)
  {
    if(thresh == false)
    {
      threshold = spp;
    }
  }
  else
  {
    getoptions();
  }
  curtree = functions.tree_new(nonodes, spp);
  alloctree(&(curtree->nodep));
  setuptree(curtree->nodep);
  allocrest();
}  /* doinit */


void inputoptions(void)
{
  /* input the information on the options */
  long i;
  if(justwts)
  {
    if(firstset)
    {
      scan_eoln(infile);
      if (ancvar)
      {
        inputancestors(anczero0, ancone0);
      }
    }
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    inputweights(chars, weight, &weights);
  }
  else
  {
    if (!firstset)
    {
      samenumsp(&chars, ith);
      reallocchars();
    }
    scan_eoln(infile);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (ancvar)
      inputancestors(anczero0, ancone0);
    if (weights)
      inputweights(chars, weight, &weights);
  }
  if ((weights || justwts) && printdata)
    printweights(outfile, 0, chars, weight, "Characters");
  for (i = 0; i < chars; i++)
  {
    if (!ancvar)
    {
      anczero[i] = true;
      ancone[i] = false;
    }
    else
    {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && printdata)
    printancestors(outfile, anczero, ancone);
  questions = false;
  for (i = 0; i < chars; i++)
  {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */


void doinput(void)
{
  /* reads the input data */
  inputoptions();
  if(!justwts || firstset)
    inputdata(curtree->nodep, dollo, printdata, outfile);
}  /* doinput */


void dollop_count(node *p, steptr numsone, steptr numszero)
{
  /* Counts the number of steps in a fork of the tree.
     The program spends much of its time in this procedure. */
  long i, j, l;

  if (dollo)
  {
    for (i = 0; i < words; i++)
      steps[i] = (curtree->nodep[p->back->index - 1]->stateone[i] & p->statezero[i] & zeroanc[i]) | (curtree->nodep[p->back->index - 1]->statezero[i] & p->stateone[i] & fullset & (~zeroanc[i]));
  }
  else
  {
    for (i = 0; i < words; i++)
      steps[i] = curtree->nodep[p->back->index - 1]->stateone[i] & curtree->nodep[p->back->index - 1]->statezero[i] & p->stateone[i] & p->statezero[i];
  }
  j = 1;
  l = 0;
  for (i = 0; i < chars; i++)
  {
    l++;
    if (l > bits)
    {
      l = 1;
      j++;
    }
    if (((1L << l) & steps[j - 1]) != 0)
    {
      assert(j <= words);   /* checking array indexing */
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
}  /* dollop_count */


void preorder(node *p, steptr numsone, steptr numszero, long words, boolean dollo, long fullset, bitptr zeroanc, pointptr treenode)
{
  /* go back up tree setting up and counting interior node
     states */

  if (!p->tip)
  {
    correct(p, fullset, dollo, zeroanc, treenode);
    preorder(p->next->back, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
    preorder(p->next->next->back, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
  }
  if (p->back != NULL)
    dollop_count(p, numsone, numszero);
}  /* preorder */


void evaluate(node *r)
{
  /* Determines the number of losses or polymorphisms needed
     for a tree. This is the minimum number needed to evolve
     chars on this tree */
  long i, stepnum, smaller;
  double sum, term;

  sum = 0.0;
  for (i = 0; i < chars; i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < words; i++)
    zeroanc[i] = fullset;
  postorder(r);
  preorder(r, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
  for (i = 0; i < words; i++)
    zeroanc[i] = 0;
  postorder(r);
  preorder(r, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
  for (i = 0; i < chars; i++) {
    smaller = 2 * spp * weight[i];
    numsteps[i] = smaller;
    if (anczero[i]) {
      numsteps[i] = numszero[i];
      smaller = numszero[i];
    }
    if (ancone[i] && (numsone[i] < smaller))
      numsteps[i] = numsone[i];
    stepnum = numsteps[i] + extras[i];
    if (stepnum <= threshwt[i])
      term = stepnum;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && (which <= maxuser))
      fsteps[which - 1][i] = term;
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && (numszero[i] < numsone[i])))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && (numsone[i] < numszero[i])))
      guess[i] = '1';
  }
  if (usertree && (which <= maxuser)) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */


void savetree(void)
{
  /* record in place where each species has to be
     added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[curtree->root->index - 1] = 1;
  for (i = 1; i <= spp; i++)
  {
    p = curtree->nodep[i - 1];
    while (place[p->index - 1] == 0)
    {
      place[p->index - 1] = i;
      p = p->back;
      if (p != NULL)
        p = curtree->nodep[p->index - 1];
    }
    if (i > 1)
    {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done)
      {
        place[p->index - 1] = spp + i - 1;
        p = curtree->nodep[p->index - 1];
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* savetree */


void dollop_addtree(long *pos)
{
  /*puts tree from ARRAY place in its proper position
    in ARRAY bestrees */
  long i;
  for (i =nextree - 1; i >= (*pos); i--)
  {
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp * sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < spp; i++)
    bestrees[(*pos) - 1].btree[i] = place[i];
  nextree++;
}  /* dollop_addtree */


void tryadd(node *p, node **item, node **nufork)
{
  /* Temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  long pos;
  boolean found;

  add(p, *item, *nufork, &(curtree->root), curtree->nodep);
  evaluate(curtree->root);
  if (lastrearr)
  {
    if (like >= bstlike2)
    {
      savetree();
      if (like > bstlike2)
      {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        dollop_addtree(&pos);
      }
      else
      {
        pos = 0;
        findtree(&found, &pos, nextree, place, bestrees);
        /* findtree calls for a bestelm* but is getting */
        /* a long**, LM                                 */
        if (!found)
        {
          if (nextree <= maxtrees)
            dollop_addtree(&pos);
        }
      }
    }
  }
  if (like > bestyet)
  {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork, &(curtree->root), curtree->nodep);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *item= item_;
  node *nufork = nufork_;

  if (p == NULL)
    return;
  tryadd(p, &item, &nufork);
  if (!p->tip)
  {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, node **r, boolean *success)
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = curtree->nodep[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  re_move(&p, &forknode, &(curtree->root), curtree->nodep);
  add(whereto, p, forknode, &(curtree->root), curtree->nodep);
  evaluate(*r);
  if (like <= oldlike)
  {
    re_move(&p, &forknode, &(curtree->root), curtree->nodep);
    add(frombelow, p, forknode, &(curtree->root), curtree->nodep);
  }
  else
  {
    (*success) = true;
    bestyet = like;
  }
}  /* tryrearr */


void repreorder(node *p, node **r, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, r, success);
  if (!p->tip)
  {
    repreorder(p->next->back, r, success);
    repreorder(p->next->next->back, r, success);
  }
}  /* repreorder */


void rearrange(node **r_)
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  node **r         = r_;
  boolean success  = true;

  while (success)
  {
    success = false;
    repreorder(*r, r, &success);
  }
}  /* rearrange */


void describe(void)
{
  /* prints ancestors, steps and table of numbers of steps in
     each character */

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", -like);
  if (stepbox)
  {
    putc('\n', outfile);
    writesteps(weights, dollo, numsteps);
  }
  if (questions)
    guesstates(guess);
  if (ancseq)
  {
    hypstates(fullset, dollo, guess, curtree->nodep, curtree->root, garbage, zeroanc, oneanc);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout)
  {
    col = 0;
    treeout(curtree->root, nextree, &col, curtree->root);
  }
}  /* describe */


void initdollopnode(tree *treep, node **p, long len, long nodei, long *ntips, long *parens, initops whichinit, pointarray treenode, Char *str, Char *ch, FILE *intree)
{
  /* initializes a node */
  /* LM 7/27  I added this function and the commented lines around  */
  /* treeread() to get the program running, but all 4 move programs */
  /* are improperly integrated into the v4.0 support files.  As is  */
  /* this is a patchwork function                                   */
  boolean minusread;
  double valyew, divisor;

  (void)len;                            // RSGnote: Parameter never used.
  (void)ntips;                          // RSGnote: Parameter never used.
  (void)treenode;                       // RSGnote: Parameter never used.

  switch (whichinit)
  {
    case bottom:
      treep->nodep[nodei - 1] = *p;
      break;
    case nonbottom:
      break;
    case tip:
      match_names_to_data (str, treep->nodep, p, spp);
      break;
    case length:         /* if there is a length, read it and discard value */
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      break;
    default:      /* cases hslength, hsnolength, treewt, unittrwt, iter, */
      break;
  }
} /* initdollopnode */


void maketree(void)
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees, nextnode;
  double gotlike;
  node *item, *nufork, *dummy, *p;

  fullset = (1L << (bits + 1)) - (1L << 1);
  if (!usertree)
  {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree->root = curtree->nodep[enterorder[0] - 1];
    add(curtree->nodep[enterorder[0] - 1], curtree->nodep[enterorder[1] - 1], curtree->nodep[spp], &curtree->root, curtree->nodep);
    if (progress)
    {
      sprintf(progbuf, "Adding species:\n");
      print_progress(progbuf);
      writename(0, 2, enterorder);
      phyFillScreenColor();
    }
    lastrearr = false;
    for (i = 3; i <= spp; i++)
    {
      bestyet = -350.0 * spp * chars;
      item = curtree->nodep[enterorder[i - 1] - 1];
      nufork = curtree->nodep[spp + i - 2];
      addpreorder(curtree->root, item, nufork);
      add(there, item, nufork, &curtree->root, curtree->nodep);
      like = bestyet;
      rearrange(&curtree->root);
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
        do {
          if (progress)
          {
            sprintf(progbuf, "   ");
            print_progress(progbuf);
          }
          gotlike = bestlike;
          for (j = 0; j < nonodes; j++)
          {
            bestyet = - 350.0 * spp * chars;
            item = curtree->nodep[j];
            if (item != curtree->root)
            {
              nufork = curtree->nodep[j]->back;
              re_move(&item, &nufork, &curtree->root, curtree->nodep);
              there = curtree->root;
              addpreorder(curtree->root, item, nufork);
              add(there, item, nufork, &curtree->root, curtree->nodep);
            }
            if (progress)
            {
              if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              {
                sprintf(progbuf, ".");
                print_progress(progbuf);
              }
            }
          }
          if (progress)
          {
            sprintf(progbuf, "\n");
            print_progress(progbuf);
            phyFillScreenColor();
          }
        } while (bestlike > gotlike);
      }
    }
    if (progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
    for (i = spp - 1; i >= 1; i--)
      re_move(&(curtree->nodep[i]), &dummy, &curtree->root, curtree->nodep);
    if (jumb == njumble)
    {
      if (treeprint)
      {
        putc('\n', outfile);
        if (nextree == 2)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1)
      {
        if (treeprint)
          fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      for (i = 0; i <= (nextree - 2); i++)
      {
        curtree->root = curtree->nodep[0];
        add(curtree->nodep[0], curtree->nodep[1], curtree->nodep[spp], &curtree->root, curtree->nodep);
        for (j = 3; j <= spp; j++)
        {
          add(curtree->nodep[bestrees[i].btree[j - 1] - 1], curtree->nodep[j - 1], curtree->nodep[spp + j - 2], &curtree->root, curtree->nodep);}
        evaluate(curtree->root);
        printree(1.0, treeprint, curtree->root);
        describe();
        for (j = 1; j < spp; j++)
          re_move(&(curtree->nodep[j]), &dummy, &curtree->root, curtree->nodep);
      }
    }
  }
  else
  {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(intree);
    if (numtrees > MAXNUMTREES)
    {
      printf("\nERROR:  Number of input trees is read incorrectly from %s.\n\n", intreename);
      exxit(-1);
    }
    if (numtrees > 2)
    {
      initseed(&inseed, &inseed0, seed);
      printf("\n");
    }
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    names = (boolean *)Malloc(spp * sizeof(boolean));
    which = 1;
    firsttree = true;                       /**/
#if 0 // RSGbugfix: Global variable never used.
    nodep = NULL;                           /**/
#endif
    nextnode = 0;                           /**/
    haslengths = 0;                         /**/
    phirst = 0;                             /**/
    while (which <= numtrees)
    {
      treeread(curtree, intree, &(curtree->root), curtree->nodep, &goteof, &firsttree, &nextnode, &haslengths, initdollopnode, false, nonodes);
      for (i = spp; i < nonodes; i++)
      {
        p = curtree->nodep[i];
        for (j = 1; j <= 3; j++)
        {
          p->stateone = (bitptr)Malloc(words * sizeof(long));
          p->statezero = (bitptr)Malloc(words * sizeof(long));
          p = p->next;
        }
      } /* debug: see comment at initdollopnode() */
      if (treeprint)
        fprintf(outfile, "\n\n");
      evaluate(curtree->root);
      printree(1.0, treeprint, curtree->root);
      describe();
      which++;
    }
    FClose(intree);
    fprintf(outfile, "\n");
    if (numtrees > 1 && chars > 1)
      standev(numtrees, minwhich, minsteps, nsteps, fsteps, seed);
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
    }
    if (ancseq)
      freegarbage(&garbage);
  }
}  /* maketree */


void dolloprun(void)
{
  // debug printout // JRMdebug
  /*
    printf("maxtrees: %li\n", maxtrees);
    printf("ancvar: %i\n", ancvar);
    printf("dollo: %i\n", dollo);
    printf("jumble: %i\n", jumble);
    printf("njumble: %li\n", njumble);
    printf("thresh: %i\n", thresh);
    printf("threshold: %f\n", threshold);
    printf("trout: %i\n", trout);
    printf("usertree: %i\n", usertree);
    printf("goteof: %i\n", goteof);
    printf("weights: %i\n", weights);
    printf("justwts: %i\n", justwts);
    printf("mulsets: %i\n", mulsets);
    printf("printdata: %i\n", printdata);
    printf("progress: %i\n", progress);
    printf("treeprint: %i\n", treeprint);
    printf("stepbox: %i\n", stepbox);
    printf("ancseq: %i\n", ancseq);
    fflush(stdout);
  */
  // do the work
  if (dollo)
    fprintf(outfile, "Dollo");
  else
    fprintf(outfile, "Polymorphism");
  fprintf(outfile, " parsimony method\n\n");
  if (printdata && justwts)
    fprintf(outfile, "%2ld species, %3ld  characters\n\n", spp, chars);

  for (ith = 1; ith <= msets; ith++) {
    if (msets > 1 && !justwts) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    if (justwts) {
      fprintf(outfile, "Weights set # %ld:\n\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nWeights set # %ld:\n\n", ith);
        print_progress(progbuf);
      }
    }
    if (printdata && !justwts)
      fprintf(outfile, "%2ld species, %3ld  characters\n\n", spp, chars);
    doinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    fflush(outfile);
    fflush(outtree);
  }
}


void dollop(
  char *infilename,
  char *intreename,
  char *weightsfilename,
  char *ancfilename,
  char *OutfileName,
  char *outfileopt,
  char *OuttreeName,
  char *outtreeopt,
  int BestTree,
  int UseDollo,
  int TreeMax,
  int InputOrder,
  int RandNum,
  int NumJumble,
  int ThreshDollop,
  double ThreshVal,
  int UseAncStates,
  int SitesWeighted,
  int AnalyzeMult,
  int MultDataset,
  int NumMult,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int PrintSteps,
  int PrintSeq,
  int WriteTree)
{
  //printf("Hello from Dollop!\n"); // JRMdebug
  //fflush(stdout);

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Dollop";
  phylipinit(argc, argv, NULL, true);
  progname = argv[0];
  /*
  //maxtrees = 10000;
  //ancvar = false;
  //dollo = true;
  //jumble = false;
  //njumble = 1;
  //thresh = false;
  //threshold = spp;
  //trout = true;
  //usertree = false;
  //goteof = false;
  //weights = false;
  //justwts = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  //stepbox = false;
  //ancseq = false;
  char *infile,
  char *intree,
  char *weightfile,
  char *ancfile,
  char *outfile,
  char *outfileopt,
  char *outtree,
  char *outtreeopt,
  //int BestTree,
  //int UseDollo,
  //int TreeMax,
  //int InputOrder,
  //int RandNum,
  //int NumJumble,
  //int ThreshDollop,
  //double ThreshVal,
  //int UseAncStates,
  //int SitesWeighted,
  //int AnalyzeMult,
  //int MultDataset,
  //int NumMult,
  //int PrintData,
  //int PrintInd,
  //int PrintTree,
  //int PrintSteps,
  //int PrintSeq,
  //int WriteTree)

  */
  if (BestTree != 0)
  {
    usertree = false;
  }
  else
  {
    usertree = true;
  }

  if (UseDollo != 0)
  {
    dollo = true;
  }
  else
  {
    dollo = false;
  }

  maxtrees = TreeMax;

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

  if (ThreshDollop != 0)
  {
    thresh = true;
    threshold = ThreshVal;
  }
  else
  {
    thresh = false;
    threshold = spp;
  }

  if (UseAncStates != 0)
  {
    ancvar = true;
  }
  else
  {
    ancvar = false;
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
    trout = true;
  }
  else
  {
    trout = false;
  }

  goteof = false;

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
    intree = fopen(intreename, "r");
  }

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  firstset = true;
  bits = 8 * sizeof(long) - 1;
  doinit();

  if (weights || justwts)
  {
    weightfile = fopen(weightsfilename, "r");
  }

  if(ancvar)
  {
    ancfile = fopen(ancfilename, "r");
  }

  if (trout)
  {
    outtree = fopen(OuttreeName, outtreeopt);
    strcpy(outtreename, OuttreeName);
  }

  if (jumble) {
    initjumble(&inseed, &inseed0, seed, &njumble);
  }
  //printf("calling dolloprun\n"); // JRMdebug
  //fflush(stdout);

  dolloprun();

  FClose(infile);
  if (weights || justwts)
  {
    FClose(weightfile);
  }
  if(ancvar)
  {
    FClose(ancfile);
  }
  FClose(outfile);
  if (trout)
  {
    FClose(outtree);
  }

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* Dollo or polymorphism parsimony by uphill search */
#ifdef MAC
  argc = 1;           /* macsetup("Dollop", "");  */
  argv[0] = "Dollop";
#endif
  phylipinit(argc, argv, NULL, false);
  /* reads in spp, chars, and the data. Then calls maketree to construct the tree */
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  mulsets = false;
  msets = 1;
  firstset = true;
  bits = 8 * sizeof(long) - 1;
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if(ancvar)
    openfile(&ancfile, ANCFILE, "ancestors file", "r", argv[0], ancfilename);

  dolloprun();

  /* this would be an appropriate place to deallocate memory, including these items:
     free(steps);
  */
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
  printf("\nDone.\n\n");

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  phyRestoreConsoleAttributes();
  return 0;
}  /* Dollo or polymorphism parsimony by uphill search */


// End.
