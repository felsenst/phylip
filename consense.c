/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Hisashi Horino,
   Akiko Fuseki, Dan Fineman, Sean Lamont, and Andrew Keeffe.
   Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cons.h"

/* The following extern's refer to things declared in cons.c */

extern int tree_pairing;

extern Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
extern node *root;

extern long numopts, outgrno, col;
extern long maxgrp;               /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets, progress, treeprint, goteof, strict, mr, mre, ml;
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern long tipy;

extern double trweight, ntrees, mlfrac;

#ifndef OLDC
/* function prototypes */
void   doinit(void);
void   getoptions(void);

#if 0                                   // RSGbugfix: This function seems to do nothing.
void   count_siblings(node **p);
#endif

void   treeout(node *);
void   consenserun(void);
void   consense(char *intreename, char *outfilename, char *outfileopt, char *outtreename, char *outtreeopt,
                char *ConsType, double Fraction, int OutRoot, int OutNum, int Rooted, int PrintData,
                int PrintInd, int PrintTree, int WriteTree);
/* function prototypes */
#endif


void doinit(void)
{
  fprintf(outfile, "\nConsensus tree");
  fprintf(outfile, " program, version %s\n\n", VERSION);

  /* Initial settings */
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  firsttree      = true;
  spp            = 0 ;
  col            = 0 ;

  /* This is needed so functions in cons.c work */
  tree_pairing   = NO_PAIRING ;

  if (!javarun)
  {
    getoptions();
  }
}


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;
  boolean done, done1;

  putchar('\n');
  strict = false;
  mr = false;
  mre = true;
  ml = false;
  mlfrac = 0.5;
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  prntsets = true;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nConsensus tree");
    printf(" program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf(" C        Consensus type (MRe, strict, MR, Ml):");
    if (strict)
      printf("  strict\n");
    else if (mr)
      printf("  Majority rule\n");
    else if (mre)
      printf("  Majority rule (extended)\n");
    else if (ml)
      printf("  Ml\n");
    else printf("  Adams\n");
    if (noroot)
    {
      printf(" O                               Outgroup root:");
      if (outgropt)
        printf("  Yes, at species number%3ld\n", outgrno);
      else
        printf("  No, use as outgroup species%3ld\n", outgrno);
    }
    printf(" R               Trees to be treated as Rooted:");
    if (noroot)
      printf("  No\n");
    else
      printf("  Yes\n");
    printf(" T          Terminal type (IBM PC, ANSI, none):");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf(" 1               Print out the sets of species:");
    if (prntsets)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf(" 2        Print indications of progress of run:  %s\n",
           (progress ? "Yes" : "No"));
    printf(" 3                              Print out tree:");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf(" 4              Write out trees onto tree file:");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");

    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done)
    {
      if ((noroot && (ch == 'O')) || strchr("CRT1234", ch) != NULL)
      {
        switch (ch)
        {
          case 'C':
            if (strict)
            {
              strict = false;
              mr = true;
            }
            else
            {
              if (ml)
              {
                ml = false;
                mre = true;
              }
              else
              {
                if (mre)
                {
                  mre = false;
                  strict = true;
                }
                else
                {
                  if (mr)
                  {
                    mr = false;
                    ml = true;
                  }
                }
              }
            }
            break;

          case 'O':
            outgropt = !outgropt;
            if (outgropt)
            {
              numopts++;
              loopcount2 = 0;
              do {
                printf("Type number of the outgroup:\n");
                phyFillScreenColor();
                if(scanf("%ld%*[^\n]", &outgrno)) {} // Read number and scan to EOL.
                (void)getchar();
                done1 = (outgrno >= 1);
                if (!done1)
                {
                  printf("\nERROR:  Bad outgroup number: %ld.\n", outgrno);
                  printf("        Must be greater than zero.\n");
                }
                countup(&loopcount2, 10);
              } while (done1 != true);
            }
            break;

          case 'R':
            noroot = !noroot;
            break;

          case 'T':
            initterminal(&ibmpc, &ansi);
            break;

          case '1':
            prntsets = !prntsets;
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
  if (ml)
  {
    do {
      printf("\nFraction (l) of times a branch must appear\n");
      if(scanf("%lf%*[^\n]", &mlfrac)) {} // Read number and scan to EOL.
      (void)getchar();
    } while ((mlfrac < 0.5) || (mlfrac > 1.0));
  }
}  /* getoptions */


#if 0                        // RSGbugfix: This function seems to do nothing.

void count_siblings(node **p)
{
  node *tmp_node;
  int i;

  if (!(*p))
  {
    /* This is a leaf, */
    return;
  }
  else
  {
    tmp_node = (*p)->next;
  }

  for (i = 0 ; i < 1000; i++)
  {
    if (tmp_node == (*p))
    {
      /* When we've gone through all the siblings, */
      break;
    }
    else if (tmp_node)
    {
      tmp_node = tmp_node->next;
    }
    else
    {
      /* Should this be executed? */
      return ;
    }
  }
} /* count_siblings */

#endif


void treeout(node *p)
{
  /* write out file with representation of final tree */
  long i, n = 0;
  Char c;
  node *q;
  double x;

#if 0                                   // RSGbugfix: This function seems to do nothing.
  count_siblings (&p);
#endif

  if (p->tip)
  {
    /* If we're at a node which is a leaf, figure out how long the
       name is and print it out. */
    for (i = 1; i <= MAXNCH; i++)
    {
      if (p->nayme[i - 1] != '\0')
        n = i;
    }
    for (i = 0; i < n; i++)
    {
      c = p->nayme[i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  }
  else
  {
    /* If we're at a furcation, print out the proper formatting, loop
       through all the children, calling the procedure recursively. */
    putc('(', outtree);
    col++;
    q = p->next;
    while (q != p)
    {
      /* This should terminate when we've gone through all the
         siblings, */
      treeout(q->back);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      col++;
      if (col > 60)
      {
        putc('\n', outtree);
        col = 0;
      }
    }
    putc(')', outtree);
    col++;
  }

  if (p->tip)
    x = ntrees;
  else
    x = (double)p->deltav;

  if (p == root)
  {
    /* When we're all done with this tree, */
    fprintf(outtree, ";\n");
    return;
  }

  /* Figure out how many characters the branch length requires: */
  else
  {
    if (!strict)
    {
      if (x >= 100.0)
      {
        fprintf(outtree, ":%5.1f", x);
        col += 4;
      }
      else if (x >= 10.0)
      {
        fprintf(outtree, ":%4.1f", x);
        col += 3;
      }
      else if (x >= 1.0)
      {
        fprintf(outtree, ":%4.2f", x);
        col += 2;
      }
    }
  }
}  /* treeout */


void consenserun(void)
{
  pattern_elm  ***pattern_array;
  long trees_in = 0;
  long i, j;

  /*
  // debug printout // JRMdebug
  printf("strict: %i\n", strict);
  printf("mr: %i\n", mr);
  printf("mre: %i\n", mre);
  printf("ml: %i\n", ml);
  printf("mlfrac: %f\n", mlfrac);
  printf("noroot: %i\n", noroot);
  printf("numopts: %li\n", numopts);
  printf("outgrno: %li\n", outgrno);
  printf("outgropt: %i\n", outgropt);
  printf("trout: %i\n", trout);
  printf("prntsets: %i\n", prntsets);
  printf("progress: %i\n", progress);
  printf("treeprint: %i\n", treeprint);
  */

  // do the actual work
  if (prntsets)
    fprintf(outfile, "Species in order: \n\n");

  ntrees = 0.0;
  maxgrp = 32767;   /* initial size of set hash table */
  lasti  = -1;

  trees_in = countsemic(intree);

  /* Read the tree file and put together grouping, order, and timesseen */
  read_groups (&pattern_array, trees_in, trees_in, intree);
  /* Compute the consensus tree. */
  putc('\n', outfile);
  tree * curtree = functions.tree_new(2 * spp, spp);
  for (i = 0; i < spp; i++)
  {
    for (j = 0; j < MAXNCH; j++)
      curtree->nodep[i]->nayme[j] = '\0';
    strncpy(curtree->nodep[i]->nayme, nayme[i], MAXNCH);
  }
  consensus(curtree, pattern_array, trees_in);
  printf("\n");
  if (trout) {
    treeout(root);
    if (progress)
    {
      sprintf(progbuf, "Consensus tree written to file \"%s\".\n\n", outtreename);
      print_progress(progbuf);
    }
  }
  if (progress)
  {
    sprintf(progbuf, "Output written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
  }

#if 0
  // probably do need to free some structures here
  for (i = 0; i < spp; i++)
    free(nodep[i]);
  for (i = spp; i < 2*(1 + spp); i++)
  {
    if (nodep[i] != NULL)
    {
      p = nodep[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != nodep[i]);
      free(p);
    }
  }
  free(nodep);
#endif
}


void consense (
  char *intreename,
  char *OutfileName,
  char *outfileopt,
  char *OuttreeName,
  char *outtreeopt,
  char *ConsType,
  double Fraction,
  int OutRoot,
  int OutNum,
  int Rooted,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int WriteTree)
{
  initdata *funcs;
  //printf("Hello from Consense!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Consense";
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = cons_node_new;
  phylipinit(argc, argv, funcs, true);

  /*
  //strict = false;
  //mr = false;
  //mre = true;
  //ml = false;
  //mlfrac = 0.5;
  //noroot = true;
  numopts = 0; // this does not seem to do anything
  //outgrno = 1;
  //outgropt = false;
  //trout = true;
  //prntsets = true;
  //progress = true;
  //treeprint = true;
  //char *intree,
  //char *outfile,
  //char *outfileopt,
  //char *outtree,
  //char *outtreeopt,
  //char *ConsType,
  //double Fraction,
  //int OutRoot,
  //int OutNum,
  //int Rooted,
  //int PrintData,
  //int PrintInd,
  //int PrintTree,
  //int WriteTree)
  */

  mlfrac = 0.5;
  numopts = 0; // this does not seem to do anything

  if (!strcmp(ConsType, "extended"))
  {
    strict = false;
    mr = false;
    mre = true;
    ml = false;
  }
  else if (!strcmp(ConsType, "strict"))
  {
    strict = true;
    mr = false;
    mre = false;
    ml = false;
  }
  else if (!strcmp(ConsType, "majority"))
  {
    strict = false;
    mr = true;
    mre = false;
    ml = false;
  }
  else // msubl
  {
    strict = false;
    mr = false;
    mre = false;
    ml = true;
    mlfrac = Fraction;
  }

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

  if (Rooted != 0)
  {
    noroot = false;
  }
  else
  {
    noroot = true;
  }

  if (PrintData != 0)
  {
    prntsets = true;
  }
  else
  {
    prntsets = false;
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

  if (WriteTree != 0)
  {
    treeprint = true;
  }
  else
  {
    treeprint = false;
  }

  // everything translated, start the run
  intree = fopen(intreename, "rb");
  outfile = fopen(OutfileName, outfileopt);
  strcpy(outfilename, OutfileName);

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
  doinit();

  consenserun();  // do the actual work

  FClose(intree);
  FClose(outfile);
  if (trout)
  {
    FClose(outtree);
  }
  if (progress)
  {
    FClose(progfile);
  }
  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{
  /* Local variables added by Dan F. */
#if 0
  pattern_elm  ***pattern_array;
  long trees_in = 0;
  long i, j;
  node *p, *q;
#endif
  initdata *funcs;

#ifdef MAC
  argc = 1;                /* macsetup("Consense", "");        */
  argv[0] = "Consense";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = cons_node_new;
  phylipinit(argc, argv, funcs, false);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree file", "rb", argv[0], intreename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  /* Initialize option-based variables, then ask for changes regarding
     their values. */
  doinit();

  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  consenserun();

  FClose(outtree);
  FClose(intree);
  FClose(outfile);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  printf("Done.\n\n");

  phyRestoreConsoleAttributes();

  free(funcs);

  return 0;
}  /* main */


// End.
