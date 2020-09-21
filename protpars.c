/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef DEBUG
#define TESTING(x) fprintf(outfile, "TEST(%c)\n", x); TEST1; TEST2
#define TEST1 printree(curtree)
#define TEST2 fprintf(outfile, "\nrequires a total of %10.3f\n", curtree->score / -10)
#define TEST3(x) fprintf(outfile, "score(%c)=%g\n", x, t->score)
#define TEST4 if(root) for(i=0;i<endsite;i++) printf("%d ", ((pars_node*)root)->numsteps[i])
#define TEST5 fprintf(outfile, "maybe = %d ", *maybe)
#define TEST6 fprintf(outfile, "nonzero = %d ", *nonzero)
#define TEST7 TEST7a; fprintf(outfile, "\n")
#define TEST7a fprintf(outfile, "site[%d] = %d ", i, ((protpars_node*)r)->siteset[i][0])
#define TEST7b fprintf(outfile, "site[%d] = %d ", i, ((protpars_node*)r)->siteset[i][1])
#define TEST7c fprintf(outfile, "site[%d] = %d ", i, ((protpars_node*)r)->siteset[i][2])
#define TEST8 fprintf(outfile, "in prothypstates\n")
#define TEST9 fprintf(outfile, "i=%d, hypset=%g, siteset1=%g, siteset2=%g, seq[i]=%c, k=%d\n", \
                      i, hypset[i], ((protpars_node*)r->next->back)->siteset[i], \
                      ((protpars_node*)r->next->next->back)->siteset[i], temparray->seq[i], k)
#else
#define TESTING(x)
#define TEST1
#define TEST2
#define TEST3(x)
#define TEST4
#define TEST5
#define TEST6
#define TEST7
#define TEST8
#define TEST9
#endif

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "Slist.h"
#include "seq.h"
#include "parsimony.h"

extern long nextree;

typedef enum {
  universal, ciliate, mito, vertmito, flymito, yeastmito
} codetype;

/* nodes will form a binary tree */

typedef struct gseq {
  seqptr seq;
  struct gseq *next;
} gseq;

typedef struct protpars_node {
  pars_node pars_node;
  aas *seq;                  /* the sequence */
  seqptr siteset;            /* temporary storage for aa's used in protpars*/
} protpars_node;


typedef struct protpars_tree {
  pars_tree pars_tree;
} protpars_tree;

#ifndef OLDC
/* function prototypes */
void   protgnu(gseq **);
void   protchuck(gseq *);
void   code(void);
void   setup(void);
void   getoptions(void);
void   protalloctree(void);
void   allocrest(void);
void   doinit(void);
void   protinputdata(void);

void   prot_makevalues(tree* t, boolean usertree);
void   doinput(void);
void   protfillin(node *, node *, node *);
void   protpreorder(node *);
void   protadd(node *, node *, node *);
void   protre_move(node **, node **);
void   evaluate(node *);
void   protreroot(node *);
void   protgetch(Char *);
void   protaddelement(node **, long *, long *, boolean *, Char *);
void   prottreeread(void);
void   protancestset(long *, long *, long *, long *, long *);

void   prothyprint(long, long, boolean *, node *, boolean *, boolean *);
void   prothyptrav(node *, sitearray *, long, long, long *, boolean *, sitearray);
void   prothypstates(long *);
void   describe(void);
void   maketree(void);
void   reallocnode(node* p);
void   reallocchars(void);
node * protpars_node_new(node_type type, long index);
tree * protpars_tree_new(long nonodes, long spp);
void   protpars_node_init(node* n, node_type type, long index);
void   protpars_tree_init(tree *t, long nonodes, long spp);
void   protpars_tree_nuview(tree* t, node* p);
double protpars_tree_evaluate(tree* t, node* p, boolean );
void protparsrun(void);
void protpars(char * infilename, char * intreename, char * outfilename, char * outfileopt, char * weightsfilename,
              char * outtreename, char * outtreeopt, int searchbest, int treesave, int inputorder, int RandNum,
              int NumJumble, int OutRoot, int OutNum, int ThreshPars, double ThreshVal, char * GeneCode, int SitesWeighted,
              int AnalyzeMult, int NumMult, int MultDataSet, int InputSeq, int PrintData, int DotDiff, int PrintInd,
              int PrintTree, int PrintSteps, int PrintSeq, int WriteTree);
/* function prototypes */
#endif

long jumb = 0, nonodes = 0;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
long chars, col, msets, ith, njumble;
/*   chars = number of sites in actual sequences */
extern long maxtrees;
long inseed, inseed0;
boolean jumble, usertree, weights, thresh, trout, progress, stepbox, justwts, ancseq, mulsets, firstset, rearrfirst = true;
codetype whichcode;
steptr oldweight; /* to make writesteps happy */
long fullset, fulldel;
double threshold;
double *threshwt;
longer seed;
long *enterorder;
sitearray translate[(long)quest - (long)ala + 1];
aas trans[4][4][4];
long **fsteps;
bestelm *bestrees;
boolean dummy;
gseq *garbage;
/* Char ch; */
aas tmpa;
char *progname;
tree* curtree;

/* Local variables for maketree, propagated globally for C version: */
long minwhich;
double like, bestyet, bestlike, minsteps, bstlike2;
boolean lastrearr, recompute;
node *there;
double nsteps[maxuser];
long *place;
boolean *names;


tree* protpars_tree_new(long nonodes, long spp)
{
  tree* t = Malloc(sizeof(protpars_tree));
  protpars_tree_init(t, nonodes, spp);

  return t;
}


void protpars_tree_init(tree *t, long nonodes, long spp)
{
  pars_tree_init((tree*)t, nonodes, spp);
  t->nuview = protpars_tree_nuview;
  t->evaluate = protpars_tree_evaluate;
}


node* protpars_node_new(node_type type, long index) // RSGbugfix
{
  node* n = Malloc(sizeof(protpars_node));
  protpars_node_init(n, type, index);
  return n;
}


void protpars_node_init(node* node, node_type type, long index)
{
  protpars_node *n = (protpars_node *)node;
  pars_node_init(node, type, index);
  node->init = protpars_node_init;
  if ( n != NULL && n->seq != NULL )
    free(n->seq);
  n->seq =  (aas*)Malloc(chars * sizeof(aas));
  if ( n->siteset)
    free(n->siteset);
  n->siteset = (seqptr)Malloc(chars * sizeof(sitearray));
}


void protgnu(gseq **p)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL)
  {
    *p = garbage;
    free((*p)->seq);
    (*p)->seq = (seqptr)Malloc(chars * sizeof(sitearray));
    garbage = garbage->next;
  }
  else
  {
    *p = (gseq *)Malloc(sizeof(gseq));
    (*p)->seq = (seqptr)Malloc(chars * sizeof(sitearray));
  }
  (*p)->next = NULL;
}  /* protgnu */


void protchuck(gseq *p)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* protchuck */


void code(void)
{
  /* make up table of the code 0 = u, 1 = c, 2 = a, 3 = g */
  trans[0][0][0] = phe;
  trans[0][0][1] = phe;
  trans[0][0][2] = leu;
  trans[0][0][3] = leu;
  trans[0][1][0] = ser1;
  trans[0][1][1] = ser1;
  trans[0][1][2] = ser1;
  trans[0][1][3] = ser1;
  trans[0][2][0] = tyr;
  trans[0][2][1] = tyr;
  trans[0][2][2] = stop;
  trans[0][2][3] = stop;
  trans[0][3][0] = cys;
  trans[0][3][1] = cys;
  trans[0][3][2] = stop;
  trans[0][3][3] = trp;
  trans[1][0][0] = leu;
  trans[1][0][1] = leu;
  trans[1][0][2] = leu;
  trans[1][0][3] = leu;
  trans[1][1][0] = pro;
  trans[1][1][1] = pro;
  trans[1][1][2] = pro;
  trans[1][1][3] = pro;
  trans[1][2][0] = his;
  trans[1][2][1] = his;
  trans[1][2][2] = gln;
  trans[1][2][3] = gln;
  trans[1][3][0] = arg;
  trans[1][3][1] = arg;
  trans[1][3][2] = arg;
  trans[1][3][3] = arg;
  trans[2][0][0] = ileu;
  trans[2][0][1] = ileu;
  trans[2][0][2] = ileu;
  trans[2][0][3] = met;
  trans[2][1][0] = thr;
  trans[2][1][1] = thr;
  trans[2][1][2] = thr;
  trans[2][1][3] = thr;
  trans[2][2][0] = asn;
  trans[2][2][1] = asn;
  trans[2][2][2] = lys;
  trans[2][2][3] = lys;
  trans[2][3][0] = ser2;
  trans[2][3][1] = ser2;
  trans[2][3][2] = arg;
  trans[2][3][3] = arg;
  trans[3][0][0] = val;
  trans[3][0][1] = val;
  trans[3][0][2] = val;
  trans[3][0][3] = val;
  trans[3][1][0] = ala;
  trans[3][1][1] = ala;
  trans[3][1][2] = ala;
  trans[3][1][3] = ala;
  trans[3][2][0] = asp;
  trans[3][2][1] = asp;
  trans[3][2][2] = glu;
  trans[3][2][3] = glu;
  trans[3][3][0] = gly;
  trans[3][3][1] = gly;
  trans[3][3][2] = gly;
  trans[3][3][3] = gly;
  if (whichcode == mito)
    trans[0][3][2] = trp;
  if (whichcode == vertmito)
  {
    trans[0][3][2] = trp;
    trans[2][3][2] = stop;
    trans[2][3][3] = stop;
    trans[2][0][2] = met;
  }
  if (whichcode == flymito)
  {
    trans[0][3][2] = trp;
    trans[2][0][2] = met;
    trans[2][3][2] = ser2;
  }
  if (whichcode == yeastmito)
  {
    trans[0][3][2] = trp;
    trans[1][0][2] = thr;
    trans[2][0][2] = met;
  }
} /* code */


void setup(void)
{
  /* set up set table to get aasets from aas */
  aas a, b;
  long i, j, k, l, s;

  for (a = ala; (long)a <= (long)stop; a = (aas)((long)a + 1))
  {
    translate[(long)a - (long)ala][0] = 1L << ((long)a);
    translate[(long)a - (long)ala][1] = 1L << ((long)a);
  }
  for (i = 0; i <= 3; i++)
  {
    for (j = 0; j <= 3; j++)
    {
      for (k = 0; k <= 3; k++)
      {
        for (l = 0; l <= 3; l++)
        {
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[l][j][k]);
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[i][l][k]);
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[i][j][l]);
        }
      }
    }
  }
  translate[(long)del - (long)ala][1] = 1L << ((long)del);
  fulldel = (1L << ((long)stop + 1)) - (1L << ((long)ala));
  fullset = fulldel & (~(1L << ((long)del)));
  translate[(long)asx - (long)ala][0]
    = (1L << ((long)asn)) | (1L << ((long)asp));
  translate[(long)glx - (long)ala][0]
    = (1L << ((long)gln)) | (1L << ((long)glu));
  translate[(long)ser - (long)ala][0]
    = (1L << ((long)ser1)) | (1L << ((long)ser2));
  translate[(long)unk - (long)ala][0] = fullset;
  translate[(long)quest - (long)ala][0] = fulldel;
  translate[(long)asx - (long)ala][1] = translate[(long)asn - (long)ala][1]
                                       | translate[(long)asp - (long)ala][1];
  translate[(long)glx - (long)ala][1] = translate[(long)gln - (long)ala][1]
                                       | translate[(long)glu - (long)ala][1];
  translate[(long)ser - (long)ala][1] = translate[(long)ser1 - (long)ala][1]
                                       | translate[(long)ser2 - (long)ala][1];
  translate[(long)unk - (long)ala][1] = fullset;
  translate[(long)quest - (long)ala][1] = fulldel;
  for (a = ala; (long)a <= (long)quest; a = (aas)((long)a + 1))
  {
    s = 0;
    for (b = ala; (long)b <= (long)stop; b = (aas)((long)b + 1))
    {
      if (((1L << ((long)b)) & translate[(long)a - (long)ala][1]) != 0)
        s |= translate[(long)b - (long)ala][1];
    }
    translate[(long)a - (long)ala][2] = s;
  }
}  /* setup */


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
  thresh = false;
  trout = true;
  usertree = false;
  weights = false;
  whichcode = universal;
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
    printf("\nProtein parsimony algorithm, version %s\n\n", VERSION);
    printf("Setting for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree)
    {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  C               Use which genetic code?  %s\n",
           (whichcode == universal) ? "Universal"                  :
           (whichcode == ciliate)   ? "Ciliate"                    :
           (whichcode == mito)      ? "Universal mitochondrial"    :
           (whichcode == vertmito)  ? "Vertebrate mitochondrial"   :
           (whichcode == flymito)   ? "Fly mitochondrial"          :
           (whichcode == yeastmito) ? "Yeast mitochondrial"        : "");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    if (ancseq || printdata)
      printf("  .  Use dot-differencing to display them  %s\n",
             dotdiff ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts)
    {
      printf("WARNING:  W option and Multiple Weights options are both on.  ");
      printf("The W menu option is unnecessary and has no additional effect. \n");
    }
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WCJOTUMI12345.60", ch) != NULL))
        || (usertree && ((strchr("WCOTUMI12345.60", ch) != NULL))))
    {
      switch (ch)
      {
        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          else outgrno = 1;
          break;

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
          break;

        case 'C':
          printf("\nWhich genetic code?\n");
          printf(" type         for\n\n");
          printf("   U           Universal\n");
          printf("   M           Mitochondrial\n");
          printf("   V           Vertebrate mitochondrial\n");
          printf("   F           Fly mitochondrial\n");
          printf("   Y           Yeast mitochondrial\n\n");
          loopcount2 = 0;
          do {
            printf("type U, M, V, F, or Y\n");
            if(scanf("%c%*[^\n]", &ch)) {} // Read char and scan to EOL.
            (void)getchar();
            if (ch == '\n')
              ch = ' ';
            uppercase(&ch);
            countup(&loopcount2, 10);
          } while (ch != 'U' && ch != 'M' && ch != 'V'
                   && ch != 'F' && ch != 'Y');
          switch (ch)
          {
            case 'U':
              whichcode = universal;
              break;

            case 'M':
              whichcode = mito;
              break;

            case 'V':
              whichcode = vertmito;
              break;

            case 'F':
              whichcode = flymito;
              break;

            case 'Y':
              whichcode = yeastmito;
              break;
          }
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

        case 'U':
          usertree = !usertree;
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
{ /* reallocates variables that are dependand on the number of chars
   * do we need to reallocate the garbage list too? */
  long i;
  node *p;

  if (usertree)
    for (i = 0; i < maxuser; i++)
    {
      free(fsteps[i]);
      fsteps[i] = (long *)Malloc(chars * sizeof(long));
    }

  for (i = 0; i < nonodes; i++)
  {
    curtree->nodep[i]->init(curtree->nodep[i], i<spp, i);
    if (i >= spp)
    {
      p = curtree->nodep[i]->next;
      while (p != curtree->nodep[i])
      {
        p->init(p, false, i);
        p = p->next;
      }
    }
  }

  free(weight);
  free(oldweight);
  free(threshwt);

  weight = (steptr)Malloc(chars * sizeof(long));
  oldweight = (steptr)Malloc(chars * sizeof(long));
  threshwt = (double*)Malloc(chars * sizeof(double));
}


void allocrest(void)
{ /* allocate remaining global arrays and variables dynamically */
  long i;

  if (usertree)
  {
    fsteps = (long **)Malloc(maxuser * sizeof(long *));
    for (i = 0; i < maxuser; i++)
      fsteps[i] = (long *)Malloc(chars * sizeof(long));
  }
  bestrees = (bestelm *)Malloc(maxtrees * sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(spp * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  enterorder = (long *)Malloc(spp * sizeof(long));
  place = (long *)Malloc(nonodes * sizeof(long));
  weight = (steptr)Malloc(chars * sizeof(long));
  threshwt = (double*)Malloc(chars * sizeof(double));
  oldweight = (long *)Malloc(chars * sizeof(long));
  /*  alias = (long *)Malloc(chars * sizeof(long)); */
  ally = (long *)Malloc(chars * sizeof(long));
  location = (long *)Malloc(chars * sizeof(long));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nProtein parsimony algorithm, version %s\n\n", VERSION);

  inputnumbers(&spp, &chars, &nonodes, 1);
  endsite = chars;

  if (!javarun)
  {
    getoptions();
  }

  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n\n", spp, chars);

  curtree = (tree*)protpars_tree_new(nonodes, spp);

  allocrest();
}  /* doinit*/


void protinputdata(void)
{
  /* input the names and sequences for each species */
  long i, j, k, l, aasread, aasnew = 0;
  Char charstate;
  boolean allread, done;

  // RSGnote: "aa" may be referenced before being assigned in original version.
  // It is initialized here merely to silence a compiler warning (in case this will never happen).
  aas aa = ala;   /* temporary amino acid for input */

  if (printdata)
    headings(chars, "Sequences", "---------");
  aasread = 0;
  allread = false;
  while (!(allread))
  {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
    {
      scan_eoln(infile);
    }
    i = 1;
    while (i <= spp)
    {
      if ((interleaved && aasread == 0) || !interleaved)
        initname(i - 1);
      j = interleaved ? aasread : 0;
      done = false;
      while (!done && !eoff(infile))
      {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile)))
        {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((!isalpha(charstate) && charstate != '?' &&
               charstate != '-' && charstate != '*') || charstate == 'J' ||
              charstate == 'O' || charstate == 'U')
          {
            printf("WARNING -- BAD AMINO ACID:%c", charstate);
            printf(" AT POSITION%5ld OF SPECIES %3ld\n", j, i);
            exxit(-1);
          }
          j++;
          aa = (charstate == 'A') ?  ala :
            (charstate == 'B') ?  asx :
            (charstate == 'C') ?  cys :
            (charstate == 'D') ?  asp :
            (charstate == 'E') ?  glu :
            (charstate == 'F') ?  phe :
            (charstate == 'G') ?  gly : aa;
          aa = (charstate == 'H') ?  his :
            (charstate == 'I') ? ileu :
            (charstate == 'K') ?  lys :
            (charstate == 'L') ?  leu :
            (charstate == 'M') ?  met :
            (charstate == 'N') ?  asn :
            (charstate == 'P') ?  pro :
            (charstate == 'Q') ?  gln :
            (charstate == 'R') ?  arg : aa;
          aa = (charstate == 'S') ?  ser :
            (charstate == 'T') ?  thr :
            (charstate == 'V') ?  val :
            (charstate == 'W') ?  trp :
            (charstate == 'X') ?  unk :
            (charstate == 'Y') ?  tyr :
            (charstate == 'Z') ?  glx :
            (charstate == '*') ? stop :
            (charstate == '?') ? quest:
            (charstate == '-') ? del  :  aa;

          // RSGnote: Variable "aa" may be referenced before being initialized here.
          ((protpars_node*)curtree->nodep[i - 1])->seq[j - 1] = aa;
          memcpy(((protpars_node*)curtree->nodep[i - 1])->siteset[j - 1], translate[(long)aa - (long)ala], sizeof(sitearray));
        }
        if (interleaved)
          continue;
        if (j < chars)
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        aasnew = j;
      scan_eoln(infile);
      if ((interleaved && j != aasnew) || ((!interleaved) && j != chars))
      {
        printf("\nERROR:  SEQUENCES OUT OF ALIGNMENT.\n");
        exxit(-1);
      }
      i++;
    }
    if (interleaved)
    {
      aasread = aasnew;
      allread = (aasread == chars);
    }
    else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (printdata)
  {
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++)
    {
      for (j = 1; j <= spp; j++)
      {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars)
          l = chars;
        for (k = (i - 1) * 60 + 1; k <= l; k++)
        {
          if (j > 1 && ((protpars_node*)curtree->nodep[j - 1])->seq[k - 1] == ((protpars_node*)curtree->nodep[0])->seq[k - 1])
            charstate = '.';
          else
          {
            tmpa = ((protpars_node*)curtree->nodep[j-1])->seq[k-1];
            charstate =  (tmpa == ala) ? 'A' :
              (tmpa == asx) ? 'B' :
              (tmpa == cys) ? 'C' :
              (tmpa == asp) ? 'D' :
              (tmpa == glu) ? 'E' :
              (tmpa == phe) ? 'F' :
              (tmpa == gly) ? 'G' :
              (tmpa == his) ? 'H' :
              (tmpa ==ileu) ? 'I' :
              (tmpa == lys) ? 'K' :
              (tmpa == leu) ? 'L' : charstate;
            charstate =  (tmpa == met) ? 'M' :
              (tmpa == asn) ? 'N' :
              (tmpa == pro) ? 'P' :
              (tmpa == gln) ? 'Q' :
              (tmpa == arg) ? 'R' :
              (tmpa == ser) ? 'S' :
              (tmpa ==ser1) ? 'S' :
              (tmpa ==ser2) ? 'S' : charstate;
            charstate =  (tmpa == thr) ? 'T' :
              (tmpa == val) ? 'V' :
              (tmpa == trp) ? 'W' :
              (tmpa == unk) ? 'X' :
              (tmpa == tyr) ? 'Y' :
              (tmpa == glx) ? 'Z' :
              (tmpa == del) ? '-' :
              (tmpa ==stop) ? '*' :
              (tmpa==quest) ? '?' : charstate;
          }
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
  putc('\n', outfile);
}  /* protinputdata */


void prot_makevalues(tree* t, boolean usertree)
{
  /* init nodes? */
  long i, j;
  node *p;

  (void)usertree;                       // RSGnote: Parameter never used.

#ifdef DEBUG
  fprintf(outfile, "protpars:protmakevalues\n");   /* debug */
#endif
  for (i = 1; i <= nonodes; i++)
  {
    ((node*)t->nodep[i - 1])->back = NULL;
    ((node*)t->nodep[i - 1])->tip = (i <= spp);
    ((node*)t->nodep[i - 1])->index = i;
    for (j = 0; j < chars; j++)
      ((pars_node*)t->nodep[i - 1])->numsteps[j] = 0;
    if (i > spp)
    {
      p = ((node*)t->nodep[i - 1])->next;
      while (p != ((node*)t->nodep[i - 1]))
      {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        for (j = 0; j < chars; j++)
          ((pars_node*)p)->numsteps[j] = 0;
        p = p->next;
      }
    }
  }
}  /* prot_makevalues */


void doinput(void)
{
  /* reads the input data */
  long i;

  if (justwts)
  {
    if (firstset)
      protinputdata();
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
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (weights)
    {
      inputweights(chars, weight, &weights);
    }
    if (weights)
      printweights(outfile, 0, chars, weight, "Sites");
    protinputdata();
  }
  if(!thresh)
    threshold = spp * 3.0;
  for(i = 0 ; i < chars ; i++)
  {
    threshwt[i] = (threshold * weight[i]);
  }
  prot_makevalues(curtree, usertree);
}  /* doinput */


double protpars_tree_evaluate(tree* t, node* p, boolean dummy)
{
  protpars_node *q = (protpars_node*)p->back;
  protpars_node *r = (protpars_node*)p;
  long m, sum = 0;
  sitearray base;

  generic_tree_evaluate(t, p, dummy);

  for ( m = 0 ; m < chars ; m++ )
  {
    sum += ((pars_node*)p)->numsteps[m];
    if ( q == NULL ) continue;
    sum += ((pars_node*)q)->numsteps[m];
    base[0] = q->siteset[m][0] & r->siteset[m][0];
    if ( base[0] == 0 )
      sum += 1;
    base[0] = (q->siteset[m][0] & r->siteset[m][1]) |
      (q->siteset[m][1] & r->siteset[m][0]);
    if ( base[0] == 0 )
      sum += 1;
    base[0] = (q->siteset[m][0] & r->siteset[m][2]) |
      (q->siteset[m][1] & r->siteset[m][1]) |
      (q->siteset[m][2] & r->siteset[m][2]);
    if ( base[0] == 0 )
      sum += 1;
  }
  t->score = -sum;
  TEST3('a');
  return -sum;
}


void protpars_tree_nuview(tree* t, node* p)
{
  generic_tree_nuview(t, p);
  protfillin(p, p->next->back, p->next->next->back);
}


void protfillin(node *p, node *left, node *rt)
{
  /* sets up for each node in the tree the aa set for site m
     at that point and counts the changes.  The program
     spends much of its time in this function */
  boolean counted, done;
  aas aa;
  long s = 0;
  sitearray ls, rs, qs;
  long i, j, m, n;

  for (m = 0; m < chars; m++)
  {
    if (left != NULL)
      memcpy(ls, ((protpars_node*)left)->siteset[m], sizeof(sitearray));
    if (rt != NULL)
      memcpy(rs, ((protpars_node*)rt)->siteset[m], sizeof(sitearray));
    if (left == NULL)
    {
      n = ((pars_node*)rt)->numsteps[m];
      memcpy(qs, rs, sizeof(sitearray));
    }
    else if (rt == NULL)
    {
      n = ((pars_node*)left)->numsteps[m];
      memcpy(qs, ls, sizeof(sitearray));
    }
    else
    {
      n = ((pars_node*)left)->numsteps[m] + ((pars_node*)rt)->numsteps[m];
      if ((ls[0] == rs[0]) && (ls[1] == rs[1]) && (ls[2] == rs[2]))
      {
        qs[0] = ls[0];
        qs[1] = ls[1];
        qs[2] = ls[2];
      }
      else
      {
        counted = false;
        for (i = 0; (!counted) && (i <= 3); i++)
        {
          switch (i)
          {
            case 0:
              s = ls[0] & rs[0];
              break;

            case 1:
              s = (ls[0] & rs[1]) | (ls[1] & rs[0]);
              break;

            case 2:
              s = (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
              break;

            case 3:
              s = ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
              break;

          }
          if (s != 0)
          {
            qs[0] = s;
            counted = true;
          }
          else
            n += weight[m];
        }
        /*        qs[1] = 0;    */
        /*        qs[2] = 0;    */
        switch (i)
        {
          case 1:
            qs[1] = qs[0] | (ls[0] & rs[1]) | (ls[1] & rs[0]);
            qs[2] = qs[1] | (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
            break;
          case 2:
            qs[1] = qs[0] | (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
            qs[2] = qs[1] | ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
            break;
          case 3:
            qs[1] = qs[0] | ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
            qs[2] = qs[1] | ls[1] | (ls[2] & rs[2]) | rs[1];
            break;
          case 4:
            qs[1] = qs[0] | ls[1] | (ls[2] & rs[2]) | rs[1];
            qs[2] = qs[1] | ls[2] | rs[2];
            break;
        }
        for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1))
        {
          done = false;
          for (i = 0; (!done) && (i <= 1); i++)
          {
            if (((1L << ((long)aa)) & qs[i]) != 0)
            {
              for (j = i+1; j <= 2; j++)
                qs[j] |= translate[(long)aa - (long)ala][j-i];
              done = true;
            }
          }
        }
      }
    }
    ((pars_node*)p)->numsteps[m] = n;
    memcpy(((protpars_node*)p)->siteset[m], qs, sizeof(sitearray));
  }
}  /* protfillin */


void protpreorder(node *p)
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account */
  if (p != NULL && !p->tip)
  {
    protfillin (p->next, p->next->next->back, p->back);
    protfillin (p->next->next, p->back, p->next->back);
    protpreorder (p->next->back);
    protpreorder (p->next->next->back);
  }
} /* protpreorder */



void protreroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q;

  if (outgroup->back->index == curtree->root->index)
    return;
  p = curtree->root->next;
  q = curtree->root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* protreroot */


void protgetch(Char *c)
{
  /* get next nonblank character */
  do {
    if (eoln(intree))
      scan_eoln(intree);
    *c = gettc(intree);
    if (*c == '\n' || *c == '\t')
      *c = ' ';
  } while (!(*c != ' ' || eoff(intree)));
}  /* protgetch */


void protaddelement(node **p, long *nextnode, long *lparens, boolean *names, Char *ch)
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  long i, n;
  boolean found;
  Char str[nmlngth];

  protgetch(ch);

  if (*ch == '(' )
  {
    if ((*lparens) >= spp - 1)
    {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exxit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = curtree->nodep[(*nextnode) - 1];
    protaddelement(&q->next->back, nextnode, lparens, names, ch);
    q->next->back->back = q->next;
    findch(',', ch, which);
    protaddelement(&q->next->next->back, nextnode, lparens, names, ch);
    q->next->next->back->back = q->next->next;
    findch(')', ch, which);
    *p = q;
    return;
  }
  for (i = 0; i < nmlngth; i++)
    str[i] = ' ';
  n = 1;
  do {
    if (*ch == '_')
      *ch = ' ';
    str[n - 1] = *ch;
    if (eoln(intree))
      scan_eoln(intree);
    *ch = gettc(intree);
    n++;
  } while (*ch != ',' && *ch != ')' && *ch != ':' && n <= nmlngth);
  n = 1;
  do {
    found = true;
    for (i = 0; i < nmlngth; i++)
      found = (found && ((str[i] == nayme[n - 1][i]) ||
                         ((nayme[n - 1][i] == '_') && (str[i] == ' '))));
    if (found)
    {
      if (names[n - 1] == false)
      {
        *p = curtree->nodep[n - 1];
        names[n - 1] = true;
      }
      else
      {
        printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
        for (i = 0; i < nmlngth; i++)
          putchar(nayme[n - 1][i]);
        putchar('\n');
        exxit(-1);
      }
    }
    else
      n++;
  } while (!(n > spp || found));
  if (n <= spp)
    return;
  printf("CANNOT FIND SPECIES: ");
  for (i = 0; i < nmlngth; i++)
    putchar(str[i]);
  putchar('\n');
}  /* protaddelement */


void prottreeread(void)
{
  /* read in user-defined tree and set it up */
  long nextnode, lparens, i;
  char ch;

  curtree->root = curtree->nodep[spp];
  nextnode = spp;
  curtree->root->back = NULL;
  names = (boolean *)Malloc(spp * sizeof(boolean));
  for (i = 0; i < spp; i++)
    names[i] = false;
  lparens = 0;
  protaddelement(&curtree->root, &nextnode, &lparens, names, &ch);
  if (ch == '[')
  {
    do
      ch = gettc(intree);
    while (ch != ']');
    ch = gettc(intree);
  }
  findch(';', &ch, which);
  if (progress)
  {
    sprintf(progbuf, "\n\n");
    print_progress(progbuf);
  }
  scan_eoln(intree);
/********  gettc(intree);
  what was this here for?  ***/
  free(names);
}  /* prottreeread */


void protancestset(long *a, long *b, long *c, long *d, long *k)
{
  /* sets up the aa set array. */
  aas aa;
  long s, sa, sb;
  long i, j, m, n;
  boolean counted;

  counted = false;
  *k = 0;
  for (i = 0; i <= 5; i++)
  {
    if (*k < 3)
    {
      s = 0;
      if (i > 3)
        n = i - 3;
      else
        n = 0;
      for (j = n; j <= (i - n); j++)
      {
        if (j < 3)
          sa = a[j];
        else
          sa = fullset;
        for (m = n; m <= (i - j - n); m++)
        {
          if (m < 3)
            sb = sa & b[m];
          else
            sb = sa;
          if (i - j - m < 3)
            sb &= c[i - j - m];
          s |= sb;
        }
      }
      if (counted || s != 0)
      {
        d[*k] = s;
        (*k)++;
        counted = true;
      }
    }
  }
  for (i = 0; i <= 1; i++)
  {
    for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1))
    {
      if (((1L << ((long)aa)) & d[i]) != 0)
      {
        for (j = i + 1; j <= 2; j++)
          d[j] |= translate[(long)aa - (long)ala][j - i];
      }
    }
  }
}  /* protancestset */


void prothyprint(long b1, long b2, boolean *bottom, node *r,
                 boolean *nonzero, boolean *maybe)
{
  /* print out states in sites b1 through b2 at node */
  long i;
  boolean dot;
  Char ch = 0;
  aas aa;

  if (*bottom)
  {
    if (!outgropt)
      fprintf(outfile, "      ");
    else
      fprintf(outfile, "root  ");
  }
  else
    fprintf(outfile, "%3ld   ", r->back->index - spp);
  if (r->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[r->index - 1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", r->index - spp);
  TEST5;
  TEST6;
  if (*bottom)
    fprintf(outfile, "          ");
  else if (*nonzero)
    fprintf(outfile, "   yes    ");
  else if (*maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1 - 1; i < b2; i++)
  {
    aa = ((protpars_node*)r)->seq[i];
    switch (aa)
    {
      case ala:
        ch = 'A';
        break;

      case asx:
        ch = 'B';
        break;

      case cys:
        ch = 'C';
        break;

      case asp:
        ch = 'D';
        break;

      case glu:
        ch = 'E';
        break;

      case phe:
        ch = 'F';
        break;

      case gly:
        ch = 'G';
        break;

      case his:
        ch = 'H';
        break;

      case ileu:
        ch = 'I';
        break;

      case lys:
        ch = 'K';
        break;

      case leu:
        ch = 'L';
        break;

      case met:
        ch = 'M';
        break;

      case asn:
        ch = 'N';
        break;

      case pro:
        ch = 'P';
        break;

      case gln:
        ch = 'Q';
        break;

      case arg:
        ch = 'R';
        break;

      case ser:
        ch = 'S';
        break;

      case ser1:
        ch = 'S';
        break;

      case ser2:
        ch = 'S';
        break;

      case thr:
        ch = 'T';
        break;

      case trp:
        ch = 'W';
        break;

      case tyr:
        ch = 'Y';
        break;

      case val:
        ch = 'V';
        break;

      case glx:
        ch = 'Z';
        break;

      case del:
        ch = '-';
        break;

      case stop:
        ch = '*';
        break;

      case unk:
        ch = 'X';
        break;

      case quest:
        ch = '?';
        break;
    }
    if (!(*bottom) && dotdiff)
      dot = (((protpars_node*)r)->siteset[i] [0] ==
             ((protpars_node*)curtree->nodep[r->back->index - 1])->siteset[i][0]
             || ((((protpars_node*)r)->siteset[i][0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                     (1L << ((long)ser))))) == 0 &&
                 (((protpars_node*)curtree->nodep[r->back->index - 1]) ->siteset[i] [0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                     (1L << ((long)ser))))) == 0));
    else
      dot = false;
    if (dot)
      putc('.', outfile);
    else
      putc(ch, outfile);
    if ((i + 1) % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* prothyprint */


void prothyptrav(node *r, sitearray *hypset, long b1, long b2, long *k,
                 boolean *bottom, sitearray nothing)
{
  boolean maybe, nonzero;
  long i;
  aas aa;
  long anc = 0, hset;
  gseq *ancset, *temparray;

  protgnu(&ancset);
  protgnu(&temparray);
  maybe = false;
  nonzero = false;
  for (i = b1 - 1; i < b2; i++)
  {
    if (!r->tip)
    {
      protancestset(hypset[i], ((protpars_node*)r->next->back)->siteset[i], ((protpars_node*)r->next->next->back)->siteset[i], temparray->seq[i], k);
      memcpy(((protpars_node*)r)->siteset[i], temparray->seq[i], sizeof(sitearray));
      TEST9;
    }
    if (!(*bottom))
      anc = ((protpars_node*)curtree->nodep[r->back->index - 1])
        ->siteset[i][0];
    if (!r->tip)
    {
      hset = ((protpars_node*)r)->siteset[i][0];
      ((protpars_node*)r)->seq[i] = quest;
      for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1))
      {
        if (hset == 1L << ((long)aa))
          ((protpars_node*)r)->seq[i] = aa;
      }
      if (hset == ((1L << ((long)asn)) | (1L << ((long)asp))))
        ((protpars_node*)r)->seq[i] = asx;
      if (hset == ((1L << ((long)gln)) | (1L << ((long)gly))))
        ((protpars_node*)r)->seq[i] = glx;
      if (hset == ((1L << ((long)ser1)) | (1L << ((long)ser2))))
        ((protpars_node*)r)->seq[i] = ser;
      if (hset == fullset)
        ((protpars_node*)r)->seq[i] = unk;
    }
    TEST7;
    nonzero = (nonzero || (((protpars_node*)r)->siteset[i][0] & anc) == 0);
    maybe = (maybe || ((protpars_node*)r)->siteset[i][0] != anc);
  }
  prothyprint(b1, b2, bottom, r, &nonzero, &maybe);
  *bottom = false;
  if (!r->tip)
  {
    memcpy(temparray->seq, ((protpars_node*)r->next->back)->siteset, chars * sizeof(sitearray));
    for (i = b1 - 1; i < b2; i++)
      protancestset(hypset[i], ((protpars_node*)r->next->next->back)->siteset[i], nothing, ancset->seq[i], k);
    prothyptrav(r->next->back, ancset->seq, b1, b2, k, bottom, nothing );
    for (i = b1 - 1; i < b2; i++)
      protancestset(hypset[i], temparray->seq[i], nothing, ancset->seq[i], k);
    prothyptrav(r->next->next->back, ancset->seq, b1, b2, k, bottom, nothing);
  }
  protchuck(temparray);
  protchuck(ancset);
}  /* prothyptrav */


void prothypstates(long *k)
{
  /* fill in and describe states at interior nodes */
  boolean bottom;
  sitearray nothing;
  long i, n;
  seqptr hypset;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                             ");
  fprintf(outfile, "( . means same as in the node below it on tree)\n\n");
  memcpy(nothing, translate[(long)quest - (long)ala], sizeof(sitearray));
  hypset = (seqptr)Malloc(chars * sizeof(sitearray));
  for (i = 0; i < chars; i++)
    memcpy(hypset[i], nothing, sizeof(sitearray));
  bottom = true;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++)
  {
    putc('\n', outfile);
    TEST8;
    n = i * 40;
    if (n > chars)
      n = chars;
    bottom = true;
    prothyptrav(curtree->root, hypset, i * 40 - 39, n, k, &bottom, nothing);
  }
  free(hypset);
}  /* prothypstates */


void describe(void)
{
  void debugtree(tree*, FILE*);
  long k;
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", -(curtree->score));
  if (stepbox)
  {
    putc('\n', outfile);
    writesteps(curtree, chars, weights, oldweight);
  }
  if (ancseq)
  {
    prothypstates(&k);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout)
  {
    col = 0;
    treeout(curtree->root, nextree, &col, curtree->root);
#ifdef DEBUG
    debugtree(curtree, outtree);
#endif
  }
}  /* describe */


void maketree(void)                     // RSGbugfix
{
  /* constructs a binary tree from the pointers in curtree->nodep.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees;
  boolean done; /* tst */

  if (!usertree)
  {
    lastrearr=false;    /* tst */
    hsbut(curtree, false, jumble, seed, progress);
    TESTING('1');        /* testing      */
    if (progress)
    {
      sprintf(progbuf, "\nDoing global rearrangements");
      print_progress(progbuf);
      if (rearrfirst)
      {
        sprintf(progbuf, " on the first of the trees tied for best\n");
      }
      else
      {
        sprintf(progbuf, " on all trees tied for best\n");
      }
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
      phyFillScreenColor();
    }

    done = false;                       /* tst */
    (void)done;                         // RSGnote: Variable set but never used.

    lastrearr = true;
    grandrearr(curtree, bestree, progress, rearrfirst);

    if (progress)
    {
      sprintf(progbuf, "\n");
      print_progress(progbuf);
    }
    if (jumb == njumble)
    {
      /*      collapsebestrees(curtree, bestrees, place, chars, progress); */
      TESTING('2');
      if (treeprint)
      {
        putc('\n', outfile);
        if (nextree == 1)
          fprintf(outfile, "One most parsimonious tree found:\n");
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
        protreroot(curtree->nodep[outgrno - 1]);
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
    if (treeprint)
    {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n\n\n");
    }
    which = 1;
    while (which <= numtrees)
    {
      prottreeread();
      if (outgropt)
        protreroot(curtree->nodep[outgrno - 1]);
      curtree->evaluate(curtree, curtree->root, 0);
      printree(curtree);
      describe();
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 )
      standev(chars, numtrees, minwhich, minsteps, nsteps, fsteps, seed);
  }

  if (jumb == njumble && progress)
  {
    sprintf(progbuf, "Output written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Trees also written onto file \"%s\".\n\n", outtreename);
      print_progress(progbuf);
    }
  }
}  /* maketree */


void protparsrun(void)
{
  // JRMdebug printout
  /*
    printf("\njumble: %i\n", jumble);
    printf("njumble: %li\n", njumble);
    printf("outgrno %li\n", outgrno);
    printf("outgropt %i\n", outgropt);
    printf("thresh %i\n", thresh);
    printf("trout %i\n", trout);
    printf("usertree %i\n", usertree);
    printf("weights %i\n", weights);
    printf("whichcode %i\n", whichcode);
    printf("printdata %i\n", printdata);
    printf("progress %i\n", progress);
    printf("treeprint %i\n", treeprint);
    printf("stepbox %i\n", stepbox);
    printf("ancseq %i\n", ancseq);
    printf("dotdiff %i\n", dotdiff);
    printf("interleaved %i\n", interleaved);
    fflush(stdout);
  */
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


void protpars(
  char * infilename,
  char * intreename,
  char * OutfileName,
  char * outfileopt,
  char * weightsfilename,
  char * OuttreeName,
  char * outtreeopt,
  int searchbest,
  int treesave,
  int inputorder,
  int RandNum,
  int NumJumble,
  int OutRoot,
  int OutNum,
  int ThreshPars,
  double ThreshVal,
  char * GeneCode,
  int SitesWeighted,
  int AnalyzeMult,
  int NumMult,
  int MultDataSet,
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

  (void)searchbest;                     // RSGnote: Parameter never used.
  //printf("Hello from ProtPars!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Protpars";

  memset(&funcs, 0, sizeof(initdata));
  funcs.node_new = protpars_node_new;
  funcs.tree_new = protpars_tree_new;

  phylipinit(argc, argv, &funcs, true);
  progname = argv[0];


  /*
  //jumble = false;
  //njumble = 1;
  //outgrno = 1;
  //outgropt = false;
  //thresh = false;
  //trout = true;
  //usertree = false;
  //weights = false;
  //whichcode = universal;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  //stepbox = false;
  //ancseq = false;
  //dotdiff = true;
  //interleaved = true;

  char * infilename,
  char * intreename,
  char * outfilename,
  char * outfileopt,
  char * weightfilename,
  char * outtreename,
  char * outtreeopt,
  int searchbest,
  //int treesave,
  //int inputorder,
  //int RandNum,
  //int NumJumble,
  //int OutRoot,
  //int OutNum,
  //int ThreshPars,
  //double ThreshVal,
  //char * GeneCode,
  //int SitesWeighted,
  //int AnalyzeMult,
  //int NumMult,
  //int MultDataSet,
  //int InputSeq,
  //int PrintData,
  //int DotDiff,
  //int PrintInd,
  //int PrintTree,
  //int PrintSteps,
  //int PrintSeq,
  //int WriteTree)

  */
  maxtrees = treesave;

  if (inputorder != 0)
  {
    jumble = true;
  }
  else
  {
    jumble = false;
  }

  inseed0 =  RandNum;

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

  if (ThreshPars != 0)
  {
    thresh = true;
  }
  else
  {
    thresh = false;
  }

  threshold = ThreshVal;

  if (!strcmp(GeneCode, "Universal"))
  {
    whichcode = universal;
  }
  else if (!strcmp(GeneCode, "Mitochondrial"))
  {
    whichcode = mito;
  }
  else if (!strcmp(GeneCode, "Vertebrate"))
  {
    whichcode = vertmito;
  }
  else if (!strcmp(GeneCode, "Fly"))
  {
    whichcode = flymito;
  }
  else //if (!strcmp(GeneCode, "Yeast"))
  {
    whichcode = yeastmito;
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

  garbage = NULL;
  firstset = true;
  code();
  setup();
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

  if (jumble)
  {
    initjumble(&inseed, &inseed0, seed, &njumble);
  }

  protparsrun();  // do the actual work

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
{  /* Protein parsimony by uphill search */
  initdata funcs;

#ifdef MAC
  argc = 1;         /* macsetup("Protpars", "");                */
  argv[0] = "Protpars";
#endif

  memset(&funcs, 0, sizeof(initdata));
  funcs.node_new = protpars_node_new;
  funcs.tree_new = protpars_tree_new;
  phylipinit(argc, argv, &funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  mulsets = false;
  msets = 1;
  firstset = true;
  code();
  setup();
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  protparsrun();

  FClose(infile);
  FClose(outfile);
  FClose(outtree);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif

  return 0;
}  /* Protein parsimony by uphill search */


#ifdef DEBUG
void debugtree(tree* t, FILE *f)
{
  int i, j, k;
  node *p, *q;
  fprintf(f, "debugtree:\n");
  fprintf(f, "spp=%ld\n", spp);
  fprintf(f, "nonodes=%ld\n", nonodes);
  fprintf(f, "\n");
  fprintf(f, "t->score=%f\n", t->score);
  if(t->root) fprintf(f, "t->root=%p (%ld)\n", t->root, t->root->index);
  fprintf(f, "t->spp=%ld\n", t->spp);
  fprintf(f, "t->nonodes=%ld\n", t->nonodes);
  fprintf(f, "\n");
  fprintf(f, "\n-------------------------------------\n");
  for(i=0;i<nonodes;i++)
  {
    p=t->nodep[i];
    if(p) fprintf(f, "nodep[%d]->index=%ld\n", i, p->index);
    if(i<spp && nayme && nayme[i]) fprintf(f, "nayme[%d]=%10.10s\n", i, nayme[i]);
    fprintf(f, "t->nodep[%d]=%p\n", i, p);
    if(p && p->nayme && p->nayme[0]) fprintf(f, "nodep[%d]->nayme=%s\n", i, p->nayme);
    if(p && p->back) fprintf(f, "nodep[%d]->back=%p (%ld)\n", i, p->back, p->back->index);
    q=p;
#if 0
    if(p && p->next) fprintf(f, "nodep[%d]->next=%p (%ld)\n", i, p->next, p->next->index);
    if(p && p->next && p->next->back)
    {
      fprintf(f, "nodep[%d]->next->back=%p (%ld)\n", i, p->next->back, p->next->back->index);
    }
    if(q) q=q->next;
    if(q) q=q->next;
    if(q)
    {
      for(j=2; q && q!=p && j<20; q=q->next, j++)
      {
        fprintf(f, "nodep[%d]", i);
        for(k=0; k<j; k++) fprintf(f, "->next");
        fprintf(f, "=%p (%ld)\n", q, q->index);
        if(!q->back) continue;
        fprintf(f, "nodep[%d]", i);
        for(k=0; k<j; k++) fprintf(f, "->next");
        fprintf(f, "->back=%p (%ld)\n", q->back, q->back->index);
      }
    }
#else
    if(p) q=p->next;
    if(q)
    {
      for(j=1; j<20; q=q->next, j++)
      {
        if(!q || q==p) break;
        if(!q->back) continue;
        fprintf(f, "nodep[%d]", i);
        for(k=0; k<j; k++) fprintf(f, "->next");
        fprintf(f, "->back=%p (%ld)\n", q->back, q->back->index);
      }
    }
    if(p) q=p->next;
    if(q)
    {
      for(j=1; j<20; q=q->next, j++)
      {
        fprintf(f, "nodep[%d]", i);
        for(k=0; k<j; k++)
          fprintf(f, "->next");
        if(q)
          fprintf(f, "=%p (%ld)\n", q, q->index);
        else
          fprintf(f, "=%p \n", q);
        if(!q || q==p)
          break;
      }
    }
#endif
    if(p)
      fprintf(f, "nodep[%d]->tip=%d\n", i, p->tip);
    fprintf(f, "\n-------------------------------------\n");
  }
}
#endif


// End.
