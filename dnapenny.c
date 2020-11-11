/* Version 4.0.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "seq.h"
#include "dnaparsimony.h"
#include <unistd.h>

extern long nextree;    /* (parsimony.c) */
extern long maxtrees;   /* (parsimony.c) max number trees to be printed out */
#define often          10000 /* how often to notify how many trees examined */
#define many           1000  /* how many multiples of howoften before stop  */

typedef node **pointptr;
typedef double *valptr;
typedef node **placeptr;

typedef struct dnapenny_tree {
  dnapars_tree dnapars_tree;
  long m;
  long n;
  valptr valyew;
  placeptr place;
}dnapenny_tree;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   makeweights(void);
void   doinput(void);
long   dnapenny_supplement(tree*, long);
void   addit(long );
void   dnapenny_reroot(node *);
void   describe(void);
void   maketree(void);
void   reallocchars(void);
tree*  dnapenny_tree_new(long nonodes, long spp);
void   dnapenny_tree_init(tree* t, long nonodes, long spp);
boolean dnapenny_tree_try_insert_(tree* t, node* item, node* p, node* where,
                                   double *bestyet, tree* dummy, boolean,
                                   boolean, boolean, double*);
long calculate_supplement(long i);
void nodeshellsort(double *a, node **b, long n);
void dnapennyrun(void);
void dnapenny(char * infilename, char * outfilename, char * outfileopt,
               char * weightsfilename, char * outtreename, char * outtreeopt,
               int TreeGroups, int ReportFreq, int BaBSimple, int OutRoot,
               int OutNum, int ThreshPars, double ThreshVal,
               int SitesWeighted, int AnalyzeMult, int NumMult,
               int MultDataSet, int InputSeq, int PrintData, int DotDiff,
               int PrintInd, int PrintTree, int PrintSteps, int PrintSeq,
               int WriteTree);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH],
      weightfilename[FNMLNGTH];
node *p;
long chars, howmanny, howoften, col, msets, ith, nonodes = 0;
boolean weights, thresh, simple, trout, progress, stepbox, ancseq, mulsets,
         firstset, justwts;
double threshold, bestfound;
steptr oldweight;
tree *curtree, *bestree;
dnapenny_tree* dcurtree;
double fracdone, fracinc;
boolean *added;
Char basechar[]="ACMGRSVTWYHKDBNO???????????????";

valptr *mbestval;
placeptr *mbestplace;
valptr *mvalyew;
placeptr *mplace;
long **suplements;

/* Variables for maketree, propagated globally for C version: */
long examined, mults;
boolean firsttime, done;
double like, bestyet;
bestelm *bestorders;
placeptr current;
long* order, *save;
double *threshwt;
long **supplements;
boolean *supplement_cache;
baseptr nothing;
long suppno[] = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                  2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4};

long suppset[] = {          /* extra steps needed if in this state at a tip */
  1 << ((long)A),
  1 << ((long)C),
  1 << ((long)G),
  1 << ((long)T),
  1 << ((long)O),
  (1 << ((long)A)) | (1 << ((long)C)),
  (1 << ((long)A)) | (1 << ((long)G)),
  (1 << ((long)A)) | (1 << ((long)T)),
  (1 << ((long)A)) | (1 << ((long)O)),
  (1 << ((long)C)) | (1 << ((long)G)),
  (1 << ((long)C)) | (1 << ((long)T)),
  (1 << ((long)C)) | (1 << ((long)O)),
  (1 << ((long)G)) | (1 << ((long)T)),
  (1 << ((long)G)) | (1 << ((long)O)),
  (1 << ((long)T)) | (1 << ((long)O)),
  (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)G)),
  (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)T)),
  (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)O)),
  (1 << ((long)A)) | (1 << ((long)G)) | (1 << ((long)T)),
  (1 << ((long)A)) | (1 << ((long)G)) | (1 << ((long)O)),
  (1 << ((long)A)) | (1 << ((long)T)) | (1 << ((long)O)),
  (1 << ((long)C)) | (1 << ((long)G)) | (1 << ((long)T)),
  (1 << ((long)C)) | (1 << ((long)G)) | (1 << ((long)O)),
  (1 << ((long)C)) | (1 << ((long)T)) | (1 << ((long)O)),
  (1 << ((long)G)) | (1 << ((long)T)) | (1 << ((long)O)),
  (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G))|(1 << ((long)T)),
  (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G))|(1 << ((long)O)),
  (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)T))|(1 << ((long)O)),
  (1 << ((long)A))|(1 << ((long)G))|(1 << ((long)T))|(1 << ((long)O)),
  (1 << ((long)C))|(1 << ((long)G))|(1 << ((long)T)) | (1 << ((long)O)),
  (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G)) | (1 << ((long)T)) | (1 << ((long)O))};


tree* dnapenny_tree_new(long nonodes, long spp)
{
  /* allocate a new tree and call function to initialize it */
  tree *t = Malloc(sizeof(dnapenny_tree));
  generic_tree_init(t, nonodes, spp);
  dnapenny_tree_init(t, nonodes, spp);
  return t;
} /* dnapenny_tree_new */


void dnapenny_tree_init(tree* t, long nonodes, long spp)
{
  /* initialize variables in a new tree */

  dnapars_tree_init(t, nonodes, spp);
  t->try_insert_ = dnapenny_tree_try_insert_;
  ((pars_tree*)t)->supplement = dnapenny_supplement;
} /* dnapennt_tree_init */


void dnapenny_tree_setup(long nonodes, long spp)
{
  /* create and initialize the necessary trees */

  curtree = dnapenny_tree_new(nonodes, spp);
} /* proml_tree_setup */


void nodeshellsort(double *a, node **b, long n)
{
  /* Shell sort keeping a, b in same order
   * used by dnapenny, dolpenny, penny, and threshml */
  long gap, i, j;
  node* itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0)
  {
    for (i = gap + 1; i <= n; i++)
    {
      j = i - gap;
      while (j > 0)
      {
        if (a[j - 1] > a[j + gap - 1])
        {
          rtemp = a[j - 1];
          a[j - 1] = a[j + gap - 1];
          a[j + gap - 1] = rtemp;
          itemp = b[j - 1];
          b[j - 1] = b[j + gap - 1];
          b[j + gap - 1] = itemp;
        }
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* nodeshellsort */


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  howoften = often;
  howmanny = many;
  outgrno = 1;
  outgropt = false;
  simple = true;
  thresh = false;
  threshold = spp;  //???
  trout = true;
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
    printf("\nPenny algorithm for DNA, version %s\n", VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  H        How many groups of %4ld trees:%6ld\n", howoften, howmanny);
    printf("  F        How often to report, in trees:  %4ld\n", howoften);
    printf("  S           Branch and bound is simple?  %s\n", (simple ?  "Yes" : "No. reconsiders order of species"));
    printf("  O                        Outgroup root?  %s%3ld\n", (outgropt ? "Yes, at sequence number" :
            "No, use as outgroup species"), outgrno); printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  W                       Sites weighted?  %s\n", (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets, (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n", (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n", (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n", (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n", (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n", (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n", (stepbox ? "Yes" : "No" ));
    printf("  5  Print sequences at all nodes of tree  %s\n", (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n", (trout ? "Yes" : "No"));
    if(weights && justwts)
    {
      printf("WARNING:  W option and Multiple Weights options are both on.  ");
      printf("The W menu option is unnecessary and has no additional effect. \n");
    }
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    uppercase(&ch);
    if ((strchr("WHMSOFTI1234560", ch)) != NULL)
    {
      switch (ch)
      {
        case 'H':
          inithowmany(&howmanny, howoften);
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
          }
          break;

        case 'F':
          inithowoften(&howoften);
          break;

        case 'S':
          simple = !simple;
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          else
            outgrno = 1;
          break;

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
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


void allocrest(void)
{
  long i;

  inputSequences = (Char **)Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));
  weight = (long *)Malloc(chars * sizeof(long));
  oldweight = (long *)Malloc(chars * sizeof(long));
  alias = (steptr)Malloc(chars * sizeof(long));
  ally = (steptr)Malloc(chars * sizeof(long));
  location = (steptr)Malloc(chars * sizeof(long));
  nayme = (naym *)Malloc(spp * sizeof(naym));
  bestorders = (bestelm *)Malloc(maxtrees * sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestorders[i - 1].btree = (long*)Malloc(spp * sizeof(long));
  current = (placeptr)Malloc(spp * sizeof(node*));
  order = (long*)Malloc(nonodes * sizeof(long));
  added = (boolean *)Malloc(nonodes * sizeof(boolean));
  save = (long*)Malloc(nonodes * sizeof(long));
  supplements = (long**)Malloc(spp * sizeof(long*));
  supplement_cache = (boolean*)Malloc(spp * sizeof(boolean));
  mvalyew = (valptr*)Malloc(spp * sizeof(valptr));
  mbestval = (valptr*)Malloc(spp * sizeof(valptr));
  mplace = (placeptr*)Malloc(spp * sizeof(placeptr));
  mbestplace = (placeptr*)Malloc(spp * sizeof(mbestplace));
  for ( i = 0 ; i < spp ; i++ )
  {
    mvalyew[i] = (valptr)Malloc(nonodes * sizeof(double));
    mbestval[i] = (valptr)Malloc(nonodes * sizeof(double));
    mplace[i] = (placeptr)Malloc(nonodes * sizeof(node*));
    mbestplace[i] = (placeptr)Malloc(nonodes * sizeof(node*));
  }
}  /* allocrest */


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
  free(alias);
  free(ally);
  free(location);

  weight = (long *)Malloc(chars * sizeof(long));
  oldweight = (long *)Malloc(chars * sizeof(long));
  alias = (steptr)Malloc(chars * sizeof(long));
  ally = (steptr)Malloc(chars * sizeof(long));
  location = (steptr)Malloc(chars * sizeof(long));
} /* reallocchars */


void doinit(void)
{
  /* initializes variables */
  fprintf(outfile, "\nPenny algorithm for DNA, version %s\n", VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");

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
      for ( i = 0 ; i < spp ; i++ )
        free(supplements[i]);
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
  for ( i = 0 ; i < spp ; i++ )
  {
    supplement_cache[i] = false;
    supplements[i] = (long*)Malloc(endsite  * sizeof(long));
  }
  dnapenny_tree_setup(nonodes, spp);
  dcurtree = (dnapenny_tree*)curtree;
  dna_makevalues(curtree, false);
}  /* doinput */


long dnapenny_supplement(tree* t, long i)
{
  long j;
  long m = ((dnapenny_tree*)t)->m;

  if ( simple && supplement_cache[m - 1] == true )
    return supplements[m - 1][i];

  if( !simple)
    return calculate_supplement(i);

  for ( j = 0 ; j < endsite ; j++ )
  {
    supplements[m-1][j] = calculate_supplement(j);
  }
  supplement_cache[m - 1] = true;
  return supplements[m - 1][i];
}


long calculate_supplement(long i)
{
  /* determine minimum number of steps at site i which will
     be added when rest of species are put in tree */
  long j, k, has, sum;
  boolean addedmayhave, nonaddedhave;

  nonaddedhave = 0;;
  addedmayhave = 0;
  for (k = 0; k < spp; k++)
  {
    has = (((dnapars_node*)curtree->nodep[k])->base[i]);
    if (has != 31)
    {
      if (added[k])
        addedmayhave |= has;
      else
      {
        if ((has == 1) || (has == 2) || (has == 4) || (has == 8) || (has == 16))
          nonaddedhave |= has;
      }
    }
  }
  sum = 0;
  j = 1;
  for (k = 1; k <= 5; k++)
  {
    if ((j & nonaddedhave) != 0)
      if ((j & addedmayhave) == 0)
        sum++;
    j += j;
  }
  return(sum * weight[i]);
}  /* supplement */


boolean dnapenny_tree_try_insert_(tree* t, node* item, node* p, node* where, double *bestyet, tree* dummy, boolean thorough, boolean dummy3, boolean multf, double* bestfound)
{
  node* qwhere;
  double like;
  boolean root = false;
  long i;

  if (done)
    return true;
  if (true ||((dnapenny_tree*)t)->m <= 2 )
  {
    if ( p->back == NULL )
    {
      hookup(item, p);
      root= true;
    }
    else
      t->insert_(t, item, p, false);
    for ( i = 0 ; i < nonodes ; i++ )
    {
      t->nodep[i]->initialized = false;
      if ( i >= spp )
      {
        t->nodep[i]->next->initialized = false;
        t->nodep[i]->next->next->initialized = false;
      }
    }
    ((dnapenny_tree*)t)->n++;
    like = t->evaluate(t, item, 0);
    if ( !root)
      t->re_move(t, item, &qwhere, false);
    else
    {
      t->root->back = NULL;
      item->back = NULL;
    }
    examined++;
    if (examined == howoften)
    {
      examined = 0;
      mults++;
      if (mults == howmanny)
        done = true;
      if (progress)
      {
        sprintf(progbuf, "%7ld", mults);
        print_progress(progbuf);
        if (*bestyet <= 0)
        {
          sprintf(progbuf, "%16.1f", -*bestyet);
        }
        else
        {
          sprintf(progbuf, "            -   ");
        }
        print_progress(progbuf);
        sprintf(progbuf, "%17ld%20.2f\n", nextree - 1, fracdone * 100);
        print_progress(progbuf);
        phyFillScreenColor();
      }
    }
    ((dnapenny_tree*)t)->valyew[(((dnapenny_tree*)t)->n) - 1] = -like;

    ((dnapenny_tree*)t)->place[(((dnapenny_tree*)t)->n) - 1] = p;
  }
  return true;
} /* dnapenny_tree_try_insert */


void addit(long m)
{
  /* adds the species one by one, recursively */
  long i, j, n1, besttoadd=0;
  double oldfrac, oldfdone, sum, bestsum;
  valptr bestval;
  placeptr bestplace;
  valptr valyew;
  placeptr place;
  boolean multf;
  boolean root=false;
  node *dum, *qwhere;
  long n = 0;                           // RSGnote: Formerly not initialized and later potentially referenced before before being set.

  valyew = mvalyew[m-1];
  bestval = mbestval[m-1];
  place = mplace[m-1];
  bestplace = mbestplace[m-1];
  if (simple && !firsttime)
  {
    n = 0;
    added[order[m-1] - 1] = true;
    dcurtree->valyew = valyew;
    dcurtree->place = place;
    dcurtree->n = n;
    curtree->addtraverse(curtree, curtree->nodep[order[m - 1] - 1],
     curtree->root, true, qwhere, &bestyet, bestree, true, true, true, &bestfound);
    n = dcurtree->n;
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes * sizeof(long));
    memcpy(bestval, valyew, nonodes * sizeof(double));
  }
  else
  {
    bestsum = UNDEFINED;
    for (i = 1; i <= spp; i++)
    {
      if (!added[i - 1])
      {
        n = 0;
        added[i - 1] = true;
        dcurtree->m = m;
        dcurtree->valyew = valyew;
        dcurtree->place = place;
        dcurtree->n = n;
        curtree->addtraverse(curtree, curtree->nodep[i - 1], curtree->root,
            true, qwhere, &bestyet, bestree, true, true, multf, &bestfound);
/* debug
boolean	pars_addtraverse(tree*, node*, node*, boolean, node*, double*,
                       bestelm*, boolean, boolean, boolean, double*);
boolean generic_tree_addtraverse(tree* t, node* p, node* q, boolean contin,
                           node* qwherein, double* bestyet,
                           tree* bestree, boolean thorough, boolean storing,
                           boolean atstart, double* bestfound)
debug  */
        n = dcurtree->n;
        added[i - 1] = false;
        sum = 0.0;
        for (j = 0; j < n; j++)
          sum += -valyew[j];
        if (sum < bestsum || bestsum == UNDEFINED)
        {
          bestsum = sum;
          besttoadd = i;
          memcpy(bestplace, place, nonodes * sizeof(long));
          memcpy(bestval, valyew, nonodes * sizeof(double));
        }
      }
    }
  }
  order[m - 1] = besttoadd;
  memcpy(place, bestplace, nonodes * sizeof(long));
  memcpy(valyew, bestval, nonodes * sizeof(double));
  nodeshellsort(valyew, place, n);
  oldfrac = fracinc;
  oldfdone = fracdone;
  n1 = 0;
  for (i = 0; i < n; i++)
  {
    if (-valyew[i] >= bestyet || bestyet == UNDEFINED)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < n; i++)
  {
    if (-valyew[i] >= bestyet || bestyet > 0.0 )
    {
      current[m - 1] = place[i];
      if ( place[i]->back == NULL )
      {
        hookup(place[i], curtree->nodep[besttoadd - 1]);
        root = true;
      }
      else
        curtree->insert_(curtree, curtree->nodep[besttoadd - 1],
                          place[i], false);
      added[besttoadd - 1] = true;
      if (m < spp)
      {
        addit(m + 1);
      }
      else
      {
        if (-valyew[i] > bestyet || bestyet > 0.0)
        {
          nextree = 1;
          bestyet = -valyew[i];
        }
        if (nextree <= maxtrees)
        {
          savetree(curtree, save);
          memcpy(bestorders[nextree - 1].btree, save, spp * sizeof(long));
        }
        nextree++;
        firsttime = false;
      }
      if ( !root)
        curtree->re_move(curtree, curtree->nodep[besttoadd - 1], &dum, true);
      else
      {
        curtree->nodep[besttoadd - 1]->back->back = NULL;
        curtree->nodep[besttoadd - 1]->back = NULL;
      }
      added[besttoadd - 1] = false;
    }
    fracdone += fracinc;
  }
  fracinc = oldfrac;
  fracdone = oldfdone;

}  /* addit */


void dnapenny_reroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom;

  if (outgroup->back->index == curtree->root->index)
    return;
  newbottom = outgroup->back;
  p = curtree->nodep[newbottom->index - 1]->back;
  while (p->index != curtree->root->index)
  {
    oldbottom = curtree->nodep[p->index - 1];
    curtree->nodep[p->index - 1] = p;
    p = oldbottom->back;
  }
  p = curtree->root->next;
  q = curtree->root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = curtree->root->next->next;
  outgroup->back = curtree->root->next;
  curtree->nodep[newbottom->index - 1] = newbottom;
}  /* dnapenny_reroot */


void describe(void)
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  long indent;

  if (treeprint)
  {
    fprintf(outfile, "\n  between      and       length\n");
    fprintf(outfile, "  -------      ---       ------\n");
    printbranchlengths(curtree->root);
  }
  if (stepbox)
    writesteps(curtree, chars, weights, oldweight);
  if (ancseq)
  {
    dna_hypstates(curtree, chars, basechar);
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


void maketree(void)                     // RSGbugfix
{
  /* tree construction recursively by branch and bound */
  long i, k;

  if (progress)
  {
    sprintf(progbuf, "\nHow many\n");
    print_progress(progbuf);
    sprintf(progbuf, "trees looked                                       Approximate\n");
    print_progress(progbuf);
    sprintf(progbuf, "at so far      Length of        How many           percentage\n");
    print_progress(progbuf);
    sprintf(progbuf, "(multiples     shortest tree    trees this short   searched\n");
    print_progress(progbuf);
    sprintf(progbuf, "of %4ld):      found so far     found so far       so far\n", howoften);
    print_progress(progbuf);
    sprintf(progbuf, "----------     ------------     ------------       ------------\n");
    print_progress(progbuf);
  }
  phyFillScreenColor();
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  curtree->root = curtree->nodep[0];
  firsttime = true;
  for (i = 0; i < spp; i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = UNDEFINED;

  addit(k);

  if (done)
  {
    if (progress)
    {
      sprintf(progbuf, "Search broken off!  Not guaranteed to\n");
      print_progress(progbuf);
      sprintf(progbuf, " have found the most parsimonious trees.\n");
      print_progress(progbuf);
      sprintf(progbuf, " To ensure finding the most parsimonious trees,\n");
      print_progress(progbuf);
      sprintf(progbuf, " try runs which have larger numbers set in the\n");
      print_progress(progbuf);
      sprintf(progbuf, " H and F menu options.\n\n");
      print_progress(progbuf);
    }
    if (treeprint)
    {
      fprintf(outfile, "Search broken off!  Not guaranteed to\n");
      fprintf(outfile, " have found the most parsimonious\n");
      fprintf(outfile, " trees, but here is what we found.\n");
      fprintf(outfile, " To ensure finding the most parsimonious trees,\n");
      fprintf(outfile, " try runs which have larger numbers set in the\n");
      fprintf(outfile, " H and F menu options.\n");
    }
  }

  long outCount = 0;
  collapsebestrees(curtree, bestorders, order, chars, progress, &outCount);
  long missedCount = nextree - 1 - maxtrees;
  assert(outCount > 0);
  if (treeprint)

  {
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", -bestyet );
    if (outCount == 1)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
    {
      if (missedCount > 0)
      {
        fprintf(outfile, "as many as %ld trees may have been found\n", missedCount + outCount);
        fprintf(outfile, "here are the first %4ld of them\n", outCount );
      }
      else
      {
        fprintf(outfile, "%6ld trees in all found\n", outCount);
      }
    }
  }

  if (treeprint)
    putc('\n', outfile);

  for (i = 0; i < spp; i++)
    added[i] = true;

  for (i = 0; i < outCount; i++)
  {
    load_tree(curtree, i, bestorders);
    like = curtree->evaluate(curtree, curtree->root, 0);
    curtree->root = root_tree(curtree, curtree->nodep[outgrno - 1]->back);
    dna_treelength(curtree->root, chars, curtree->nodep);
    printree(curtree);
    describe();
    reroot_tree(curtree, curtree->root); // RSGbugfix: Name change.
  }

  if (progress)
  {
    sprintf(progbuf, "\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Trees also written onto file \"%s\".\n\n", outtreename);
      print_progress(progbuf);
    }
  }
}  /* maketree */


void dnapennyrun(void)
{
  // JRM debug
  /*
    printf("\nhowoften %li\n", howoften);
    printf("howmanny %li\n", howmanny);
    printf("outgrno %li\n", outgrno);
    printf("outgropt %i\n", outgropt);
    printf("simple %i\n", simple);
    printf("thresh  %i\n", thresh);
    printf("threshold %f\n", threshold);
    printf("trout %i\n", trout);
    printf("weights %i\n", weights);
    printf("justwts %i\n", justwts);
    printf("printdata %i\n", printdata);
    printf("dotdiff %i\n", dotdiff);
    printf("progress %i\n", progress);
    printf("treeprint %i\n", treeprint);
    printf("stepbox %i\n", stepbox);
    printf("ancseq %i\n", ancseq);
    printf("interleaved %i\n", interleaved);
  */
  //do the work
  for (ith = 1; ith <= msets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (msets > 1 && !justwts) {
      fprintf(outfile, "\nData set # %ld:\n", ith);
      if (progress)
      {
        sprintf(progbuf, "\nData set # %ld:\n", ith);
        print_progress(progbuf);
      }
    }
    maketree();
    fflush(outfile);
    fflush(outtree);
    free(threshwt);
  }
}


void dnapenny(
  char * infilename,
  char * OutfileName,
  char * outfileopt,
  char * weightsfilename,
  char * OuttreeName,
  char * outtreeopt,
  int TreeGroups,
  int ReportFreq,
  int BaBSimple,
  int OutRoot,
  int OutNum,
  int ThreshPars,
  double ThreshVal,
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

  //printf("Hello from DnaPenny!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Dnapenny";
  maxtrees = 10000;

  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = dnapars_node_new;
  funcs.tree_new = dnapenny_tree_new;
  phylipinit(argc, argv, &funcs, true);

  /*
  //howoften = often;
  //howmanny = many;
  //outgrno = 1;
  //outgropt = false;
  //simple = true;
  //thresh = false;
  //threshold = spp;
  //trout = true;
  //weights = false;
  //justwts = false;
  //printdata = false;
  //dotdiff = true;
  //progress = true;
  //treeprint = true;
  //stepbox = false;
  //ancseq = false;
  //interleaved = true;

  //char * infile,
  //char * outfile,
  //char * outfileopt,
  //char * weightfile,
  //char * outtree,
  //char * outtreeopt,
  //int TreeGroups,
  //int ReportFreq,
  //int BaBSimple,
  //int OutRoot,
  //int OutNum,
  //int ThreshPars,
  //double ThreshVal,
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
  howoften = ReportFreq;
  howmanny = TreeGroups;

  if (BaBSimple != 0)
  {
    simple = true;
  }
  else
  {
    simple = false;
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

  if (ThreshPars != 0)
  {
    thresh = true;
  }
  else
  {
    thresh = false;
  }

  threshold = ThreshVal;

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

  dnapennyrun();

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

  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* Penny's branch-and-bound method for DNA sequences */
  initdata funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Dnapenny", "");        */
  argv[0] = "Dnapenny";
#endif
  maxtrees = 10000;

  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = dnapars_node_new;
  funcs.tree_new = dnapenny_tree_new;
  phylipinit(argc, argv, &funcs, false);

  /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  msets = 1;
  firstset = true;
  doinit();
  if (weights || justwts)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);

  dnapennyrun();

  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  phyRestoreConsoleAttributes();
  return 0;
}  /* Penny's branch-and-bound method for DNA sequences */


// End.
