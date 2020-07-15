/* Version 4.0. (c) Copyright 2009-2013 by the University of Washington.
   Written by Joseph Felsenstein, Lucas Mix, Elizabeth Walkup, Eric Rynes,
   Akiko Fuseki, Sean Lamont, Andrew Keeffe, Dan Fineman,
   and Patrick Colacurcio.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "seq.h"
#include "ml.h"

typedef long vall[maxcategs];
typedef double contribarr[maxcategs];

#if 0
typedef double eigCalcT;
eigCalcT QR_accuracy    = 1e-6;
#endif

typedef long double eigCalcT;
eigCalcT QR_accuracy    = 1e-24;

double LIKE_ACCURACY  =  1e-8;

typedef enum {
  universal, ciliate, mito, vertmito, flymito, yeastmito
} codetype;

typedef struct codon_tree {
  ml_tree ml_tree;
} codon_tree;

#define NUM_ALL_CODONS          64
#define NUM_SENSE_CODONS        61
#define NUM_NUCLEOTIDE_STATES    4
#define NUM_NUCS_IN_CODON        3

//#define JAVADEBUG

#ifndef OLDC
/* function prototypes */
void   debugtree(tree*);
tree * codon_tree_new(long, long);
void   codon_tree_init(tree*, long, long);
void   init_nucSubRates(double);
void   fill_codon64_to(long, long, long, double);
void   init_codon64matrix(double);
void   matrix64to61(double **, double ***);
void   line64to61(double*, double*);
void   codml_empiricalfreqs(double *, double *, double *, double *, steptr, pointarray);
void   initomega(double*);
void   getoptions(void);
void   free_pmatrix(long);
void   alloc_pmatrix(long);
void   codon_freetable(void);
void   codon_inittable(void);
void   reallocsites(long, long);
void   allocrest(long, long);
void   code(void);
void   doinit(void);
void   inputoptions(void);
void   input_codondata(long);
void   makeweights(boolean);
void   aa_codon_makevalues(long, pointarray, long, long, sequence, steptr);
void   nuc_codon_makevalues(long, pointarray, long, long, sequence, steptr);
void   getcodonbasefreqs (void);
void   init_SymMatrix(void);
boolean is_symmetric(double **, long);
boolean codonfreqs_sum_to_one(void);
void   makecodonfreqs(void);
void   freelrsaves(void);
void   codml_alloclrsaves(void);
boolean is_f_times_m_symmetric(double **);
void   givens(eigCalcT **, long, long, long, eigCalcT, eigCalcT, boolean);
void   coeffs(eigCalcT, eigCalcT , eigCalcT *, eigCalcT *);
void   tridiag(eigCalcT **, eigCalcT **, long);
void   shiftqr(eigCalcT **, eigCalcT **, long );
void   qreigen(double **, long, long);
void   maketables(void);
void   getinput(void);
void   codon_tree_calc_nuview(tree *, node *, cphenotype, double *);
void   codon_tree_nuview(tree *, node *);
void   codon_slopecurv(node *, double, double *, double *, double *);
void   codon_tree_makenewv(tree *, node *);
void   codml_restore_traverses(tree *t, node *p, node * q);
boolean codml_node_good(tree *t, node *n);
void   make_pmatrix(double **, double **, double **, long, double lz, double, double *, double **);
double codon_tree_evaluate(tree*, node *, boolean);
void   codmlcopy(tree *, tree *, long, long);
void   codml_coordinates(node *, double, long *, double *);
void   codml_printree(void);
void   sigma(node *, double *, double *, double *);
void   describe(node *);
void   codon_reconstr(node *, long);
void   rectrav(node *, long, long);
void   summarize(void);
void   dnaml_treeout(node *);
void   free_all_codonx (long, pointarray);
void   codml_reroot(tree *);            // RSGbugfix: Name change.
void   initcodonmlnode(tree *, node **, long, long, long *, long *, initops, pointarray, Char *, Char *, FILE *);
void   maketree (void);
void   clean_up(void);
void   codmlrun(void);
void   codml(char * Infilename, char * Intreename, char * Wgtsfilename, char * Outfilename,
             char * outfileopt, char * Outtreename, char * outtreeopt, char * TreeUseMethod,
             int UseLengths, int InputNuc, char * NucSubModel, double TTratio, double NSratio,
             int useEmpBF, double BaseFreqA, double BaseFreqC, double BaseFreqG, double BaseFreqTU,
             char * GeneCode, int OneOmega, int AdjOmegasCor, int OmegaBlockLen, int CodonsWgted,
             int NumOmegas, double OmegaVal1, double OmegaVal2, double OmegaVal3, double OmegaVal4,
             double OmegaVal5, double OmegaVal6, double OmegaVal7, double OmegaVal8, double OmegaVal9,
             double OmegaProb1, double OmegaProb2, double OmegaProb3, double OmegaProb4, double OmegaProb5,
             double OmegaProb6, double OmegaProb7, double OmegaProb8, double OmegaProb9, int SpeedAn,
             int GlobalRe, int RandInput, int RandNum, int Njumble, int OutRoot, int OutNum,
             int MultData, int MultDSet, int NumSeqs, int InputSeq, int PrintData, int PrintInd,
             int PrintTree, int DotDiff, int WriteTree, int RecHypo);
/* function prototypes */
#endif

boolean haslengths;

Char infilename[100], outfilename[100], intreename[100], outtreename[100], catfilename[100], weightfilename[100];
double *omegavalue, *probcat;
long nonodes2, numsensecodons, dataSites, workingSites, weightsum, rcategs, datasets, ith, njumble, jumb = 0;
long inseed, inseed0, parens;
boolean f84, hky, k2p, jc, global, jumble, weights, trout, usertree, reusertree, rctgry, auto_, freqsfrom, hypstate,
  progress, mulsets, justwts, firstset, improve, smoothit, polishing, lngths, gama, invar, inputDataIsNucleotideData;
tree *curtree, *bestree, *bestree2, *priortree;
double probsum, cv, alpha, lambda, invarfrac, ttratio, omega;
long *enterorder;
steptr aliasweight, weight;
contribarr *contribution, like, nulike, clai;
double **basesincodon, **term, **slopeterm, **curveterm;
longer seed;
char *progname;
codetype whichcode = universal;

double ** codon64matrix;    /* [NUM_ALL_CODONS][NUM_ALL_CODONS] */
double    nucSubRate[NUM_NUCLEOTIDE_STATES][NUM_NUCLEOTIDE_STATES];
double    nucFreq[NUM_NUCLEOTIDE_STATES];


/* Standard Codon Table based on indices       0123 <=> ACGT */
/*     which can be modified in the code for other models    */
/* Note that protdist uses a trans[][][] with  0123 <=> TCAG */

Char codonStrings[NUM_SENSE_CODONS][NUM_NUCS_IN_CODON+1] =
{"AAA", "AAC", "AAG", "AAT",
 "ACA", "ACC", "ACG", "ACT",
 "AGA", "AGC", "AGG", "AGT",
 "ATA", "ATC", "ATG", "ATT",
 "CAA", "CAC", "CAG", "CAT",
 "CCA", "CCC", "CCG", "CCT",
 "CGA", "CGC", "CGG", "CGT",
 "CTA", "CTC", "CTG", "CTT",
 "GAA", "GAC", "GAG", "GAT",
 "GCA", "GCC", "GCG", "GCT",
 "GGA", "GGC", "GGG", "GGT",
 "GTA", "GTC", "GTG", "GTT",
 "TAC", "TAT",
 "TCA", "TCC", "TCG", "TCT",
 "TGC", "TGG", "TGT",
 "TTA", "TTC", "TTG", "TTT"};

aas transx[4][4][4] = {{{lys, asn , lys , asn}, {thr , thr , thr , thr },
                        {arg , ser , arg , ser}, {ileu, ileu, met , ileu} },
                       { {gln , his , gln , his}, {pro , pro , pro , pro },
                         {arg , arg , arg , arg}, {leu , leu , leu , leu } },
                       { {glu , asp , glu , asp}, {ala , ala , ala , ala },
                         {gly , gly , gly , gly}, {val , val , val , val } },
                       { {stop, tyr , stop, tyr}, {ser , ser , ser , ser },
                         {stop, cys , trp , cys}, {leu , phe , leu , phe } } };


/* Local variables for maketree, propagated globally for C version: */
long k, nextsp, numtrees, maxwhich, mx, mx0, mx1, shimotrees;
double dummy, maxlogl;
boolean smoothed = false;
double **l0gf;
double *l0gl;
double *tbl;
Char ch, ch2;
long col;
vall *mp;


/* Variables introduced to allow for protein probability calculations   */
long   max_num_sibs;            /* maximum number of siblings used in a */
                                /* nuview calculation.  determines size */
                                /* final size of pmatrices              */
double **eigmat;                /* eig matrix variable                  */
double ***probmat;              /* prob matrix variable                 */
double ***dpmatrix;            /* derivative of pmatrix                 */
double ***ddpmatrix;           /* derivative of pmatrix                 */
double ****pmatrices;           /* matrix of probabilities of codon     */
                                /* conversion.  The 4 subscripts refer  */
                                /* to sibs, rcategs, final and  */
                                /* initial states, respectively.        */
double codonfreq[NUM_SENSE_CODONS];    /* codon frequencies         */
double **  InstMatrix;   /* Instantaneous Matrix */
double **  SymMatrix;    /* Symmetric Matrix     */


void debugtree(tree* t)
{
  int i;
  node *p;

  printf("DEBUGTREE spp = %ld nonodes2 = %ld ", spp, nonodes2);
  if(t->root)
  {
    printf("root = %p (%ld)\n", (void *)t->root, t->root->index);
  }
  else
  {
    printf("UNROOTED\n");
  }

  // could add here t->spp ; t->nonodes2 ; t->score

  for(i=0;i<nonodes2;i++)
  {
    printf("-------------------------------------\n");
    p=t->nodep[i];
    if(!p)
    {
      printf("%3s nodep[%d] == NULL\n", " ", i);
      continue;
    }

    printf("%3s nodep[%d] == %p ->index = %ld ", p->tip > 0 ? "TIP" : "", i, (void *)p, p->index);
    if(i<spp && nayme && nayme[i]) printf("nayme[%d]=%10.10s", i, nayme[i]);
    printf("\n");

    if(p != NULL)
    {
      boolean firstTime=true;
      boolean nulledOut=false;
      node * firstNode = p;

      while((firstTime || (p != firstNode)) && !nulledOut)
      {
        if (p == NULL)
        {
          if(! firstNode->tip) printf("%3s p:NULL\n", "");
          nulledOut = true;
        }
        else
        {
          printf("%3s p:%p b:%p ", "", (void *)p, (void *)p->back);
          if(p->back != NULL)
          {
            printf("(%ld)", p->back->index);
          }
          printf(" branch is %e", p->v);
          printf(" init'd is %d", p->initialized);
          printf(" iter is %d", p->iter);
          long siteIndex;

          if (!(p->back == NULL))
            //if (!p->tip && !(p->back == NULL))
          {
            while(false)
              // for(siteIndex=0;siteIndex < 1; siteIndex++)
              // for(siteIndex=0;siteIndex < endsite; siteIndex++)
            {
              long codonIndex;
              double highestVal = 0.0;
              double sumVal = 0.0;
              for(codonIndex=0;codonIndex < NUM_SENSE_CODONS;codonIndex++)
              {
                double val = ((codon_node*)p)->codonx[siteIndex][0][codonIndex];
                sumVal += val;
                if(val > highestVal) highestVal = val;
              }
              if(highestVal > 0.0)
              {
                printf("\n        codonx for %2ld: ", alias[siteIndex]);
                printf("sum: %lf ", sumVal);
                for(codonIndex=0;codonIndex < NUM_SENSE_CODONS;codonIndex++)
                {
                  double val = ((codon_node*)p)->codonx[siteIndex][0][codonIndex];
                  if(val > highestVal*0.01)
                    //if(val > 0.0)
                    printf("%2ld(%8lf) ", codonIndex, val);
                }
              }
            }
          }
          printf("\n");
          p = p->next;
        }
        firstTime = false;
      }
    }
  }
} /* debugtree */


tree* codon_tree_new(long nonodes, long spp)
{
  /* set up variables and then set up identities of functions */
  tree* t = Malloc(sizeof(codon_tree));

  t = generic_tree_new(nonodes, spp);
  ml_tree_init(t, nonodes, spp);
  codon_tree_init(t, nonodes, spp);
  return t;
} /* codon_tree_new */



void codon_tree_init(tree* t, long nonodes, long spp)
{
  /* set up functions for a proml_tree */

  t->evaluate = codon_tree_evaluate;
  t->try_insert_ = ml_tree_try_insert_;
  t->nuview = codon_tree_nuview;
  t->makenewv = codon_tree_makenewv;
  t->tree_print_f = debugtree;
  t->restore_traverses = codml_restore_traverses;
  t->node_good_f = codml_node_good;
} /* codon_tree_init */


void codon_tree_setup(long nonodes, long spp)
{
  /* create and initialize the necessary trees */

  curtree = codon_tree_new(nonodes, spp);
  bestree = codon_tree_new(nonodes, spp);
  bestree2 = codon_tree_new(nonodes, spp);
  priortree = codon_tree_new(nonodes, spp);
} /* codon_tree_setup */


void init_nucSubRates(double ttratio)
{
  /* establishes instantaneous rates of change for nucleotides */
  /* on the basis of the F84, HKY, K2P, or JC models.                            */
  long i, j;
  double freqr;                         /* purine freq         */
  double freqy;                         /* pyrimidine freq     */
  double freqag, freqct;                /* freqa*freqg...      */
  double bbb, aa1, aa2, aaa, aaar, rho;

  freqr = nucFreq[0] + nucFreq[2];  /* from pages 201-202 of ... */
  freqy = nucFreq[1] + nucFreq[3];  /* ... "Inferring Phylogenies" */
  freqag = nucFreq[0] * nucFreq[2];
  freqct = nucFreq[1] * nucFreq[3];
  if (f84 || k2p || jc)  /* doesn't matter for K2P and JC anyway */
    rho = 1.0;
  else
    rho = freqr /freqy;
  if (jc)
  {
    aaar = 0.0;
    aaa = 0.0;
    bbb =  4.0/3.0;
  }
  else
  {
    aa1 = (freqr * freqy * ttratio) - freqag - freqct;
    aa2 = 2.0 * (1.0 + ttratio) * ((freqy*freqag)*rho + (freqr*freqct));
    aaa = aa1 / aa2;
    aaar = rho * aaa;
    bbb = 1.0 / (2.0*freqr*freqy*(1.0 + ttratio));
  }

  // BUG.969
  // how to give error message for this and restart properly
  // if the ttratio/freqs combination is bad? see other programs.
  // assert(aa1 >= 0.0);

  nucSubRate[0][1] = bbb*nucFreq[1];
  nucSubRate[0][2] = bbb*nucFreq[2] + aaar*nucFreq[2]/freqr;
  nucSubRate[0][3] = bbb*nucFreq[3];

  nucSubRate[1][0] = bbb*nucFreq[0];
  nucSubRate[1][2] = bbb*nucFreq[2];
  nucSubRate[1][3] = bbb*nucFreq[3] + aaa*nucFreq[3]/freqy;

  nucSubRate[2][0] = bbb*nucFreq[0] + aaar*nucFreq[0]/freqr;
  nucSubRate[2][1] = bbb*nucFreq[1];
  nucSubRate[2][3] = bbb*nucFreq[3];

  nucSubRate[3][0] = bbb*nucFreq[0];
  nucSubRate[3][1] = bbb*nucFreq[1] + aaa*nucFreq[1]/freqy;
  nucSubRate[3][2] = bbb*nucFreq[2];

  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)
  {
    nucSubRate[i][i] = 0.0;
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)
      if (j != i)
        nucSubRate[i][i] -= nucSubRate[i][j];
  }
} /* init_nucSubRates */


void fill_codon64_to(long i, long j, long k, double omeg)
{                                    /* for each third position */

  long m;
  long toIndex = i*16 + j*4 + k;
  long fromIndex;
  double subrat;      /* substitutionRate */

  if (transx[i][j][k] == stop)
    /* skip initial stop codons, leaving matrix entries 0.0 */
  {
    return;
  }


  for (m = 0; m < NUM_NUCLEOTIDE_STATES; m++)
  {
    if(transx[m][j][k] != stop && m != i)
      /* going from m to i in first  position */
    {
      subrat = (transx[m][j][k] == transx[i][j][k] ? 1.0 : omeg);
      fromIndex = m*16 + j*4 + k;
      codon64matrix[fromIndex][toIndex] = subrat * nucSubRate[m][i];
    }

    if(transx[i][m][k] != stop && m != j)
      /* going from m to j in middle position */
    {
      subrat = (transx[i][m][k] == transx[i][j][k] ? 1.0 : omeg);
      fromIndex = i*16 + m*4 + k;
      codon64matrix[fromIndex][toIndex] = subrat * nucSubRate[m][j];
    }

    if(transx[i][j][m] != stop && m != k)
      /* going from m to k in last   position */
    {
      subrat = (transx[i][j][m] == transx[i][j][k] ? 1.0 : omeg);
      fromIndex = i*16 + j*4 + m;
      codon64matrix[fromIndex][toIndex] = subrat * nucSubRate[m][k];
    }
  }
} /* fill_codon64_to */


void init_codon64matrix(double omeg)
{
  /* Make a 64 by 64 matrix of transition probabilities such that       */
  /* based on the probability of nucleotide mutation rate and omega,    */
  /* the nonsynonymous substitution rate factor.                        */

  long row, column;
  long i, j, k;                              /* counters                   */

  /* Most of the matrix elements are zero.        */
  /* Start by setting each element to zero.       */
  codon64matrix = (double **) Malloc(NUM_ALL_CODONS * sizeof(double*));
  for (row = 0; row < NUM_ALL_CODONS; row++)
  {
    codon64matrix[row] = (double *) Malloc(NUM_ALL_CODONS * sizeof(double));
    for (column = 0; column < NUM_ALL_CODONS; column++)
      codon64matrix[row][column] = 0.0;
  }

  /* Fill in the non-zero elements of the matrix. */
  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)       /* first position  */
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)     /* second position */
      for (k = 0; k < NUM_NUCLEOTIDE_STATES; k++)   /* third position  */
        fill_codon64_to(i, j, k, omeg);

  /* rates on the diagonal are these are inverse (negative) of the sum   */
  /* of all items in each row, because the row is the "from" index */
  for (row = 0; row < NUM_ALL_CODONS; row++)
  {
    double rowSum = 0.0;
    for (column = 0; column < NUM_ALL_CODONS; column++)
    {
      rowSum += codon64matrix[row][column];
    }
    codon64matrix[row][row] = (-1.0) * rowSum;
  }
} /* init_codon64matrix */


void matrix64to61(double** full_matrix, double*** no_stop_matrix)
{
  int i, j, k, row, column;
  boolean transcut[NUM_ALL_CODONS];
  int stopcounteri, stopcounterj;

  *no_stop_matrix = (double **) Malloc(NUM_SENSE_CODONS * sizeof(double*));
  for (row = 0; row < NUM_SENSE_CODONS; row++)
  {
    (*no_stop_matrix)[row] = (double *) Malloc(NUM_SENSE_CODONS * sizeof(double));
    for (column = 0; column < NUM_SENSE_CODONS; column++)
    {
      (*no_stop_matrix)[row][column] = 0.0;
    }
  }

  for (i = 0; i < NUM_ALL_CODONS; i++)       /* set arrays to 1, non-stop */
    transcut[i] = 1;

  /* identify stop codons for elimination from matrix */
  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)
  {                                            /* for each first position */
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)
    {                                       /* for each second position */
      for (k = 0; k < NUM_NUCLEOTIDE_STATES; k++)
      {                                    /* for each third position */
        if (transx[i][j][k] == stop)
          transcut[(i*16)+(j*4)+(k)] = 0;
      }
    }
  }

  /* create matrix61 from non-stop regions of matrix64 */
  for (i = 0, stopcounteri = 0; i < NUM_ALL_CODONS; i++)
  {
    if (transcut[i] == 0)
    {
      stopcounteri++;
    }
    else
    {
      for (j = 0, stopcounterj = 0; j < NUM_ALL_CODONS; j++)
      {
        if (transcut[j] == 0)
        {
          stopcounterj++;
        }
        else
        {
          assert(i-stopcounteri < NUM_SENSE_CODONS);
          assert(j-stopcounterj < NUM_SENSE_CODONS);
          (*no_stop_matrix)[i-stopcounteri][j-stopcounterj] =
            full_matrix[i][j];
        }
      }
      assert(stopcounterj == NUM_ALL_CODONS - NUM_SENSE_CODONS);
    }
  }
  assert(stopcounteri == NUM_ALL_CODONS - NUM_SENSE_CODONS);

} /* matrix64to61 */


void line64to61(double * full_line, double * no_stop_line)
{
  int i, j, k;
  boolean transcut[NUM_ALL_CODONS];
  int stopcounter;

  for (i = 0; i < NUM_ALL_CODONS; i++)                     /* set arrays to 1, non-stop */
    transcut[i] = true;

  /* identify stop codons for elimination from matrix */
  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)
  {                                            /* for each first position */
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)
    {                                       /* for each second position */
      for (k = 0; k < NUM_NUCLEOTIDE_STATES; k++)
      {                                    /* for each third position */
        if (transx[i][j][k] == stop)
          transcut[(i*16)+(j*4)+(k)] = false;
      }
    }
  }

  /* create matrix61 from non-stop regions of matrix64 */
  for (i = 0, stopcounter = 0; i < NUM_ALL_CODONS; i++)
  {
    if (transcut[i] == false)
    {
      stopcounter++;
    }
    else
    {
      no_stop_line[i-stopcounter] = full_line[i];
    }
  }
  assert(full_line[NUM_ALL_CODONS-1] == no_stop_line[NUM_SENSE_CODONS-1]);
  assert(stopcounter == NUM_ALL_CODONS - NUM_SENSE_CODONS);
} /* line64to61 */


void initomega(double *omega)
{ /* input omega, nonsynonymous/synonymous ratio */
  long loopcount;

  loopcount = 0;
  do {
    printf("Nonsynonymous/synonymous substitution ratio (omega)?\n");
    fflush(stdout);
    if(scanf("%lf%*[^\n]", omega)) {}   // Read number and scan to EOL.
    (void)getchar();
    countup(&loopcount, 10);
  } while (*omega < 0.0);
}  /* initomega */


void codml_empiricalfreqs(double *freqa, double *freqc, double *freqg, double *freqt, steptr weight, pointarray treenode)
{
  /* Get empirical base frequencies from the data
     this is a straightforward EM method.  If done from the tips,
     with the same code we can handle both amino acid and DNA cases */
  long i, j, k, m;
  double sum, suma, sumc, sumg, sumt, tempsum, w;

  *freqa = 0.25;
  *freqc = 0.25;
  *freqg = 0.25;
  *freqt = 0.25;
  for (k = 1; k <= 8; k++)
  {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;
    for (i = 0; i < spp; i++)
    {
      for (j = 0; j < endsite; j++)
      {
        w = weight[j];
        sum = 0.0;
        for (k = 0; k < NUM_SENSE_CODONS; k++) /* denominator */
        {
          tempsum = 1.0;
          for (m = 0; m < 4; m++)
            tempsum *= pow(nucFreq[m], basesincodon[m][k]);
          sum +=  tempsum * ((codon_node*)treenode[i])->codonx[j][0][k];
        }
        for (k = 0; k < NUM_SENSE_CODONS; k++) /* numerators */
        {
          tempsum = 1.0;
          for (m = 0; m < 4; m++)
            tempsum *= pow(nucFreq[m], basesincodon[m][k]);
          tempsum *= w * ((codon_node*)treenode[i])->codonx[j][0][k] / sum;
          suma += w * tempsum * basesincodon[0][k];
          sumc += w * tempsum * basesincodon[1][k];
          sumg += w * tempsum * basesincodon[2][k];
          sumt += w * tempsum * basesincodon[3][k];
        }
      }
    }
    sum = suma + sumc + sumg + sumt;
    *freqa = suma / sum;
    *freqc = sumc / sum;
    *freqg = sumg / sum;
    *freqt = sumt / sum;
  }
  if (*freqa <= 0.0)
    *freqa = 0.000001;
  if (*freqc <= 0.0)
    *freqc = 0.000001;
  if (*freqg <= 0.0)
    *freqg = 0.000001;
  if (*freqt <= 0.0)
    *freqt = 0.000001;
}  /* empiricalfreqs */


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;
  boolean didchangercat;
  char* string;

  putchar('\n');
  f84 = true;
  hky = false;
  k2p = false;
  jc = false;
  rctgry = false;
  didchangercat = false;
  rcategs = 1;
  auto_ = false;
  gama = false;
  freqsfrom = true;
  global = false;
  inputDataIsNucleotideData = true;
  hypstate = false;
  improve = false;
  invar = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  lambda = 1.0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  ttratio = 2.0;
  usertree = false;
  reusertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  omega = 1;
  for (;;)
  {
    cleerhome();
    printf("Codon Maximum Likelihood");
    printf(" method, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    if ( reusertree ) string = "Yes, rearrange on user tree";
    else if ( usertree ) string = "No, use user trees";
    else string = "Yes";
    printf("  U  Search for best tree (yes, no, rearrange)?  %s\n", string);
    if (usertree && !reusertree)
    {
      printf("  L                Use lengths from user trees?  %s\n", (lngths ? "Yes" : "No"));
    }
    printf("  D                               Input data as  %s\n", (inputDataIsNucleotideData ? "nucleotides" : "amino acids"));
    printf("  K              Nucleotide substitution model:  ");
    if (f84)
      printf("F84\n");
    else
    {
      if (hky)
        printf("HKY\n");
      else
      {
        if (k2p)
          printf("Kimura 2-parameter\n");
        else
          printf("Jukes-Cantor\n");
      }
    }
    if (!jc)
      printf("  T              Transition/transversion ratio:%8.4f\n",
             ttratio);
    if (!rctgry)
      printf("  Z     Nonsynonymous/synonymous ratio (omega): %7.4f\n",
             omega);

    if (!jc && !k2p)
      printf("  F             Use empirical base frequencies?  %s\n",
             (freqsfrom ? "Yes" : "No"));
    printf("  Q                                Codon table:  %s\n",
           (whichcode == universal) ? "Universal"                  :
           (whichcode == ciliate)   ? "Ciliate"                    :
           (whichcode == mito)      ? "Universal mitochondrial"    :
           (whichcode == vertmito)  ? "Vertebrate mitochondrial"   :
           (whichcode == flymito)   ? "Fly mitochondrial"          :
           (whichcode == yeastmito) ? "Yeast mitochondrial"        : "");
    printf("  N                         One value of omega?");
    if (!rctgry)
      printf("  Yes\n");
    else
    {
      printf("  multiple omegas\n");
      printf("  A        Omegas at adjacent sites correlated?");
      if (!auto_)
        printf("  No, they are independent\n");
      else
        printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    printf("  W                            Codons weighted?  %s\n", (weights ? "Yes" : "No"));
    if (!usertree || reusertree)
    {
      printf("  S              Speedier but rougher analysis?  %s\n", (improve ? "No, not rough" : "Yes"));
      printf("  G                      Global rearrangements?  %s\n", (global ? "Yes" : "No"));
    }
    if (!usertree)
    {
      printf("  J         Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                              Outgroup root?  %s%3ld\n", (outgropt ? "Yes, at sequence number" : "No, use species"), outgrno);
    printf("  M                 Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets, (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I                Input sequences interleaved?  %s\n", (interleaved ? "Yes" : "No, sequential"));
    printf("  0         Terminal type (IBM PC, ANSI, none)?  %s\n", (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1          Print out the data at start of run  %s\n", (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n", (progress ? "Yes" : "No"));
    printf("  3                              Print out tree  %s\n", (treeprint ? "Yes" : "No"));
    printf("  4             Write out trees onto tree file?  %s\n", (trout ? "Yes" : "No"));
    printf("  5         Reconstruct hypothetical sequences?  %s\n", (hypstate ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && !rctgry && (strchr("UKTZDFQNWSGJOMI012345", ch) != NULL))
        || (usertree && !rctgry && !reusertree && ((strchr("ULKTZDFQNWSGJOMI012345", ch) != NULL)))
        || (usertree && !rctgry && reusertree && ((strchr("UKTZDFQNWSGOMI012345", ch) != NULL)))
        || ((!usertree) && rctgry && (strchr("UKTZDFQNAWSGJOMI012345", ch) != NULL))
        || (usertree && !reusertree && rctgry && ((strchr("ULKTZDFQNAWSGJOMI012345", ch) != NULL)))
        || (usertree &&  reusertree && rctgry && ((strchr("UKTZDFQNAWSGOMI012345", ch) != NULL))))
    {
      // all items (lowercase == not yet)
      // !usertree                U TZDfqcWSGJOMI012345
      //  usertree && !reusertree ULTZDFQCWSGJOMI012345
      //  usertree &&  reusertree U TZDfqcWSG OMI012345

      switch (ch)
      {
        case 'U':
          if ( !usertree && !reusertree )
          {
            usertree = true;
            reusertree = false;
          }
          else if (!reusertree && usertree)
          {
            reusertree = true;
          }
          else
          {
            usertree = false;
            reusertree = false;
          }
          break;

        case 'K':
          if (f84)                      /* cycles among F84, HKY, K2P and JC */
          {
            f84 = false;
            hky = true;
          }
          else
          {
            if (hky)
            {
              hky = false;
              k2p = true;
            }
            else
            {
              if (k2p)
              {
                k2p = false;
                jc = true;
              }
              else
              {
                jc = false;
                f84 = true;
              }
            }
          }
          break;

        case 'T':
          initratio(&ttratio);
          break;

        case 'N':
          rctgry = !rctgry;
          break;

        case 'A':
          auto_ = !auto_;
          if (auto_)
            initlambda(&lambda);
          break;

        case 'Z':
          initomega(&omega);
          break;

        case 'D':
          inputDataIsNucleotideData = !inputDataIsNucleotideData;
          break;

        case 'F':
          freqsfrom = !freqsfrom;
          break;

        case 'Q':
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
          } while (ch != 'U' && ch != 'M' && ch != 'V' && ch != 'F' && ch != 'Y');

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

        case 'W':
          weights = !weights;
          break;

        case 'S':
          improve = !improve;
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
              justweights(&datasets);
            else
              initdatasets(&datasets);
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

        case '5':
          hypstate = !hypstate;
          break;
      }
    }
    else
      // BUG.969
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }

  if (!rctgry)
    auto_ = false;
  if (rctgry)
  {
    printf("\nOmega values in HMM\n");
    initcatn(&rcategs);
    printf("\n");
    if (probcat)
    {
      free(probcat);
      free(omegavalue);
    }
    probcat = (double *) Malloc(rcategs * sizeof(double));
    omegavalue   = (double *) Malloc(rcategs * sizeof(double));
    didchangercat = true;
    initcategs(rcategs, omegavalue);
    initprobcat(rcategs, &probsum, probcat);
  }
  if (!didchangercat)
  {
    omegavalue      = (double *) Malloc(rcategs * sizeof(double));
    probcat    = (double *) Malloc(rcategs * sizeof(double));
    omegavalue[0]   = omega;
    probcat[0] = 1.0;
  }

  if (f84)
    fprintf(outfile, "\nF84 base substitution model\n");
  else
  {
    if (hky)
      fprintf(outfile, "\nHKY base substitution model\n");
    else
    {
      if (k2p)
        fprintf(outfile, "\nKimura 2-parameter base substitution model\n");
      else
        fprintf(outfile, "\nJukes-Cantor base substitution model\n");
    }
  }

  if (!jc)
  {
    if (!justwts || firstset)
      fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n", ttratio);
  }
  fprintf(outfile, "\n");
}  /* getoptions */


void reallocsites(long dataSites, long workingSites)
{
  long i;
  for (i = 0; i < spp; i++)
  {
    free(inputSequences[i]);
    inputSequences[i] = (Char *) Malloc(dataSites * sizeof(Char));
  }

  free(category);
  free(weight);
  free(alias);
  free(ally);
  free(location);
  free(aliasweight);

  category    = (long *) Malloc(workingSites * sizeof(long));
  weight      = (long *) Malloc(workingSites * sizeof(long));
  alias       = (long *) Malloc(workingSites * sizeof(long));
  ally        = (long *) Malloc(workingSites * sizeof(long));
  location    = (long *) Malloc(workingSites * sizeof(long));
  aliasweight = (long *) Malloc(workingSites * sizeof(long));
} /* reallocsites */


void allocrest(long dataSites, long workingSites)
{
  long i;

  inputSequences = (Char **) Malloc(spp * sizeof(Char *));
  for (i = 0; i < spp; i++)
    inputSequences[i] = (Char *) Malloc(dataSites * sizeof(Char));

  nayme       = (naym *) Malloc(spp * sizeof(naym));
  enterorder  = (long *) Malloc(spp * sizeof(long));
  category    = (long *) Malloc(workingSites * sizeof(long));
  weight      = (long *) Malloc(workingSites * sizeof(long));
  alias       = (long *) Malloc(workingSites * sizeof(long));
  ally        = (long *) Malloc(workingSites * sizeof(long));
  location    = (long *) Malloc(workingSites * sizeof(long));
  aliasweight = (long *) Malloc(workingSites * sizeof(long));
}  /* allocrest */


void code(void)
{
  /* make up table of the code 0 = u, 1 = c, 2 = a, 3 = g */
  /* the universal code values are initialized when  transx  is declared */
  if (whichcode == mito)
  {
    transx[0][3][2] = trp;
    numsensecodons--;
  }
  if (whichcode == vertmito)
  {
    transx[0][3][2] = trp;
    transx[2][3][2] = stop;
    transx[2][3][3] = stop;
    transx[2][0][2] = met;
    numsensecodons--;
  }
  if (whichcode == flymito)
  {
    transx[0][3][2] = trp;
    transx[2][0][2] = met;
    transx[2][3][2] = ser2;
    numsensecodons++;
  }
  if (whichcode == yeastmito)
  {
    transx[0][3][2] = trp;
    transx[1][0][2] = thr;
    transx[2][0][2] = met;
    numsensecodons++;
  }
} /* code */


void doinit(void)
{ /* initializes variables */
  inputnumbers(&spp, &dataSites, &nonodes2, 1);
  numsensecodons = NUM_SENSE_CODONS;

  fprintf(outfile, "\nCodon model Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n", VERSION);

  if (!javarun)
  {
    getoptions();
  }

  if (f84)
    fprintf(outfile, "\nF84 base substitution model\n");
  else
  {
    if (hky)
      fprintf(outfile, "\nHKY base substitution model\n");
    else
    {
      if (k2p)
        fprintf(outfile, "\nKimura 2-parameter base substitution model\n");
      else
        fprintf(outfile, "\nJukes-Cantor base substitution model\n");
    }
  }

  if (!jc)
  {
    if (!justwts || firstset)
      fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n", ttratio);
  }
  fprintf(outfile, "\n");

  if(inputDataIsNucleotideData)
  {
    if (dataSites % NUM_NUCS_IN_CODON != 0)
    {
      printf("ERROR:  Number of sites (%ld) is not a multiple of (%ld); cannot\n", dataSites, (long)NUM_NUCS_IN_CODON);
      printf("        interpret these sites as codons.\n");
      exxit(-1);
    }
    workingSites = dataSites / NUM_NUCS_IN_CODON;
  }
  else
    workingSites = dataSites;

  // BUG.969
#if 0
  if (!usertree)
    nonodes2--;
#endif

  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, dataSites);

  // BUG.969
#if 0
  alloctree(&curtree.nodep, nonodes2, usertree);
#endif

  code();
  allocrest(dataSites, workingSites);

  // BUG.969
#if 0
  if (usertree)
    return;
  alloctree(&bestree.nodep, nonodes2, 0);
  alloctree(&priortree.nodep, nonodes2, 0);
  if (njumble <= 1)
    return;
  alloctree(&bestree2.nodep, nonodes2, 0);
#endif

  codon_inittable();
}  /* doinit */


void inputoptions(void)
{
  long i;

  if (!firstset && !justwts)
  {
    samenumsp(&dataSites, ith);
    reallocsites(dataSites, workingSites);
  }

  if (firstset || !justwts)
  {
    for (i = 0; i < workingSites; i++)
      category[i] = 1;
  }
  for (i = 0; i < workingSites; i++)
    weight[i] = 1;

  if (justwts || weights)
    inputweights(workingSites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < workingSites; i++)
    weightsum += weight[i];
  if (weights && printdata)
    printweights(outfile, 0, workingSites, weight, "Sites");
  if (!rctgry)
    fprintf(outfile, "omega = %8.4f\n", omega);
}  /* inputoptions */


void free_pmatrix(long sib)
{
  long j, l;

  for (j = 0; j < rcategs; j++)
  {
    for (l = 0; l < numsensecodons; l++)
      free(pmatrices[sib][j][l]);
    free(pmatrices[sib][j]);
  }
  free(pmatrices[sib]);
  pmatrices[sib] = NULL;
} /* free_pmatrix */


void alloc_pmatrix(long sib)
{
  /* Allocate memory for a new pmatrix.  Called iff num_sibs>max_num_sibs */
  long j, l;
  double ***temp_matrix;

  temp_matrix = (double ***) Malloc (rcategs * sizeof(double **));
  for (j = 0; j < rcategs; j++)
  {
    temp_matrix[j] = (double **) Malloc(numsensecodons * sizeof (double *));
    for (l = 0; l < numsensecodons; l++)
      temp_matrix[j][l] = (double *) Malloc(numsensecodons * sizeof(double));
  }
  pmatrices[sib] = temp_matrix;
  max_num_sibs++;
}  /* alloc_pmatrix */


void codon_freetable(void)
{
  long i, j, l;
  for (j = 0; j < rcategs; j++)
  {
    for (l = 0; l < numsensecodons; l++)
      free(ddpmatrix[j][l]);
    free(ddpmatrix[j]);
  }
  free(ddpmatrix);

  for (j = 0; j < rcategs; j++)
  {
    for (l = 0; l < numsensecodons; l++)
      free(dpmatrix[j][l]);
    free(dpmatrix[j]);
  }
  free(dpmatrix);

  free(tbl);

  for ( i = 0 ; i < 4 ; i++)
  {
    free(basesincodon[i]);
  }
  free(basesincodon);

  for ( i = 0 ; i < max_num_sibs ; i++ )
    free_pmatrix(i);
  free(pmatrices);
} /* codon_freetable */


void codon_inittable(void)
{
  /* Define a lookup table. Precompute values and print them out in tables */
  /* Allocate memory for the pmatrices, dpmatices and ddpmatrices          */
  long i, j, l, m;

  /* Allocate memory for pmatrices, the array of pointers to pmatrices     */

  pmatrices = (double ****) Malloc (spp * sizeof(double ***));

  /* Allocate memory for the first 2 pmatrices, the matrix of conversion   */
  /* probabilities, but only once per run (aka not on the second jumble.   */

  alloc_pmatrix(0);
  alloc_pmatrix(1);

  /*  Allocate memory for one dpmatrix for each omega value.
      It will be the first derivative matrix        */

  dpmatrix = (double ***) Malloc(rcategs * sizeof(double **));
  for (j = 0; j < rcategs; j++)
  {
    dpmatrix[j] = (double **) Malloc(numsensecodons * sizeof(double *));
    for (l = 0; l < numsensecodons; l++)
      dpmatrix[j][l] = (double *) Malloc(numsensecodons * sizeof(double));
  }

  /*  Allocate memory for one ddpmatrix, the second derivative matrix      */
  ddpmatrix = (double ***) Malloc(rcategs * sizeof(double **));
  for (j = 0; j < rcategs; j++)
  {
    ddpmatrix[j] = (double **) Malloc(numsensecodons * sizeof(double *));
    for (l = 0; l < numsensecodons; l++)
      ddpmatrix[j][l] = (double *) Malloc(numsensecodons * sizeof(double));
  }

  /* Allocate memory and assign values to tbl, the matrix of possible rates*/

  tbl = (double *) Malloc(rcategs * sizeof(double));
  for (j = 0; j < rcategs; j++)
    tbl[j] = omegavalue[j];

  if(jumb > 1)
    return;

  if (rcategs > 1)
  {
    fprintf(outfile, "\nState in HMM    Omega for state   Probability\n\n");
    for (i = 0; i < rcategs; i++)
      if (probcat[i] < 0.0001)
        fprintf(outfile, "%9ld%16.3f%20.6f\n", i+1, omegavalue[i], probcat[i]);
      else if (probcat[i] < 0.001)
        fprintf(outfile, "%9ld%16.3f%19.5f\n", i+1, omegavalue[i], probcat[i]);
      else if (probcat[i] < 0.01)
        fprintf(outfile, "%9ld%16.3f%18.4f\n", i+1, omegavalue[i], probcat[i]);
      else
        fprintf(outfile, "%9ld%16.3f%17.3f\n", i+1, omegavalue[i], probcat[i]);
    putc('\n', outfile);
    if (auto_)
      fprintf(outfile,
              "Expected length of a patch of sites having the same value of omega = %8.3f\n",
              1/lambda);
    putc('\n', outfile);
  }

  /* allocate table used to empirically estimate base frequencies */
  basesincodon = (double **) Malloc(4 * sizeof(double *));
  for (i = 0; i < 4; i++)
    basesincodon[i] = (double *) Malloc(numsensecodons * sizeof(double));
  for (i = 0; i < 4; i++)   /* and fill it in */
    for (j = 0; j < numsensecodons; j++)
      basesincodon[i][j] = 0.0;
  m = 0;    /* counter for which of the codons we are at */
  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)
      for (k = 0; k < NUM_NUCLEOTIDE_STATES; k++)
      {
        if (transx[i][j][k] != stop)
        {
          basesincodon[i][m] += 1.0;
          basesincodon[j][m] += 1.0;
          basesincodon[k][m] += 1.0;
          m++;
        }
      }
}  /* codon_inittable */


void input_codondata(long chars)
{
  /* input the names and sequences for each species */
  /* used by codml */
  long i, j, k, l, basesread, basesnew;
  Char charstate;
  /*
    Char charstate, codonStrings[NUM_NUCS_IN_CODON + 1] = {'\0','\0','\0','\0'};
  */
  boolean allread, done;

  if (inputDataIsNucleotideData && chars % NUM_NUCS_IN_CODON != 0)
  {
    printf("ERROR:  Number of sites (%ld) is not a multiple of (%ld); cannot ", chars, (long)NUM_NUCS_IN_CODON);
    printf("interpret as codons.\n");
    exxit(-1);
  }

  if (printdata)
    headings(chars, "Sequences", "---------");
  basesread = 0;
  basesnew = 0;
  allread = false;
  while (!(allread))
  {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp)
    {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i - 1); /* reads "nmlength" chars from "infile" into nayme[i-1] */
      j = (interleaved) ? basesread : 0;
      done = false; /* "done" is only meaningful for sequential format */
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
          if (inputDataIsNucleotideData)
          {
            if ((strchr("ABCDGHKMNRSTUVWXY?O-", charstate)) == NULL)
            {
              printf("ERROR:  Bad nucleotide symbol: %c at position %ld of species %ld.\n", charstate, j+1, i);
              if (charstate == '.')
              {
                printf("        Periods (.) may not be used as gap characters.\n");
                printf("        The correct gap character is (-).\n");
              }
              exxit(-1);
            }
          }
          else
          {
            if ((strchr("ABCDEFGHIKLMNPQRSTVWXYZ*?-", charstate)) == NULL)
            {
              printf("ERROR:  Bad amino acid symbol: %c at position %ld of species %ld.\n", charstate, j+1, i);
              if (charstate == '.')
              {
                printf("        Periods (.) may not be used as gap characters.\n");
                printf("        The correct gap character is (-).\n");
              }
              exxit(-1);
            }
          }
          inputSequences[i-1][j] = charstate;
          j++;
        }
        if (interleaved)
          continue; /* same effect as "break", because done==true */
        if (j < chars) /* then eol or eof was detected */
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j; /* running tally of num. chars appearing in species 1 */

      scan_eoln(infile);

      if ((interleaved && j != basesnew) ||
          (!interleaved && j != chars))
      {
        printf("ERROR:  SEQUENCES OUT OF ALIGNMENT AT OR BEFORE POSITION %ld.\n", j);
        exxit(-1);
      }
      i++;
    }

    if (interleaved)
    {
      basesread = basesnew;
      allread = (basesread == chars);
    }
    else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
  if (!printdata)
    return;
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
        if (j > 1 && inputSequences[j - 1][k - 1] == inputSequences[0][k - 1])
          charstate = '.';
        else
          charstate = inputSequences[j - 1][k - 1];
        putc(charstate, outfile);
        if (k % 10 == 0 && k % 60 != 0)
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* input_codondata */


void makeweights(boolean inputDataIsNucleotideData)
{
  /* make up bookkeeping vectors to avoid duplicate computations.
     for explanation of variables  alias, aliasweight, and ally see seq.c
     routines  sitesort2, sitecombine2, and sitescrunch2
     location[i-1] ends up as the position in the lexicographic order
     (lexicographic by site rate category and site pattern) of site i
     if site i is one whose category and site pattern represents the
     others that are tied with it.  Otherwise location[i-1] is 0
     endsite is the number (less 1) of such representative sites  */
  long i;

  for (i = 1; i <= workingSites; i++)
  {
    alias[i - 1] = i;
    ally[i - 1] = i;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }

  if(!inputDataIsNucleotideData)
  {
    sitesort2   (workingSites, aliasweight);
    sitecombine2(workingSites, aliasweight);
    sitescrunch2(workingSites, 1, 2, aliasweight);
  }

  endsite = 0;
  for (i = 1; i <= workingSites; i++)
  {
    if (ally[i - 1] == i)
      endsite++;
  }

  for (i = 1; i <= endsite; i++)
  {
    location[alias[i - 1] - 1] = i;
  }

  /* these are being allocated here because now we know how big "endsite" is*/
  term = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    term[i] = (double *) Malloc(rcategs * sizeof(double));
  slopeterm = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    slopeterm[i] = (double *) Malloc(rcategs * sizeof(double));
  curveterm = (double **) Malloc(endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
    curveterm[i] = (double *) Malloc(rcategs * sizeof(double));
  mp = (vall *) Malloc(workingSites * sizeof(vall));
  contribution = (contribarr *) Malloc(endsite * sizeof(contribarr));
}  /* makeweights */


void aa_codon_makevalues(long rcategs, pointarray treenode, long endsite, long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips         */
  /*      when data is input as amino acids        */
  /* a version of prot_makevalues found in proml.c */
  /* used by codml                                 */

  long i, j, k, l, m, n, o, p;
  double temp_line[64];

  for (k = 0; k < endsite; k++)
  {
    j = alias[k];
    for (i = 0; i < spp; i++)
    {
      for (l = 0; l < rcategs; l++)
      {
        memset(temp_line, 0, sizeof(double)*NUM_ALL_CODONS);
        switch (y[i][j - 1])
        {
          case 'A':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == ala)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'B':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == asp || transx[m][n][o] == asn)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'C':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == cys)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'D':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == asp)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'E':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == glu)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'F':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == phe)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'G':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == gly)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'H':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == his)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'I':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == ileu)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'K':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == lys)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'L':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == leu)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'M':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == met)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'N':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == asn)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'P':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == pro)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'Q':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == gln)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'R':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == arg)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'S':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == ser)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'T':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == thr)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'V':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == val)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'W':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == trp)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'X':
          case '?':
          case '-':
            for (p = 0; p < 64; p++)
              temp_line[p] = 1.0;
            break;
          case 'Y':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == tyr)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          case 'Z':
            for (m = 0; m < 4; m++)
              for (n = 0; n < 4; n++)
                for (o = 0; o < 4; o++)
                  if (transx[m][n][o] == gln || transx[m][n][o] == glu)
                    temp_line[(m*16) + (n*4) + o] = 1.0;
            break;
          default:
            printf("ERROR:  Unexpected amino acid symbol \'%c\' created internally; program error.", y[i][j - 1]);
            exxit(-1);
        }

        memset(((codon_node*)(treenode[i]))->codonx[k][l], 0, sizeof(double)*numsensecodons);
        line64to61(temp_line, ((codon_node*)(treenode[i]))->codonx[k][l]);
        treenode[i]->initialized = 1; // BUG.969 -- why?
      }
    }
  }
}  /* aa_codon_makevalues */


void nuc_codon_makevalues(long rcategs, pointarray treenode, long endsite, long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips         */
  /*      when data is input as nucleotides        */
  /* a version of prot_makevalues found in proml.c */
  /* used by codml                                 */

  long i, k, l, m, n, o, p, q;
  boolean temp_index[NUM_NUCS_IN_CODON][NUM_NUCLEOTIDE_STATES];
  double temp_line[64];
  (void)alias;                          // RSGnote: Never used.

  for (k = 0; k < endsite; k++ )
  {
    for (i = 0; i < spp; i++)
    {
      for (l = 0; l < rcategs; l++)
      {
        memset(temp_line, 0, sizeof(double)*NUM_ALL_CODONS);

        for (p = 0; p < NUM_NUCS_IN_CODON; p++)
          for (q = 0; q < NUM_NUCLEOTIDE_STATES; q++)
            temp_index[p][q] = 0;

        for (p = 0; p < NUM_NUCS_IN_CODON; p++)
          switch (y[i][3*k+p])
          {
            case 'A':
              temp_index[p][0] = 1;
              break;
            case 'C':
              temp_index[p][1] = 1;
              break;
            case 'G':
              temp_index[p][2] = 1;
              break;
            case 'T':
              temp_index[p][3] = 1;
              break;
            case 'R':                         /* purine        */
              temp_index[p][0] = 1;
              temp_index[p][2] = 1;
              break;
            case 'Y':                         /* pyrimidine    */
              temp_index[p][1] = 1;
              temp_index[p][3] = 1;
              break;
            case 'K':                         /* keto          */
              temp_index[p][2] = 1;
              temp_index[p][3] = 1;
              break;
            case 'M':                         /* amino         */
              temp_index[p][0] = 1;
              temp_index[p][1] = 1;
              break;
            case 'S':                         /* strong        */
              temp_index[p][1] = 1;
              temp_index[p][2] = 1;
              break;
            case 'W':                         /* weak          */
              temp_index[p][0] = 1;
              temp_index[p][3] = 1;
              break;
            case 'B':                         /* not A         */
              temp_index[p][1] = 1;
              temp_index[p][2] = 1;
              temp_index[p][3] = 1;
              break;
            case 'D':                         /* not C         */
              temp_index[p][0] = 1;
              temp_index[p][2] = 1;
              temp_index[p][3] = 1;
              break;
            case 'H':                         /* not G         */
              temp_index[p][0] = 1;
              temp_index[p][1] = 1;
              temp_index[p][3] = 1;
              break;
            case 'V':                         /* not T (or U)  */
              temp_index[p][1] = 1;
              temp_index[p][2] = 1;
              temp_index[p][3] = 1;
              break;
            case 'N':                         /* anything      */
            case 'X':
            case '?':
            case 'O':
            case '-':
              temp_index[p][1] = 1;
              temp_index[p][2] = 1;
              temp_index[p][3] = 1;
              break;
            default:
              printf("ERROR:  Unexpected nucleotide symbol \'%c\' created internally; program error.", y[i][3*k + p]);
              exxit(-1);
          }
        for (m = 0; m < NUM_NUCLEOTIDE_STATES; m++)
          for (n = 0; n < NUM_NUCLEOTIDE_STATES; n++)
            for (o = 0; o < NUM_NUCLEOTIDE_STATES; o++)
              if (temp_index[0][m] == 1 && temp_index[1][n] == 1 && temp_index[2][o] == 1)
                temp_line[(m*16) + (n*4) + o] = 1.0;
        memset(((codon_node*)(treenode[i]))->codonx[k][l], 0, sizeof(double)*numsensecodons);
        line64to61(temp_line, ((codon_node*)(treenode[i]))->codonx[k][l]);
        treenode[i]->initialized = 1; // BUG.969 -- why?
      }
    }
  }
} /* nuc_codon_makevalues */


void getcodonbasefreqs (void)
{
  /* read or compute the base frequencies */

  if (hky || f84)
  {
    if (!freqsfrom)
      initfreqs(&nucFreq[0], &nucFreq[1], &nucFreq[2], &nucFreq[3]);
    else
      codml_empiricalfreqs(&nucFreq[0], &nucFreq[1], &nucFreq[2], &nucFreq[3], weight, curtree->nodep);
  }
  else
  {
    nucFreq[0] = 0.25;
    nucFreq[1] = 0.25;
    nucFreq[2] = 0.25;
    nucFreq[3] = 0.25;
  };
  putc('\n', outfile);
  fprintf(outfile, "Base Frequencies:\n\n");
  fprintf(outfile, "   A    %10.5f\n", nucFreq[0]);
  fprintf(outfile, "   C    %10.5f\n", nucFreq[1]);
  fprintf(outfile, "   G    %10.5f\n", nucFreq[2]);
  fprintf(outfile, "  T(U)  %10.5f\n", nucFreq[3]);
  putc('\n', outfile);
} /* getcodonbasefreqs */


void init_SymMatrix(void)
{
  /* Turn InstMatrix into a symmetrical matrix for eigenvalue calcs */
  long i, j;

  SymMatrix = (double **) Malloc(numsensecodons * sizeof(double*));

  for (i = 0; i < numsensecodons; i++)
  {
    SymMatrix[i] = (double*) Malloc(numsensecodons * sizeof(double));
  }
  for (i = 0; i < numsensecodons; i++)
  {
    for (j = 0; j < numsensecodons; j++)
    {
      if (i == j)
      {
        SymMatrix[i][i] = InstMatrix[i][i];
      }
      else
      {
        SymMatrix[i][j] = InstMatrix[i][j] * sqrt(codonfreq[i]) / sqrt(codonfreq[j]);
      }
    }
  }
  assert(is_symmetric(SymMatrix, numsensecodons));
}  /* init_SymMatrix */


boolean codonfreqs_sum_to_one(void)
{
  double fdiff, tmpFreq = 0.0;
  long i;
  for (i = 0; i < numsensecodons; i++)
  {
    tmpFreq += codonfreq[i];
  }
  fdiff = fabs(1.0 - tmpFreq);
  assert(fdiff < LIKE_ACCURACY);
  return(fdiff < LIKE_ACCURACY);
}


void makecodonfreqs(void)
{
  /* calculate codon frequencies based on nucleotide frequencies */
  long i, j, k;
  double stopfreq = 0;
  double codonfreq64[64];

  for (i = 0; i < NUM_NUCLEOTIDE_STATES; i++)
    for (j = 0; j < NUM_NUCLEOTIDE_STATES; j++)
      for (k = 0; k < NUM_NUCLEOTIDE_STATES; k++)
      {
        long index = i*16 + j*4 + k;
        codonfreq64[index] = nucFreq[i] * nucFreq[j] * nucFreq[k];
        if (transx[i][j][k] == stop)
          stopfreq += codonfreq64[index];
      }
  for (i = 0; i < NUM_ALL_CODONS; i++)
    codonfreq64[i] *= 1.0 / (1.0 - stopfreq);
  line64to61(codonfreq64, codonfreq);
  assert(codonfreqs_sum_to_one());
} /* makecodonfreqs */


boolean is_f_times_m_symmetric(double ** mymatrix)
{
  long i, j;
  for(i=0;i<numsensecodons;i++)
  {
    for(j=i+1;j<numsensecodons;j++)
    {
      double ij = codonfreq[i] * mymatrix[i][j];
      double ji = codonfreq[j] * mymatrix[j][i];
      assert(fabs(ij-ji) < LIKE_ACCURACY);
      if(fabs(ij-ji) >= LIKE_ACCURACY) return false;
    }
  }
  return true;
}


void givens(eigCalcT **prob, long i, long j, long n, eigCalcT ctheta,
            eigCalcT stheta, boolean left)
{ /* Givens transform at i, j for 1..n with angle theta */
  long k;
  eigCalcT d;

  for (k = 0; k < n; k++)
  {
    if (left)
    {
      d = ctheta * prob[i - 1][k] + stheta * prob[j - 1][k];
      prob[j - 1][k] = ctheta * prob[j - 1][k] - stheta * prob[i - 1][k];
      prob[i - 1][k] = d;
    }
    else
    {
      d = ctheta * prob[k][i - 1] + stheta * prob[k][j - 1];
      prob[k][j - 1] = ctheta * prob[k][j - 1] - stheta * prob[k][i - 1];
      prob[k][i - 1] = d;
    }
  }
}  /* givens */


void coeffs(eigCalcT x, eigCalcT y, eigCalcT *c, eigCalcT *s )
{ /* compute cosine and sine of theta */
  eigCalcT root;

  root = sqrt(x * x + y * y);
  if (root < QR_accuracy)
  {
    *c = 1.0;
    *s = 0.0;
  }
  else
  {
    *c = x / root;
    *s = y / root;
  }
}  /* coeffs */


void tridiag(eigCalcT **prob, eigCalcT **eigvecs, long n)
{
  /* Givens tridiagonalization */
  long i, j;
  eigCalcT s, c;

  for (i = 2; i < n; i++)
  {
    for (j = i + 1; j <= n; j++)
    {
      coeffs(prob[i - 2][i - 1], prob[i - 2][j - 1], &c, &s);
      givens(prob, i, j, n, c, s, true);
      givens(prob, i, j, n, c, s, false);
      givens(eigvecs, i, j, n, c, s, true);
    }
  }
}  /* tridiag */


void shiftqr(eigCalcT **prob, eigCalcT **eigvecs, long n)
{ /* QR eigenvalue-finder */
  long i, j;
  eigCalcT approx, s, c, d, TEMP, TEMP1;

  for (i = n; i >= 2; i--)
  {
    do {
      TEMP = prob[i - 2][i - 2] - prob[i - 1][i - 1];
      TEMP1 = prob[i - 1][i - 2];
      d = sqrt(TEMP * TEMP + TEMP1 * TEMP1);
      approx = prob[i - 2][i - 2] + prob[i - 1][i - 1];
      if (prob[i - 1][i - 1] < prob[i - 2][i - 2])
        approx = (approx - d) / 2.0;
      else
        approx = (approx + d) / 2.0;
      for (j = 0; j < i; j++)
        prob[j][j] -= approx;
      for (j = 1; j < i; j++)
      {
        coeffs(prob[j - 1][j - 1], prob[j][j - 1], &c, &s );
        givens(prob, j, j + 1, i, c, s, true);
        givens(prob, j, j + 1, i, c, s, false);
        givens(eigvecs, j, j + 1, n, c, s, true);
      }
      for (j = 0; j < i; j++)
        prob[j][j] += approx;
    } while (fabs(prob[i - 1][i - 2]) > QR_accuracy);
  }
}  /* shiftqr */


boolean is_symmetric(double ** matrix, long n)
{
  long i, j;

  for (i = 0; i < n; i++)
  {
    for (j = i+1; j < n; j++)
    {
      double diff = matrix[i][j] - matrix[j][i];
      assert (fabs(diff) < 1e-6);
      if (fabs(diff) >= 1e-6) return false;
    }
  }
  return true;
}


void qreigen(double ** inMatrix, long n, long categ)
{ /* QR eigenvector/eigenvalue method for symmetric matrix       */
  /* taken from protdist.c and only changed for last few lines.  */
  /* eigenvectors left in inMatrix, eigenvalues in eigmat[categ] */

  long i, j;
  eigCalcT **prob, **eigvecs;

  assert(is_symmetric(inMatrix, n));

  prob = (eigCalcT **)Malloc(n * sizeof(eigCalcT *));
  eigvecs = Malloc(n * sizeof(eigCalcT *));
  for (i = 0; i < n; i++)
  {
    prob[i] = (eigCalcT *)Malloc(n * sizeof(eigCalcT));
    eigvecs[i] = (eigCalcT *)Malloc(n * sizeof(eigCalcT));
    for (j = 0; j < n; j++)
    {
      prob[i][j] = (eigCalcT)(inMatrix[i][j]);
      if (i == j)
        eigvecs[i][j] = 1.0;
      else
        eigvecs[i][j] = 0.0;
    }
  }

  tridiag(prob, eigvecs, n);
  shiftqr(prob, eigvecs, n);

  /* added to conform to codml */
  /* BUG.969 -- change so we use eig instead of eigmat ?? */
  eigmat[categ] = (double *)Malloc(n * sizeof(double));
  for (i = 0; i < n; i++)
    eigmat[categ][i] = (double)prob[i][i];

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      /* converting back to double from long double */
      inMatrix[i][j] = (double)eigvecs[i][j];
      /* BUG.969 -- original had conversion using codonfreqs
       * here, but unless we do that in all our qreigen
       * users, I'd rather not do it here -- it's not part
       * of the eigenvector calculation proper, just a
       * shortcut which is bound to be confusing
       */
    }
  }

  for (i = 0; i < n; i++)
  {
    free(prob[i]);
    free(eigvecs[i]);
  }
  free(prob);
  free(eigvecs);
}  /* qreigen */


void freelrsaves(void)
{
  long i, j;
  for ( i = 0 ; i < NLRSAVES ; i++ )
  {
    for (j = 0; j < endsite; j++)
      free(((codon_node*)(lrsaves[i]))->codonx[j]);
    free(((codon_node*)(lrsaves[i]))->codonx);
    free(((ml_node*)(lrsaves[i]))->underflows);
    free((codon_node*)(lrsaves[i]));
  }
  free(lrsaves);
} /* freelrsaves */


void codml_alloclrsaves(void)
{
  long i, j;
  lrsaves = Malloc(NLRSAVES * sizeof(codon_node*));
  //oldendsite = endsite;
  for ( i = 0 ; i < NLRSAVES ; i++ )
  {
    lrsaves[i] = Malloc(sizeof(codon_node));
    ((codon_node*)(lrsaves[i]))->codonx = Malloc(endsite * sizeof(ratelike));
    ((ml_node*)(lrsaves[i]))->underflows = Malloc(endsite * sizeof (double));
    for (j = 0; j < endsite; j++)
      ((codon_node*)(lrsaves[i]))->codonx[j]  = (cratelike)Malloc(rcategs * sizeof(csitelike));
  }
} /* codml_alloclrsaves */


void maketables(void)
{
  /* sets up all the necessary frequencies, probabilities, eigenvectors, */
  /* and eigenvalues                                                     */
  long i, j;

  /* making copy so we can compare */
  eigmat = (double **)Malloc(rcategs * sizeof(double));
  probmat = (double ***) Malloc(rcategs * sizeof(double **));
  for (i = 0; i < rcategs; i++)
  {
    init_nucSubRates(ttratio);
    makecodonfreqs();

    init_codon64matrix(omegavalue[i]);
    matrix64to61(codon64matrix, &InstMatrix);
    assert(is_f_times_m_symmetric(InstMatrix));

    init_SymMatrix();

    probmat[i] = (double **) Malloc(numsensecodons * sizeof(double *));
    for(j=0; j<numsensecodons; j++)
    {
      probmat[i][j] = (double *) Malloc(numsensecodons * sizeof(double));
      for(k=0; k<numsensecodons; k++)
        probmat[i][j][k] = SymMatrix[j][k];
    }
    qreigen(probmat[i], numsensecodons, i);
  }
} /* maketables */


void getinput(void)
{
  /* reads the input data */
  // BUG.969 -- took this away as other tests are inside inputoptions -- if (!justwts || firstset)
  inputoptions();
  if (!justwts || firstset)
    input_codondata(dataSites);
  if ( !firstset ) freelrsaves();
  makeweights(inputDataIsNucleotideData);
  codon_tree_setup(nonodes2, spp);
  codml_alloclrsaves();

  // BUG.969 -- probably OK already -- setuptree2(&curtree);
#if 0
  if (!usertree)
  {
    setuptree2(&bestree);
    setuptree2(&priortree);
    if (njumble > 1)
      setuptree2(&bestree2);
  }

  // BUG.969 -- these are probably being handled by the "new" phylip c++-like scheme.
  codon_allocx(nonodes2, rcategs, curtree.nodep, usertree);
  if (!usertree)
  {
    codon_allocx(nonodes2, rcategs, bestree.nodep, 0);
    codon_allocx(nonodes2, rcategs, priortree.nodep, 0);
    if (njumble > 1)
      codon_allocx(nonodes2, rcategs, bestree2.nodep, 0);
  }
#endif

  if (inputDataIsNucleotideData)
  {
    nuc_codon_makevalues(rcategs, curtree->nodep, endsite, spp, inputSequences, alias);
  }
  else
  {
    aa_codon_makevalues(rcategs, curtree->nodep, endsite, spp, inputSequences, alias);
  }
  getcodonbasefreqs();
  maketables();
}  /* getinput */


void codon_tree_calc_nuview(tree * t, node * p, cphenotype cph, double *  underflows)
{
  long i, j, l, m, num_sibs, sib_index; // RSGnote: "k" set but never used.
  node *sib_ptr, *sib_back_ptr;
  csitelike codon_xx, x2;
  double lw, prod7;
  double **pmat;
  double maxx;
  double correction;
  (void)t;                              // RSGnote: Parameter "t" never used.

  /* Figure out how many siblings the current node has  */
  /* and be sure that pmatrices is large enough         */
  num_sibs = count_sibs(p);
  for (i = 0; i < num_sibs; i++)
    if (pmatrices[i] == NULL)
      alloc_pmatrix(i);

  /* Make pmatrices for all possible combinations of category, rcateg, and sib. */
  sib_ptr = p;                                /* return to p */
  for (sib_index=0; sib_index < num_sibs; sib_index++)
  {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    lw = sib_back_ptr->v;

    for (j = 0; j < rcategs; j++)
      make_pmatrix(pmatrices[sib_index][j], NULL, NULL, 0, lw, tbl[j], eigmat[j], probmat[j]);
  }

  for (i = 0; i < endsite; i++)
  {
    maxx = 0;
    correction = 0;
    double origCorrection = ((ml_node*)p)->underflows[i]; /* just used for debugging */

    // RSGnote: "k" set but never used.
    // k = category[alias[i]-1] - 1;

    for (j = 0; j < rcategs; j++)
    {

      /* initialize to 1 all values of codon_xx */
      for (m = 0; m < numsensecodons; m++)
        codon_xx[m] = 1.0;

      sib_ptr = p;                        /* return to p */
      /* loop through all sibs and calculate likelihoods for all possible*/
      /* amino acid combinations                                         */
      for (sib_index=0; sib_index < num_sibs; sib_index++)
      {
        sib_ptr      = sib_ptr->next;
        sib_back_ptr = sib_ptr->back;

        if ( j == 0)
          correction += ((ml_node*)sib_back_ptr)->underflows[i];

        memcpy(x2, ((codon_node*)sib_back_ptr)->codonx[i][j], sizeof(csitelike));
        assert(sib_back_ptr->initialized || sib_back_ptr->tip);

        pmat = pmatrices[sib_index][j];
        for (m = 0; m < numsensecodons; m++)
        {

          prod7 = 0;
          for (l = 0; l < numsensecodons; l++)
            prod7 += x2[l] * pmat[m][l];

          codon_xx[m] *= prod7;
          if ( codon_xx[m] > maxx && sib_index == (num_sibs - 1))
            maxx = codon_xx[m];

        }
      }

      /* And the final point of this whole function: */
      memcpy(cph[i][j], codon_xx, sizeof(csitelike));
    }

    underflows[i] = 0;
    if ( maxx < MIN_DOUBLE && maxx > 0.0)
    {
      underflows[i] += log(maxx);
    }
    underflows[i] += correction;
    assert(maxx < MIN_DOUBLE || origCorrection == correction);
  }
} /* codon_tree_calc_nuview */


void codon_tree_nuview(tree* t, node *p)
{
  long siteIndex, rcategIndex;
  cphenotype temp_codonx;

  temp_codonx = (cphenotype)Malloc(endsite * sizeof(cratelike));
  for ( siteIndex = 0 ; siteIndex < endsite ; siteIndex++ )
  {
    temp_codonx[siteIndex] = (cratelike)Malloc(rcategs * sizeof(csitelike));
  }

  double * temp_underflows;
  temp_underflows = (double*)Malloc(endsite * sizeof(double));

  generic_tree_nuview(t, p);

  /* passing into separate calculation function. this is so
   * we can use it to calculate nuviews without updating them,
   * should we need to do that for debugging
   */
  codon_tree_calc_nuview(t, p, temp_codonx, temp_underflows);

  for ( siteIndex = 0 ; siteIndex < endsite ; siteIndex++ )
    for ( rcategIndex = 0 ; rcategIndex < rcategs ; rcategIndex++ )
    {
      memcpy(((codon_node*)p)->codonx[siteIndex][rcategIndex], temp_codonx[siteIndex][rcategIndex], sizeof(csitelike));
    }

  memcpy(((ml_node*)p)->underflows, temp_underflows, endsite * sizeof(double));

  p->initialized = true;
  inittrav(t, p->back);

  for ( siteIndex = 0 ; siteIndex < endsite ; siteIndex++ )
    free(temp_codonx[siteIndex]);
  free(temp_codonx);
  free(temp_underflows);

}  /* codon_tree_nuview */


void codon_slopecurv(node *p, double y, double *like, double *slope, double *curve)
{
  /* compute log likelihood, slope and curvature at node p */
  long i, j, l, m, lai;
  double sum, sumc, sumterm, lterm, sumcs, sumcc, sum2, slope2, curve2;
  double frexm = 0;                     /* frexm = codonfreq*x1[m] */
                                        /* frexml = frexm*x2[l]    */
  double prod4m, prod5m, prod6m;        /* elements of prod4-5 for */
                                        /* each m                  */
  double **pmat, **dpmat, **ddpmat;     /* local pointers to global*/
                                        /* matrices                */
  double prod4, prod5, prod6;
  contribarr thelike, nulike, nuslope, nucurve,
    theslope, thecurve, clai, cslai, cclai;
  node *q;
  csitelike x1, x2;

  q = p->back;
  sum = 0.0;
  for (j = 0; j < rcategs; j++)
    make_pmatrix(pmatrices[0][j], dpmatrix[j], ddpmatrix[j],
                 2, y, tbl[j], eigmat[j], probmat[j]);
  for (i = 0; i < endsite; i++)
  {
    for (j = 0; j < rcategs; j++)
    {
      memcpy(x1, ((codon_node*)p)->codonx[i][j], sizeof(csitelike));
      memcpy(x2, ((codon_node*)q)->codonx[i][j], sizeof(csitelike));
      assert(p->initialized || p->tip);
      assert(q->initialized || q->tip);

      pmat = pmatrices[0][j];
      dpmat = dpmatrix[j];
      ddpmat = ddpmatrix[j];
      prod4 = 0.0;
      prod5 = 0.0;
      prod6 = 0.0;
      for (m = 0; m < numsensecodons; m++)
      {
        prod4m = 0.0;
        prod5m = 0.0;
        prod6m = 0.0;
        frexm = x1[m] * codonfreq[m];
        for (l = 0; l < numsensecodons; l++)
        {
          /* BUG.969 --backwards ? */
          prod4m += x2[l] * pmat[m][l];
          prod5m += x2[l] * dpmat[m][l];
          prod6m += x2[l] * ddpmat[m][l];
        }
        prod4 += frexm * prod4m;
        prod5 += frexm * prod5m;
        prod6 += frexm * prod6m;
      }

      term[i][j] = prod4;
      slopeterm[i][j] = prod5;
      curveterm[i][j] = prod6;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * term[i][j];
    for (j = 0; j < rcategs; j++)
    {
      term[i][j] = term[i][j] / sumterm;
      slopeterm[i][j] = slopeterm[i][j] / sumterm;
      curveterm[i][j] = curveterm[i][j] / sumterm;
    }
    /* munging this here to avoid taking a bad log */
    /* BUG.969 -- don't want to use this new sumterm above */
    /* debug   was:     if (sumterm <= LIKE_ACCURACY)  */
    if (sumterm <= 0.0)
      sumterm = LIKE_ACCURACY;
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    sum += (aliasweight[i] * lterm);
  }
  for (i = 0; i < rcategs; i++)
  {
    thelike[i] = 1.0;
    theslope[i] = 0.0;
    thecurve[i] = 0.0;
  }
  for (i = 0; i < workingSites; i++)
  {
    sumc = 0.0;
    sumcs = 0.0;
    sumcc = 0.0;
    for (k = 0; k < rcategs; k++)
    {
      sumc += probcat[k] * thelike[k];
      sumcs += probcat[k] * theslope[k];
      sumcc += probcat[k] * thecurve[k];
    }
    sumc *= lambda;
    sumcs *= lambda;
    sumcc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0))
    {
      lai = location[ally[i] - 1];
      memcpy(clai, term[lai - 1], rcategs * sizeof(double));
      memcpy(cslai, slopeterm[lai - 1], rcategs * sizeof(double));
      memcpy(cclai, curveterm[lai - 1], rcategs * sizeof(double));
      if (weight[i] > 1)
      {
        for (j = 0; j < rcategs; j++)
        {
          if (clai[j] > 0.0)
            clai[j] = exp(weight[i]*log(clai[j]));
          else clai[j] = 0.0;
          if (cslai[j] > 0.0)
            cslai[j] = exp(weight[i]*log(cslai[j]));
          else cslai[j] = 0.0;
          if (cclai[j] > 0.0)
            cclai[j] = exp(weight[i]*log(cclai[j]));
          else cclai[j] = 0.0;
        }
      }
      for (j = 0; j < rcategs; j++)
      {
        nulike[j]  = ((1.0 - lambda) * thelike[j]  + sumc) *  clai[j];
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs) * clai[j]
          + ((1.0 - lambda) * thelike[j]  + sumc) *  cslai[j];
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc) * clai[j]
          + 2.0 * ((1.0 - lambda) * theslope[j] + sumcs) * cslai[j]
          + ((1.0 - lambda) * thelike[j]  + sumc) *  cclai[j];
      }
    }
    else
    {
      for (j = 0; j < rcategs; j++)
      {
        nulike[j]  = ((1.0 - lambda) * thelike[j]  + sumc);
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs);
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc);
      }
    }
    memcpy(thelike, nulike, rcategs * sizeof(double));
    memcpy(theslope, nuslope, rcategs * sizeof(double));
    memcpy(thecurve, nucurve, rcategs * sizeof(double));
  }
  sum2 = 0.0;
  slope2 = 0.0;
  curve2 = 0.0;
  for (i = 0; i < rcategs; i++)
  {
    sum2 += probcat[i] * thelike[i];
    slope2 += probcat[i] * theslope[i];
    curve2 += probcat[i] * thecurve[i];
  }
  if(sum2 <= 0.0) sum2 = LIKE_ACCURACY;
  sum += log(sum2);
  (*like) = sum;
  (*slope) = slope2 / sum2;

  /* Expressed in terms of *slope to prevent overflow */
  (*curve) = curve2 / sum2 - *slope * *slope;

  assert(*curve >= 0.0 || *curve < 0.0); /* to catch NAN */

} /* codon_slopecurv */


void codon_tree_makenewv(tree* t, node *p)
{
  /* Newton-Raphson algorithm improvement of a branch length */
  long it, ite;
  double y, yold=0, yorig, like, slope, curve, oldlike=0;
  boolean done, firsttime, better;
  node *q;

  q = p->back;
  y = p->v;
  assert(y >= 0.0);
  assert(y >  0.0);
  yorig = y;
  done = false;
  firsttime = true;
  it = 1;
  ite = 0;
  while ((it < iterations) && (ite < 20) && (!done))
  {
    codon_slopecurv(p, y, &like, &slope, &curve);
    better = false;
    if (firsttime)                      /* if no older value of y to compare with */
    {
      yold = y;
      oldlike = like;
      firsttime = false;
      better = true;
    }
    else
    {
      if (like > oldlike)               /* update the value of yold if it was better */
      {
        yold = y;
        oldlike = like;
        better = true;
        it++;
      }
    }
    if (better)
    {
      assert(curve != 0.0);
      y = y + slope/fabs(curve);        /* Newton-Raphson, forced uphill-wards */
      if (y < epsilon)
        y = epsilon;
    }
    else
    {
      if (fabs(y - yold) < epsilon)
        ite = 20;
      y = (y + (7 * yold)) / 8;         /* retract 87% of way back */
    }
    ite++;
    done = fabs(y-yold) < 0.1*epsilon;
  }
  smoothed = (fabs(yold-yorig) < epsilon) && (yorig > 1000.0*epsilon);
  p->v = yold;                          /* the last one that had better likelihood */
  q->v = yold;

  /* BUG.969 -- is there a better place to make sure it's always done? */
  inittrav(t, p);
  assert(q == p->back);
  inittrav(t, p->back);

  ((tree*)t)->score = oldlike;

}  /* codon_tree_makenewv */


void codml_restore_traverses(tree *t, node *p, node * q)
{
  /* BUG.969 -- should this be ml's version (which doesn't
     exist, but could some day? */
  generic_tree_restore_traverses(t, p, q);

  /* BUG.969 -- moved updating of {p, q}->back-> to the generic version */
}


boolean codml_node_good(tree *t, node *n)
{
  long siteIndex, rateCatIndex, codonIndex;
  cphenotype temp_codonx;

  boolean basicGood = generic_node_good(t, n);
  assert(basicGood);                    // RSGdebug

  if (! basicGood) return false;
  if (n->tip) return true;
  if (n->back == NULL) return true;
  return true;                          // BUG.969

  // RSGnote: If we return unconditionally above, what is all the rest of this for?
  temp_codonx = (cphenotype)Malloc(endsite * sizeof(cratelike));
  for ( siteIndex = 0 ; siteIndex < endsite ; siteIndex++ )
    temp_codonx[siteIndex] = (cratelike)Malloc(rcategs * sizeof(csitelike));

  double * temp_underflows;
  temp_underflows = (double*)Malloc(endsite * sizeof(double));

  codon_tree_calc_nuview(t, n, temp_codonx, temp_underflows);

  for(siteIndex=0; siteIndex < endsite; siteIndex++)
  {
    boolean underMatch =  (temp_underflows[siteIndex] == ((ml_node*)n)->underflows[siteIndex]);

    // if (! underMatch ) printf("BUG.969 %e %e\n", temp_underflows[siteIndex], ((ml_node*)n)->underflows[siteIndex]);
    assert(underMatch);
    if ( !underMatch ) return false;
  }

  for(siteIndex=0; siteIndex < endsite; siteIndex++)
  {
    for(rateCatIndex=0; rateCatIndex < rcategs; rateCatIndex++)
    {
      for(codonIndex=0; codonIndex < numsensecodons; codonIndex++)
      {
        double orig = ((codon_node*)n)->codonx[siteIndex][rateCatIndex][codonIndex];
        double newV = temp_codonx[siteIndex][rateCatIndex][codonIndex];
        double diff = fabs(newV - orig);
        assert(diff < 1e-5);
        if (diff > 1e-5) return false;
      }

    }
  }

  return true;
}


void make_pmatrix(double **mymatrix, double **dmat, double **ddmat, long derivative, double lz, double rat, double *eigm, double **probm)
{
  /* Computes the R matrix such that matrix[m][l] is the joint probability */
  /* of m and l.                                                           */
  /* Computes a P matrix such that matrix[m][l] is the conditional         */
  /* probability of m given l.  This is accomplished by dividing all terms */
  /* in the R matrix by freqaa[m], the frequency of l.                     */

  long k, l, m;                 /* (l) original character state */
                                /* (m) final    character state */
                                /* (k) lambda counter           */
  double p0, p1, p2, q;
  double elambdat[NUM_ALL_CODONS];
  double delambdat[NUM_ALL_CODONS];
  double ddelambdat[NUM_ALL_CODONS];
  /* exponential term for matrix  */
  /* and both derivative matrices */
  assert(lz > 0.0);
  for (k = 0; k < numsensecodons; k++)
  {
    elambdat[k] = exp(lz * rat * eigm[k]);

    /* BUG.969-- can this be 0.0? any upper limit?*/

    if(derivative != 0)
    {
      delambdat[k] = (elambdat[k] * rat * eigm[k]);
      ddelambdat[k] = (delambdat[k] * rat * eigm[k]);
    }
  }
  for (m = 0; m < numsensecodons; m++)
  {
    for (l = 0; l < numsensecodons; l++)
    {
      p0 = 0.0;
      p1 = 0.0;
      p2 = 0.0;
      for (k = 0; k < numsensecodons; k++)
      {
        q = probm[k][m] * probm[k][l];
        p0 += (q * elambdat[k]);
        if(derivative !=0)
        {
          p1 += (q * delambdat[k]);
          p2 += (q * ddelambdat[k]);
        }
      }

      mymatrix[m][l] = p0 * sqrt(codonfreq[l]) / sqrt(codonfreq[m]);
      if(derivative != 0)
      {
        dmat[m][l] = p1 * sqrt(codonfreq[l]) / sqrt(codonfreq[m]);
        ddmat[m][l] = p2 * sqrt(codonfreq[l]) / sqrt(codonfreq[m]);
      }

    }
  }

  assert(is_f_times_m_symmetric(mymatrix));
}  /* make_pmatrix */


double codon_tree_evaluate(tree*t, node *p, boolean saveit)
{
  contribarr tterm;
  double sum, sum2, sumc, y, prod4, prodl, frexm, sumterm, lterm;
  double **pmat;
  long i, j, k, l, m, lai;
  node *q;
  csitelike x1, x2;

  generic_tree_evaluate(t, p, saveit);

  sum = 0.0;
  q = p->back;
  y = p->v;
  for (j = 0; j < rcategs; j++)
    make_pmatrix(pmatrices[0][j], NULL, NULL, 0, y, tbl[j], eigmat[j], probmat[j]);

  for (i = 0; i < endsite; i++)
  {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++)
    {
      memcpy(x1, ((codon_node*)p)->codonx[i][j], sizeof(csitelike));
      memcpy(x2, ((codon_node*)q)->codonx[i][j], sizeof(csitelike));
      assert(p->initialized || p->tip);
      assert(q->initialized || q->tip);

      prod4 = 0.0;
      pmat = pmatrices[0][j];
      for (m = 0; m < numsensecodons; m++)
      {
        prodl = 0.0;
        for (l = 0; l < numsensecodons; l++)
          prodl += x2[l] * pmat[m][l];
        frexm = x1[m] * codonfreq[m];
        prod4 += prodl * frexm;
      }
      tterm[j] = prod4;

    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
    {
      sumterm += probcat[j] * tterm[j];
    }
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, rcategs * sizeof(double));
    /* debug    if (sumterm < LIKE_ACCURACY) sumterm = LIKE_ACCURACY; */
    if (sumterm <= 0.0) sumterm = LIKE_ACCURACY;
    lterm = log(sumterm) + ((ml_node*)p)->underflows[i] +
      ((ml_node*)q)->underflows[i];
    if (saveit && !auto_ && usertree && (which <= shimotrees) && !reusertree)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }

  for (j = 0; j < rcategs; j++)
    like[j] = 1.0;
  for (i = 0; i < workingSites; i++)
  {
    sumc = 0.0;
    for (k = 0; k < rcategs; k++)
      sumc += probcat[k] * like[k];
    sumc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0))
    {
      lai = location[ally[i] - 1];
      memcpy(clai, contribution[lai - 1], rcategs * sizeof(double));
      for (j = 0; j < rcategs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
    }
    else
    {
      for (j = 0; j < rcategs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc);
    }
    memcpy(like, nulike, rcategs * sizeof(double));
  }
  sum2 = 0.0;
  for (i = 0; i < rcategs; i++)
    sum2 += probcat[i] * like[i];
  /* debug   if(sum2 <= LIKE_ACCURACY) sum2 = LIKE_ACCURACY;   */
  if(sum2 <= 0.0) sum2 = LIKE_ACCURACY;
  sum += log(sum2);
  curtree->score = sum;
  if (!saveit || auto_ || !usertree || reusertree)
    return sum;
  if(which <= shimotrees)
    l0gl[which - 1] = sum;
  if (which == 1)
  {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl)
  {
    maxwhich = which;
    maxlogl = sum;
  }

  return sum;
}  /* codon_tree_evaluate */


void codmlcopy(tree *a, tree *b, long nonodes, long rcategs)
{
  /* copy tree a to tree b */
  long i, j=0;
  node *p, *q;
  (void)rcategs;                        // RSGnote: Parameter never used.

  for (i = 0; i < spp; i++)
  {
    codon_node_copy(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back)
    {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++)
  {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++)            /* BUGBUG what's 3? what's this loop doing? ask Joe */
    {
      codon_node_copy(p, q);
      if (p->back)
      {
        if (p->back == a->nodep[p->back->index - 1])
          q->back = b->nodep[p->back->index - 1];
        else if (p->back == a->nodep[p->back->index - 1]->next)
          q->back = b->nodep[p->back->index - 1]->next;
        else
          q->back = b->nodep[p->back->index - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->score = a->score;
  b->root = a->root;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* codmlcopy */


void codml_coordinates(node *p, double lengthsum, long *tipy, double *tipmax)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

  if (p->tip)
  {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    xx = q->v;
    if (xx > 100.0)
      xx = 100.0;
    codml_coordinates(q->back, lengthsum + xx, tipy, tipmax);
    q = q->next;
  } while ((p == curtree->root || p != q) &&
           (p != curtree->root || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == curtree->root)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* codml_coordinates */


void codml_printree(void)
{
  /* prints out diagram of the tree2 */
  long tipy;
  double scale, tipmax;
  long i;

  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  codml_coordinates(curtree->root, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline2(i, scale, curtree);
  putc('\n', outfile);
}  /* codml_printree */


void sigma(node *p, double *sumlr, double *s1, double *s2)
{
  /* compute standard deviation */
  double tt, aa, like, slope, curv;

  codon_slopecurv(p, p->v, &like, &slope, &curv);
  tt = p->v;
  p->v = epsilon;
  p->back->v = epsilon;
  aa = curtree->evaluate(curtree, p, false);
  p->v = tt;
  p->back->v = tt;
  (*sumlr) = curtree->evaluate(curtree, p, false) - aa;
  if (curv < -epsilon)
  {
    (*s1) = p->v + (-slope - sqrt(slope * slope -  3.841 * curv)) / curv;
    (*s2) = p->v + (-slope + sqrt(slope * slope -  3.841 * curv)) / curv;
  }
  else
  {
    (*s1) = -1.0;
    (*s2) = -1.0;
  }
}  /* sigma */


void describe(node *p)
{
  /* print out information for one branch */
  long i, num_sibs;
  node *q, *sib_ptr;
  double sumlr, sigma1, sigma2;

  if (!p->tip && !p->initialized)
    curtree->nuview(curtree, p);
  if (!p->back->tip && !p->back->initialized)
    curtree->nuview(curtree, p->back);
  q = p->back;
  if (q->tip)
  {
    fprintf(outfile, " ");
    for (i = 0; i < nmlngth; i++)
      putc(nayme[q->index-1][i], outfile);
    fprintf(outfile, "    ");
  }
  else
    fprintf(outfile, "  %4ld          ", q->index - spp);
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.5f", q->v);
  if (!usertree || (usertree && !lngths) || p->iter || reusertree)
  {
    sigma(q, &sumlr, &sigma1, &sigma2);
    if (sigma1 <= sigma2)
      fprintf(outfile, "     (     zero,    infinity)");
    else
    {
      fprintf(outfile, "     (");
      if (sigma2 <= 0.0)
        fprintf(outfile, "     zero");
      else
        fprintf(outfile, "%9.5f", sigma2);
      fprintf(outfile, ",%12.5f", sigma1);
      putc(')', outfile);
    }
    if (sumlr > 1.9205)
      fprintf(outfile, " *");
    if (sumlr > 2.995)
      putc('*', outfile);
  }
  putc('\n', outfile);
  if (!p->tip)
  {
    num_sibs = count_sibs(p);
    sib_ptr  = p;
    for (i=0; i < num_sibs; i++)
    {
      sib_ptr = sib_ptr->next;
      describe(sib_ptr->back);
    }
  }
}  /* describe */


void codon_reconstr(node *p, long n)
{
  /* reconstruct and print out acid at site n+1 at node p */
  long i, j, k, first, num_sibs = 0;
  double f, sum, xx[NUM_ALL_CODONS];
  node *q = NULL;

  if (p->tip)
  {
    if(workingSites == dataSites)
    {
      putc(inputSequences[p->index-1][n], outfile);
      putc(' ', outfile);
      putc(' ', outfile);
      putc(' ', outfile);
    }
    else
    {
      putc(inputSequences[p->index-1][n*3], outfile);
      putc(inputSequences[p->index-1][n*3+1], outfile);
      putc(inputSequences[p->index-1][n*3+2], outfile);
      putc(' ', outfile);
    }
  }
  else
  {
    num_sibs = count_sibs(p);
    if ((ally[n] == 0) || (location[ally[n]-1] == 0))
      putc('.', outfile);
    else
    {
      assert(p->initialized || p->tip);
      j = location[ally[n]-1] - 1;
      sum = 0;
      for (i = 0; i < numsensecodons; i++)
      {
        f = ((codon_node*)p)->codonx[j][mx-1][i];
        if (!p->tip)
        {
          q = p;
          for (k = 0; k < num_sibs; k++)
          {
            q = q->next;
            f *= ((codon_node*)q)->codonx[j][mx-1][i];
          }
        }
        if (f > 0.0)
          f = exp(log(f)/(num_sibs-1.0));
        xx[i] = f * codonfreq[i];
        sum += xx[i];
      }
      for (i = 0; i <= 60; i++)
        xx[i] /= sum;
      first = 0;
      for (i = 0; i <= 60; i++)
        if (xx[i] > xx[first])
          first = i;
      if (xx[first] > 0.95)
        fprintf(outfile, "%s ", codonStrings[first]); // BUG.969 -- check
      else
      {
        char lowCodon[NUM_NUCS_IN_CODON+1];
        int i;
        for(i=0; i < NUM_NUCS_IN_CODON+1; i++)
        {
          lowCodon[i] = tolower(codonStrings[first][i]);
        }
        fprintf(outfile, "%s ", lowCodon); // BUG.969 -- check
      }
      if (rctgry && rcategs > 1)
        mx = mp[n][mx - 1];
      else
        mx = 1;
    }
  }
} /* codon_reconstr */


void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  node *q;
  long i;

  putc(' ', outfile);
  if (p->tip)
  {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  }
  else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "  ");
  mx = mx0;
  for (i = m; i <= n; i++)
  {
    if ((i % 5 == 0) && (i != m))
      putc(' ', outfile);
    codon_reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip)
  {
    for (q = p->next; q != p; q = q->next)
    {
      rectrav(q->back, m, n);
    }
  }
  mx1 = mx;
}  /* rectrav */


void summarize(void)
{
  /* print out branch length information and node numbers */
  long i, j, num_sibs;
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;
  node   *sib_ptr;
  long mm=0;                            // RSGnote: Variable used possibly uninitialized.

  if (!treeprint)
    return;
  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree->score);
  fprintf(outfile, "\n Between        And            Length");
  if (!(usertree && lngths && haslengths ) || reusertree)
    fprintf(outfile, "      Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, " -------        ---            ------");
  if (!(usertree && lngths && haslengths) || reusertree)
    fprintf(outfile, "      ------- ---------- ------");
  fprintf(outfile, "\n\n");
  for (i = spp; i < nonodes2; i++)
  {
    /* So this works with arbitrary multifurcations */
    if (curtree->nodep[i])
    {
      num_sibs = count_sibs (curtree->nodep[i]);
      sib_ptr  = curtree->nodep[i];
      for (j = 0; j < num_sibs; j++)
      {
        sib_ptr->initialized = false;
        sib_ptr              = sib_ptr->next;
      }
    }
  }

  describe(curtree->root->back);

  /* So this works with arbitrary multifurcations */
  num_sibs = count_sibs(curtree->root);
  sib_ptr  = curtree->root;
  for (i=0; i < num_sibs; i++)
  {
    sib_ptr = sib_ptr->next;
    describe(sib_ptr->back);
  }

  fprintf(outfile, "\n");
  if (!(usertree && lngths && haslengths)|| reusertree)
  {
    fprintf(outfile, "     *  = significantly positive, P < 0.05\n");
    fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n");
  }
  curtree->evaluate(curtree, curtree->root, false);
  if (rctgry && rcategs > 1)
  {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = workingSites - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
        mp[i][j] = j + 1;
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
          {
            if (lambda * probcat[k - 1] * like[k - 1] > nulike[j])
            {
              nulike[j] = lambda * probcat[k - 1] * like[k - 1];
              mp[i][j] = k;
            }
          }
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= rcategs; i++)
    {
      if (probcat[i - 1] * like[i - 1] > mode)
      {
        mx = i;
        mode = probcat[i - 1] * like[i - 1];
      }
    }
    mx0 = mx;
    fprintf(outfile, "Combination of omega values that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 1; i <= workingSites; i++)
    {
      fprintf(outfile, "%ld", mx);
      if (i % 10 == 0)
        putc(' ', outfile);
      if (i % 60 == 0 && i != workingSites)
      {
        putc('\n', outfile);
        for (j = 1; j <= nmlngth + 3; j++)
          putc(' ', outfile);
      }
      mx = mp[i - 1][mx - 1];
    }
    fprintf(outfile, "\n\n");
    marginal = (double **) Malloc(workingSites * sizeof(double *));
    for (i = 0; i < workingSites; i++)
      marginal[i] = (double *) Malloc(rcategs * sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = workingSites - 1; i >= 0; i--)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
            nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] /= sum;
        marginal[i][j] = nulike[j];
      }
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = 0; i < workingSites; i++)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
      {
        nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++)
        {
          if (k != j + 1)
            nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        marginal[i][j] *= like[j] * probcat[j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        sum += marginal[i][j];
      for (j = 0; j < rcategs; j++)
        marginal[i][j] /= sum;
    }
    fprintf(outfile, "Most probable omega value at each codon\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < workingSites; i++)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum)
        {
          sum = marginal[i][j];
          mm = j;
        }
      fprintf(outfile, "%ld", mm+1);
      if ((i+1) % 60 == 0)
      {
        if (i != 0)
        {
          putc('\n', outfile);
          for (j = 1; j <= nmlngth + 3; j++)
            putc(' ', outfile);
        }
      }
      else if ((i+1) % 10 == 0)
        putc(' ', outfile);
    }
    putc('\n', outfile);
    putc('\n', outfile);
    fprintf(outfile, "Most probable omega value at each codon if > 0.95");
    fprintf(outfile, " probability (\".\" otherwise)\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < workingSites; i++)
    {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum)
        {
          sum = marginal[i][j];
          mm = j;
        }
      if (sum >= 0.95)
        fprintf(outfile, "%ld", mm+1);
      else
        putc('.', outfile);
      if ((i+1) % 60 == 0)
      {
        if (i != 0)
        {
          putc('\n', outfile);
          for (j = 1; j <= nmlngth + 3; j++)
            putc(' ', outfile);
        }
      }
      else if ((i+1) % 10 == 0)
        putc(' ', outfile);
    }
    putc('\n', outfile);
    for (i = 0; i < workingSites; i++)
      free(marginal[i]);
    free(marginal);
  }
  putc('\n', outfile);
  if (hypstate)
  {
    fprintf(outfile, "Probable sequences at interior nodes:\n\n");
    fprintf(outfile, "  node       ");
    for (i = 0; (i < 13) && (i < ((workingSites + (workingSites-1)/10 - 39) / 2)); i++)
      putc(' ', outfile);
    fprintf(outfile, "Reconstructed sequence (caps if > 0.95)\n\n");
    if (!rctgry || (rcategs == 1))
      mx0 = 1;
    for (i = 0; i < workingSites; i += 15)
    {
      k = i + 14;
      if (k >= workingSites)
        k = workingSites - 1;
      rectrav(curtree->root, i, k);
      rectrav(curtree->root->back, i, k);
      putc('\n', outfile);
      mx0 = mx1;
    }
  }
}  /* summarize */


void dnaml_treeout(node *p)
{
  /* write out file with representation of final tree2 */
  /* Only works for bifurcations! */
  long i, n, w;
  Char c;
  double x;
  node *q;
  boolean inloop;

  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)
    {
      if (nayme[p->index-1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)
    {
      c = nayme[p->index-1][i];
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

    inloop = false;
    q = p->next;
    do  {
      if (inloop)
      {
        putc(',', outtree);
        col++;
        if (col > 45)
        {
          putc('\n', outtree);
          col = 0;
        }
      }
      inloop = true;
      dnaml_treeout(q->back);
      q = q->next;
    } while ((p == curtree->root || p != q) && (p != curtree->root || p->next != q));

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
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* dnaml_treeout */


void free_all_codonx (long nonodes, pointarray treenode)
{
  /* used in proml */
  long i, j, k;
  node *p;

  /* Zero to (but not including) spp are tips. */
  for (i = 0; i < spp; i++)
  {
    for (j = 0; j < endsite; j++)
      free(((codon_node*)treenode[i])->codonx[j]);
    free(((codon_node*)treenode[i])->codonx);
  }

  /* The rest are FORKRINGS (i.e. triads). */
  for (i = spp; i < nonodes; i++)
  {
    if (treenode[i] != NULL)
    {
      p = treenode[i];
      do {
        for (k = 0; k < endsite; k++)
          free(((codon_node*)p)->codonx[k]);
        free(((codon_node*)p)->codonx);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
}  /* free_all_codonx */


void codml_reroot(tree* t)              // RSGbugfix: Name change.
{
  node *q;
  double newl;
  node *r = t->root;
  long numsibs = count_sibs(r);

  if ( numsibs > 2)
  {
    q = r;
    while ( q->next != r )
      q = q->next;
    q->next = r->next;
    t->release_forknode(t, r);
    t->nodep[spp] = q;
  }
  else
  {
    assert(r->back == NULL);            // RSGnote: This assumes the FORKRING being manipulated has the ROOT FORKNODE pointing to NULL.

    newl = r->next->oldlen + r->next->next->oldlen;
    r->next->back->oldlen = newl;
    r->next->next->back->oldlen = newl;

    newl = r->next->v + r->next->next->v;
    r->next->back->v = newl;
    r->next->next->back->v = newl;

    r->next->back->back = r->next->next->back;
    r->next->next->back->back = r->next->back;

    t->release_fork(t, r);
  }

  t->root = t->nodep[0]->back;          // Reset ROOT; moved from line just after call to CODML_REROOT.

} /* codml_reroot */


void initcodonmlnode(tree *treep, node **p, long len, long nodei, long *ntips, long *parens, initops whichinit, pointarray nodep, Char *str, Char *ch, FILE *intree)
{
  boolean minusread;
  double valyew, divisor;
  (void)len;                            // RSGnote: Parameter never used.
  (void)ntips;                          // RSGnote: Parameter never used.

  switch (whichinit)
  {
    case bottom:
      *p = treep->get_forknode(treep, nodei);
      codon_node_allocx((ml_node*)(*p), endsite, rcategs);
      nodep[(*p)->index - 1] = (*p);
      break;
    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      codon_node_allocx((ml_node*)*p, endsite, rcategs);
      (*p)->index = nodei;
      break;
    case tip:
      match_names_to_data(str, nodep, p, spp);
      break;
    case iter:
      (*p)->initialized = false;
      (*p)->v = initialv;
      (*p)->iter = true;
      if ((*p)->back != NULL)
      {
        (*p)->back->iter = true;
        (*p)->back->v = initialv;
        (*p)->back->initialized = false;
      }
      break;
    case length:
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      (*p)->v = valyew / divisor;
      (*p)->iter = false;
      if ((*p)->back != NULL)
      {
        (*p)->back->v = (*p)->v;
        (*p)->back->iter = false;
      }
      break;
    case hsnolength:
      haslengths = false;
      break;
    default:      /* cases hslength, treewt, unittrwt */
      break;      /* should never occur               */
  }
} /* initcodonmlnode */


void maketree(void)
{
  long i;
  boolean dummy_first, goteof;
  long nextnode;
  node * qwhere = NULL;
  double bestyet;

  if (usertree)
  {
    openfile(&intree, INTREE, "input tree file", "r", progname, intreename);
    numtrees = countsemic(intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *) Malloc(shimotrees * sizeof(double));
    l0gf = (double **) Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *) Malloc(endsite * sizeof(double));
    if (treeprint)
    {
      fprintf(outfile, "\nUser-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;

    /* This taken out of tree read, used to be [spp-1], but referring to [0] produces output identical to what the pre-modified dnaml produced. */
    for (which = 1 ; which <= numtrees; which++)
    {
      /* These initializations required each time through the loop since multiple trees require re-initialization */
      haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      preparetree(curtree);
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof, &dummy_first, &nextnode, &haslengths, initcodonmlnode, false, nonodes2);
      fixtree(curtree);
      codml_reroot(curtree);                              // RSGbugfix: Name change.

      if (goteof && (which <= numtrees))
      {
        /* if we hit the end of the file prematurely */
        printf ("\nERROR:  Trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit (-1);
      }

      // RSGbugfix: Reset of ROOT moved to inside CODML_REROOT.

      if ( outgropt )
        curtree->root = curtree->nodep[outgrno - 1]->back;
      ml_treevaluate(curtree, improve, reusertree, global, progress, priortree, bestree, ml_initialvtrav);
      // BUG.969 -- proml has reusertree stuff in here
      if (treeprint)
      {
        codml_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        dnaml_treeout(curtree->root);
      }
      if(which < numtrees)
      {
        codon_freex_notip(nextnode, curtree->nodep);
      }
      else
        nonodes2 = nextnode;
    }
    FClose(intree);
    putc('\n', outfile);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite-1, maxlogl,
               l0gl, l0gf, aliasweight, seed);
  }
  else                                  /* no user tree */
  {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    if (progress)
    {
      sprintf(progbuf, "\nAdding species:\n");
      print_progress(progbuf);
      writename(0, 3, enterorder);
      phyFillScreenColor();
    }
    nextsp = 3;
    polishing = false;
    destruct_tree(curtree);
    buildsimpletree(curtree, enterorder);
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    smoothit = improve;
    nextsp = 4;

    while (nextsp <= spp)
    {
      assert(curtree->tree_good_f(curtree));    /* debug:  tree_good functions unneeded */

      bestyet = UNDEFINED;
      if (smoothit)
      {
        curtree->copy(curtree, priortree);
      }

      curtree->addtraverse(curtree,
                            curtree->nodep[enterorder[nextsp - 1] - 1],
                            curtree->root, true, qwhere, &bestyet,
                            bestree, smoothit);
      if (smoothit)
      {
        bestree->copy(bestree, curtree);
      }
      else
      {
        smoothit = true;
        curtree->insert_(curtree, curtree->nodep[enterorder[nextsp - 1] - 1],
                          qwhere, false);
        smoothit = false;
        curtree->copy(curtree, bestree);
        bestyet = curtree->score;
      }
      if (progress)
      {
        writename(nextsp - 1, 1, enterorder);
        phyFillScreenColor();
      }
      assert(curtree->tree_good_f(curtree)); /* BUG.969 */

      if (global && nextsp == spp)
      {
        curtree->globrearrange(curtree, progress, smoothit);
      }
      else
      {
        curtree->locrearrange(curtree, curtree->nodep[enterorder[0] - 1], smoothit, priortree, bestree) ;
      }

      if(!smoothit)
      {
        curtree->smoothall(curtree, curtree->root);
        bestyet = curtree->score;
      }

      assert(curtree->tree_good_f(curtree)); /* BUG.969 */
      nextsp++;
    }

    if (global && progress)
    {
      putchar('\n');
      fflush(stdout);
      phyFillScreenColor();
    }
    curtree->copy(curtree, bestree);

    if (njumble > 1)
    {
      if (jumb == 1)
        bestree->copy(bestree, bestree2);
      else
        if (bestree2->score < bestree->score)
          bestree->copy(bestree, bestree2);
    }
    if (jumb == njumble)
    {
      if (njumble > 1)
        bestree2->copy(bestree2, curtree);
      curtree->root = curtree->nodep[outgrno - 1]->back;
      ml_treevaluate(curtree, improve, reusertree, global, progress, priortree, bestree, ml_initialvtrav);
      if (treeprint)
      {
        codml_printree();
        summarize();
      }
      if (trout)
      {
        col = 0;
        dnaml_treeout(curtree->root);
      }
    }
  }

  if (usertree)
  {
    free(l0gl);
    for (i=0; i < shimotrees; i++)
      free(l0gf[i]);
    free(l0gf);
  }
  // codon_freetable();
  if (jumb < njumble)
    return;
  free(contribution);
  free(mp);
  for (i=0; i < endsite; i++)
    free(term[i]);
  free(term);
  for (i=0; i < endsite; i++)
    free(slopeterm[i]);
  free(slopeterm);
  for (i=0; i < endsite; i++)
    free(curveterm[i]);
  free(curveterm);

#if 0 //  BUG.969
  free_all_codonx(nonodes2, curtree->nodep);
  if (!usertree)
  {
    free_all_codonx(nonodes2, bestree->nodep);
    free_all_codonx(nonodes2, priortree->nodep);
    if (njumble > 1)
      free_all_codonx(nonodes2, bestree2->nodep);
  }
#endif

  if (progress)
  {
    sprintf(progbuf, "\n\nOutput written to file \"%s\".\n\n", outfilename);
    print_progress(progbuf);
    if (trout)
    {
      sprintf(progbuf, "Tree also written onto file \"%s\".\n", outtreename);
      print_progress(progbuf);
    }
    sprintf(progbuf, "\n");
    print_progress(progbuf);
  }
}  /* maketree */


void clean_up(void)
{
  /* Free and/or close stuff */
  long i, j;

  codon_freetable();
  freelrsaves();

  free (omegavalue);
  free (probcat);
  /* Seems to require freeing every time... */
  for (i = 0; i < spp; i++)
    free (inputSequences[i]);
  free (inputSequences);
  free (nayme);
  free (enterorder);
  free (category);
  free (weight);
  free (alias);
  free (ally);
  free (location);
  free (aliasweight);
  for (j = 0; j < rcategs; j++)
  {
    for (i = 0; i < numsensecodons; i++)
      free (probmat[j][i]);
    free (probmat[j]);
  }
  for (j = 0; j < numsensecodons; j++)
  {
    free (InstMatrix[j]);
    free (SymMatrix[j]);
  }
  free (probmat);
  free (InstMatrix);
  free (SymMatrix);

  for (j=0; j < rcategs; j++)
  {
    free(eigmat[j]);
  }
  free (eigmat);

  for (i = 0; i < NUM_ALL_CODONS; i++)
  {
    free (codon64matrix[i]);
  }
  free (codon64matrix);

  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
}   /* clean_up */


void codmlrun(void)
{
  // debug printout // JRMdebug
#ifdef JAVADEBUG
  printf("\nf84: %i\n", f84);
  printf("hky: %i\n", hky);
  printf("k2p: %i\n", k2p);
  printf("jc: %i\n", jc);
  printf("rctgry: %i\n", rctgry);
  printf("rcategs: %li\n", rcategs);
  printf("auto_: %i\n", auto_);
  printf("freqsfrom: %i\n", freqsfrom);
  printf("global: %i\n", global);
  printf("inputDataIsNucleotideData: %i\n", inputDataIsNucleotideData);
  printf("hypstate: %i\n", hypstate);
  printf("improve: %i\n", improve);
  printf("invar: %i\n", invar);
  printf("jumble: %i\n", jumble);
  printf("njumble: %li\n", njumble);
  printf("lngths: %i\n", lngths);
  printf("lambda: %f\n", lambda);
  printf("outgrno: %li\n", outgrno);
  printf("outgropt: %i\n", outgropt);
  printf("trout: %i\n", trout);
  printf("ttratio: %f\n", ttratio);
  printf("usertree: %i\n", usertree);
  printf("reusertree: %i\n", reusertree);
  printf("weights: %i\n", weights);
  printf("printdata: %i\n", printdata);
  printf("progress: %i\n", progress);
  printf("treeprint: %i\n", treeprint);
  printf("interleaved: %i\n", interleaved);
  printf("whichcode: %i\n", whichcode);
  printf("mulsets: %i\n", mulsets);
  printf("datasets: %li\n", datasets);
  printf("justwts: %i\n", justwts);
  printf("omega: %f\n", omega);
#endif

  for (ith = 1; ith <= datasets; ith++)
  {
    if (datasets > 1)
    {
      fprintf(outfile, "Data set # %ld:\n", ith);
      printf("\nData set # %ld:\n", ith);
    }
    getinput();

    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      long myi;
      for(myi = 0; myi < max_num_sibs; myi++)
      {
        free_pmatrix(myi);
      }
      max_num_sibs = 0;
      maketree();
    }
    curtree->free(curtree);
    bestree->free(bestree);
    priortree->free(priortree);
    bestree2->free(bestree2);
  }

  clean_up();
}


void codml(
  char * Infilename,
  char * Intreename,
  char * Wgtsfilename,
  char * Outfilename,
  char * outfileopt,
  char * Outtreename,
  char * outtreeopt,
  char * TreeUseMethod,
  int UseLengths,
  int InputNuc,
  char * NucSubModel,
  double TTratio,
  double NSratio,
  int useEmpBF,
  double BaseFreqA,
  double BaseFreqC,
  double BaseFreqG,
  double BaseFreqTU,
  char * GeneCode,
  int OneOmega,
  int AdjOmegasCor,
  int OmegaBlockLen,
  int CodonsWgted,
  int NumOmegas,
  // these are explicitly named because JNA doesn't pass arrays gracefully
  double OmegaVal1,
  double OmegaVal2,
  double OmegaVal3,
  double OmegaVal4,
  double OmegaVal5,
  double OmegaVal6,
  double OmegaVal7,
  double OmegaVal8,
  double OmegaVal9,
  // these are explicitly named because JNA doesn't pass arrays gracefully
  double OmegaProb1,
  double OmegaProb2,
  double OmegaProb3,
  double OmegaProb4,
  double OmegaProb5,
  double OmegaProb6,
  double OmegaProb7,
  double OmegaProb8,
  double OmegaProb9,
  int SpeedAn,
  int GlobalRe,
  int RandInput,
  int RandNum,
  int Njumble,
  int OutRoot,
  int OutNum,
  int MultData,
  int MultDSet,
  int NumSeqs,
  int InputSeq,
  int PrintData,
  int PrintInd,
  int PrintTree,
  int DotDiff,
  int WriteTree,
  int RecHypo)
{
  initdata funcs;

  //printf("Hello from CodML!\n"); // JRMdebug
  javarun = true;

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "CodML";

  (void)BaseFreqA;
  (void)BaseFreqC;
  (void) BaseFreqG;
  (void)BaseFreqTU;

  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = codon_node_new;
  funcs.tree_new = codon_tree_new;
  progname = argv[0];

  phylipinit(argc, argv, &funcs, true);

  //printf("init done\n"); // JRMdebug
  /* // JRMdebug
  // from java
  String infile,
  String intree,
  String wgtsfile,
  String outfile,
  String outfileopt,
  String outtree,
  String outtreeopt,
  //String TreeUseMethod,
  //boolean UseLengths,
  //boolean InputNuc,
  //String NucSubModel,
  //double TTratio,
  //double NSratio,
  //boolean useEmpBF,
  double BaseFreqA,
  double BaseFreqC,
  double BaseFreqG,
  double BaseFreqTU,
  //String GeneCode,
  //boolean OneOmega,
  //boolean AdjOmegasCor,
  //int OmegaBlockLen,
  //boolean CodonsWgted,
  //int NumOmegas,

  //these are explicitly named because JNA doesn't pass arrays gracefully
  //double OmegaVal1,
  //double OmegaVal2,
  //double OmegaVal3,
  //double OmegaVal4,
  //double OmegaVal5,
  //double OmegaVal6,
  //double OmegaVal7,
  //double OmegaVal8,
  //ouble OmegaVal9,
  //
  // these are explicitly named because JNA doesn't pass arrays gracefully
  //double OmegaProb1,
  //double OmegaProb2,
  //double OmegaProb3,
  //double OmegaProb4,
  //double OmegaProb5,
  //double OmegaProb6,
  //double OmegaProb7,
  //double OmegaProb8,
  //double OmegaProb9,

  //boolean SpeedAn,
  //boolean GlobalRe,
  //boolean RandInput,
  //int RandNum,
  //int Njumble,
  //boolean OutRoot,
  //int OutNum,
  //boolean MultData,
  //boolean MultDSet,
  //int NumSeqs,
  //boolean InputSeq,
  //boolean PrintData,
  //boolean PrintInd,
  //boolean PrintTree,
  //boolean DotDiff,
  //boolean WriteTree,
  //boolean RecHypo

  // internal variables
  //f84 = true;
  //hky = false;
  //k2p = false;
  //jc = false;
  //rctgry = false;
  ##didchangercat = false;
  //rcategs = 1;
  //auto_ = false;
  ??gama = false;
  //freqsfrom = true;
  //global = false;
  //inputDataIsNucleotideData = true;
  //hypstate = false;
  //improve = false;
  ??invar = false;
  //jumble = false;
  //njumble = 1;
  //lngths = false;
  //lambda = 1.0;
  //outgrno = 1;
  //outgropt = false;
  //trout = true;
  //ttratio = 2.0;
  //usertree = false;
  //reusertree = false;
  //weights = false;
  //printdata = false;
  //progress = true;
  //treeprint = true;
  //interleaved = true;
  ##loopcount = 0;
  //omega = 1;
  //whichcode = universal;
  //mulsets = false;
  //datasets = 1;
  //??justwts
  */
  // transfer java data to local variables
  if (UseLengths != 0)
  {
    lngths = true;
  }
  else
  {
    lngths = false;
  }

  if (!strcmp(TreeUseMethod, "No"))
  {
    // No, use user trees in input file
    reusertree = false;
    usertree   = true;
  }
  else if (!strcmp(TreeUseMethod, "rearrange"))
  {
    // Yes, rearrange on user tree
    reusertree = true;
    usertree   = true;
  }
  else
  {
    // Yes
    reusertree = false;
    usertree   = false;
    lngths     = false; // shut it off no matter what the user entered
  }

  if (useEmpBF != 0)
  {
    freqsfrom = true;
  }
  else
  {
    freqsfrom = false;
  }

  if (!strcmp(NucSubModel, "F84"))
  {
    // F84
    f84 = true;
    hky = false;
    k2p = false;
    jc = false;
  }
  else if (!strcmp(NucSubModel, "HKY"))
  {
    // HKY
    f84 = false;
    hky = true;
    k2p = false;
    jc = false;
  }
  else if (!strcmp(NucSubModel, "Kimura"))
  {
    // Kimura 2-parameter
    f84 = false;
    hky = false;
    k2p = true;
    jc = false;
  }
  else
  {
    // Jukes-Cantor
    f84 = false;
    hky = false;
    k2p = false;
    jc = true;
  }

  if (!strcmp(GeneCode, "Universal"))
  {
    whichcode = universal;
  }
  else if (!strcmp(GeneCode, "Ciliate"))
  {
    whichcode = ciliate;
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

  if (InputNuc != 0)
  {
    inputDataIsNucleotideData = true;
  }
  else
  {
    inputDataIsNucleotideData = false;
  }

  if (AdjOmegasCor != 0)
  {
    auto_ = true;
    lambda = 1.0 / OmegaBlockLen;
  }
  else
  {
    auto_ = false;
    lambda = 1.0;
  }

  if (OneOmega != 0)
  {
    rctgry = false;
  }
  else
  {
    rctgry = true;
  }

  if (CodonsWgted != 0)
  {
    weights = true;
  }
  else
  {
    weights = false;
  }

  if (SpeedAn != 0)
  {
    improve = false;
  }
  else
  {
    improve = true;
  }

  if (GlobalRe != 0)
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
    datasets = NumSeqs;
  }
  else
  {
    mulsets = false;
    datasets = 1;
  }

  if (MultDSet != 0)
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

  if (WriteTree != 0)
  {
    trout = true;
  }
  else
  {
    trout = false;
  }

  if (DotDiff != 0)
  {
    dotdiff = true;
  }
  else
  {
    dotdiff = false;
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

  if (RecHypo != 0)
  {
    hypstate = true;
  }
  else
  {
    hypstate = false;
  }

  // transfer values
  ttratio   = TTratio;
  omega     = NSratio;
  rcategs   = NumOmegas;

  omegavalue = (double *) Malloc(maxcategs * sizeof(double));
  omegavalue[0]   = OmegaVal1;
  omegavalue[1]   = OmegaVal2;
  omegavalue[2]   = OmegaVal3;
  omegavalue[3]   = OmegaVal4;
  omegavalue[4]   = OmegaVal5;
  omegavalue[5]   = OmegaVal6;
  omegavalue[6]   = OmegaVal7;
  omegavalue[7]   = OmegaVal8;
  omegavalue[8]   = OmegaVal9;

  probcat = (double *) Malloc(maxcategs * sizeof(double));
  probcat[0]   = OmegaProb1;
  probcat[1]   = OmegaProb2;
  probcat[2]   = OmegaProb3;
  probcat[3]   = OmegaProb4;
  probcat[4]   = OmegaProb5;
  probcat[5]   = OmegaProb6;
  probcat[6]   = OmegaProb7;
  probcat[7]   = OmegaProb8;
  probcat[8]   = OmegaProb9;

  // everything translated, start the run
  infile = fopen(Infilename, "r");
  outfile = fopen(Outfilename, outfileopt);
  strcpy(outfilename, Outfilename);

  if (progress)
  {
    progfile = fopen("progress.txt", "w");
    fclose(progfile); // make sure it is there for the Java code to detect
    progfile = fopen("progress.txt", "w");
  }

  if (usertree)
  {
    intree = fopen(Intreename, "r");
  }

  if (weights || justwts)
  {
    weightfile = fopen(Wgtsfilename, "r");
  }

  if (trout)
  {
    outtree = fopen(Outtreename, outtreeopt);
    strcpy(outtreename, Outtreename);
  }

  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinit();

  //printf("calling codmlrun\n");
  codmlrun();                           // do the actual work

  FClose(infile);
  FClose(outfile);

  if (weights || justwts)
  {
    FClose(weightfile);
  }

  if (usertree)
  {
    FClose(intree);
  }

  if (trout)
  {
    FClose(outtree);
  }
  //printf("\ndone\n"); // JRMdebug
}


int main(int argc, Char *argv[])
{  /* Codon Maximum Likelihood */
  initdata *funcs;
  javarun = false;

#ifdef MAC
  argc = 1;             /* macsetup("CodML", "");        */
  argv[0] = "CodML";
#endif
  funcs = Malloc(sizeof(initdata));
  funcs->node_new = codon_node_new;
  funcs->tree_new = codon_tree_new;
  phylipinit(argc, argv, funcs, false);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);
  mulsets = false;
  datasets = 1;
  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinit();
  if (weights || justwts)
  {
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  }
  if (trout)
  {
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  }
  codmlrun();  // do the actual work

#if 0
  for (ith = 1; ith <= datasets; ith++)
  {
    if (datasets > 1)
    {
      fprintf(outfile, "Data set # %ld:\n", ith);
      printf("\nData set # %ld:\n", ith);
    }
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
    {
      long myi;
      for(myi = 0; myi < max_num_sibs; myi++)
      {
        free_pmatrix(myi);
      }
      max_num_sibs = 0;
      maketree();
    }
    curtree->free(curtree);
    bestree->free(bestree);
    priortree->free(priortree);
    bestree2->free(bestree2);
  }

  clean_up();
#endif

  free(funcs);
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* Codon Maximum Likelihood */


// End.
