/* Version 4.0. Copyright 2022
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   */


/* debug:  Need to
(1) Put warnings if too many variables (in Within case) to invert matrix --
    use generalized inverses.
(2) complete adding fossil machinery.  Can we use morphometrics in that case?
(3) add ability to compute correlations and regressions from the CovA and CovE
    matrices in the within-species case
(4) add ancestor estimation (need to use variables transformed to independence?
(5) A phylogenetic Procrustes superposition?  An iterative weighted Procrustes??
debug)   */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylip.h"
#include "ml.h"
#include "cont.h"

typedef double** matrix;

#ifndef OLDC
/* function prototypes */
void   contrast_tree_new(struct tree**, long, long, long);
void   contrast_tree_init(struct tree* t, long nonodes, long spp);
struct node* contrast_node_new(node_type, long, long);
void   contrast_node_init(struct node*, node_type, long);
void   getoptions(void);
void   bookstein(void);
double linesearch (double*, double, long, long, boolean*);
double linesearchsz (double*, double, long, long, boolean*);
void   rotateform(long, long);
void   resizeform(long, long);
/* void   felsie(void);   debug:  maybe use later */
void   boasfit (void);
void   nmalterspecimen(long, long, double*);
void   initnmiterate(void);
void   nmiterate(void);
void   felsmorph(void);
void   morph(void);
void   getdata(void);
void   allocrest(void);
void   doinit(void);
void   rotate(long, long, double);
void   resize(long, long, double);
void   copyztotemp(long, long);
void   copytemptoz(long, long);
void   copyztox(long, long);
void   copyztow(long, long);
void   copywtoz(long, long);
void   copyallztox(void);
void   copyallztow(void);
void   copyallwtoz(void);
void   contwithin(void);
void   contbetween(node *);
void   makecontrasts(node *);
void   getcovariances (void);
void   getmeans (void);
void   writesuper(void);
void   writecontrasts(void);
void   getsizevar(void);
void   writemethods(void);
double logdet(double **);
double glogdet(double *);
void   invert(double **);
void   initcovars(boolean);
double normdiff(boolean);
void   matcopy(double **, double **, long);
void   newcovars(boolean, boolean);
void   printcovariances(boolean, boolean);
void   calculateregressions(matrix, long);
void   reportregressions(matrix, matrix);
void   reportlrttests(void);
void   emiterate(boolean, boolean);
void   initcontrastnode(tree *, node **, long, long, long *,
                        long *, initops, pointarray, Char *, Char *, FILE *);
void   readthetree (void);
void   givens(matrix, long, long, long, double, double, boolean);
void   coeffs(double, double, double *, double *, double);
void   tridiag(matrix, long, double);
void   shiftqr(matrix, long, double);
void   qreigen(matrix, long);
void   reportmatrix(matrix, long);
void   reportpca(long);
void   reportlogL (double, long);
void   writescales(void);
void   writevarsz(void);
void   writeallom(void);
void   getshapecovars (double **);
void   reportshapes(void);
void   writemeans (void);
void   contrast_node_copy(node *, node *);
void   contrast_node_init(node*, node_type, long);
void   contrast_node_reinit(node*);
node*  contrast_node_make(tree *, node_type, long);
void   makenewbranches(void);
void   updatebounds(node *);
double evaluate(node *);
#if 0    /* fossil stuff commented out for now */
double fevaluate(double, node *, boolean);
void   locatefossilonbranch (node *, node *, node *, boolean);
boolean placefossilonbranch(node *, node *);
void   bltimetraverse (node *, double *);
void   placetraverse(long, node *);
void   replacefossil(void);
void   placeonefossil(long);
void   makefossilcovars(node *);
void   placeallfossils(void);
#endif        /* end commenting-out of fossil stuff */
void   writereportforonecovar(boolean, matrix, double*);
void   writereports(void);
void   calculatecovsetc(void);
void   maketree(void);
/* function prototypes */
#endif

typedef struct contrast_tree {
  struct ml_tree ml_tree;
} contrast_tree;

typedef struct contrast_node {
  struct ml_node ml_node;
} contrast_node;

/* debug:
typedef struct contrast_node {
  cont_node_type cont_node_var;
} contrast_node;
debig */


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH];
long nonodes, chars, startmchar, endmchar, dropchars=0,
     charsd, charsd2, charsp, charspp, df, bestplace,*sample, *nset,
     contnum, i, *eigorder, *fossilsp, contno, n1, n2, numtrees,
     ndatas, morphchars, dimensions, numfossils, ndiv, nmult, nparams,
     nmpoints, nmworstone, nmsecondworstone, nmbestone;
long ith, ithwas, jth, jthwas, kth;
phenotype3 **x, **y, **z, **w, **specsize, *rotation, **cntrast, *ssqcont, meanz;
double **vara, **vare, **oldvara, **oldvare, **Bax, **Bex, **temp1, **temp2,
  **temp3, **temp4, **temp5, **temp6, **temp7, *temp8, **parameters,
  *nmmeanparams, *reflectedparameters, *expandedparameters, *slopes,
  *contractedparameters, *shrunkparameters, *loglike, **sumprod, *mean, *umean,
  **regressions, **correlations, *allom, **eigvecs, *eig, *fossiltime;
double logL, bestlogL, logLvara, logLnocorr, logLnovara, multiplier,
  multiplier0, jacobian, pi=3.141592653, minmultmult, maxmultmult, bestwhere,
  bestmult, howworse, tolerance, tollimit, tolstart, unorm, varz, zsum,
  incparam0, startlogL, endlogL, nmepsilon, nmupsilon, worstlogL, alpha,
  gammma, rho, sigma;
boolean nomorph, bookmorph, mlrots, justprocrust, centroidsize, mlsizes,
  linearsize, morphall, sizes, sizes, sizechar, shapes, nophylo, printdata,
  progress, reg, multrees, muldata, treeswithin,
  datawithin, cross, fossil, inferscale, varywithin,
  nocorr, writecont, bifurcating, pca, firsttime,
  firstplace, reportplacefossils, changed, somechanged,
  superposition, weightsuper, supshapesonly, omitheaders, includesppchars,
  toladjusted, *isfossil;
Char ch;

/* Local variables for maketree, propagated globally for C version: */
tree *curtree;

/* Variables declared just to make treeread happy */
boolean haslengths, goteof, first;
double trweight;


void contrast_tree_new(struct tree** treep, long nonodes, long spp, long treesize)
{
  /* set up variables and then set up identities of functions */
  tree* t;

  generic_tree_new(treep, nonodes, spp, treesize);
  t = *treep;
  t->setupfunctions = generic_tree_setupfunctions;
  generic_tree_init(t, nonodes, spp);
} /* contrast_tree_new */


void contrast_tree_init(struct tree* t, long nonodes, long spp)
{
  /* set up functions for a contrast_tree */

  t->try_insert_ = ml_tree_try_insert_;
  t->get_fork = generic_tree_get_fork;
} /* contrast_tree_init */


struct node* contrast_node_new(node_type type, long index, long nodesize)
{
  /* make new contrast_node */
  struct node *n;

  nodesize = (long)sizeof(contrast_node);
  n = ml_node_new(type, index, nodesize);
  contrast_node_init(n, type, index);
  return n;
} /* contrast_node_new */


void contrast_node_init(node* n, node_type type, long index)
{
  /* initialize a contrast_node */
  contrast_node *cn = (contrast_node*)n;

  generic_node_init(&(cn->cont_node_var.node_var), type, index);
  ((cont_node_type*)(cn))->view = (phenotype3)Malloc((long)charspp * sizeof(double));
  cn->cont_node_var.node_var.copy = contrast_node_copy;
  cn->cont_node_var.node_var.reinit = contrast_node_reinit;
} /* contrast_node_init */



void getoptions(void)
{
  /* interactively set options */
  long loopcount, i;
  Char ch;
  boolean done;

  nomorph = true;
  bookmorph = false;
  mlrots = false;
  justprocrust = false;
  centroidsize = false;;
  sizes = false;
  mlsizes = true;
  linearsize = false;
  dimensions = 2;
  numtrees = 1;
  ndatas = 1;
  firsttime = true;
  muldata = false;
  multrees = false;
  cross = false;
  treeswithin = false;
  datawithin = false;
  nophylo = false;
  printdata = false;
  progress = true;
  fossil = false;
  numfossils = 0;
  ndiv = 100;
  inferscale = true;
  multiplier = 1.0;
  nmult = 10.0;
  reportplacefossils = true;
  howworse = 3.0;
  varywithin = false;
  nocorr = false;
  writecont = false;
  pca = true;
  superposition = false;
  supshapesonly = false;
  weightsuper = false;
  omitheaders = false;
  includesppchars = false;
  dropchars = 0;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nContinuous character comparative analysis, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  B               Morphometric transformations?");
    if (nomorph)
      printf("  No, not morphometric data\n");
    else {
      if (bookmorph)
        printf("  Bookstein transform\n");
      else {
        if (mlrots)
          printf("  ML optimum rotations\n");
        else {
          if (justprocrust && centroidsize)
            printf("  Procrustes transform\n");
          else
            printf("  Boas transform\n");
        }
      }
    }
    if (!nomorph) {
      printf("  Z                                Infer sizes?");
      if (sizes && !bookmorph) {
        printf("  Yes\n");
        printf("  I       Method of inferring sizes and shapes?");
        if (mlsizes)
          printf("  by ML\n");
        else if (linearsize)
            printf("  linear inference of sizes\n");
          else {
            if (centroidsize)
              printf("  centroid size\n");
          }
        printf("  E                     Report size variation?");
        if (sizes)
          printf("  Yes, in addition to shape\n");
        else
          printf("  No, just covariation of coordinates\n");
      } else {
        printf("  No\n");
      }
    }
    printf("  W        Within-population variation in data?");
    if (varywithin)
      printf("  Yes, multiple individuals\n");
    else {
      printf("  No, species values are means\n");
    }
    printf("  R     Print out correlations and regressions?  %s\n",
           (reg ? "Yes" : "No"));
    if (varywithin) {
      printf("  A      LRT test of no phylogenetic component?");
      if (nophylo)
        printf("  Yes, with and without VarA\n");
      else
        printf("  No, just assume it is there\n");
    }
    printf("  T                 LRT test of no correlation?");
    if (nocorr)
      printf("  Yes, test correlations\n");
    else
      printf("  No\n");
    if (!varywithin)
      printf("  C                        Print out contrasts?  %s\n",
             (writecont? "Yes" : "No"));
    if (!varywithin) {
      printf("  X               Get Principal Component Axes?");
      if (pca)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
#if 0            /* comment this out for now as undeveloped */
    printf("  F                Infer position of fossil(s)? ");
    if (fossil)
      printf(" Yes\n");
    else
      printf(" No\n");
    if (fossil) {
      printf("  G       Report progress of placing fossil(s)? ");
      if (reportplacefossils)
        printf(" Yes\n");
      else
        printf(" No\n");
      printf("  S    Infer time scale per unit branch length? " );
      if (inferscale) {
        printf(" Yes\n");
        printf("  V    Number of values of time scaling to try: ");
        printf(" %ld\n", nmult);
      }
      else {
        printf(" No, use multiplier provided\n");
        printf("  L  Multiplier to get time from branch length: ");
        if (multiplier < 10.0)
          printf(" %6.5f\n", multiplier);
        else if (multiplier < 100.00)
          printf(" %7.5f\n", multiplier);
        else if (multiplier < 1000.00)
          printf(" %8.5f\n", multiplier);
        else if (multiplier < 10000.00)
          printf(" %9.5f\n", multiplier);
        else printf(" %10.5f\n", multiplier);
      }
      printf("  N     Number of places to try on each branch: ");
      printf(" %ld\n", ndiv);
      printf("  H   How close to best-yet logL to show value: ");
      printf(" %4.2f\n", howworse);
    }
#endif
    printf("  M        Multiple trees?  Multiple data sets?");
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
    printf("  P           Write out inferred superposition?  %s\n",
           superposition ? "Yes" : "No");
    if (superposition) {
      printf("  U                      Type of superposition?  %s\n",
             weightsuper? "Weighted" : "Centroid");
      if (sizes)
        printf("  Q               Superposition of shapes only?  %s\n",
               supshapesonly? "Shapes, they are resized" : "No, forms, not resized");
      printf("  O           Superpositions computer readable?  ");
      if (!omitheaders)
        printf("No, include headers\n");
      else if (includesppchars)
          printf("Yes, PHYLIP format\n");
        else
          printf("Yes, R-compatible format\n");
    }
    printf("  0         Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC"  :
           ansi  ? "ANSI"    : "(none)");
    printf("  1          Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("BDEIZRATMFGSVLNHWCXPUQOG120", ch) != NULL) {
        switch (ch)
        {
          case 'B':
            if (nomorph) {
              nomorph = false;
              mlrots = true;
              }
            else {
              if (mlrots) {
              mlrots = false;
              justprocrust = true;
              }
              else {
                if (justprocrust) {
                  justprocrust = false;
                  bookmorph = true;
                  }
                else
                  if (bookmorph) {
                    bookmorph = false;
                    nomorph = true;
                  }
                }
              }
            break;

          case 'D':
            if (dimensions == 2)
              dimensions = 3;
            else
              dimensions = 2;
            break;

          case 'Z':
            sizes = !sizes;
            if (sizes) {
              centroidsize = justprocrust;
              mlsizes = !justprocrust;
            }
            break;

          case 'I':
            if (sizes) {
              if (mlsizes) {
                mlsizes = false;
                linearsize = true;
              }
              else if (linearsize) {
                  centroidsize = true;
                  linearsize = false;
                } else if (centroidsize) {
                    mlsizes = true;
                    centroidsize = false;
                  }
            }
            break;

          case 'E':
            sizechar = !sizechar;
            break;

          case 'R':
            reg = !reg;
            break;

          case 'A':
            nophylo = !nophylo;
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

          case 'C':
            writecont = !writecont;
            break;

          case 'W':
            varywithin = !varywithin;
            break;

          case 'T':
            nocorr = !nocorr;
            break;

          case 'X':
            pca = !pca;
            break;

          case 'F':
            fossil = !fossil;
            break;

          case 'G':
            reportplacefossils = !reportplacefossils;
            break;

          case 'S':
            inferscale = !inferscale;
            break;

          case 'L':
            loopcount = 0;
            do {
              printf("Multiplier to convert branch length to time scale?\n");
              fflush(stdout);
              if(scanf("%lf%*[^\n]", &multiplier)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (multiplier < 0.0);
            break;

          case 'N':
            loopcount = 0;
            do {
              printf("Number of points on each branch to try?\n");
              fflush(stdout);
              if(scanf("%ld%*[^\n]", &ndiv)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (ndiv <= 0.0);
            break;

          case 'V':
            loopcount = 0;
            do {
              printf("Number of values of time scaling to try on each branch?\n");
              fflush(stdout);
              if(scanf("%ld%*[^\n]", &nmult)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (nmult <= 0.0);
            break;

          case 'H':
            loopcount = 0;
            do {
              printf("How many units of logL worse than best to show?\n");
              fflush(stdout);
              if(scanf("%lf%*[^\n]", &howworse)) {} // Read number and scan to EOL.
              (void)getchar();
              countup(&loopcount, 10);
            } while (howworse <= 0.0);
            break;

          case 'P':
            superposition = !superposition;
            break;

          case 'U':
            weightsuper = !weightsuper;
            break;

          case 'Q':
            supshapesonly = !supshapesonly;
            break;

          case 'O':
            if (!omitheaders)
              omitheaders = true;
            else {
              if (!includesppchars)
                includesppchars = true;
              else {
                includesppchars = false;
                omitheaders = false;
              }
            }
            break;

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
  printf("\n");
  if (bookmorph || mlrots || justprocrust) {
    loopcount = 0;
    done = false;
    while (!done) {
      printf("Are all characters morphometric coordinates (type Y or N)?\n");
      if(scanf("%c%*[^\n]", &ch)) {}    // Read char and scan to EOL.
      (void)getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      done = (ch == 'Y') || (ch == 'N');
      countup(&loopcount, 100);
    };
    morphall = (ch == 'Y');
    if (!morphall) {
      loopcount = 0;
      done = false;
      while (!done) {
        printf("Number of character which is first morphometric coordinate?\n");
        if(scanf("%ld%*[^\n]", &startmchar)) {} // Read number and scan to EOL.
        (void)getchar();
        done = (startmchar > 0);
        countup(&loopcount, 100);
      }
      loopcount = 0;
      done = false;
      while (!done) {
        printf("Number of character which is last morphometric coordinate?\n");
        if(scanf("%ld%*[^\n]", &endmchar)) {} // Read number and scan to EOL.
        (void)getchar();
        done = (endmchar > startmchar);
        countup(&loopcount, 100);
      }
    } else {
      startmchar = 1;
      endmchar = chars;
    }
  }
  else {
    startmchar = 0;
    endmchar = 0;
  };
  if (bookmorph) { /* Bookstein transformation without log(size) character */
    dropchars = 4;
    charsd = chars - dropchars;       /* rank of superposed forms */
    charsd2 = chars - (dropchars-2);  /* number of non-location dimensions */
  } else {
    if (sizes && (!sizechar))
      dropchars = 4;
    else if (mlrots || justprocrust)
      dropchars = 3;
      else 
        dropchars = 0;
    charsd = chars - dropchars;       /* rank of superposed forms */
    charsd2 = chars - (dropchars-1);  /* number of non-location dimensions */
  }
  charsp = chars;                         /* an inferred size character */
  if (sizes && sizechar)  /* make room for an inferred size character */ 
    charspp = chars + 1;
  else
    charspp = chars;
  if (fossil) {
    printf("How many fossil species?\n");
    if(scanf("%ld%*[^\n]", &numfossils)) {} // Read number and scan to EOL.
    (void)getchar();
  }
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
  if (nocorr) {
    printf("\nTo test whether two sets of characters are uncorrelated,\n");
    printf(" type how many characters are in the first set:\n");
    printf(" (all other characters will be assumed to be in the other set).\n");
    if(scanf("%ld%*[^\n]", &n1)) {}     // Read number and scan to EOL.
    nset = (long *)Malloc((long)(charspp) * sizeof(long));
    for (i = 0; i < chars; i++)  /* set up indicators of character sets */
      nset[i] = 0;
    printf("\nType the numbers of the individual characters in the set,\n");
    printf(" on one line and separated by blanks:\n");
    for (i = 0; i < n1; i++) {
      if(scanf("%ld", &n2)) {}          // Read number and scan to EOL.
      nset[n2-1] = 1;
    }
  }
}  /* getoptions */


void getdata(void)
{
  /* read species data */
  long i, j, k, l, spp2, chars2, nonodes2;

  if (printdata && (reg || pca || nophylo || nocorr)) {
    fprintf(outfile, "Name");
    fprintf(outfile, "                       Phenotypes\n");
    fprintf(outfile, "----");
    fprintf(outfile, "                       ----------\n\n");
  }
  x = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3 *));
  y = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3 *));
  z = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3 *));
  w = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3 *));
  specsize = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3));
  cntrast = (phenotype3 **)Malloc((long)spp * sizeof(phenotype3 *));
  ssqcont = (phenotype3 *)Malloc((long)spp * sizeof(phenotype3 *));
  rotation = (phenotype3 *)Malloc((long)spp * sizeof(phenotype3 *));
  contnum = spp-1;
  if (fossil) {    /* set aside space for list of fossil species, times */
    fossilsp = (long *)Malloc((long)numfossils * sizeof(long));
    fossiltime = (double *)Malloc((long)numfossils * sizeof(double));
    isfossil = (boolean *)Malloc((spp+1) * sizeof(boolean)); /* is it a fossil? */
    for (i = 1; i <= spp; i++)  /* set up array showing which are fossils */
      isfossil[i] = false;
  }
  if (ith > 1) {
    inputnumbers(&spp2, &chars2, &nonodes2, 1);
    if (spp2 != spp) {
      printf("\nERROR:  Number of species in data set %ld is not the same\n", ith);
      printf("        as in all the previous data sets.\n");
    }
    if (chars2 != chars) {
      printf("\nERROR:  Number of characters in data set %ld is not the same\n", ith);
      printf("        as in all the previous data sets.\n");
    }
    if ((spp2 != spp) || (chars2 != chars)) {
      printf("\n");
      exxit(-1);
    }
  }
  scan_eoln(infile);
  for (i = 0; i < spp; i++) {
    initname(i);
    if (varywithin) {
      if (fscanf(infile, "%ld", &sample[i]) == 1) {
        contnum += sample[i]-1;
        scan_eoln(infile);
      }
    }
    else sample[i] = 1;
    if (printdata)
      for(j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
    x[i] = (phenotype3 *)Malloc((long)sample[i] * sizeof(phenotype3));
    y[i] = (phenotype3 *)Malloc((long)sample[i] * sizeof(phenotype3));
    z[i] = (phenotype3 *)Malloc((long)sample[i] * sizeof(phenotype3));
    w[i] = (phenotype3 *)Malloc((long)sample[i] * sizeof(phenotype3));
    specsize[i] = (phenotype3 *)Malloc((long)(sample[i] * sizeof(double)));
    cntrast[i] = (phenotype3 *)Malloc((long)(sample[i] * sizeof(phenotype3)));
    ssqcont[i] = (double *)Malloc((long)(sample[i] * sizeof(double)));
    rotation[i] = (double *)Malloc((long)(sample[i] * sizeof(double)));
    for (k = 0; k <= sample[i]-1; k++) {
      x[i][k] = (phenotype3)Malloc((long)charspp * sizeof(double));
      y[i][k] = (phenotype3)Malloc((long)charspp * sizeof(double));
      z[i][k] = (phenotype3)Malloc((long)charspp * sizeof(double));
      w[i][k] = (phenotype3)Malloc((long)charspp * sizeof(double));
      cntrast[i][k] = (phenotype3)Malloc((long)charspp * sizeof(double));
      specsize[i][k] = (phenotype3)Malloc((long)((long)charspp * sizeof(double)));
      for (j = 1; j <= chars; j++) {
        if (eoln(infile))
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][k][j - 1]) != 1) {
          printf("Error in input file at species %ld, sample %ld, character %ld.\n", i+1, k+1, j);
          exxit(-1);
        }
printf("x[%ld][%ld][%ld] = %12.9f\n", i, k, j-1, x[i][k][j-1]);
        if (printdata) {
          fprintf(outfile, " %14.9f", x[i][k][j - 1]);
          if (j % 6 == 0) {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth; l++)
              putc(' ', outfile);
          }
        }
        y[i][k][j-1] = x[i][k][j-1];    /* copy it into  y */
      }
    }
    /* got all the data we need for that member, */
    /* read the next line if we still have more members to define */
    if (k != sample[i]-1)
      scan_eoln(infile);
    putc('\n', outfile);
  }
  if (printdata)
    putc('\n', outfile);
  checknames(spp);                      // Check NAYME array for duplicates.
  if (fossil) { /* debug: ultimately this will be replaced by reading from a file */
    printf("Type the numbers of the fossil species whose connections will be inferred\n");
    printf(" (separated by blanks)\n");
    for (i = 0; i < numfossils; i++) {
      if (scanf("%ld", &fossilsp[i]) != 1) {
        if (i == 0)
          printf("Error in list of fossils at first one\n");
        if (i == 1)
          printf("Error in list of fossils at second one\n");
        if (i == 2)
          printf("Error in list of fossils at third one\n");
        if (i > 2)
          printf("Error in list of fossils at %ldth one\n", i+1);
        exxit(-1);
      }
      if (fossilsp[i] > spp) {
        printf("Error: number of a fossil (%ld) is greater than total number of species (%ld)\n",fossilsp[i], spp);
        exxit(-1);
      }
    }
    (void)getchar();
    printf("Type the ages of these fossil species (separated by blanks)\n");
    for (i = 0; i < numfossils; i++) {
      if (scanf("%lf", &fossiltime[i]) != 1) {
        if (i == 0)
          printf("Error in age of fossil at first one\n");
        if (i == 1)
          printf("Error in age of fossil at second one\n");
        if (i == 2)
          printf("Error in age of fossil at third one\n");
        if (i > 2)
          printf("Error in age of fossil at %ldth one\n", i+1);
        exxit(-1);
      }
    }
    (void)getchar();
    for (i = 0; i < numfossils; i++)
      isfossil[fossilsp[i]] = true;
  }
}  /* getdata */


void allocrest(void)
{
  long i, j;

  /* otherwise if individual variation, these are allocated in getdata  */
  sample = (long *)Malloc((long)spp * sizeof(long));
  nayme = (naym *)Malloc((long)spp * sizeof(naym));
  vara = (double **)Malloc((long)charspp * sizeof(double *));
  oldvara = (double **)Malloc((long)charspp * sizeof(double *));
  vare = (double **)Malloc((long)charspp * sizeof(double *));
  oldvare = (double **)Malloc((long)charspp * sizeof(double *));
  Bax = (double **)Malloc((long)charspp * sizeof(double *));
  Bex = (double **)Malloc((long)charspp * sizeof(double *));
  temp1 = (double **)Malloc((long)charspp * sizeof(double *));
  temp2 = (double **)Malloc((long)charspp * sizeof(double *));
  temp3 = (double **)Malloc((long)charspp * sizeof(double *));
  temp4 = (double **)Malloc((long)charspp * sizeof(double *));
  temp5 = (double **)Malloc((long)charspp * sizeof(double *));
  temp6 = (double **)Malloc((long)charspp * sizeof(double *));
  temp7 = (double **)Malloc((long)charspp * sizeof(double *));
  for (i = 0; i < charspp; i++) {
    vara[i] = (double *)Malloc((long)charspp * sizeof(double));
    oldvara[i] = (double *)Malloc((long)charspp * sizeof(double));
    vare[i] = (double *)Malloc((long)charspp * sizeof(double));
    oldvare[i] = (double *)Malloc((long)charspp * sizeof(double));
    Bax[i] = (double *)Malloc((long)charspp * sizeof(double));
    Bex[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp1[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp2[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp3[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp4[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp5[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp6[i] = (double *)Malloc((long)charspp * sizeof(double));
    temp7[i] = (double *)Malloc((long)charspp * sizeof(double));
  }
  sumprod = (double **)Malloc((long)charspp * sizeof(double *));
  for (i = 0; i < charspp; i++) {
    sumprod[i] = (double *)Malloc((long)charspp * sizeof(double));
    for (j = 0; j < charspp; j++)
      sumprod[i][j] = 0.0;
  }
  mean = (double *)Malloc((long)charsp * sizeof(double));
  umean = (double *)Malloc((long)charsp * sizeof(double));
  regressions = (double **)Malloc((long)charspp * sizeof(double *));
  for (i = 0; i < charspp; i++) {
    regressions[i] = (double *)Malloc((long)charspp * sizeof(double));
  }
  correlations = (double **)Malloc((long)charspp * sizeof(double *));
  for (i = 0; i < charspp; i++) {
    correlations[i] = (double *)Malloc((long)charspp * sizeof(double));
  }
  mean = (double *)Malloc((long)charsp * sizeof(double));
  allom = (double *)Malloc((long)charsp * sizeof(double));
  eigvecs =  (double **)Malloc(charspp * sizeof(double *));   /* eigenvectors */
  for (i = 0; i < charspp; i++)
    eigvecs[i] = (double *)Malloc(charspp * sizeof(double));
  eig = (double *)Malloc(charspp * sizeof(double));      /* eigenvalues */
  eigorder = (long *)Malloc(charspp * sizeof(long));    /* sort tags for them */
  temp8 = (double *)Malloc((long)charspp * sizeof(double));  /* temporary storage for x[i][] */
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  allocrest();
}  /* doinit */


void rotate(long i, long j, double sintheta) {
  /* rotate the j-th sample of the i-th form in z by theta */
  long k;
  double costheta, tempx, tempy;

  costheta = sqrt(1.0 - sintheta*sintheta);
  for (k = startmchar-1; k <= endmchar-2; k += 2) {
    tempx = costheta * z[i][j][k] - sintheta * z[i][j][k+1];
    tempy = sintheta * z[i][j][k] + costheta * z[i][j][k+1];
    z[i][j][k] = tempx;
    z[i][j][k+1] = tempy;
  }
  rotation[i][j] += asin(sintheta);
  copyztox(i, j);
} /* rotate */


void resize(long i, long j, double logsz) {
  /* resize the j-th sample of the i-th form in z by logsz,
     alter the log-size in x  */
  long k;

  for (k = startmchar-1; k <= endmchar-1; k += 1) {
    z[i][j][k] *= exp(-logsz);     /* resize the centered specimens */
  }
  if (sizechar)
    z[i][j][charspp-1] += logsz;  /* change the log-size accordingly */
/* debug:   to allow compiling   specsize[i][j] = z[i][j][charspp-1];   also separate log specimen size */
  copyztox(i, j);              /* synch  x  with  z  */
} /* resize */


void copyztotemp(long i, long j) {
  /* copy j-th sample of the i-th morphometric form from x to temporary
     storage in array temp8 */
  long k, m;

  m = 0;
  for (k = 0; k < startmchar-1; k++) { /* copy those before the morphometrics */
    temp8[m] = z[i][j][k];
    m++;
  }
  for (k = startmchar-1; k < endmchar; k++) { /* copy, dropping as needed */
    temp8[m] = z[i][j][k];
    m++;
  }
  for (k = endmchar; k < charspp; k++) { /* copy the rest */
    temp8[m] = z[i][j][k];
    m++;
  }
} /* copyztotemp */


void copytemptoz(long i, long j) {
  /* retrieve j-th sample of the i-th morphometric form from temporary
     storage in array temp8 */
  long k, m;

  m = 0;
  for (k = 0; k < startmchar-1; k++) { /* copy those before the morphometrics */
    z[i][j][k] = temp8[m];
    m++;
  }
  for (k = startmchar-1; k < endmchar; k++) { /* copy, dropping as needed */
    z[i][j][k] = temp8[m];
    m++;
  }
  for (k = endmchar; k < charspp; k++) { /* copy the rest */
    z[i][j][k] = temp8[m];
    m++;
  }
} /* copytemptoz */


void copyztox(long i, long j)
{
  /* copy j-th sample of the i-th morphometric form from z to x
   * -- note there are  charspp  characters, including maybe log(size) 
   */
  long k, m;

  m = 0;
  if (!nomorph) {
    for (k = 0; k < startmchar-1; k++) { /* copy those before morphometrics */
      x[i][j][m] = z[i][j][k];
      m++;
    }
    for (k = startmchar-1; k < endmchar; k++) { /* copy, dropping as needed */
      x[i][j][m] = z[i][j][k];
      m++;
    }
    for (k = endmchar; k < charspp; k++) { /* copy the rest */
      x[i][j][m] = z[i][j][k];
      m++;
    }
  } else {
    for (k = 0; k < charspp; k++) { /* copy them all */
      x[i][j][m] = z[i][j][k];
      m++;
    }
  }
} /* copyztox */


void copyztow(long i, long j) {
  /* copy j-th sample of the i-th morphometric form from z to w
    -- note there are  charspp  characters copied, including maybe log(size) */
  long k, m;

  m = 0;
  if (!nomorph) {
    for (k = 0; k < startmchar-1; k++) { /* copy those before morphometrics */
      w[i][j][m] = z[i][j][k];
      m++;
    }
    for (k = startmchar-1; k < endmchar; k++) { /* copy, dropping as needed */
      w[i][j][m] = z[i][j][k];
      m++;
    }
    for (k = endmchar; k < charspp; k++) { /* copy the rest */
      w[i][j][m] = z[i][j][k];
      m++;
    }
  } else {
    for (k = 0; k < charspp; k++) { /* copy them all */
      w[i][j][m] = z[i][j][k];
      m++;
    }
  }
} /* copyztow */


void copywtoz(long i, long j) {
  /* copy j-th sample of the i-th morphometric form from w to z
    -- note there are  charspp  characters copied, including maybe log(size) */
  long k, m;

  m = 0;
  if (!nomorph) {
    for (k = 0; k < startmchar-1; k++) { /* copy those before morphometrics */
      z[i][j][m] = w[i][j][k];
      m++;
    }
    for (k = startmchar-1; k < endmchar; k++) { /* copy, dropping as needed */
      z[i][j][m] = w[i][j][k];
      m++;
    }
    for (k = endmchar; k < charspp; k++) { /* copy the rest */
      z[i][j][m] = w[i][j][k];
      m++;
    }
  } else {
    for (k = 0; k < charspp; k++) { /* copy them all */
      z[i][j][m] = w[i][j][k];
      m++;
    }
  }
} /* copywtoz */


void copyallztox(void) {
/* copy all forms from z to x */
  long i;

  for (i = 0; i < spp; i++)
    copyztox(i, 0);
} /* copyallztox */


void copyallztow(void) {
/* copy all forms from z to w */
  long i;

  for (i = 0; i < spp; i++)
    copyztow(i, 0);
} /* copyallztow */


void copyallwtoz(void) {
/* copy all forms from w to z */
  long i;

  for (i = 0; i < spp; i++)
    copywtoz(i, 0);
} /* copyallwtoz */


void bookstein(void)
{
  /* carry out the Bookstein Transformation to rotate forms */
  long i, j, k, l, n;
  double sumprodmuxmuy, sumsqmux, sumsqmuy, tempx, tempy,
    theta, costheta, sintheta, shrinkto, sum;
  phenotype3 meanform = NULL;
  double *J3 = NULL, *J4 = NULL, *tempsum = NULL;
  double **IJTJ = NULL;

  if (firsttime)
    meanform = (phenotype3)Malloc((long)morphchars * sizeof(double));
  for (i = startmchar-1; i <= endmchar-1; i++) /* initialize form mean to 0 */
    meanform[i] = 0.0;
  n = 0;
  for (i = 0; i < spp; i++)              /* compute means of forms */
    for (j = 0; j < sample[i]; j++) {
      n++;
      for (k = startmchar-1; k <= endmchar-1; k++)
        meanform[k] += z[i][j][k];       /* sum over individuals */
    }
  for (k = startmchar-1; k <= endmchar-1; k ++)
    meanform[k] /= n;                    /* make into means */
  /* standardize and rotate mean forms */
  sumprodmuxmuy = 0.0;
  sumsqmux = 0.0;
  sumsqmuy = 0.0;
  for (k = startmchar-1; k <= endmchar-2; k += 2) {
    tempx = meanform[k];
    tempy = meanform[k+1];
    sumprodmuxmuy += tempx * tempy;
    sumsqmux += tempx * tempx;
    sumsqmuy += tempy * tempy;
  }
  if (fabs(sumsqmux - sumsqmuy) < 1.0e-10)  /* rotation to make axes OK */
    theta = pi/4;
  else theta = 0.5 * atan(2.0*sumprodmuxmuy/(sumsqmuy-sumsqmux));
  costheta = cos(theta);
  sintheta = sin(theta);
  shrinkto = 1.0/sqrt(sumsqmux+sumsqmuy);
  for (k = startmchar-1; k <= endmchar-2; k += 2) { /* rotate, scale means */
    tempx = costheta * meanform[k] - sintheta * meanform[k+1];
    tempy = sintheta * meanform[k] + costheta * meanform[k+1];
    meanz[k] = tempx;                 /* save an unrescaled copy too */
    meanz[k+1] = tempy;
    meanform[k] = tempx * shrinkto;
    meanform[k+1] = tempy * shrinkto;
  }
  for (i = 0; i < spp; i++)         /* rotate the z's by the same amount */
    for (j = 0; j < sample[i]; j++)
      for (k = startmchar-1; k <= endmchar-2; k += 2) {
        tempx = costheta * z[i][j][k] - sintheta * z[i][j][k+1];
        tempy = sintheta * z[i][j][k] + costheta * z[i][j][k+1];
        z[i][j][k] = tempx;
        z[i][j][k+1] = tempy;
      }
  if (firsttime) {
    IJTJ = (double **)Malloc((long)morphchars * sizeof(double *));
    for (i = 0; i < morphchars; i++)
      IJTJ[i] = (double *)Malloc((long)morphchars * sizeof(double));
    J3 = (double *)Malloc((long)morphchars * sizeof(double));
    J4 = (double *)Malloc((long)morphchars * sizeof(double));
    tempsum = (double *)Malloc((long)morphchars * sizeof(double));
  }
  for (i = 0; i < morphchars; i++) {   /* initialize transformation matrix */
    for (j = 0; j < morphchars; j++) {
      IJTJ[i][j] = 0.0;
    }
    IJTJ[i][i] = 1.0;
  }
  for (i = 0; i < morphchars; i++)
    J3[i] = meanform[startmchar-1+i];   /* mu_x1, mu_y1, mu_x2, ... */
  for (i = 0; i < morphchars; i += 2) {
    J4[i] = -meanform[startmchar+i];  /* -mu_y1, mu_x1, -mu_y2, ... */
    J4[i+1] = meanform[startmchar-1+i];
  }
  for (i = 0; i < morphchars; i++)      /* set up  I-J^t J  matrix */
    for (j = 0; j < morphchars; j++)
      IJTJ[i][j] -= J3[i]*J3[j] + J4[i]*J4[j];
  for (i = 0; i < spp; i++) {          /* do the Bookstein I-J^t J transform */
    for (j = 0; j < sample[i]; j++) {
      for (k = 0; k < morphchars; k++) {
        sum = 0.0;
        for (l = 0; l < morphchars; l++)
          sum += IJTJ[k][l] * z[i][j][startmchar-1+l];
        tempsum[k] = sum;
      }
      for (k = 0; k < morphchars; k++)   /* ... and also add mean back in */
        z[i][j][startmchar-1+k] = tempsum[k] + meanz[k];
    }
  }
  for (i = 0; i < spp; i++)         /* unrotate the z's by theta */
    for (j = 0; j < sample[i]; j++) {
      for (k = startmchar-1; k <= endmchar-2; k += 2) {
        tempx =  costheta * z[i][j][k] + sintheta * z[i][j][k+1];
        tempy = -sintheta * z[i][j][k] + costheta * z[i][j][k+1];
        z[i][j][k] = tempx;
        z[i][j][k+1] = tempy;
      }
      copyztox(i, j);
    }
} /* bookstein */


double linesearch (double* theta, double tol, long n, long m, boolean* changed) {
  /* line search in one direction in steps of  tol  of species n, individual m.
     return theta and whether changed  */
  boolean better;
  double logL;

  copyztotemp(n, m);  /* set aside best-yet value */
  copyztox(n, m);
  do {
    rotate(n, m, tol);
    *theta += tol;
    copyztox(n, m);              /* synch  x  with  z  */
    logL = evaluate(curtree->root);
    better = (logL > bestlogL);
    if (better) {
      bestlogL = logL;
      *changed = true;
      copyztotemp(n, m);  /* set aside best-yet value */
    }
  } while (better);
  *theta -= tol;
  copytemptoz(n, m);  /* restore best-yet value */
  copyztox (n, m);    /* and put it in x too  */
  bestlogL = evaluate(curtree->root);
  logL = bestlogL;
  return(*theta);
} /* linesearch */


double linesearchsz (double* logsz, double tol, long n, long m, boolean* changed) {
  /* do line search in one direction in steps of tol. return delta-log(size) and
   * whether changed  */
  boolean better;
  double logszold, logszcurrent, logsznew;

  logszold = *logsz;
  logszcurrent = logszold;
  copyztotemp(n, m);  /* set aside best-yet value */
  do {
    logsznew = logszcurrent + tol;
    resize(n, m, tol);
    logL = evaluate(curtree->root);
    better = (logL > bestlogL);
    if (better) {
      bestlogL = logL;
      logszcurrent = logsznew;
      copyztotemp(n, m);  /* set aside best-yet value */
      *changed = true;
    }
  } while (better);
  copytemptoz(n, m);  /* set aside best-yet value */
  copyztox (n, m);    /* and put it in x too  */
  bestlogL = evaluate(curtree->root);
  return(logszcurrent);
} /* linesearchsz */


// void rotateform(long i, long j) {
//  /* rotate one tip form in tip i in steps of tolerance to optimal angle */
//  double theta, tol2;
//
//  rotation[i][j] = 0.0;
//  tol2 = tolerance;
//  changed = false;
//  do {
//    theta = 0.0;  /* start with no rotation? */
//    theta = linesearch(&theta, tol2, i, j, &changed);  /* search upwards */
//    if (!changed) {
//      theta = linesearch(&theta, -tol2, i, j, &changed); /* and downwards */
//      }
//    tol2 = tol2/10.0;
//    somechanged = somechanged || changed;
//    rotation[i][j] += theta;   /* store rotation */
//  }
//  while (tol2 > tollimit);
//  if (changed)
//    printf(" angle of %ld  changes to: %15.12f,      ln L now = %12.8f\n", i,
//rotation[i][j], bestlogL);
//} /* rotateform */
//
//
//void resizeform(long i, long j) {
//  /* change size of tip form  i  in steps of tolerance to optimal size */
//  double logsize;
//  double tol2;
//
//  logsize = x[i][j][charspp-1];    /* start size after Procrustes */
//  changed = false;
//  tol2 = tolerance;
//  do {
//    logsize = linesearchsz(&logsize, tol2, i, j, &changed);  /* search up */
//    if (!changed) {
//      logsize = linesearchsz(&logsize, -tol2, i, j, &changed); /* and down */
//    }
//    tol2 = tol2/2.0;
//    somechanged = somechanged || changed;
//  }
//  while (tol2 > tollimit);
//  if (changed)
//    printf(" size  of %ld  changes to: %15.12f,      ln L now = %12.8f\n", i,
//logsize, bestlogL);
//} /* resizeform */


void boasfit(void)
{
  /* make a Boas least squares fit of the forms in z to each other for a rough start */
  long i, j, k, num, specimens;
  double meanx, meany, eps, maxrotation, A, B, D, oldx, oldy,
         absrot, suma, targetsize, suma2, lnsize;
  boolean firstcycle;

  nmpoints = morphchars/2;               /* ... as these are (x,y) pairs */
  if (firsttime)                        /* allocate only the first time */
    meanz = (phenotype3)Malloc((long)charsp * sizeof(double));
  for (i = 0; i < spp; i++) {   /* start with zero rotations */
    for (j = 0; j < sample[i]; j++) {
      rotation[i][j] = 0.0;
    }
  }
  for (i = 0; i < spp; i++) {
    for (j = 0; j < sample[i]; j++) {   /* center each individual on (0,0) */
      meanx = 0.0;                   /* obtain the mean for the form */
      for (k = startmchar-1; k <= endmchar-2; k += 2)
        meanx += z[i][j][k];
      meany = 0.0;
      for (k = startmchar; k <= endmchar-1; k += 2)
        meany += z[i][j][k];
      meanx /= nmpoints;
      meany /= nmpoints;
      for (k = startmchar-1; k <= endmchar-2; k += 2)
        z[i][j][k] -= meanx;        /* subtract out the mean for the form */
      for (k = startmchar; k <= endmchar-1; k += 2)
        z[i][j][k] -= meany;
    }
  }
  firstcycle = true;
  eps = 0.0000001;
  specimens = 0;
  targetsize = 0.0;    /* compute initial mean centroid size */
  for (i = 0; i < spp; i++) {
    for (j = 0; j < sample[i]; j++) {
      suma = 0.0;
      for (k = startmchar-1; k <= endmchar-1; k++)
        suma += z[i][j][k]*z[i][j][k];     
      suma = sqrt(suma/(morphchars/2.0));  /* centroid size (2D only) */
      targetsize += suma;   /* centroid size for this specimen */
      specimens++;
    }
  } 
  targetsize /= specimens;  /* mean centroid size of all individuals */
  if (sizes && (!linearsize)) {     /* if sizes being computed */
    for (i = 0; i < spp; i++) {
      for (j = 0; j < sample[i]; j++) {
        suma2 = 0.0;
        for (k = startmchar-1; k <= endmchar-1; k += 1) {
          suma2 += z[i][j][k]*z[i][j][k];
        }
        lnsize = log(sqrt(suma2/(morphchars/2.0))/targetsize); 
        resize(i, j, lnsize);
      }
    }
  }
  do {
    maxrotation = 0.0;
    if (firstcycle) {
      for (j = startmchar-1; j <= endmchar-1; j++)
        meanz[j] = z[0][0][j]; /* fitting all to the first one */
      firstcycle = false;
    } else {    /* after that, each time fit them all to the mean form */
      for (k = startmchar-1; k <= endmchar-1; k++)
        meanz[k] = 0.0;
      num = 0;
      for (i = 0; i < spp; i++)
        for (j = 0; j < sample[i]; j++) {
          num++;
          for (k = startmchar-1; k <= endmchar-1; k++) {
            meanz[k] += z[i][j][k];
          }
        }
      for (k = startmchar-1; k <= endmchar-1; k++)
        meanz[k] /= num;
    }
    for (i = 0; i < spp; i++) {
      for (j = 0; j < sample[i]; j++) {  /* find rotation needed */
        A = 0.0;
        B = 0.0;
        for (k = startmchar-1; k <= endmchar-2; k += 2) {
          A += z[i][j][k]*meanz[k+1] - z[i][j][k+1]*meanz[k];
          B += z[i][j][k]*meanz[k] + z[i][j][k+1]*meanz[k+1];
        }
        D = sqrt(A*A+B*B);
        for (k = startmchar-1; k <= endmchar-2; k += 2) { /* rotate it */
          oldx = z[i][j][k];
          oldy = z[i][j][k+1];
          z[i][j][k] = (B/D)*oldx - (A/D)*oldy;
          z[i][j][k+1] = (A/D)*oldx + (B/D)*oldy;
        }
        absrot = fabs(asin(A/D));   /* keep track of max of rotation */
        rotation[i][j] += asin(A/D);  /* save the rotation */
        if (absrot > maxrotation)
          maxrotation = absrot;
        copyztox(i, j);
      }
    }
  firstcycle = false;
  } while (maxrotation > eps);  /* repeatedly Boas-fit until little change */
} /* boasfit */


// void felsie(void)
// {
//   /* do optimal rotations using likelihood on the tree (which amounts to
//      minimizing the determinant of the estimated covariance matrix,
//      which amounts to minimizing  x'C^{-1}x where C is the covariances
//      of the contrasts of the other species) */
//   /* note -- this function is not being used right now at all but parts of
//      it may be used later */
//   node *storedroot, *wherefrom, *nearestinternalnode;
//   double leftbranchlength, rightbranchlength, maxangle, theta, cos2theta,
//     A1122, A12, alpha11, alpha12, alpha22, costheta, sintheta, oldx, oldy;
//   long i, j, k;
//   boolean wasrootedatinteriornode;
// 
//   // RSGnote: "maxtheta" was originally not initialized and is referenced as such
//   // below.  This is definately a bug.  Initialized here only to silence compiler warning.
//   double maxtheta = 0.0;
// 
//   storedroot = curtree->root; /* store the root of the tree to restore later */
//   (void)storedroot;          // RSGnote: Variable set but never used.
// 
//   wasrootedatinteriornode = (curtree->root->back == NULL);  /* debug: works? */
//   if (!wasrootedatinteriornode)
//   {
//     leftbranchlength = curtree->root->next->v;
//     rightbranchlength = curtree->root->next->next->v;
//     (void)leftbranchlength;             // RSGnote: Variable set but never used.
//     (void)rightbranchlength;            // RSGnote: Variable set but never used.
//   }
//   do {   /* repeated rounds of improvement until angles all small enough */
//     maxangle = 0.0;
//     for (i = 0; i < spp; i++) { /* for each species reroot where connects */
//       wherefrom = curtree->root;
//       nearestinternalnode = wherefrom;
//       if (wherefrom->tip)
//         nearestinternalnode = wherefrom->back;
//       (void)nearestinternalnode;        // RSGnote: Variable set but never used.
//       /* code to remove root, fuse branches, put it near species i here */
//       /* compute covariance matrix for contrasts down to there */
//       /* double check the next statement: is it calling contrasting on the
//          correct nearby node?  Must have inserted root, done this call
//          in a correctly coordinated way    -- debug */
//       makecontrasts(curtree->root->next->back);
//       getcovariances();
//       /* now compute inverse of contrasts, plus get JC^(-1)I etc. */
//       /*   question -- what to do if C is not of full rank? */
//       /* find optimal rotation of that species to minimize det(Cov)*/
//       /* debug  dummy*/
//       alpha11 = 0.0;
//       alpha12 = 0.0;
//       alpha22 = 0.0;
//       /* put here code to compute quadratic forms.  Helmertize first? */
//       A1122 = alpha11-alpha22;
//       A12 = 2.0*alpha12;
//       cos2theta = A1122/sqrt(A1122*A1122+A12*A12);
//       costheta = sqrt((1.0+cos2theta)/2.0);    /* use these to rotate */
//       sintheta = sqrt((1.0-cos2theta)/2.0);
//       theta = acos(costheta);   /* just to test convergence */
//       if (fabs(theta) > maxtheta)
//         maxangle = fabs(theta);
//       for (j = 0; j < sample[i]; j++) {            /* find rotation needed */
//         for (k = startmchar-1; k <= endmchar-2; k += 2) {     /* rotate it */
//           oldx = z[i][j][k];
//           oldy = z[i][j][k+1];
//           z[i][j][k] =   costheta*oldx - sintheta*oldy;
//           z[i][j][k+1] = sintheta*oldx + costheta*oldy;
//         }
//       }
//     }  /* end of loop over species */
//     /* this next is just to be used while we are debugging  -- debug */
//     printf("largest absolute angle changed = %10.8f\n", maxangle);
//   } while (maxangle > 0.00001);
//   /* then finally restore root to its original location */
// } /* felsie */
//

void nmalterspecimen (long i, long j, double* newparameters) {
/* rotate or resize  z  depending on  k  */

    if (i < spp)
      rotate(i, j, newparameters[i]-rotation[i][0]);
    if (i >= spp)
      resize(i, j, newparameters[i-spp]);  /* debug  make sure is OK */
} /* nmalterspecimen */


void findbestandworst (void) {
/* for Nelder-Mead algorithm, find points that are best, worst */
  long i;

  if (loglike[0] > loglike[1]) {  /* initialize best, worst, secondworst */
    nmsecondworstone = 1;
    nmbestone = 0;
  } else {
    nmsecondworstone = 0;
    nmbestone = 1;
  }
  if (loglike[2] > loglike[nmbestone]) {  /* sort-of assumes no ties */
    nmworstone = nmsecondworstone;
    nmsecondworstone = nmbestone;
    nmbestone = 2;
  } else {
    if (loglike[2] > loglike[nmsecondworstone]) {
      nmworstone = nmsecondworstone;
      nmsecondworstone = 2;
    } else {
      nmworstone = 2;
    }
  }
  for (i = 3; i < nmpoints; i++) {  /* update these. Assumes nmpoints>3 */
    if (loglike[i] > loglike[nmbestone]) {
      nmbestone = i;
    } else {
      if (loglike[i] < loglike[nmworstone]) {
        nmsecondworstone = nmworstone;
        nmworstone = i;
      } else {
        if (loglike[i] < loglike[nmsecondworstone]) {
          nmsecondworstone = i;
        }
      }
    }
  }
  bestlogL = loglike[nmbestone];
  worstlogL = loglike[nmworstone];
}; /* findbestandworst */


void getloglikefromparameters (double* newparameters) {
/* gets the log-like for a new suggested set of parameters */
  long i;

  copyallwtoz();
  for (i = 0; i < spp; i++) {  /* update all specimens */
    nmalterspecimen(i, 0, newparameters);
    }
  logL = evaluate(curtree->root);
  copyallwtoz();
} /* getloglikefromparameters */


void initnmiterate() {
/* initialize a full set of points for the Nelder-Mead algorithm */
  long i, j;

/* debug   needs extension to deal with multiple specimens per species */
  if (sizes)
    nparams = 2*spp;
  else 
    nparams = spp;
  nmpoints = nparams + 1;
  parameters = (double **)Malloc((long)nmpoints * sizeof(double *));
  for (i = 0; i < nmpoints; i++)
    parameters[i] = (double *)Malloc((long)nparams * sizeof(double));
  nmmeanparams = (double *)Malloc((long)nparams * sizeof(double));
  reflectedparameters = (double *)Malloc((long)nparams * sizeof(double));
  expandedparameters = (double *)Malloc((long)nparams * sizeof(double));
  contractedparameters = (double *)Malloc((long)nparams * sizeof(double));
  shrunkparameters = (double *)Malloc((long)nparams * sizeof(double));
  loglike = (double *)Malloc((long)nmpoints * sizeof(double));
  copyallztow();  /* set aside the reference forms */
  for (i = 0; i < nmpoints; i++) {  /* set up with initial rotations, sizes */
    for (j = 0; j < nparams; j++) {
      if (j < spp)
        parameters[i][j] = rotation[j][0];
      else
        parameters[i][j] = z[j-spp][0][charspp];
    }
  }
  for (i = 0; i < nmpoints; i++) {    /* all but first point changed some */
    if (i > 0) {
      parameters[i][i-1] += nmupsilon;
      }
    getloglikefromparameters(parameters[i]);
    loglike[i] = logL;
    }
} /* initnmiterate */


void copyparameters(double* newparameters, long where) {
/* for Nelder-Mead: replace a point in parameter space with new values */
  long i;

  for (i = 0; i < nparams; i++)
    parameters[where][i] = newparameters[i];
} /* copyparameters */


void nmiterateonepoint() {
/* do one step of the Nelder-Mead algorithm */
  long i, j;
  double logLreflected;

  findbestandworst();
printf("replacing point %3ld ", nmworstone);  /* debug */
  nmepsilon = 0.0;
  for (i = 0; i < nparams; i++) {
    nmepsilon += fabs(parameters[nmworstone][i]-parameters[nmbestone][i]);
  }
  printf("nmepsilon: %15.12f ", nmepsilon);
  for (j = 0; j < nparams; j++) {        /* get mean of all parameter sets */
    nmmeanparams[j] = 0.0;
    for (i = 0; i < nmpoints; i++) {
      nmmeanparams[j] += parameters[i][j];
    }
    nmmeanparams[j] /= nmpoints;
  }
  for (i = 0; i < nparams; i++) {       /* mean of all but worst */
    nmmeanparams[i] = (nmpoints*nmmeanparams[i] - parameters[nmworstone][i])
                   / (nmpoints - 1);
  }
  for (i = 0; i < nparams; i++) {  /* Step 3 in NM Wikipedia page algorithm */
    reflectedparameters[i] = nmmeanparams[i]
                        + alpha*(nmmeanparams[i] - parameters[nmworstone][i]);
    }
  getloglikefromparameters(reflectedparameters);
  logLreflected = logL;
  if ((logL > loglike[nmsecondworstone]) && (!(logL > loglike[nmbestone]))) {
    copyparameters(reflectedparameters, nmworstone);
printf("reflect point  %3ld", nmworstone);  /* debug */
    loglike[nmworstone] = logL;
  } else { /* step 4 in NM Wikpedia page algorithm */
    getloglikefromparameters(parameters[nmbestone]);
    bestlogL = logL;
    if (logLreflected > bestlogL) {
      for (i = 0; i < nparams; i++) {  /* compute expanded point */
        expandedparameters[i] = nmmeanparams[i]
                           + gammma*(reflectedparameters[i] - nmmeanparams[i]);
      }
      getloglikefromparameters(expandedparameters);
      if (logL > logLreflected) {
        copyparameters(expandedparameters, nmworstone);
printf("expand point   %3ld", nmworstone);  /* debug */
      } else {
        copyparameters(reflectedparameters, nmworstone);
        getloglikefromparameters(reflectedparameters);
printf("reflect point  %3ld", nmworstone);  /* debug */
      }
      loglike[nmworstone] = logL;
    } else {  /* Step 5 in NM Wikipedia page algorithm */
      for (i = 0; i < nparams; i++) {  /* compute contracted point */
        contractedparameters[i] = nmmeanparams[i]
                         + rho*(parameters[nmworstone][i] - nmmeanparams[i]);
      }
      getloglikefromparameters(contractedparameters);
      if (logL > loglike[nmworstone]) {
        copyparameters(contractedparameters, nmworstone);
printf("contract point %3ld", nmworstone);  /* debug */
        loglike[nmworstone] = logL;
      } else { /* Step 6 in NM Wikipedia page algorithm */
printf("shrink all but %3ld", nmbestone);  /* debug */
        for (i = 0; i < nmpoints; i++) {  /* for all points but best one ... */
          if (i != nmbestone) {
            for (j = 0; j < nparams; j++) {
              shrunkparameters[j] = parameters[nmbestone][j]
                        + sigma*(parameters[i][j] - parameters[nmbestone][j]);
            }
            copyparameters(shrunkparameters, i);
            getloglikefromparameters(parameters[i]);
            loglike[i] = logL;
          }
        }
      }
    }
  }
  findbestandworst();
printf(" best (%3ld) is: %15.10f, worst (%3ld) is %15.10f\n", nmbestone, bestlogL, nmworstone, worstlogL); /* debug */
} /* nmiterateonepoint */


void felsmorph() {
 /* do log-likelihood rotation, resizing (if called for) */

   alpha = 1.0;  /* constants used in the Nelder-Mead algorithm */
   gammma = 2.0;
   rho = 0.3;
   sigma = 0.5;
   nmupsilon = 0.1;
   initnmiterate();   /*  initialize, including first set of steps */
   do {
      nmiterateonepoint(); /* make a move in the Nelder-Mead algorithm */
   } while (fabs(nmepsilon) > 0.00000001);
   /* make sure that here there is  x  set up to have best point */
} /* felsmorph */


void felsmorph2() { /* do log-likelihood rotation, resizing (if called for)
                       using simple uphill search one variable at a time */
 double incparam, sintheta, changeofparam, logLnow, oldlogL,
   bestrot, bestsize;
 long i, j, k, kmax;

 for (i = 0; i < spp; i++) {
   for (j = 0; j < sample[i]; j++) {
     copyztox(i,j);
   }
 }
 do {
   oldlogL = bestlogL;
   if (sizes)
     kmax = 1;
   else
     kmax = 0;
   for (k = 0; k <= kmax; k++) {  /* loop angles, k=0, and sizes, k=1 if relevant */
     for (i = 0; i < spp; i++) {   /* do for all specimens except the first */
       for (j = 0; j < sample[i]; j++) {
         if (!(((i == 0) && (j == 0)) || ((i == spp) && (j == 0)))) {
             /* skip over the rotation and resizing of first one */
             copyztotemp(i, j);
             changeofparam = 0.0;
             copyztox(i, j);
             bestlogL = evaluate(curtree->root);
             logLnow = bestlogL;
             bestrot = rotation[i][j];
             bestsize = specsize[i][j][charspp-1];
             incparam = incparam0;
             do {           /* line search of parameter (angle or size) */
               if (sizes || (k == 0)) {  /* rotate or resize, depending */
                 if (k == 0) {   /* change angle */
                   sintheta = sin(incparam);
                   rotate(i, j, sintheta);
                 }
                 if (k == 1) {  /* change size */
                   resize(i, j, incparam);  /* change log-size by amount */
                 }
                 changeofparam += incparam;
                 copyztox(i, j);
                 logLnow = evaluate(curtree->root);
                 if (logLnow > bestlogL) { /* if the log-likelihood is better */
                   bestlogL = logLnow;
                   copyztotemp(i, j);
                   bestrot = rotation[i][j];
                   bestsize = specsize[i][j][charspp-1];
                   incparam = 2.0*incparam;
                 }
                 else  {
                   changeofparam -= incparam;
                    if (k == 0)
                     rotation[i][j] = bestrot;
                   else
                     specsize[i][j][charspp-1] = bestsize;
                   copytemptoz(i, j);
                   incparam = -incparam*0.3;
                 }
               }
               else incparam = 0.0;      /* so loop will terminate right away */
             } while (abs(incparam) > 0.00000001);  /* end loop for k */
             copytemptoz(i, j);
             copyztox(i, j);
           }
         }
       }
if (k == 0)
  printf("after rotations:  best Ln: %15.10f\n", bestlogL);  /* debug */
if ((k == 1) && sizes)
  printf("after resizings:  best Ln: %15.10f\n", bestlogL);  /* debug */
    }
  } while (bestlogL - oldlogL > 0.0000000001);
} /* felsmorph2 */


void felsmorph3() { /* successive uphill line searches using slopes */
 double slopeinc, incparam, sintheta, logLnow, logLnow2, oldlogL, olderlogL; 
   
 long nparams, i, which;

 nparams = spp;
 if (sizes)
   nparams += spp;
 slopes = (double *)Malloc((long)nparams * sizeof(double));
 slopeinc = 0.000001;       /* may need to tune this */
 copyallztow();  /* set aside forms */
 do {
   copyallwtoz();  /* start from same forms */
   copyallztox();
   oldlogL = evaluate(curtree->root);
   olderlogL = oldlogL;
   for (i = 0; i < nparams; i++) {      /* compute approximate slope vector */
     which = 0;
     if (!((i == 0) || (i == spp))) {
       if (i >= spp)
         which = i - spp;
       else
         which = i;
       if (i < spp) {   /* change angle of all but first specimen */
         sintheta = sin(slopeinc);
         rotate(which, 0, sintheta);
       } else { /* change size */
         resize(which, 0, slopeinc);  /* ... or change log-size by */
       }
     }
     copyztox(which, 0);
     logLnow = evaluate(curtree->root);
     if (i > 0) {
       if (i < spp) {   /* change angle of all but first specimen */
         sintheta = sin(-slopeinc);
         rotate(which, 0, sintheta);
       } else { /* change size */
         resize(which, 0, -slopeinc);  /* ... or change log-size by */
       }
     }
     copyztox(which, 0);
     logLnow2 = evaluate(curtree->root);
     slopes[i] = logLnow - logLnow2;
     copywtoz(which, 0);  /* restore original form */
   }
   incparam = 0.0001;     /* may need to tune this */
   copyallztox();
   bestlogL = evaluate(curtree->root);
   logLnow = bestlogL;
   do {           /* line search of parameter (angle or size) */
     copyallwtoz();
     for (i = 0; i < nparams; i++) {   /* do for all parameters */
       if (i >= spp)
         which = i - spp;
       else
         which = i;
       if (i < spp) {   /* change angle */
         sintheta = sin(incparam*slopes[i]);
         rotate(which, 0, sintheta);
       } else { /* change size */
         resize(which, 0, incparam*slopes[i]);  /* change log-size by amount */
       }
     }
     copyallztox();
     logLnow = evaluate(curtree->root);
/* printf("incparam = %15.12f, bestlogL = %20.12f, logL = %20.12f\n", incparam, bestlogL, logLnow);A   debug */
     if (logLnow > bestlogL) { /* if the log-likelihood is better */
       bestlogL = logLnow;
       incparam += incparam;    /* double it */
       copyallztow();     /* update the best form */
     }
     else  {
       incparam = -incparam*0.3;
     }
   } while (fabs(incparam) > 0.00000000001);  /* end loop for line search */
 if (sizes)
   printf("after rotations and resizings:  best Ln: %15.10f, slopeinc = %15.10f\n", bestlogL, slopeinc);  /* debug */
else
  printf("after rotations:  best Ln: %15.10f, slopeinc = %15.10f\n", bestlogL, slopeinc);  /* debug */
 slopeinc *= 1.0;
 } while (bestlogL - olderlogL > 0.0000000001);
  copyallwtoz();
} /* felsmorph3 */


void morph(void)
{
  /* do morphometric transforms if needed, then put result into x */
  long i, j, k;

  for (i = 0; i < spp; i++)      /* copy the original data, y, into array z */
    for (j = 0; j < sample[i]; j++) {
      for (k = 0; k < charsp; k++)
        z[i][j][k] = y[i][j][k];
      if (sizes) {
        z[i][j][charspp] = 0.0;   /* initialize size character  debug: this will be deleted */
      }
      specsize[i][j][charspp] = 0.0;  /* initialize sizes array */
    }
  if (bookmorph || mlrots || justprocrust) {
    if (morphall) {
      startmchar = 1;
      endmchar = chars;
    }
    if (((endmchar-startmchar+1)%2) > 0) {  /* debug:   change for 3D */
      printf("\nERROR:   Number of characters which are morphometric coordinates\n");
      printf("         must be even, but it is %ld.\n\n", endmchar-startmchar+1);
      exxit(-1);
    }
    morphchars = endmchar - startmchar + 1;
    boasfit(); /*  for a rough start, Boas-fit z's to each other */
  }

  if (mlrots) {
    do {
      startlogL = bestlogL;
      incparam0 = 0.0006;
      felsmorph2(); /*  debug */
      felsmorph3(); /*  debug */
      incparam0 = 0.0002;
      felsmorph2(); /*  debug */
      felsmorph3(); /*  debug */
      incparam0 = 0.0001;
      felsmorph2(); /*  debug */
      felsmorph3(); /*  debug */
      endlogL = bestlogL;
    } while (endlogL - startlogL > 0.00001);
  }
  else {
    logL = evaluate(curtree->root);
    bestlogL = logL;
  }
  if (bookmorph)
    bookstein();       /* the Bookstein Transformation */

  for (i = 0; i < spp; i++)      /* copy the possibly transformed data to x */
    for (j = 0; j < sample[i]; j++) {
      copyztox(i, j);
      }
} /* morph */


void contwithin(void)
{
  /* compute all the within-species contrasts, if any */
  long i, j, k;
  double *sumphen;

  sumphen = (double *)Malloc((long)charspp * sizeof(double));
  for (i = 0; i <= spp-1 ; i++) {
    for (j = 0; j < charspp; j++)
      sumphen[j] = 0.0;
    for (k = 0; k <= (sample[i]-1); k++) {
      for (j = 0; j < charspp; j++) {
        if (k > 0)
          cntrast[i][k][j]
            = (sumphen[j] - k*x[i][k][j])/sqrt((double)(k*(k+1)));
        sumphen[j] += x[i][k][j];
        if (k == (sample[i]-1))
          ((cont_node_type*)curtree->nodep[i])->view[j] = sumphen[j]/sample[i];
        x[i][0][j] = sumphen[j]/sample[i];
      }
      if (k == 0)
        curtree->nodep[i]->ssq = 1.0/sample[i]; /* sum of squares for sp. i */
      else
        ssqcont[i][k] = 1.0;   /* if a within contrast */
    }
  }
  free(sumphen);
  contno = 0;
}  /* contwithin */


void contbetween(node *p)
{
  /* compute contrasts and views at a node, multifurcations allowed */
  long j;
  node *q, *r, *pp;
  double v0, v1, vtot, f0, f1;
  boolean atbase;

  atbase = (p == curtree->root);
  pp = p->next;
  q = pp->back;    /* starting with first two descendants ... */
  v0 = multiplier * (q->v + q->deltav);
  if (pp->v < 0.0) {
    printf("\nERROR:  Input tree has a negative branch length,");
    printf(" which is not allowed.\n\n");
    exxit(-1);
  }
  do {
    r = pp->next->back;
    v1 = multiplier * (r->v + r->deltav);
/* debug printf(" %10ld  %10ld  %15.10lf  %15.10lf\n", q->index, r->index, v0, v1);*/  
    if (r->v < 0.0) {
      printf("\nERROR:  Input tree has a negative branch length,");
      printf(" which is not allowed.\n\n");
      exxit(-1);
    }
    vtot = v0 + v1;
    if (vtot > 0.0)
      f0 = v1 / vtot;
    else
      f0 = 0.5;
    f1 = 1.0 - f0;
    if (pp == p->next) {
      p->ssq = q->ssq;
      p->deltav = v1 * f1 / multiplier;
    } else
      p->deltav = p->deltav * f1 / multiplier;
    p->ssq = f0*f0*p->ssq + f1*f1*r->ssq;
    if (pp == p->next) {  /* view using first two descendants */
      for (j = 0; j < charspp; j++) {
        ((cont_node_type*)p)->view[j] = f0 * ((cont_node_type*)q)->view[j]
          + f1 * ((cont_node_type*)r)->view[j];
        cntrast[contno][0][j] = (((cont_node_type*)q)->view[j] -
                                 ((cont_node_type*)r)->view[j])/sqrt(vtot);
/* debug printf(" %10ld  %10ld  %15.12f %15.10f \n", q->index, r->index, vtot, cntrast[contno][0][j]);  */  
      }
      ssqcont[contno][0] = q->ssq + r->ssq;
    } else {              /* use next one to update the view, the delta */
      for (j = 0; j < charspp; j++) {
        ((cont_node_type*)p)->view[j] = f0 * ((cont_node_type*)p)->view[j]
          + f1 * ((cont_node_type*)r)->view[j];
        cntrast[contno][0][j] = (((cont_node_type*)p)->view[j] -
                                 ((cont_node_type*)r)->view[j])/sqrt(vtot);
      }
      ssqcont[contno][0] = p->ssq + r->ssq;    /* debug: q? */
      /* the variance of the contrast */
    }
/* jacobian -= 0.5 * charsd * log((v0+v1));   debug */
   jacobian -= 0.5 * charsd * log(ssqcont[contno][0]);
    contno++;       /* added one more between-species contrast */
    v0 = multiplier * p->deltav;
    pp = pp->next;
  } while (((!atbase) && ((pp->next) != p))
           || (atbase && (((curtree->root->back == NULL)
                           && (pp->next != curtree->root))
                          || ((curtree->root->back != NULL)
                              && (pp != p)))));
  df = charsd;           /* reduce df if too few characters */
  if (charsd > contno)
    df = contno;      /* debug   correct?  */
}  /* contbetween */


void makecontrasts(node *p)
{
  /* compute the contrasts, recursively */
  node *pp;
  boolean atbase;

  if (p->tip) {
    if (sizes && (!linearsize))
      zsum += ((cont_node_type*)p)->view[charspp-1];
    return;
  }
  atbase = (p == curtree->root);
  pp = p->next;
  do {   /* go around the ring making contrasts on descendants */
    makecontrasts(pp->back);
    pp = pp->next;
  } while (((!atbase) && (pp != p))
           || (atbase && (pp != p->next) && (pp->back != NULL)));
  contbetween(p);  /* make contrasts, view at this node */
}  /* makecontrasts */


void getcovariances (void) {
  /* infer the covariances from the contrasts */
  long i, j, k;

  /* debug -- check that the covariances are not done already if fossil case */
// printf("contrast:   %5ld\n", contno);/* debug */  
// printf("how many:   %5ld\n", charspp);/* debug */  
  for (i = 0; i < contno; i++) {
// printf("contrast:   %5ld\n", i);/* debug  */  
    for (j = 0; j < charspp; j++) {
      for (k = 0; k < charspp; k++) {
        if (i == 0)
          sumprod[j][k] = 0.0;
        sumprod[j][k] += cntrast[i][0][j] * cntrast[i][0][k];
/* debug       sumprod[j][k] += cntrast[i][0][j] * cntrast[i][0][k] / ssqcont[i][0]; */
      }
    }
  }
  for (i = 0; i < charspp; i++) {  /* compute covariance from sum or products */
    for (j = 0; j < charspp; j++) {
      sumprod[i][j] /= contno;
// printf(" %20.12lf", sumprod[i][j]);  /* debug */
    }
// printf("\n");  /* debug */
  }
} /* getcovariances */


void getmeans (void) {
  /* infer the mean from the phenotype pruned to the root
   * Note: could be inferred from pruning to other places too,
   * which means these will only be the same in the asymptote
   * of small changes */
  long i;

  for (i = 0; i < charspp; i++)
    mean[i] = ((cont_node_type*)curtree->root)->view[i];
  unorm = 0.0;
  for (i = 0; i < charsp; i++)  /* compute the norm of the mean vector */
    unorm += mean[i]*mean[i];
  unorm = sqrt(unorm);
} /* getmeans */


void writecontrasts(void)
{
  /* write out the contrasts */
  long i, j;
  double xx;

  if (printdata || reg) {  /* debug   ??? correct condition? */
    fprintf(outfile, "\nContrasts (columns are different characters)\n");
    fprintf(outfile, "--------- -------- --- --------- -----------\n\n");
  }
  for (i = 0; i <= contno - 1; i++) {
    for (j = 0; j < charspp; j++) {
      xx = cntrast[i][0][j];
/* debug       xx = cntrast[i][0][j]/sqrt(ssqcont[i][0]);   */
      if (xx < 0.0001)
        fprintf(outfile, " %14.10g ", xx);
      else
        fprintf(outfile, " %12.10f ", xx);
      }
    putc('\n', outfile);
  }
}  /* writecontrasts */


void writesuper(void)
{
  /* write out the superposition with centroids set to the origin */
  long i, j, k, npoints;
  double meanx, meany, weightsumx, weightsumy, relprod, xx;

  if (!omitheaders) {
    fprintf(outfile, "\nSuperposition of");
    if (supshapesonly)
      fprintf(outfile, " shapes of");
    fprintf(outfile, " species (columns are coordinates)\n");
    fprintf(outfile, "------------- --");
    if (supshapesonly)
      fprintf(outfile, " ------ --");
    fprintf(outfile, "  ------ -------- --- ------------\n\n");
  }
  npoints = morphchars / 2;  /* debug  what to do if not morph? */
  if (includesppchars)
 fprintf(outfile, " %ld %ld\n", spp, endmchar-startmchar+1); /* debug: if W? */
  for (i = 0; i < spp; i++) {
    for (j = 0; j < sample[i]; j++) {   /* center each individual on (0,0) */
      meanx = 0.0;                   /* obtain the mean for the form */
      meany = 0.0;
      weightsumx = 0.0;
      weightsumy = 0.0;
      for (k = startmchar-1; k <= endmchar-2; k += 2)
        if (weightsuper) {
          meanx += x[i][j][k]/sumprod[k][k];
          weightsumx += 1.0/sumprod[k][k];
          meany += x[i][j][k+1]/sumprod[k+1][k+1];
          weightsumy += 1.0/sumprod[k+1][k+1];
        } else {
          meanx += x[i][j][k];
          meany += x[i][j][k+1];
        }
      if (weightsuper) {
        meanx /= weightsumx;
        meany /= weightsumy;
      } else {
        meanx /= npoints;
        meany /= npoints;
      }
      for (k = startmchar-1; k <= endmchar-2; k += 2)
        x[i][j][k] -= meanx;        /* subtract out the mean for the form */
      for (k = startmchar; k <= endmchar-1; k += 2)
        x[i][j][k] -= meany;
    }
  }
  if (linearsize && supshapesonly) {
    unorm = 0.0;
    for (i = 0; i < charsp; i++)  /* compute the norm of the mean vector */
      unorm += mean[i]*mean[i];
    unorm = sqrt(unorm);
  }
  for (i = 0; i < spp; i++) {   /* debug: make it write all specimens */
    for(k = 0; k < nmlngth; k++)
      putc(nayme[i][k], outfile);
    relprod = 0.0;
    if (linearsize && supshapesonly )
      for (j = startmchar-1; j <= endmchar-1; j++)
        relprod += x[i][0][j]*mean[j]/(unorm*unorm);
    for (j = startmchar-1; j <= endmchar-1; j++) {
      xx = x[i][0][j];
      if (sizes && !supshapesonly) {
        xx *= exp(x[i][0][endmchar]);
      }
      if (linearsize && supshapesonly)
        xx /= relprod;
      if (abs(xx) < 10000.0)
        fprintf(outfile, " %12.8f", xx);
      else
        fprintf(outfile, " %13.3f", xx);
      }
    putc('\n', outfile);
  }
}  /* writesuper */


void writesizes (void)
{
  /* write out a table of inferred sizes (scales) */
  long i, j;

  fprintf(outfile, "\nInferred scales (sizes)\n");
  fprintf(outfile, "-------- ----- ------\n\n");
  for (i = 0; i < spp; i++)
    for (j = 0; j < sample[i]; j++)
      fprintf(outfile, "%10.6f\n", 1.0 + specsize[i][j][charspp-1]);
} /* writesizes */


void writemeans (void)
{
  /* write out the estimated phenotype means from the root view */
  fprintf(outfile, "\nEstimated means\n");
  fprintf(outfile, "--------- -----\n\n");
  for (i = 0; i < charspp; i++)
    fprintf(outfile, " %11.7f", mean[i]);
  fprintf(outfile, "\n\n");
} /* writemeans */


void getscales(void) {
   /* compute scales (sizes) */
  long i, j, k;

  for (i = 0; i < spp; i++) {  /* debug  ignores within-species case */
    for (j = 0; j < sample[i]; j++) {
      if (linearsize) {
        specsize[i][j][charspp-1] = 0.0;
        for (k = 0; k < charsp; j++)
          specsize[i][j][charspp-1] += x[i][j][j]*mean[k]/(unorm*unorm);
      }
    else 
      specsize[i][0][charspp-1] = exp(specsize[i][0][charspp-1]);
    }
  }
} /* getscales */


void writescales(void) {
   /* report estimated scales (sizes) from ML, Procrustes, or linear case */
  long i, j;

  fprintf(outfile, "Inferred scales (sizes) of the specimens\n");
  fprintf(outfile, "-------- ------ ------- -- --- ---------\n\n");
  for (i = 0; i < spp; i++) { 
    for (j = 0; j < sample[i]; j++) {
      fprintf(outfile, "%12.8f\n", specsize[i][j][charspp-1]);
    }
  }
} /* writescales */


void writevarsz (void)
{
  /* write out variance of inferred scale */
   
  fprintf(outfile, "\nVariance of scale\n");
  fprintf(outfile, "-------- -- -----\n");
  fprintf(outfile, "%12.10f\n", varz);
} /* writevarsz */


void writerotations (void)
{
  /* write out a table of angles by which specimens have been rotated */
  long i, j;

  fprintf(outfile, "\nInferred specimen rotations)\n");
  fprintf(outfile, "-------- -------- ---------\n\n");
  for (i = 0; i < spp; i++)
    for (j = 0; j < sample[i]; j++)
      fprintf(outfile, "%10.6f\n", rotation[i][j]);
} /* writerotations */


void writeallom (void)
{
  /* write out allometry coefficients in log-likelihood and linear cases
      as well as in the Procrustes case where size character is included */
   
  fprintf(outfile, "\nAllometric coefficients: slope above 1 regressed on scale\n");
  fprintf(outfile, "---------- ------------- ----- ----- - --------- -- -----\n"); 
  fprintf(outfile, "\n");
  for (i = 0; i < charsp; i++) {
    allom[i] = sumprod[i][charsp]/varz;
    fprintf(outfile, " %12.6f", allom[i]);
  }
  fprintf(outfile, "\n\n");
} /* writeallom */


void getsizevar(void)
{
  /* compute, in ML-plus-linearsize case, the size variance,
   * plus the mean vector's norm, plus the means
   */
  long i, j;

  getmeans();
  unorm = 0.0;
  for (i = 0; i < charsp; i++)  /* compute the norm of the mean vector */
    unorm += mean[i]*mean[i];
  unorm = sqrt(unorm);
  varz = 0.0;    /*  infer the variance of the size 1+z  */
  for (i = 0; i < charsp; i++)
    for (j = 0; j < charsp; j++)
      varz += mean[i]*sumprod[i][j]*mean[j]/(unorm*unorm*unorm*unorm);
} /* getsizevar */


void getshapecovars (double **sumprods)
{
  /* make character covariances into shape covariances in linear case */
  long i, j, k;
  double sum;

  for (i = 0; i < charsp; i++)
    umean[i] = mean[i]/unorm;  /* unit mean vector */
  for (i = 0; i < charsp; i++) {  /* Cov(x, z):  V mean/unorm^2  */
    sum = 0.0;
    for (j = 0; j < charsp; j++)
      sum += sumprod[i][j]*umean[j]/unorm;
    temp8[i] = sum;
  }
  if (sizechar) {
    sumprods[charspp-1][charspp-1] = varz;  /* put in new part of diagonal */
    for (i = 0; i < charsp; i++) {  /* make that into Cov(y, z) */
      temp8[i] -= mean[i]*varz;
      sumprods[i][charspp-1] = temp8[i]; /* put it in extra row, column */
      sumprods[charspp-1][i] = temp8[i];
    }
  }
  for (i = 0; i < charsp; i++) {
    for (j = 0; j < charsp; j++) {
      temp1[i][j] = 0.0;
      for (k = 0; k < charsp; k++) {
        temp1[i][j] += umean[i]*umean[k]*sumprods[k][j]; /* s s^T V */
      }
    }
  }
  for (i = 0; i < charsp; i++)
    for (j = 0; j < charsp; j++)
      sumprods[i][j] = sumprods[i][j] - temp1[i][j] - temp1[j][i]
                        + varz*unorm*unorm*umean[i]*umean[j];
} /* getshapecovars */


void reportcovars (double **sumprods, long charspp)
{
  /* print out table of covariances among characters */
  long i, j, chars2;

  fprintf(outfile, "\nCovariance matrix\n");
  fprintf(outfile, "---------- ------\n\n");
  if (sizechar)
    chars2 = charspp;
  else
    chars2 = charsp;
  for (i = 0; i < chars2; i++) {
    for (j = 0; j < chars2; j++) {
      if (fabs(sumprod[i][j]) < 0.00001)
        fprintf(outfile, " %10.4g ", sumprods[i][j]);
      else
        fprintf(outfile, " %10.6f ", sumprods[i][j]);
      }
    putc('\n', outfile);
  }
} /* reportcovars */


void writemethods (void) {
/* write out headers indicating what was done */

  if (justprocrust && !sizes)
    fprintf(outfile, "\nBoas transform is used to align specimens\n\n");
  if (justprocrust && sizes)
    fprintf(outfile, "\nProcrustes transform is used to align specimens\n\n");
  if (bookmorph)
    fprintf(outfile, "\nBookstein Transform is used to align specimens\n\n");
  if (mlrots)
    fprintf(outfile, "\nMaximum Likelihood rotations are used to align specimens\n\n");
  if (sizes && linearsize) {
    fprintf(outfile, "\nScale, shape and allometry are inferred by a linear\n");
    fprintf(outfile, "approximate model, after covariances are inferred\n\n");
  }
  if (linearsize && (!sizes)) {
    fprintf(outfile, "\nAllometry, size variance, covariances are inferred by a linear\n");
    fprintf(outfile, "approximate model, after covariances are inferred\n\n");
  }
  if (!linearsize) {
    if (sizes) {
      if (mlsizes) {
        fprintf(outfile, "\nSizes are being inferred by ML:\n variables 1 - %ld are shape variables\n", chars);
        if (sizechar)
          fprintf(outfile, "  and variable %ld is the inferred log(size)\n\n", charspp);
      } else {
        fprintf(outfile,
          "\nSizes are being inferred by Generalized Procrustes:\n");
        fprintf(outfile, "   variables 1 - %ld are shape variables\n", chars);
      if (sizechar)
        fprintf(outfile, "    and variable %ld is the inferred log(size)\n\n", charspp);
      }
    }
  }
  if (pca) {
    fprintf(outfile, "\nPrincipal components are being calculated");
    if ((!linearsize) && sizes && sizechar)
      fprintf(outfile, " based on shape and size\n\n");
    else {
      if (sizes) {
        fprintf(outfile, " based on shape characters only\n\n");
      } else {
        fprintf(outfile, " based on all characters\n\n");
      }
    }
  }
}; /* writemethods */


void covsandshapes (void) {
  /* write out sizes, allometry, and covariances of shape */
  long i, j;
  double tempx;

  if(!superposition) {
    if (sizes) {
      fprintf(outfile, "\nCovariance matrix of shapes\n");
      fprintf(outfile, "---------- ------ -- ------\n\n");
      for (i = 0; i < charspp; i++) {
        for (j = 0; j < charspp; j++) {
          tempx = sumprod[i][j];
          if (fabs(tempx) < 0.00001)
            fprintf(outfile, " %10.4g ", tempx);
          else
            fprintf(outfile, " %10.6f ", tempx);
        }
        putc('\n', outfile);
      }
    }
  }
}; /* covsandshapes */


void calculateregressions (matrix sumprod, long charspp)
{
  /* compute regressions and correlations among contrasts */
  long i, j;

  for (i = 0; i < charspp; i++) {
    for (j = 0; j < charspp; j++) {
      regressions[i][j] = sumprod[i][j] / sumprod[i][i];
    }
  }
  for (i = 0; i < charspp; i++) {
    for (j = 0; j < charspp; j++) {
      correlations[i][j] = sumprod[i][j] / sqrt(sumprod[i][i] * sumprod[j][j]);
    }
  }
}  /* calculateregressions */


void reportregressions (matrix regressions, matrix correlations)
{
  /* compute regressions and correlations among contrasts */
  long i, j;

  if (reg) {
    fprintf(outfile, "\nRegressions (columns on rows)\n");
    fprintf(outfile, "----------- -------- -- -----\n");
    for (i = 0; i < charspp; i++) {          /* print out column numbers */
      if (i == 0)
        fprintf(outfile, "\nrow\\col ");
      fprintf(outfile, "   %3ld    ", i+1);
      if ((i > 0) && (i%10 == 10) && (charspp > i))
        fprintf(outfile, "\n        ");
      if ((i+1) == charspp)
        fprintf(outfile, "\n");
    }
    for (i = 0; i < charspp; i++) {
      fprintf(outfile, " %2ld   ", i+1);
      for (j = 0; j < charspp; j++)
        fprintf(outfile, " %9.4f", regressions[i][j]);
      putc('\n', outfile);
    }
    fprintf(outfile, "\nCorrelations\n");
    fprintf(outfile, "------------\n");
    for (i = 0; i < charspp; i++) {          /* print out column numbers */
      if (i == 0)
        fprintf(outfile, "\nrow\\col ");
      fprintf(outfile, "   %3ld    ", i+1);
      if ((i > 0) && (i%10 == 0) && (charspp > i))
        fprintf(outfile, "\n        ");
      if ((i+1) == charspp)
        fprintf(outfile, "\n");
    }
    for (i = 0; i < charspp; i++) {
      fprintf(outfile, " %2ld   ", i+1);
      for (j = 0; j < charspp; j++)
        fprintf(outfile, " %9.4f", correlations[i][j]);
      putc('\n', outfile);
    }
  }
  if (nocorr) {
    fprintf(outfile, "\nLikelihood Ratio Test\n\n");
    for (i = 0; i < charspp; i++)
      for (j = 0; j < charspp; j++)
        temp1[i][j] = sumprod[i][j] / sqrt(sumprod[i][i]*sumprod[j][j]);
    if (2*spp >= charspp+3) {
      logLvara = -0.5*contno*logdet(temp1);
    } else {
      qreigen(temp1, charspp);
      logLvara = -0.5 * glogdet(eig);  /* generalized logdet */
    }
    for (i = 0; i < charspp; i++)
      for (j = 0; j < charspp; j++) {
        if (nset[i] != nset[j])
          temp1[i][j] = 0.0;
        else
          temp1[i][j] = sumprod[i][j] / sqrt(sumprod[i][i]*sumprod[j][j]);
      }
    logLnocorr = -0.5*contno*logdet(temp1);
    fprintf(outfile,
            "    Log likelihood with correlation       = %13.5f,",
            logLvara);
    fprintf(outfile, "  %ld parameters\n\n", charsd*(charsd+1)/2);
    fprintf(outfile,
            "    Log likelihood without correlation    = %13.5f,",
            logLnocorr);
    fprintf(outfile, "  %ld parameters\n\n",
            n1*(n1+1)/2+(charsd-n1+1)*(charsd-n1)/2);
    fprintf(outfile, "                     difference    = %13.5f\n\n",
            logLvara-logLnocorr);
    fprintf(outfile, "                Chi-square value  = %12.5f,",
            2.0*(logLvara-logLnocorr));
    if (n1*(charsd-n1) == 1)
      fprintf(outfile, "  %ld  degree of freedom\n\n",
              n1*(charsd-n1));
    else
      fprintf(outfile, "  %ld  degrees of freedom\n\n",
              n1*(charsd-n1));
  }
  putc('\n', outfile);
}  /* reportregressions */


double logdet(double **a)
{
  /* Gauss-Jordan log determinant calculation.
   * in place, overwriting previous contents of a.  On exit,
   * matrix a contains the inverse. Works only for positive definite A */
  long i, j, k;
  double temp, sum;

  sum = 0.0;
  for (i = 0; i < charspp; i++) {
    if (fabs(a[i][i]) < 1.0E-37) {
      printf("ERROR:  Tried to invert singular matrix.\n");
      exxit(-1);
    }
    sum += log(a[i][i]);
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < charspp; j++)
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


double glogdet (double* eig)
{
  /* log generalized determinant, from non-nearly-zero postitive
   * eigenvalues */
  long i;
  double x;

  x = 0.0;
  for (i = 0; i < charspp; i++) {
    if (fabs(eig[i]) > 1.0e-14)         /* if this dimension has variation */
      x += log(eig[i]);
  }
  return x;
} /* glogdet */


void invert(double **a)
{
  /* Gauss-Jordan reduction -- invert chars x chars matrix a
     in place, overwriting previous contents of a.  On exit,
     matrix a contains the inverse.*/
  long i, j, k;
  double temp;

  for (i = 0; i < charspp; i++) {
    if (fabs(a[i][i]) < 1.0E-37) {
      printf("ERROR:  Tried to invert singular matrix.\n");
      exxit(-1);
    }
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < charspp; j++)
      a[i][j] *= temp;
    for (j = 0; j < charspp; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < charspp; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
}  /*invert*/


void ginverse (double** a)
{
  /* Moore-Penrose inverse, using spectral decomposition.  On exit  a
   * contains the generalized inverse.  This is obtained by taking the
   * eigenvalues and inverting them, except for ones near zero which
   * are instead zeroed, then reconstituting the matrix.  The function
   * qreigen is used, which leaves column eigenvectors in array  eigvecs
   * and eigenvalues in array  eig.   Note that eig is sorted in
   * decreasing order, and eigorder[i] tells which element of eigvecs
   * corresponds to that eigenvalue.
   */
  boolean *nearlyzero;
  long i, j, k;
  
  nearlyzero = (boolean*)Malloc(charspp*sizeof(boolean));
  qreigen(a, charspp);                /* obtains the spectral decomposition */
  for (i = 0; i < charspp; i++)      /* indicate which eigenvalues */
    if (fabs(eig[i]) < 1.0e-14)      /* are near enough to zero */
      nearlyzero[i] = true;
    else
      nearlyzero[i] = false;
  for (i = 0; i < charspp; i++) {                  /* reconstitute to get */
    for (j = 0; j < charspp; j++) {              /* M-P generalized inverse */
      a[i][j] = 0.0;
      for (k = 0; k < charspp; k++) {
        if (!nearlyzero[k])
          a[i][j] += eigvecs[i][k]*eigvecs[j][k] / eig[k];
      }
    }
  }
} /* ginverse */


void initcovars(boolean novara)
{
  /* Initialize covariance estimates */
  long i, j, k, l, contswithin;
  boolean nonzero;

  /* zero the matrices */
  for (i = 0; i < charspp; i++)
    for (j = 0; j < charspp; j++) {
      vara[i][j] = 0.0;
      vare[i][j] = 0.0;
    }
  /* estimate VE from within contrasts -- unbiasedly */
  contswithin = 0;
  for (i = 0; i < spp; i++) {
    for (j = 1; j < sample[i]; j++) {
      contswithin++;
      for (k = 0; k < charspp; k++)
        for (l = 0; l < charspp; l++)
          vare[k][l] += cntrast[i][j][k]*cntrast[i][j][l];
    }
  }
  /* estimate VA from between contrasts -- biasedly: does not take out VE */
  if (!novara) {   /* leave VarA = 0 if no A component assumed present */
    for (i = 0; i < charspp; i++)
      for (j = 0; j < charspp; j++) {
        nonzero = true;
        if (nocorr)
          nonzero = (nset[i] == nset[j]);
        for (k = 0; k < spp-1; k++) {
          if (nonzero) {
            if (ssqcont[k][0] <= 0.0)
              vara[i][j] += cntrast[k][0][i]*cntrast[k][0][j];
            else
              vara[i][j] += cntrast[k][0][i]*cntrast[k][0][j]
                / ((long)(spp-1)*ssqcont[k][0]);
          }
        }
      }
  }
  for (k = 0; k < charspp; k++)
    for (l = 0; l < charspp; l++)
      if (contswithin > 0)
        vare[k][l] /= contswithin;
      else {
        if (!novara) {
          vara[k][l] = 0.5 * vara[k][l];
          vare[k][l] = vara[k][l];
        }
      }
}  /* initcovars */


double normdiff(boolean novara)
{
  /* Get relative norm of difference between old, new covariances */
  double s;
  long i, j;

  s = 0.0;
  for (i = 0; i < charspp; i++)
    for (j = 0; j < charspp; j++) {
      if (!(novara || (nocorr && (i != j)))) {
        if (fabs(oldvara[i][j]) <= 0.00000001)
          s += vara[i][j];
        else
          s += fabs(vara[i][j]/oldvara[i][j]-1.0);
      }
      if (fabs(oldvare[i][j]) <= 0.00000001)
        s += vare[i][j];
      else
        s += fabs(vare[i][j]/oldvare[i][j]-1.0);
    }
  return s/((double)(charspp*charspp));
}  /* normdiff */


void matcopy(double **a, double **b, long m)
{
  /* Copy matrices m x m: a to b */
  long i;

  for (i = 0; i < m; i++) {
    memcpy(b[i], a[i], charspp * sizeof(double));
  }
}  /* matcopy */


void newcovars(boolean nocorr, boolean novara)
{
  /* one EM update of covariances, compute old likelihood too */
  long i, j, k, l, m;
  double sum, sum2, sum3, sqssq;

  if (!novara)
    matcopy(vara, oldvara, charspp);
  matcopy(vare, oldvare, charspp);
  sum2 = 0.0;         /* log likelihood of old parameters accumulates here */
  for (i = 0; i < charspp; i++)                    /* zero out vara and vare */
    for (j = 0; j < charspp; j++) {
      if ((!novara) || nocorr)
        vara[i][j] = 0.0;
      vare[i][j] = 0.0;
    }
  for (i = 0; i < spp-1; i++) {            /* accumulate over contrasts ... */
    if (i <= spp-2) {      /* E(aa'|x) and E(ee'|x) for "between" contrasts */
      sqssq = sqrt(ssqcont[i][0]);      /* sqrt(d) */
      for (k = 0; k < charspp; k++)     /* compute (dA+E) for this contrast */
        for (l = 0; l < charspp; l++)
          if (!novara)
            temp1[k][l] = ssqcont[i][0] * oldvara[k][l] + oldvare[k][l];
          else
            temp1[k][l] = oldvare[k][l];
      matcopy(temp1, temp2, charspp);
      if (2*spp >= charspp+3) {
        invert(temp2);                               /* compute (dA+E)^(-1) */
      } else {
        ginverse(temp2);                      /* or its generalized inverse */
      }
      matcopy(temp2, temp4, charspp);
      /* sum of - x (dA+E)^(-1) x'/2 for old A, E */
      for (k = 0; k < charspp; k++)
        for (l = 0; l < charspp; l++)
          sum2 -= cntrast[i][0][k]*temp2[k][l]*cntrast[i][0][l]/2.0;
      matcopy(temp1, temp3, charspp);
      if (2*spp >= charspp+3) {
        sum2 -= 0.5 * logdet(temp3);            /* log determinant term too */
      } else {
        qreigen(temp3, charspp);
        sum2 -= 0.5 * glogdet(eig);         /* generalized log  determinant */
      }
      if (!novara) {
        for (k = 0; k < charspp; k++)
          for (l = 0; l < charspp; l++) {
            sum = 0.0;
            for (j = 0; j < charspp; j++)
              sum +=  sqssq * oldvara[k][j] * temp4[j][l];
            Bax[k][l] = sum;            /*  Bax  = sqrt(d) * A *(dA+E)^(-1) */
          }
      }
      for (k = 0; k < charspp; k++)
        for (l = 0; l < charspp; l++) {
          sum = 0.0;
          for (j = 0; j < charspp; j++)
            sum += oldvare[k][j] * temp4[j][l];
          Bex[k][l] = sum;                        /*  Bex = (dA+E)^(-1) * E */
        }
      if (!novara) {
        for (k = 0; k < charspp; k++)
          for (l = 0; l < charspp; l++) {
            sum = 0.0;
            for (m = 0; m < charspp; m++)
              sum += Bax[k][m]
                * (cntrast[i][0][m]*cntrast[i][0][l]
                   - temp1[m][l]);
            temp2[k][l] = sum;                  /*  Bax * (xx'- (dA+E)) ... */
          }
        for (k = 0; k < charspp; k++)
          for (l = 0; l < charspp; l++) {
            sum = 0.0;
            for (m = 0; m < charspp; m++)
              sum += temp2[k][m] * Bax[l][m];
            vara[k][l] += sum;   /*   ... * Bax' */
          }
      }
      for (k = 0; k < charspp; k++)
        for (l = 0; l < charspp; l++) {
          sum = 0.0;
          for (m = 0; m < charspp; m++)
            sum += Bex[k][m] * (cntrast[i][0][m]*cntrast[i][0][l]
                                - temp1[m][l]);
          temp2[k][l] = sum;                     /*  Bex * (xx'-(dA+E)) ... */
        }
      for (k = 0; k < charspp; k++)
        for (l = 0; l < charspp; l++) {
          sum = 0.0;
          for (m = 0; m < charspp; m++)
            sum += temp2[k][m] * Bex[l][m];
          vare[k][l] += sum;                                /*   ... * Bex' */
        }
    }
  }
  matcopy(oldvare, temp2, charspp);
  if (2*spp >= charspp+3) {
    invert(temp2);                                   /* compute (dA+E)^(-1) */
  } else {
    ginverse(temp2);                          /* or its generalized inverse */
  }
  matcopy(oldvare, temp3, charspp);
  if (!novara) {
    if (2*spp >= charspp+3) {
      sum3 = 0.5 * logdet(temp3);                     /* get 1/2 log det(E) */
    } else {
      qreigen(temp3, charspp);
      sum3 = 0.5 * glogdet(eig);           /* generalized log  determinant */
    }
  }
  for (i = 0; i < spp; i++) {
    if (sample[i] > 1) {
      for (j = 1; j < sample[i]; j++) {       /* E(aa'|x) (invisibly) and
                                              E(ee'|x) for within contrasts */
        for (k = 0; k < charspp; k++)
          for (l = 0; l < charspp; l++) {
            vare[k][l] += cntrast[i][j][k] * cntrast[i][j][l] - oldvare[k][l];
            sum2 -= cntrast[i][j][k] * temp2[k][l] * cntrast[i][j][l] / 2.0;
            /* accumulate - x*E^(-1)*x'/2 for old E */
          }
        sum2 -= sum3;                           /* log determinant term too */
      }
    }
  }
  for (i = 0; i < charspp; i++)     /* complete EM by dividing by denom ... */
    for (j = 0; j < charspp; j++) {            /* ... and adding old VA, VE */
      if (!novara) {
        if (nocorr) {
          if (nset[i] != nset[j])
            vara[i][j] = 0.0;
          else {
            vara[i][j] /= (double)contnum;
            vara[i][j] += oldvara[i][j];
          }
        }
        else {
          vara[i][j] /= (double)contnum;
          vara[i][j] += oldvara[i][j];
        }
      }
      vare[i][j] /= (double)contnum;
      vare[i][j] += oldvare[i][j];
    }
  logL = sum2;                             /* log likelihood for old values */
}  /* newcovars */


void printcovariances(boolean nocorr, boolean novara)
{
  /* print out ML covariances and regressions in the error-covariance case */
  long i, j;

  fprintf(outfile, "\n\n");
  if (novara)
    fprintf(outfile, "Estimates when VarA is not in the model\n\n");
  else
    if (nocorr) {
      fprintf(outfile, "Estimates when the two sets of characters are");
      fprintf(outfile, " phylogenetically uncorrelated\n\n");
    }
    else
      fprintf(outfile, "Estimates when VarA is in the model\n\n");
  if (!novara) {
    fprintf(outfile, "Estimate of VarA\n");
    fprintf(outfile, "-------- -- ----\n");
    fprintf(outfile, "\n");
    for (i = 0; i < charspp; i++) {
      for (j = 0; j < charspp; j++)
        if (fabs(vara[i][j]) < 0.0001)
          fprintf(outfile, " %12.3g ", vara[i][j]);
        else
          fprintf(outfile, " %12.6f ", vara[i][j]);
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "Estimate of VarE\n");
  fprintf(outfile, "-------- -- ----\n");
  fprintf(outfile, "\n");
  for (i = 0; i < charspp; i++) {
    for (j = 0; j < charspp; j++)
      if (fabs(vare[i][j]) < 0.0001)
        fprintf(outfile, " %12.3g ", vare[i][j]);
      else
        fprintf(outfile, " %12.6f ", vare[i][j]);
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
  if (!novara) {
    fprintf(outfile, "VarA Regressions (columns on rows)\n");
    fprintf(outfile, "---- ----------- -------- -- -----\n\n");
    for (i = 0; i < charspp; i++) {
      for (j = 0; j < charspp; j++)
        fprintf(outfile, " %9.4f", vara[i][j] / vara[i][i]);
      putc('\n', outfile);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "VarA Correlations\n");
    fprintf(outfile, "---- ------------\n\n");
    for (i = 0; i < charspp; i++) {
      for (j = 0; j < charspp; j++)
        fprintf(outfile, " %9.4f",
                vara[i][j] / sqrt(vara[i][i] * vara[j][j]));
      putc('\n', outfile);
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "VarE Regressions (columns on rows)\n");
  fprintf(outfile, "---- ----------- -------- -- -----\n\n");
  for (i = 0; i < charspp; i++) {
    for (j = 0; j < charspp; j++)
      fprintf(outfile, " %9.4f", vare[i][j] / vare[i][i]);
    putc('\n', outfile);
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "\nVarE Correlations\n");
  fprintf(outfile, "---- ------------\n\n");
  for (i = 0; i < charspp; i++) {
    for (j = 0; j < charspp; j++)
      fprintf(outfile, " %9.4f",
              vare[i][j] / sqrt(vare[i][i] * vare[j][j]));
    putc('\n', outfile);
  }
  fprintf(outfile, "\n\n");
} /* printcovariances */


void emiterate(boolean nocorr, boolean novara)
{
  /* EM iteration of error and phylogenetic covariances */
  /* for future: How to handle missing values? */
  long its;
  double relnorm;

  initcovars(novara);
  its = 1;
  do {
    newcovars(nocorr, novara);
    relnorm = normdiff(novara);
    printf("Iteration no. %ld:  ln L = %10.5f, Norm = %10.5f\n", its, logL, relnorm);
    its++;
  } while ((relnorm > 0.00001) && (its < 100000));
  if (its == 10000) {
    printf("\nWARNING: Iterations did not converge.");
    printf("  Results may be unreliable.\n");
  }
} /* emiterate */


void givens(double **a, long i, long j, long n, double ctheta,
            double stheta, boolean left)
{
  /* Givens transform at i,j for 1..n with angle theta */
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
{
  /* compute cosine and sine of theta */
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
{
  /* Givens tridiagonalization */
  long i, j;
  double s, c;

  for (i = 2; i < n; i++) {
    for (j = i + 1; j <= n; j++) {
      coeffs(a[i - 2][i - 1], a[i - 2][j - 1], &c, &s,accuracy);
      givens(a, i, j, n, c, s, true);
      givens(a, i, j, n, c, s, false);
      givens(eigvecs, i, j, n, c, s, true);
    }
  }
}  /* tridiag */


void shiftqr(double **a, long n, double accuracy)
{
  /* QR eigenvalue-finder */
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
{
  /* QR eigenvector/eigenvalue method for symmetric matrix. Leaves
   * right eigenvectors as columns in eigvecs, and their eigenvalues in prob */
  double accuracy;
  long i, j;

  accuracy = 1.0e-15;
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


void reportmatrix (matrix QQ, long m)
{
  /* print out the axes that are independent */
  long i, j, n;

  fprintf(outfile, "row\\col ");
  for (i = 0; i < m; i++) {          /* print out column numbers */
    if ((i > 0) && (i%6 == 0))
      fprintf(outfile, "\n        ");
    fprintf(outfile, "    %3ld    ", i+1);
    if ((i+1) == m)
      fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
  for (i = 0; i < df; i++) {          /* print out the matrix */
    fprintf(outfile, " %2ld   ", i+1);
    n = eigorder[m-i-1];
    for (j = 0; j < m; j++) {
      if (fabs(QQ[n][j]) < 0.001)
        fprintf(outfile, " %10.4g ", QQ[n][j]);
      else
        fprintf(outfile, " %10.6f ", QQ[n][j]);
      if ((j > 0) && (j%6 == 5) && ((j+1) != charspp)) /* in 6-column blocks */
        fprintf(outfile, "\n       ");
      if ((j+1) == m)
        fprintf(outfile, "\n");
    }
  }
  fprintf(outfile, "\n\n");
} /* reportmatrix */


void reportpca(long m)
{
  /* print out principal component axes for a PCA of  m  characters */
  long i, j, n;
  double eigtotal, corr;

  for (i = 0; i < m; i++)   /* original order of eigenvalues */
    eigorder[i] = i;
  shellsort(eig, eigorder, m);
  fprintf(outfile, "\nPrincipal Component Analysis\n\n");
  eigtotal = 0.0;
  for (i = 0; i < m; i++)   /* sum the eigenvalues */
    eigtotal += eig[i];
  fprintf(outfile,
  "   Principal Components (rows) as linear combinations of characters (columns)\n");
  fprintf(outfile,
  "   --------- ---------- ------ -- ------ ------------ -- ---------- ---------\n\n");
  reportmatrix (eigvecs, m);
  fprintf(outfile, "  For principal     Variance of      Fraction of");
  if (linearsize) 
    fprintf(outfile, "      Correlation\n");
  else
    fprintf(outfile, "\n");
  fprintf(outfile, "   component        its change       variance");
  if (linearsize) 
    fprintf(outfile, "         with size\n");
  else
    fprintf(outfile, "\n");
  fprintf(outfile, "   --------         -----------      -----------");
  if (linearsize) 
    fprintf(outfile, "      -----------\n");
  else
    fprintf(outfile, "\n");
  for (i = 0; i < df; i++) {       /* print out the eigenvalues (variances) */
    n = m-i-1;
    corr = 0.0;
    for (j=0; j < charspp; j++)         /* correlation with the mean vector */
      corr += mean[j]*eigvecs[eigorder[n]][j];
    corr /= unorm;
    if (eig[n] < 0.0001)
      fprintf(outfile, "  %6ld             %10.5g       %8.5f",
              i+1, eig[n], eig[n]/eigtotal);
    else
      fprintf(outfile, "  %6ld           %12.7f       %8.5f",
              i+1, eig[n], eig[n]/eigtotal);
    if (linearsize) 
      fprintf(outfile, "         %8.5f\n", fabs(corr));
    else
      fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
} /* reportpca */


void reportlogL (double lnL, long dff)
{
  /* write out log likelihood and degrees of freedom */

  fprintf(outfile, "\nLog L  =  %12.5f\n", lnL);
  fprintf(outfile, "\ndf  =  %5ld\n\n", dff);
} /* reportlogL */


void contrast_node_copy(node *src, node *dst)
{
  /* make a copy of a node */
  contrast_node *c = (contrast_node *)src;
  contrast_node *d = (contrast_node *)dst;

  cont_node_copy((node*)c, (node*)d);
}  /* contrast_node_copy */


void contrast_node_reinit(node* n)
{
  /* re-init a contrast_node */
  contrast_node *cn = (contrast_node *)n;

  generic_node_reinit(&(cn->cont_node_var.node_var));
} /* contrast_node_reinit */


node* contrast_node_make(tree * treep, node_type type, long index)
{
  /* either create a new contrast_node or obtain it from the garbage */
  node* n;

  n = treep->get_forknode(treep, index);
  contrast_node_init(n, type, index);
  return (node *)n;
} /* contrast_node_make */


void initcontrastnode(tree * treep, node **p, long len, long nodei,
                       long *ntips, long *parens, initops whichinit,
                       pointarray nodep, Char *str,
                       Char *ch, FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit)
  {
    case bottom:
      (*p) = contrast_node_make(treep, FORK_NODE, nodei);
      (*p)->tip = false;
      contrast_node_init(*p, FORK_NODE, nodei);
      nodep[(*p)->index - 1] = (*p);
      break;
    case nonbottom:
      (*p) = contrast_node_make(treep, FORK_NODE, nodei);
      (*p)->tip = false;
      contrast_node_init(*p, FORK_NODE, nodei);
      (*p)->tip = (nodei <= spp);
      break;
    case tip:
      match_names_to_data (str, nodep, p, spp);
      nodei = (*p)->index; /* debug  Why not also necessary in initdnamlnode?  */
      contrast_node_init(*p, TIP_NODE, nodei);
      (*p)->deltav = 0.0;
      (*p)->tip = (nodei <= spp);
      break;
    case length:
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      (*p)->v = valyew / divisor;
      (*p)->iter = false;
      if ((*p)->back != NULL) {
        (*p)->back->v = (*p)->v;
        (*p)->back->iter = false;
      }
      break;
    default:        /* cases of hslength,iter,hsnolength,treewt,unittrwt*/
      break;        /* not handled                                            */
  }
} /* initcontrastnode */


void readthetree (void)
{
  /* read in the tree */
  long nextnode;

  alloctree(&curtree->nodep, nonodes);
  setuptree(curtree, nonodes);
  nextnode = 0;
  goteof = false;
  first = true;
  treeread (curtree, intree, &curtree->root, curtree->nodep, &goteof, &first,
            &nextnode, &haslengths,
            initcontrastnode, false, nonodes);
} /* readthetree */


void makefossilcovars(node *p)
{
  /* use contrasts computed from middle of branch  p  and fossils to
     compute the covariances for that placement of the fossil lineage */
  /* debug -- this will need expansion when within-species variation is also
     allowed (including in the fossils) */
  long i, j, k;

  for (i = 0; i <= contno - 2; i++) {
    for (j = 0; j < charspp; j++) {
      for (k = 0; k < charspp; k++) {
        if (i == 0)
          sumprod[j][k] = 0.0;
        sumprod[j][k] += cntrast[i][0][j] * cntrast[i][0][k] / ssqcont[i][0];
      }
    }
  }
  for (i = 0; i < charspp; i++) /* now compute covariance from sum or products */
    for (j = 0; j < charspp; j++)
      sumprod[i][j] /= contno;
  /* from here add on to products the term for the fossil to connection */
  /* debug:  not there yet */
} /* makefossilcovars */


void bltimetraverse (node *p, double *timedowntohere)
{
  /* from tips downward propagate longest time from tips to node,
     store in node variable  tyme */
  node *pp;
  boolean atbase;

  if (!(p->tip)) {
    atbase = (p == curtree->root);
    pp = p->next;
    do {
      bltimetraverse (pp->back, timedowntohere);
      if ((*timedowntohere) + multiplier * pp->v > p->tyme)
        p->tyme = (*timedowntohere) + multiplier * pp->v;
      *timedowntohere = p->tyme;
      pp = pp->next;
    } while (((!atbase) && (pp != p)) || (atbase && (pp != curtree->root)));
  } else {
    if (isfossil[p->index])
      p->tyme = fossiltime[p->index - 1];
    else
      p->tyme = 0.0;
  }
  (*timedowntohere) = p->tyme;
} /* bltimetraverse */


void updatebounds(node *p)
{
  /* recurse through the tree calculating bounds of multiplier wherever a
     lineage bearing only fossils joins one that bears some nonfossil tips */
  boolean atbase, allfossilsabove, somefossilsabove, allofonedescfossils;
  node *pp;
  double lowestfossil, lowestfossilabovehere;

  atbase = (p == curtree->root);
  if (p->tip)
  {
    p->onlyfossilsabove = isfossil[p->index];
    p->fossilsabove = isfossil[p->index];
    if (isfossil[p->index])
      p->lowestfossilabove = p->tyme;
  }
  else       /* traverse possibly multifurcating tree */
  {
    somefossilsabove = false;
    allfossilsabove = true;
    lowestfossilabovehere = 0.0;
    allofonedescfossils = false;
    pp = p->next;
    do {     /* figure out whether this node ancestral to only fossils
                and also whether some of its descendants are fossils
                and what age the oldest fossil above it is */
      updatebounds(pp->back);
      allofonedescfossils = allofonedescfossils || pp->back->onlyfossilsabove;
      allfossilsabove = allfossilsabove && pp->back->onlyfossilsabove;
      somefossilsabove = somefossilsabove || pp->back->fossilsabove;
      if (pp->back->lowestfossilabove > lowestfossilabovehere)
        lowestfossilabovehere = pp->back->lowestfossilabove;
      pp = pp->next;
    } while (((!atbase) && (pp != p)) || (atbase && (pp != curtree->root)));
    p->onlyfossilsabove = allfossilsabove;
    p->fossilsabove = somefossilsabove;
    if (somefossilsabove  && !allfossilsabove)   /* then need to get scaling */
    {
      lowestfossil = 0.0;
      pp = p->next;     /* start around the ring again */
      do {        /* find the oldest of the ancestral-to-fossils-only nodes */
        if (pp->back->onlyfossilsabove) {
          if (pp->back->tyme > lowestfossil)
            lowestfossil = pp->back->tyme;
        }
        pp = pp->next;
      } while (((!atbase) && (pp != p)) || (atbase && (pp != curtree->root)));
      if (!atbase)
      {
        if (lowestfossil > minmultmult * multiplier * p->tyme)
          minmultmult = lowestfossil / (p->tyme * multiplier);
      }
    }
  }
} /* updatebounds */


double evaluate (node *q) {
  /* evaluate tree and compute log of likelihood, logL */
  double ldet;

  contno = 0;
/* debug -- should we correct for jacobian of size correction here ?? */
  jacobian = 0.0;
    if (sizes && (!linearsize))
      zsum = 0.0;
  contwithin();
  makecontrasts(curtree->root);
  getcovariances();
  if (fossil)
    makefossilcovars(curtree->root);
  matcopy(sumprod, temp5, charspp);
/* debug  if ((chars == df) && !(justprocrust || bookmorph || mlrots))     if full rank */
/* debug     ldet = logdet(temp5);   */           /* can just take Log Det */
/* debug   else { */
    qreigen(temp5, charspp);      /* all eigenvalues including (near-)zeros */
    for (i = 0; i < charspp; i++)  /* set up sort tags for eigenvalue order */
      eigorder[i] = i;
    shellsort(eig, eigorder, charspp);       /* sort eig in ascending order */
// for (i = 0; i < df; i++)   /* set up sort tags for eigenvalue order */
//   printf(" %12.8f", -0.5*contno*log(eig[charspp-i-1])); /*debug */
// printf("\n");  /* debug */
    ldet = 0.0;
  for (i = 0; i < df; i++)    /* log det by dropping missing dimensions */
    for (i = 0; i < df; i++)    /* log det by dropping missing dimensions */
      ldet += log(eig[charspp-i-1]);
/* debug  }   */
  if (sizes && (!linearsize)) {
    logL = -contno*charsd*log(2*pi)/2.0
            -contno*charsd/2.0-0.5*contno*ldet+jacobian - charsd*zsum;
  }
  else {
    logL = -contno*charsd*log(2*pi)/2.0
             -contno*charsd/2.0-0.5*contno*ldet+jacobian;
  }
// printf(" %20.12f\n",-contno*charsd*log(2*pi)/2.0);
// printf(" %20.12f\n",-contno*charsd/2.0);
// printf(" %20.12f\n",-0.5*contno*ldet);
// printf(" %20.12f\n",jacobian);
// printf(" %20.12f\n",-charsd*zsum);
  return (logL);
} /* evaluate */

#if 0          /* for now, comment out all the fossil stuff */
double fevaluate (double twhere, node *q, boolean atroot)
{
  /* calculate log(L) by calling evaluate, for a given placement
     of the fossil(s) and a given scaling of branch length versus time */
  long k;
  node* r = NULL;                       // RSGnote: Initialized to silence warning.  See note below.
  node* rr;

  (void)q;                              // RSGnote: Variable never referenced.

  assert(numfossils > 0);               // RSGnote: Otherwise "r" fails to be initialized.
  for (k = 0; k < numfossils; k++)      /* branch length near fossils */
  {
    r = curtree->nodep[fossilsp[k]-1];
    rr = r->back->next->next;
    r->v = (rr->tyme - r->tyme)/multiplier;
    if (!atroot)
      rr->back->v = r->v;
  }

  logL = evaluate (curtree->root);
  if (firstplace || (logL > bestlogL))
  {
    bestlogL = logL;
    bestplace = r->back->next->back->index; // RSGnote: If "numfossils" == 0, "r" will be uninitialized here.
    bestwhere = twhere/multiplier;
    bestmult = multiplier;
  }

  firstplace = false;
  if (reportplacefossils) {
    if (logL > (bestlogL - howworse))
      fprintf(outfile, " %10.4f %10.4f %10.4f %16.8f\n",
              twhere/multiplier, multiplier, twhere, logL);
  }
  return(logL);
} /* fevaluate */


void locatefossilonbranch (node *p, node *q, node *qq, boolean atroot)
{
  /* try different placements along branch and different values of
     the scaling if branch length versus time */
  long i;
  double twhere, upperlimit;

  upperlimit = p->tyme;  /* p can connect from its own tyme ... */
  if (upperlimit < q->tyme) /* ... or from q's tyme, on back */
    upperlimit = q->tyme;
  for (i = 1; i <= ndiv; i++) {/* try ndiv equally-spaced places in branch */
    twhere = ((2*i-1)*q->back->tyme + (2*ndiv+1-2*i)*upperlimit)/(2.0*ndiv);
    p->v = (twhere - p->tyme)/multiplier;
    p->back->v = p->v;
    p->back->next->v = (twhere - q->tyme)/multiplier;
    p->back->next->back->v = p->back->next->v;
    if (p->back->next->next->back != NULL) {
      p->back->next->next->v =
        (curtree->nodep[p->back->next->next->back->index - 1]->tyme
         - twhere)/multiplier;
      p->back->next->next->back->v = p->back->next->next->v;
    }
    p->back->next->next->tyme = p->tyme + p->v;
    logL = fevaluate(twhere, curtree->root, atroot);
  }
} /* locatefossilonbranch */


boolean placefossilonbranch (node *p, node *q)
{
  /* try different placements of node  p  on branch ancestral to node  q,
     remember best one; if cannot place here, do nothing and return false */
  /* debug: do we want to remember them all?? */
  node *qq, *qtemp, *q1;
  double tup, tnew, timedowntohere;
  boolean canplacehere = false;      /* set to silence compiler warning */

  bltimetraverse (curtree->root, &timedowntohere); /* calculate  tyme  values */
  if (!(q == curtree->root))
  {
    qq = curtree->nodep[q->back->index - 1];   /* qq  is node ancestral to q */
    if (q->tyme > qq->tyme)                   /* switch them if out of time order */
    {
      qtemp = q;
      q = qq;
      qq = qtemp;
    };
    if (qq->tyme < p->tyme) {/* check if can't connect fossil to this branch */
      canplacehere = false;
    } else {
      if (reportplacefossils)
        fprintf(outfile,
                "\n from node %ld (length ago: %6.3f) to node %ld (length ago: %6.3f)\n",
                q->index, curtree->nodep[q->index-1]->tyme/multiplier,
                q->back->index, curtree->nodep[q->back->index-1]->tyme/multiplier);
      generic_tree_insert_(&curtree, p, q, false);
      tup = q->tyme;
      if (tup < p->tyme)
        tup = p->tyme;
      tnew = (tup + qq->tyme)/2.0;
      p->back->tyme = tnew;   /* place new node halfway between nearby ones */
      q->back->tyme = tnew;
      q->back->next->tyme = tnew;
      if (reportplacefossils) {
        fprintf(outfile, "     length  multiplier     time      LogDet\n");
        fprintf(outfile, "      (ago)  (to time)      (ago)           \n");
        fprintf(outfile, "     ------  -----------    -----     ------\n");
        locatefossilonbranch (p, q, qq, false);   /* debug */
      }
      generic_tree_re_move(&curtree, p, &q, false);
      q1 = curtree->nodep[q->back->index-1];
      q->v = (q1->tyme - q->tyme)/multiplier;
      q->back->v = q->v;
      canplacehere = true;  /* debug   temporary */
    }
  }
  else
  {
    generic_tree_insert_(&curtree, p, q, false);
    if (p->tyme > q->tyme)
      p->back->tyme = p->tyme + 0.1;
    else
      p->back->tyme = q->tyme + 0.1;
    p->back->next->tyme = p->back->tyme;  /* set tymes of new root node */
    p->back->next->next->tyme = p->back->tyme;
    curtree->root = p->back->next->next;
    tup = q->tyme;
    if (tup < p->tyme)
      tup = p->tyme;
    if (reportplacefossils) {
      fprintf(outfile, " (at root) minmultmult = %10.6f\n",minmultmult);  /* debug */
      fprintf(outfile, " place species %ld between node %ld (length ago: %6.3f) and -infinity\n",
              p->index, q->index, curtree->nodep[q->index-1]->tyme);
      fprintf(outfile, "     length  multiplier     time\n");
      fprintf(outfile, "      (ago)  (to time)      (ago)     Log(L)\n");
      fprintf(outfile, "     ------  -----------    -----     ------\n");
    }
    locatefossilonbranch (p, q, qq, true);
    maxmultmult = 10.0;   /*  why this particular arbitrary value?? */
    generic_tree_re_move(&curtree, p, &q, false);
  }
  return(canplacehere);
} /* placefossilonbranch */


void placetraverse (long n, node *p)
{
  /* traverse to place the n-th fossil, scaling tree each time.  Be careful not
     to dereference a null pointer at root */
  node *r, *pp;
  boolean atbase, worked;
  double timedowntohere;

  bltimetraverse (curtree->root, &timedowntohere); /* calculate  tymes */
  atbase = (p == curtree->root);
  worked = placefossilonbranch(curtree->nodep[fossilsp[n-1]-1], p);
  if (!(p->tip)) {
    pp = p->next;
    do {
      r = pp->back;
      /* find best place for  n  in branch p being careful if it is the root */
      if (worked)   /* bail out if bottom end of this branch is too far up */
        placetraverse (n, r);
      pp = pp->next;
    } while (((!atbase) && (pp != p))
             || (atbase && (((curtree->root->back == NULL) && (pp != curtree->root))
                            || ((curtree->root->back != NULL) && (pp != curtree->root->next))
                   )));
  }
} /* placetraverse */


void replacefossil(void)
{
  /* dummy for now -- this finally connects the fossil in the best place */
} /* replacefossil */


void makenewbranches (void)
{
  /* for each fossil, make a new interior node and hook it to the fossil */
  long i, emptynode;
  node *p, *q;
  tree* t = NULL;                       // RSGnote: Initialized only to silence compiler warning.

  p = NULL;                             // RSGnote: "p" initialized merely to silence compiler warning.

/* #if 0                                    debug; may not need this  
  curtree->nodep = realloc(curtree->nodep, (nonodes + numfossils) * sizeof(node *  *));
  #endif      debug  */

  *t = curtree;                         // RSGnote: Write from uniititialized pointer.

  emptynode = nonodes - numfossils;
  for (i = 0; i < numfossils; i++) {
    q = curtree->nodep[fossilsp[i]-1];   /* a pointer to that fossil */
    q->tip = true;

/* #if 0    debug    probably don't want this: */
    initcontrastnode(&p, &grbg, NULL, 0.0, emptynode, &spp, 0, bottom, NULL, curtree->nodep, NULL, NULL, NULL);
    for (j = 1; j <= 2; j++)            // complete the triangle of nodes
    {
      initcontrastnode(&p->next, &grbg, NULL, 0.0, emptynode, &spp, 0, nonbottom, NULL, curtree->nodep, NULL, NULL, NULL);
      p = p->next;
    }
    p->next = p;                        // connect last part of triangle
/* debug  #endif  */

    t->get_forknode(t, emptynode);      /* get a fork */
    curtree->nodep[emptynode] = p;       // RSGnote: "p" referenced before being initialized.
    p->index = emptynode+1;             /* set up number of new interior node */
    p->next->index = p->index;          // RSGnote: Using a non-initialized pointer for a memory write.
    p->next->next->index = p->index;
    emptynode++;
    hookup(p, q);                       /* hook up the fossil to the triangle of new interior-type nodes */
  }
} /* makenewbranches */


void placeonefossil(long n)
{
  long j;
  /* traverse over tree placing fossil in various branches,
     Then leave it where it fits best.  */
  updatebounds(curtree->root);
  curtree->nodep[fossilsp[n-1]-1]->tyme = fossiltime[n-1];
  if (reportplacefossils)
    printf("placing fossil number %ld, species %ld, on tree\n",
           n, fossilsp[n-1]);
  firstplace = true;      /* used to set up saving of best log-likelihood */
  minmultmult = 0.005; /* go over all fossils: find bounds of multiplier */
  maxmultmult = 200.0;   /*  why this particular arbitrary value?? */
  multiplier0 = multiplier; /* set aside the current multiplier */
  if (!inferscale) {
    nmult = 1;
    minmultmult = 1.0;
    maxmultmult = 1.0;
  }
  for (j = 1; j < 2*nmult; j = j+2) {  /* try different multipliers */
    multiplier = (minmultmult*(2*nmult-j)+maxmultmult*j)*multiplier0/(2*nmult);
    placetraverse (n, curtree->root);  /* start at root, traverse */
    fprintf(outfile, "\n\nBest location:  branch %ld, at branch length %lf down the branch,\n scaling multiplier  %lf, best time %lf,  log(L) = %lf\n",
             bestplace, bestwhere, bestmult, bestwhere*bestmult, bestlogL);
  }
  replacefossil();  /* debug -- what arguments to use? */
} /* placeonefossil */


void placeallfossils(void) 
{
  /*  debug -- the iterative machinery will be placed here but for now we
      just place the first fossil */
  long i;

  /*  makenewbranches ();    debug   (probably don't need that) */
  printf("\nPlacing fossils on tree:\n");
  for (i = 1; i <= numfossils; i++) { /* create interior nodes for   */
    placeonefossil(i);                 /*   all fossils and attach them */
  }
} /* placeallfossils */
#endif  /* end commenting-out of fossil stuff */


void reportlrttests()
{
  /* write out the results of the LRT tests of no phylogenetic covariances */
/* debug: make sure in other functions to save and pass in the results of the test */

  if (nocorr || nophylo) {
    fprintf(outfile, "\n\n\n    Likelihood Ratio Test");
    if (nocorr)
      fprintf(outfile,  " of no correlation");
    if (nophylo && nocorr)
      fprintf(outfile, " and");
    if (nophylo)
      fprintf(outfile,  " of no VarA component");
    fprintf(outfile, "\n");
    fprintf(outfile, "    ---------- ----- ----");
    if (nocorr)
      fprintf(outfile,  " -- -- -----------");
    if (nophylo && nocorr)
      fprintf(outfile, " ---");
    if (nophylo)
      fprintf(outfile,  " -- -- ---- ---------");
    fprintf(outfile, "\n\n");
    if (nophylo) {
      if (nocorr)
        fprintf(outfile, "    Log likelihood with VarA, correlation = %13.5f,",
                logLvara);
      else
        fprintf(outfile,
                 "    Log likelihood with VarA              = %13.5f,",
                 logLvara);
    } else {
      fprintf(outfile, "    Log likelihood with correlation       = %13.5f,",
              logLvara);
    }
    if (nocorr && nophylo)
      fprintf(outfile, "  %ld parameters\n\n\n", charsd*(charsd+1));
    else
      fprintf(outfile, "  %ld parameters\n\n", charsd*(charsd+1));
    if (nocorr) {
      fprintf(outfile,
              "    Log likelihood without correlation    = %13.5f,",
              logLnocorr);
      fprintf(outfile, "  %ld parameters\n\n",
                 n1*(n1+1)/2+(charsd-n1+1)*(charsd-n1)/2+charsd*(charsd+1)/2);
      fprintf(outfile, "                     difference    = %13.5f\n\n",
              logLvara-logLnocorr);
      fprintf(outfile, "                Chi-square value = %13.5f,",
              2.0*(logLvara-logLnocorr));
      if (n1*(charsd-n1) == 1)
        fprintf(outfile, "  %ld  degree of freedom\n\n",
                n1*(charsd-n1));
      else
        fprintf(outfile, "  %ld  degrees of freedom\n\n",
                n1*(charsd-n1));
      if (nophylo)
        fprintf(outfile, "\n");
    }
    if (nophylo) {
      fprintf(outfile,
              "    Log likelihood without varA           = %13.5f,", logLnovara);
      fprintf(outfile, "  %ld parameters\n\n", charsd*(charsd+1)/2);
      fprintf(outfile, "                     difference    = %13.5f\n\n",
              logLvara-logLnovara);
      fprintf(outfile, "                Chi-square value = %13.5f,",
              2.0*(logLvara-logLnovara));
      if (charsd*(charsd+1)/2 == 1)
        fprintf(outfile, "  %ld  degree of freedom\n\n",
                charsd*(charsd+1)/2);
      else
        fprintf(outfile, "  %ld  degrees of freedom\n\n",
                charsd*(charsd+1)/2);
    }
  }
} /* reportlrttests */


void writereportforonecovar(boolean within, matrix var, phenotype3 meanz)
{
  /* for one of the inferred covariance matrices, write out all the desired
   * reports on it (covariances, correlations, regressions, means, PCs etc
   * before these, write a string describing which covariation */

  if (superposition)
    writesuper();
  if (writecont)
    writecontrasts();
  if (!omitheaders)
    writemeans();
  if (sizes) {
    writereports();
    writerotations();
    writescales();
    writevarsz();
    writeallom();
    }
  reportcovars(sumprod, charspp);
  if (reg) {
    calculateregressions(sumprod, charspp);
    reportregressions(regressions, correlations);
    }
  else {
    printcovariances(nocorr, nophylo);
  }
  if (pca) {
    reportpca(charspp);
    }
  if(reg || pca || (sizes && !superposition)) {
    reportlogL (logL, contno*df);
    }
  if (!omitheaders)
    putc('\n', outfile);
} /* writereportforonecovar */


void writereports(void)
{
  /* for one combination of data set and tree, write out means, covariances
   * (within and between species as needed), principal components,
   * superpositions, etc. Menu settings control which are written out,
   * and can turn off headings so as to write out in computer-readable form */

/* debug: later on may move some of the principal component stuff back here */
    if (!omitheaders) {
      writemethods();
      fprintf(outfile, " Evolutionary covariation (between species)\n\n");
      }
    if (varywithin) {
      writereportforonecovar(false, vara, meanz); /* the between variation */
      writereportforonecovar(true, vare, meanz);  /* the within variation */
      }
    else {
      writereportforonecovar(false, vara, meanz); /* the between variation */
      }
    if (nocorr || nophylo)
      logLvara = logL;
    if (nocorr) {
      logLnocorr = logL;
    }
    if (nophylo) {
      logLnovara = logL;
    }
  if (nocorr || nophylo) {
    reportlrttests();
    }
} /* writereports */


void calculatecovsetc(void)
{
  /* for one combination of data set and tree, infer the covariances and
   * means and call the functions that report on them and other stuff */

  multiplier = 1.0;
  makecontrasts(curtree->root);
  getcovariances();
  getmeans();
  if (!varywithin) {
    if (sizes) {
      getscales();
      if (linearsize)  {
        getsizevar();
        getshapecovars(sumprod);
      }
      else
        varz = sumprod[charspp-1][charspp-1]; 
    } else {
      getsizevar();
      getshapecovars(sumprod);
    }
    if (reg || pca)
      calculateregressions(sumprod, charspp);
    if (pca) {
      matcopy(sumprod, temp5, charspp);
      qreigen (temp5, charspp);
    }
    if(reg || pca || (sizes && !superposition)) {
      logL = evaluate(curtree->root);
      bestlogL = logL;
    }
  }
  else {
    emiterate(false, false);
    if (nocorr || nophylo)
      logLvara = logL;
    if (nocorr) {
      emiterate(nocorr, false);
      logLnocorr = logL;
    }
    if (nophylo) {
      emiterate(nocorr, nophylo);
      logLnovara = logL;
    }
  }
} /* calculatecovsetc */


void maketree(void)
{
  /* process one data set and one tree */

  morph();              /* put y into x after doing morphometric transforms */
  firsttime = false;      /* make sure next data set isn't considered first */
  bifurcating = (curtree->root->next->next == curtree->root);
  contwithin();   /* debug: does this need to be moved? */
#if 0
  if (fossil)
  {
    placeallfossils();
  }
  else {
#endif
    calculatecovsetc();
    writereports();
#if 0
  }
#endif
} /* maketree */


int main(int argc, Char *argv[])
{
  /* main program */
  initdata *funcs;
  long ncases, datasper, treesper;
  boolean datafirst, treesfirst;

#ifdef MAC
  argc = 1;                /* macsetup("Contrast","Contrast");                */
  argv[0] = "Contrast";
#endif
  funcs->tree_new = (tree_new_t)contrast_tree_new;
  funcs->tree_init = (tree_init_t)contrast_tree_init;
  funcs->node_new = (node_new_t)contrast_node_new;
  funcs->node_init = (node_init_t)contrast_node_init;
  phylipinit(argc, argv, funcs, false);
  openfile(&infile, INFILE, "input data", "r", argv[0], infilename);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree", "rb", argv[0], intreename);
  openfile(&outfile, OUTFILE, "output", "w", argv[0], outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  reg = true;
  numtrees = 1;
  doinit();
  contrast_tree_new(&curtree, nonodes, spp, 0);        /* default tree size */
  ncases = numtrees;    /* ncases will be how many tree-data pairs are done */
  if (ndatas > numtrees)
    ncases = ndatas;
  if (cross)
    ncases = ndatas * numtrees;
  datafirst = muldata && !treeswithin;    /* which one always gets advanced */
  treesfirst = multrees && !datawithin;
  treesper = ncases / ndatas;                /* how many trees per data set */
  datasper = ncases / numtrees;              /* how many data sets per tree */
  ith = 0;                                  /* which data set was just done */
  jth = 0;                                      /* which tree was just done */
  if (progress)
    if (muldata || multrees)
      printf("\nProcessed:\n\n");
  for (kth = 1; kth <= ncases; kth++) {
    ithwas = ith;
    jthwas = jth;
    if (datafirst || ((jth%treesper) == 0)) {
      ith++;
      getdata ();                                    /* read a data set ... */
    }
    if (((printdata && (pca || reg || nocorr || nophylo))
         || pca || reg || (!writecont)) && (!superposition) &&
        (ith == 1) && (jth == 0))
    {
      fprintf(outfile,
              "\nContinuous character contrasts analysis, version %s\n\n",
              VERSION);
      fprintf(outfile, "%4ld Species, %4ld Characters\n\n", spp, chars);
    }
    if (treesfirst || ((ithwas%datasper) == 0)) {
      jth++;
      readthetree ();                                    /* read a tree ... */
    }
    if (treesfirst && (ith > ithwas)) {
      if (ndatas > 1) {
        printf(" Data set # %ld\n", ith);
      }
    }
    if ((numtrees > 1) && (jth > jthwas)) {
      if (numtrees > 1) {
        if (!omitheaders) {
          if (treesfirst)
            fprintf(outfile, "\nTree # %ld:\n\n", jth);
          else
            fprintf(outfile, "\nTree # %ld:\n\n", jth);
          }
        printf(" Tree # %ld\n", jth);
      }
    }
    if (!treesfirst) {
      if (ndatas > 1) {
        if (!omitheaders) {
          if (datafirst)
            fprintf(outfile, "\nData set # %ld:\n\n", ith);
          else
            fprintf(outfile, "  \nData set # %ld:\n\n", ith);
        }
        printf("  Data set # %ld\n", ith);
      }
    }
    maketree();                       /* process one data set with one tree */
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
  }                                     /* end loop over trees and data sets */
  for (i = 0; i < charspp; i++)                      /* probably unnecessary */
    free(sumprod[i]);
  free(sumprod);
  FClose(infile);
  FClose(outfile);
  FClose(intree);
  if (progress)
    printf("\nOutput written to file \"%s\".\n\n", outfilename);
  printf("Done.\n\n");
  phyRestoreConsoleAttributes();
  free(funcs);
  return 0;
} /* main */

/* End. */
