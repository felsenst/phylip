/* Version 4.0a.  Copyright 1993-2023.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Mike Palczewski, Doug Buxton, Dan Fineman and Bob Giansiracusa. */

/* Define VERSION string if config.h has not done that already.  This is
 * used throughout the PHYLIP package instead of having version strings
 * that have to be kept up-to-date in other places */

#ifndef PHYLIP_H
#define PHYLIP_H

#ifndef VERSION
#define VERSION "4.0a"
#endif


/* this is only for configure/make compiles, which we do not use these days */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#define true    1                             /* messing with truth itself! */
#define false   0                              /* and now nothing is false! */

/* machine-specific stuff:
   based on a number of factors in the library stdlib.h, we will try
   to determine what kind of machine/compiler this program is being
   built on.  However, it doesn't always succeed.  However, if you have
   ANSI conforming C, it will probably work.

   We will try to figure out machine type
   based on defines in stdio, and compiler-defined things as well.:

   Quite a bit of this is being paranoid about being on antique operating
   systems or antique compilers */

#ifdef WIN32                               /* if we're in Microsoft Windows */
#include <windows.h>

#else                                    /* If not, use null macros instead */
#define NULL_EXPR                       ((void)(0))

#define phySaveConsoleAttributes()      NULL_EXPR
#define phySetConsoleAttributes()       NULL_EXPR
#define phyRestoreConsoleAttributes()   NULL_EXPR
#define phyFillScreenColor()            NULL_EXPR
#define phyClearScreen()                NULL_EXPR

#endif /* WIN32 */

#ifdef  GNUDOS                     /* GNU functions to make it act like DOS */
#define DJGPP
#define DOS
#endif

#ifdef THINK_C                   /* Think C = Lightspeed C for MacOS 8 or 9 */
#define MAC
#endif

#ifdef __MWERKS__             /* the defunct Metrowerks C for Apple PowerPC */
#ifndef WIN32
#define MAC
#endif
#endif

#ifdef __CMS_OPEN  /* default file names for Java Content Management System */
#define CMS
#define EBCDIC true
#define INFILE "infile data"
#define OUTFILE "outfile data"
#define FONTFILE "fontfile data"
#define PLOTFILE "plotfile data"
#define INTREE "intree data"
#define INTREE2 "intree data 2"
#define OUTTREE "outtree data"
#define CATFILE "categories data"
#define WEIGHTFILE "weights data"
#define ANCFILE "ancestors data"
#define MIXFILE "mixture data"
#define FACTFILE "factors data"

#else               /* default file names for everything else */

#define EBCDIC false
#define INFILE "infile"
#define OUTFILE "outfile"
#define FONTFILE "fontfile"     /* on Unix maybe in /usr/local/lib/fontfile */
#define PLOTFILE "plotfile"
#define INTREE "intree"
#define INTREE2 "intree2"
#define OUTTREE "outtree"
#define CATFILE "categories"
#define WEIGHTFILE "weights"
#define ANCFILE "ancestors"
#define MIXFILE "mixture"
#define FACTFILE "factors"
#endif

#ifdef L_ctermid                     /* try and detect for Unix sysV or V7. */
#define SYSTEM_FIVE
#endif

#ifdef sequent                   /* if on a Sequent multiprocessing machine */
#define SYSTEM_FIVE                      /* (which went extinct about 2002) */
#endif

#ifndef SYSTEM_FIVE                       /* diagnosing whether on BSD Unix */
#include <stdlib.h>
#if defined(_STDLIB_H_) || defined(_H_STDLIB) || defined(H_SCCSID) || defined(unix)
#define UNIX
#define MACHINE_TYPE "BSD Unix C"
#endif
#endif

#ifdef __STDIO_LOADED                    /* diagnosing whether on a DEC VAX */
#define VMS                           /* running their VMS operating system */
#define MACHINE_TYPE "VAX/VMS C"                   /* also, now almost gone */
#endif

#ifdef __WATCOMC__               /* diagnosing whether compiler is Watcom C */
#define QUICKC                     /* which was multiplatform, now obsolete */
#define WATCOM
#define DOS
#include "graph.h"
#endif
/* Watcom-C has graphics library calls that are almost identical to
 * Microsoft Quick-C, so the "QUICKC" symbol name stays                 */

#ifdef _QC                      /* is the compiler Microsoft's old Quick C? */
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#include "graph.h"
#define DOS
#endif

#ifdef _DOS_MODE
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS                    /* DOS is always defined if on a DOS machine */
#define MSC                               /* MSC is defined for Microsoft C */
#endif

#ifdef __MSDOS__       /* TURBO C compiler, ONLY (no other DOS C compilers) */
#define DOS
#define TURBOC
#include <stdlib.h>
#include <graphics.h>
#endif

#ifdef DJGPP                        /* DJ Delorie's original Gnu C/C++ port */
#include <graphics.h>
#endif

#ifndef MACHINE_TYPE                 /* if none of the above, assume ANSI C */
#define MACHINE_TYPE "ANSI C"
#endif

#ifdef DOS              /* if running under MSDOS or something imitating it */
#define MALLOCRETURN void
#else
#define MALLOCRETURN void   /* debug:  what ... ?  then why #ifdef ? */
#endif

#ifdef VMS
#define signed                               /* signed doesn't exist in VMS */
#endif

/* default screen types */

/*  if on a DOS but not a Windows system can use IBM PC screen controls */
#ifdef DOS
#ifndef WIN32
#define IBMCRT true
#define ANSICRT false
#endif
#endif
/*  if on a Mac cannot use screen controls */
#ifdef MAC
#define IBMCRT false
#define ANSICRT false
#endif
/*  if on a Windows system can use IBM PC screen controls */
#ifdef WIN32
#define IBMCRT true
#define ANSICRT false
#endif
/* otherwise, let's assume we are on a Linux or Unix system
   with ANSI terminal controls */
#ifndef MAC
#ifndef DOS
#ifndef WIN32
#define IBMCRT false
#define ANSICRT true
#endif
#endif
#endif

#ifdef DJGPP
#undef MALLOCRETURN
#define MALLOCRETURN void
#endif


/* includes: */
#ifdef UNIX
#include <strings.h>
#else
#include <string.h>
#endif

#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "Slist.h"  /* debug: why quotes and not angle-braces here? */

#ifdef MAC
#ifdef DRAW
#include "interface.h"
#else
#include "macface.h"
#endif
#define getch gettch
#endif

/* delimiters in paths to files */
#ifdef MAC
#define DELIMITER ':'                           /* In MacOS, use the colon? */
#else
#ifdef WIN32
#define DELIMITER '\\'                             /* a backslash character */
#else                           /* otherwise the Unix/Linux style delimiter */
#define DELIMITER '/'
#endif
#endif


#define FClose(file) if (file) fclose(file) ; file=NULL
#define Malloc(x) mymalloc((long)x)   /* mymalloc is our wrapper for malloc */

typedef void *Anyptr;

#define Signed    signed
#define Const     const
#define Volatile  volatile
#define Char      char                            /* Characters (not bytes) */
#define Static    static                   /* Private global funcs and vars */
#define Local     static                                /* Nested functions */

#ifndef WIN32
typedef unsigned int boolean;
#endif

/* Number of items per machine word in when stores a binary set of bits.
 * Used in Consense, Treedist, and Clique, because can't use the sign bit */
#define SETBITS 31

/* Static memory parameters */

#define FNMLNGTH        200         /* length of array to store a file name */
#define MAXNCH          20     /* extra-long names for Drawtree or Drawgram */
#define nmlngth         10          /* number of characters in species name */
#define maxcategs       9                   /* maximum number of site types */
#define maxcategs2     11               /* maximum number of site types + 2 */
#define point           "."
#define pointe          '.'
#define down            2
#define MAXNUMTREES    10000000 /* greater than number of user trees can be */
#define MAXSHIMOTREES 100  /* SHT test. (Yes, he uses this as his nickname) */

/* Maximum likelihood parameters */
/* debug: should these go into  bl.h ? */

#define smoothings      8   /* number of passes through smoothing algorithm */
#define iterations      8             /* number of iterates for each branch */
#define epsilon         0.0001             /* small number used in makenewv */
#define EPSILON         0.00001 /* small number used in Hermite rootfinding */
#define initialv        0.1      /* starting branch length unless otherwise */
#define over            60  /* maximum width all branches of tree on screen */
#define LIKE_EPSILON    1e-10   /* Estimate of round-off error in likelihood
                                 * calculations. */

/* Math constants */

#define SQRTPI 1.7724538509055160273       /* square root of Pi, for Normal */
#define SQRT2  1.4142135623730950488                   /* square root of 2. */
#define purset ((1 << (long)A) + (1 << (long)G))        /* the purine bases */
#define pyrset ((1 << (long)C) + (1 << (long)T))         /* the pyrimidines */
/* debug: should these go into something like  dna.h ? */
#define NLRSAVES 5    /* number of views that need to be saved during local
                       * rearrangement                                      */

/* Used in Proml, Promlk, Dnaml, Dnamlk and in the parsimony programs
 * for undefined bestyet value */
#define UNDEFINED -99.99999

/* a basic stack */
typedef struct stack {
  struct stack* next;
  void *data;
} stack;

typedef long *steptr;
typedef long longer[6];
typedef char naym[MAXNCH];                             /* for species names */
typedef long *bitptr;
typedef double raterootarray[maxcategs2][maxcategs2];

typedef struct bestelm {                                    /* stores trees */
  long *btree;
  boolean gloreange;
  boolean locreange;
  boolean collapse;
} bestelm;

FILE *infile, *outfile, *intree, *intree2, *outtree, *workingplot;
FILE *weightfile, *catfile, *ancfile, *mixfile, *factfile;
FILE *progfile;

long spp;                                      /* global: number of species */
long chars;                        /* global: number of characters or sites */
long words, bits;    /* binary words, bit length for binary sets of species */
boolean ibmpc, ansi, tranvsp;       /* screen types, transversion parsimony */
naym *nayme;                                   /* array of names of species */
char progbuf[256];              /* string to display in the progress output */

#define ebcdic EBCDIC                     /* IBM character set pre-ANSI/ISO */

typedef Char plotstring[MAXNCH];

/* Approx. 1GB, used to test for memory request errors */
#define TOO_MUCH_MEMORY 1000000000      /* debug: maybe should make bigger? */


#include <limits.h>

/* minimum double we feel safe with, anything less will be considered
   underflow */
#define MIN_DOUBLE 10e-100

/* K&R says that there should be a plus in front of the number, but no
   machine we've seen actually uses one; we'll include it just in
   case. */
#define MAX_32BITS        2147483647    /* max integer if 32-bit arithmetic */
#define MAX_32BITS_PLUS  +2147483647

/* The below pre-processor commands define the type used to store
   group arrays.  We can't use #elif for Metrowerks C, so we use
   cascaded if statements */
#if INT_MAX == MAX_32BITS                  /* If ints are 4 bytes, use them */
typedef int  group_type;
#else
  #if INT_MAX == MAX_32BITS_PLUS
  typedef int group_type;
  #else                             /* else, if longs are 4 bytes, use them */
    #if LONG_MAX == MAX_32BITS
    typedef long group_type;
    #else
      #if LONG_MAX == MAX_32BITS_PLUS
      typedef long group_type;                          /* Default to longs */
      #else
        typedef long group_type;
      #endif
    #endif
  #endif
#endif


/* for many programs */

#define maxuser        10000        /* maximum number of user-defined trees */
     /* (mostly used to set up user-trees x sites arrays for KHT, SH tests) */

typedef enum {  /* for local vs. not, how much further to go in addtraverse */
  nofurther,
  onestep,
  further
} traversetype;

typedef Char** sequence;                        /* a set of arrays of bases */

typedef enum {                 /* the four DNA bases plus unknown or absent */
  A, C, G, T, O
} bases;

typedef enum {                               /* the amino acids in proteins */
  alanine, arginine, asparagine, aspartic, cysteine,
  glutamine, glutamic, glycine, histidine, isoleucine,
  leucine, lysine, methionine, phenylalanine, proline,
  serine, threonine, tryptophan, tyrosine, valine
} acids;

/* names of discrete character states for Pars */
typedef enum {
  zero = 0, one, two, three, four, five, six, seven
} discbases;

/* used by Protpars and Protdist */
typedef enum {          /* the three-letter amino acid codes, extended */
  ala = 0, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro,
  ser1, ser2, thr, trp, tyr, val, del, stop, asx, glx, ser, unk, quest
} aas;

/* arrays for likelihoods in Dnaml, Dnamlk, Dnadist ... */
typedef double sitelike[(long)T - (long)A + 1];
typedef sitelike *ratelike;
typedef ratelike *phenotype;

/* arrays for likelihoods in Proml, Prokmlk ... */
typedef double psitelike[(long)valine - (long)alanine + 1];
typedef psitelike *pratelike;
typedef pratelike *pphenotype;

/* arrays for likelihoods in Codml */
typedef double csitelike[61];
typedef csitelike *cratelike;
typedef cratelike *cphenotype;

typedef long *baseptr;        /* baseptr used in Dnapars, Dnacomp, Dnapenny */
typedef unsigned char *discbaseptr;     /* discbaseptr used in Pars         */
typedef double *phenotype3;                 /* for continuous char programs */

typedef double *vector;                     /* used in distance programs    */

typedef long nucarray[(long)O - (long)A + 1];
typedef long discnucarray[(long)seven - (long)zero + 1];

typedef enum { bottom, nonbottom, hslength, tip, iter, length,
                 hsnolength, treewt, unittrwt } initops ;


typedef double **transmatrix;
typedef transmatrix *transptr;                   /* transptr used in Restml */

typedef long sitearray[3];
typedef sitearray *seqptr;                       /* seqptr used in Protpars */

/* datastructure typedefs */
enum node_type { FORK_NODE = 0, TIP_NODE, FREE_NODE };
typedef enum node_type node_type;

/* Macros for calling dynamic functions */
/* Might be better as actual functions if performance hit is not severe */
//#define node_init(n,b,l)        (((node*)(n))->init((node*)(n),(b),(l)))
//#define node_free(np)           (((node**)(np))->free((node*)(np)))
//#define node_copy(src,dst)      (((node*)(src))->copy((node*)(src),(node*)(dst)))

/*
 * debug:  TODO: Call them like this eventually:
 *
 * #define node_init(n,b,l)        (((node*)(n))->vtable->node_init_f((node*)(n),(b),(l)))
 * #define node_free(np)           (((node**)(np))->vtable->node_free_f((node*)(np)))
 * #define node_copy(src,dst)      (((node*)(src))->vtable->node_copy_f((node*)(src),(node*)(dst)))
 */

typedef enum nodetype {                                /* what kind of data */
  NODE_T_UNKNOWN,      /* debug:  maybe rename this type "nodedatatype"? */
  NODE_T_GENERIC,
  NODE_T_ML,
  NODE_T_DNA,
  NODE_T_PROT
} nodetype;

typedef struct tree tree;                            /* forward declaration */
typedef struct node node;                            /* forward declaration */

/* prototypes of types of functions */
typedef void (*tree_new_t)(tree**, long, long, long); /* tree_new fn */
typedef void (*tree_init_t)(struct tree*, long, long);      /* tree_init fn */
typedef struct node* (*node_new_t)(node_type, long, long); /* node_new type */
typedef void (*node_init_t)(struct node*, node_type, long);  /* n_init type */
typedef void (*tree_copy_t)(struct tree*, struct tree*);
typedef void (*tree_setupfunctions_t)(struct tree*);   /* sets up functions */
typedef void (*node_reinit_t)(struct node*);
typedef void (*node_free_t)(struct node**);
typedef void (*node_copy_t)(struct node*, struct node*);
typedef void (*fork_print_t)(struct node*);
typedef void (*node_print_t)(struct node*);
typedef void (*do_branchl_on_insert_t)(struct tree*, struct node*, 
                                         struct node*);
typedef void (*do_branchl_on_re_move_t)(struct tree*, struct node*, 
                                          struct node*);
typedef boolean (*fork_good_t)(struct tree*, struct node*);   /* debug: needed for debugging */


struct node {  /* a basic node: space for "everything but the kitchen sink" */
           /* debug: in future could use polymorphism to defer some of these
            * variables to the declarations of subclasses */
  nodetype type;                                         /* Runtime type id */
  struct node *next, *back;           /* points around fork circle, outward */
  long index;                                  /* number of the fork or tip */
  boolean tip;                                /* true if node is a tip node */
  plotstring nayme;                      /* name of the tip, if it is a tip */
  double xcoord, ycoord;                                /* used by printree */
  double oldlen, naymlength;                            /*  "   "     "     */
  long ymin, ymax;                                      /*  "   "     "     */
  boolean haslength;               /* haslength used in dnamlk (and fitch?) */
  double length;                                          /* used in retree */
  boolean iter;                       /* iter used in dnaml, fitch & restml */
  boolean do_newbl;                           /* new branch lengths needed? */
  boolean initialized;              /* initialized used in dnamlk & restml  */
  boolean deleted;                      /* true if node is deleted (retree) */
  boolean hasname;                       /* true if tip has a name (retree) */
  double beyond;       /* in retree: distance beyond it to most distant tip */
  boolean deadend;     /* in retree: true if no undeleted nodes beyond this */
  boolean onebranch;     /* in retree: true if is one undeleted node beyond */
  struct node *onebranchnode;     /* if is, a pointer to that node (Retree) */
  double onebranchlength;/* if is, the distance from here to there (Retree) */
  boolean onebranchhaslength;   /* true if there is a valid combined length */
                                             /* from here to there (Retree) */
  boolean bottom;                   /* used in Dnapars & Dnacomp, disc char */
  boolean visited;                  /* used in Dnapars & Dnacomp, disc char */
  bitptr stateone, statezero;                     /* discrete char programs */
  Char state;                     /* state used in Dnamove, Dolmove & Move  */
  boolean onlyfossilsabove;        /* used in Contrast for fossil machinery */
  boolean fossilsabove;            /* used in Contrast for fossil machinery */
  double lowestfossilabove;        /* used in Contrast for fossil machinery */

/* debug: these function variables should be in vtable? */
  node_copy_t copy;                      /* functions defined for this node */
  node_free_t free;
  node_init_t node_init;    /* debug: use this or one in node_vtable? */
  node_reinit_t reinit;
  fork_print_t fork_print_f;
  node_print_t node_print_f;

  struct node_vtable *vtable;                     /* Pointer to node vtable */  /* debug: what is it? */
};                                /* end of the basic node type declaration */

struct node_vtable {
/* debug: needed here?    node_init_t node_init_f; */
  node_free_t node_free_f;
  node_copy_t node_copy_f;
} vtable;

/* debug:  extern struct node_vtable node_vtable;  */


typedef struct node **pointarray; /* type is an array of pointers to nodes
                                  * and is the type of array nodep */
typedef void (*tree_re_move_t)(struct tree*, struct node*, struct node**, boolean);
typedef boolean (*tree_addtraverse_t)(struct tree*, struct node*, struct node*, 
                           traversetype, struct node*, double*, struct tree*, 
                           boolean, boolean, boolean, double*);
typedef boolean (*tree_addtraverse_1way_t)(struct tree*, struct node*, struct node*, 
                   traversetype, struct node**, double*, struct tree*, boolean, 
                   boolean, boolean*, double*);
typedef void (*tree_insert_t)(struct tree*, struct node*, struct node*, boolean);
typedef boolean (*tree_try_insert_t)(struct tree*, struct node*, struct node*, 
                   struct node*, double*, struct tree*, boolean, 
		   boolean, boolean, double*);
typedef void (*tree_free_t)(struct tree*);
typedef void (*tree_globrearrange_t)(struct tree*, struct tree*, boolean, 
                                       boolean, double*);
typedef void (*tree_locrearrange_t)(struct tree*, struct node*, boolean, double*,
                                    struct tree*,struct tree*, boolean, double*);
typedef void (*tree_smoothall_t)(struct tree*, struct node*);
typedef double (*tree_evaluate_t)(struct tree*, struct node*, boolean);
typedef void (*tree_save_lr_nodes_t)(struct tree*, struct node*, struct node*);
typedef void (*tree_restore_lr_nodes_t)(struct tree*, struct node*, struct node*);
typedef void (*tree_save_traverses_t)(struct tree*, struct node*, struct node*);
typedef void (*tree_restore_traverses_t)(struct tree*, struct node*, struct node*);
typedef void (*tree_release_fork_t)(struct tree*, struct node*);
typedef struct node* (*tree_get_fork_t)(struct tree*,  long);
typedef struct node* (*tree_get_forknode_t)(struct tree*, long);
typedef void (*tree_release_forknode_t)(struct tree*, struct node*);
typedef void (*tree_reinit_forknode_t)(struct tree*, struct node*);
typedef void (*tree_nuview_t)(struct tree*, struct node*);
typedef void (*tree_makenewv_t)(struct tree*, struct node*);
typedef void (*tree_print_t)(struct tree*);

typedef boolean (*tree_good_t)(struct tree*);
typedef boolean (*node_good_t)(struct tree*, struct node*);   // check the individual node

typedef struct tree_vtable tree_vtable;              /* forward declaration */

struct tree_vtable { /* this is a table of tree functions to reassign as
                      * needed in subclasses */
  tree_copy_t copy;
  tree_re_move_t re_move;
  tree_addtraverse_t addtraverse;
  tree_addtraverse_t addtraverse_1way;
  tree_insert_t insert_;
  tree_insert_t tree_insert_;
  tree_try_insert_t tree_try_insert_;
  tree_free_t free;
  tree_globrearrange_t globrearrange;
  tree_smoothall_t smoothall;
  tree_evaluate_t evaluate;
  tree_locrearrange_t locrearrange;
  tree_save_lr_nodes_t save_lr_nodes;
  tree_restore_lr_nodes_t restore_lr_nodes;
  tree_save_traverses_t save_traverses;
  tree_restore_traverses_t restore_traverses;
  tree_release_fork_t release_fork;
  tree_get_fork_t get_fork;
  tree_get_forknode_t get_forknode;
  tree_release_forknode_t release_forknode;
  tree_reinit_forknode_t reinit_forknode;
  tree_nuview_t nuview;
  tree_makenewv_t makenewv;
};

typedef enum {                           /* enum type which is type of tree */
  TREE_T_UNKNOWN,
  TREE_T_GENERIC,
  TREE_T_ROOTED,
  TREE_T_UNROOTED,
  TREE_T_ML
} treetype;


struct tree {                                         /* the tree structure */
  treetype type;                                                /* its type */
  pointarray nodep;    /* the array of pointers to tips and to fork circles */
  double score;                             /* the quantity being maximized */
  struct node *root;  /* rootmost node in rootmost circle, null if unrooted */
  long outgrno;                                /* the index of the outgroup */
  long nonodes;     /* the number of nodes needed for tips and fork circles */
  long spp;                                  /* the number of tip "species" */

  struct node *temp_p, * temp_q;   /* generic temporary nodes, used to save */

  /* for local rearrangement */
  struct node **lrsaves;
  struct node *rb, *rnb, *rnnb;
  boolean mulf;
  boolean onleft;
  boolean do_newbl;

  /* fork management bookeeping stacks */
  Slist_ptr free_forks;   /* debug: I think this is no longer used */
  Slist_ptr free_fork_nodes;

  tree_setupfunctions_t setupfunctions;     /* sets up functions */
  tree_copy_t copy;
  tree_re_move_t re_move;
  tree_addtraverse_t addtraverse;
  tree_addtraverse_1way_t addtraverse_1way;
  tree_insert_t insert_;
  tree_try_insert_t try_insert_;
  tree_globrearrange_t globrearrange;
  tree_smoothall_t smoothall;
  tree_evaluate_t evaluate;
  tree_locrearrange_t locrearrange;
  tree_save_lr_nodes_t save_lr_nodes;
  tree_restore_lr_nodes_t restore_lr_nodes;
  tree_save_traverses_t save_traverses;
  tree_restore_traverses_t restore_traverses;
  tree_free_t free;
  tree_release_fork_t release_fork;
  tree_get_fork_t get_fork;
  tree_get_forknode_t get_forknode;
  tree_release_forknode_t release_forknode;
  tree_reinit_forknode_t reinit_forknode;
  tree_nuview_t nuview;
  tree_makenewv_t makenewv;
  tree_print_t tree_print_f;
  do_branchl_on_insert_t do_branchl_on_insert_f;
  do_branchl_on_re_move_t do_branchl_on_re_move_f;

  tree_good_t   tree_good_f;    /* debug: what are these for? */
  node_good_t node_good_f;
  fork_good_t   fork_good_f;

  tree_vtable *vtable;     /* debug:  is this needed?  used? */
};

typedef void (*initptr)(struct tree *, struct node **, long, long,
                         long *, long *, initops, pointarray,
                         Char *, Char *, FILE *);

/* some pointers to functions needed in the node and tree class hierarchies */

typedef struct initdata {
  tree_new_t tree_new;                                  /* makes a new tree */
  tree_init_t tree_init;                     /* initiates stuff in the tree */
  node_new_t node_new;                      /* makes a new node in the tree */
  node_init_t node_init;                     /* initiates stuff in the node */
} initdata;

initdata funcs;    /* declaration of the  funcs  function pointer structure */

boolean javarun;               /* boolean for when Java front-end is in use */

#ifndef OLDC /* need if not the old original Kernighan & Rtichie C compiler */
/* function prototypes */
void            generic_tree_new(struct tree**, long, long, long);
void            generic_tree_init(struct tree*, long, long);
struct node*    generic_node_new(node_type, long, long);
void            generic_node_init(struct node*, node_type, long);
void            no_op(void);
void            phylipinit(int, char**, initdata*, boolean);
struct node*    where_in_dest (struct tree*, struct tree*, struct node*);
void            generic_tree_copy(struct tree*, struct tree*);
void            generic_node_copy(struct node*, struct node*);
void            generic_fork_print(struct node*);
void            generic_node_print(struct node*);
void            generic_node_free(struct node**);
void            generic_node_reinit(struct node*);
void            setupnode(struct node*, long);
long            count_sibs(struct node*);
boolean         isemptyroot(struct node*);
struct node*    findroot(struct node*, boolean*);
struct node*    findrootmostandroot(struct tree*, struct node*, boolean*);
void            generic_insertroot(struct tree*, struct node*, struct node*);
void            generic_root_insert(struct tree*, struct node*);
void            generic_tree_re_move(struct tree*, struct node*,
                                      struct node**, boolean);
void            put_root_near_outgroup(struct tree*, long, boolean);
void            rooted_tree_insert_(struct tree*, struct node*, 
                                     struct node*, boolean);
void            generic_do_branchl_on_re_move(struct tree*, struct node*, 
                                               struct node*);
void            generic_tree_release_forknode(struct tree*, struct node*);
boolean         generic_tree_try_insert_(struct tree*, struct node*, 
                                  struct node*, struct node*, double*, 
                                  struct tree*, boolean, boolean, boolean, 
                                  double*);
void            buildsimpletree(struct tree*, long*);
struct node*    generic_newrootfork(struct tree*);
void            rooted_tree_re_move(struct tree*, struct node*, 
                                     struct node**, boolean);
void            hsbut(struct tree*, struct tree*, struct tree*, boolean, 
                        boolean, long, longer, boolean, double*);
void            preparetree(struct tree*);  /* debug: need this here? */
void            fixtree(struct tree*);
void            arbitrary_resolve(struct tree*) ;
void            writename(long, long, long*);
void            print_progress(char*);
void 		seetree(struct tree*);
void 		dumpnodelinks(struct node *, pointarray, long nonodes);

/* if following not in phylip.c. best to demote them downwards unless shared
   by two branches of hierarchy that split below this */

void            verify_nuview(struct node*);
void            invalidate_nuview(struct node*);
void            invalidate_traverse(struct node*);
void            inittrav_all(struct tree*);
void            initializetrav (struct tree*, struct node*);
void            inittrav (struct tree*, struct node*);
void            EOF_error(void);
void            crash_handler(int);
void            scan_eoln(FILE*);
boolean         eoff(FILE*);
boolean         eoln(FILE*);
boolean         filexists(const char*);
void            openfile(FILE**, const char*, const char*, const char*,
                          const char*, char*);
const char*     get_command_name (const char*);
void		_fgetline_finalize(void);
char*		fgetline(FILE*);
char            menu_getchar(void);
void            getstryng(char*);
void            countup(long*, long);
void            cleerhome(void);
long            readlong(const char*);
void            uppercase(Char*);
double          randum(longer);
void            randumize(longer, long*);
double          normrand(longer);
void            initseed(long*, long*, longer);
void            initjumble(long*, long*, longer, long*);
void            initoutgroup(long*, long);
void            initthreshold(double*);
void            initcatn(long*);
void            initcategs(long, double*);
void            initprobcat(long, double*, double*);
void            lgr(long, double, raterootarray);
double          logfac (long);
double          glaguerre(long, double, double);
void            initlaguerrecat(long, double, double*, double*);
double          hermite(long, double);
void            root_hermite(long, double*);
double          halfroot(double (*func)(long , double), long, double, double);
void            hermite_weight(long, double*, double*);
void            inithermitcat(long, double, double*, double*);
void            initgammacat(long, double, double*, double*);
void            inithowmany(long*, long);
void            inithowoften(long*);
void            initlambda(double*);
void            initfreqs(double*, double*, double*, double*);
void            initratio(double*);
void            initpower(double*);
void            initdatasets(long*);
void            justweights(long*);
void            initterminal(boolean*, boolean*);
void            initnumlines(long*);
void            newline(FILE*, long, long, long);
void            recursiveTreeRead( Char*, long*, FILE*, boolean*, boolean*,
                                   long*, long*, boolean*, boolean);
void            inputNumbersFromTreeFile(FILE*, long* spp, long*);
void            inputnumbers(long*, long*, long*, long);
void            inputnumbers2(long*, long*, long);
void            samenumsp(long*, long);
void            samenumsp2(long);
void            readoptions(long*, const char*);
void            matchoptions(Char*, const char*);
void            headings(long, const char*, const char*);
void            initname(long);
void            checknames(long int);
void            inputweights(long, steptr, boolean*);
void            inputweights2(long, long, long*, steptr, boolean*, const char*);
void            printweights(FILE*, long, long, steptr, const char*);
void            inputcategs(long, long, steptr, long, const char*);
void            printcategs(FILE*, long, steptr, const char*);
void            inputfactors(long, Char*, boolean*);
void            printfactors(FILE*, long, Char*, const char*);
void            findtree(boolean*, long*, long, long*, bestelm*);
void            addtree(long, long*, boolean, long*, bestelm*);
long            findunrearranged(bestelm*, long, boolean);
void            shellsort(double*, long*, long);
void            getch(Char*, long*, FILE*);
void            findch(Char, Char*, long);
void            processlength(double*,double*, Char*, boolean*, FILE*, long*);
void            commentskipper(FILE*, long*);
long            countcomma(FILE*, long*);
long            countsemic(FILE*);
void            memerror(void);
void            odd_malloc(long);
MALLOCRETURN    *mymalloc(long);

void            hookup(struct node*, struct node*);
struct node*    precursor(struct node*);
void            link_trees(long, long, long, pointarray);
void            allocate_nodep(pointarray*, FILE*, long*);
long            take_name_from_tree (Char*, Char*, FILE*);
void            match_names_to_data (Char*, pointarray, struct node**, long);
void            addelement(struct tree*, struct node**, struct node*, Char*, 
                            long*, FILE*, pointarray, boolean*, boolean*, 
                            long*, long*, boolean*, initops, boolean, long);
void            treeread (struct tree*, FILE*, struct node**, pointarray, 
                            boolean*, boolean*, long*, boolean*, initops, 
                            boolean, long);
void            addelement2(struct tree*, struct node*, Char*, long*, FILE*, 
                             boolean, double*, boolean*, long*, long*, long,
                             boolean*, boolean, long);
void            treeread2 (struct tree*, FILE*, struct node**, boolean, 
                            double*, boolean*, boolean*, 
                            long*, boolean, long);
void            exxit (int);
char            gettc(FILE*);
void            unroot(struct tree*, long);
void            unroot_here(struct tree*, struct node*, long);
void            unroot_r(struct tree*, struct node*, long);
void            release_all_forks(struct tree*);
void            destruct_tree(struct tree*);
void            rooted_tree_init(struct tree*, long, long);
void		generic_tree_setupfunctions(struct tree*);
void            generic_tree_free(struct tree*);
void            generic_tree_print(struct tree*);
boolean         generic_tree_good(struct tree*);
boolean         generic_fork_good(struct tree*, struct node*);
boolean         generic_node_good(struct tree*, struct node*);
boolean         oktoinsertthere(struct tree*, struct node*);
boolean         oktorearrangethere(struct tree*, struct node*);
void            rooted_globrearrange(struct tree*, struct tree*, boolean, 
                                       boolean, double*);
void            generic_globrearrange(struct tree*, struct tree*, boolean, 
                                       boolean, double*);
boolean         oktoputthere(struct tree*, struct node*);
boolean         generic_tree_addtraverse(struct tree*, struct node*, 
                                   struct node*, traversetype, struct node*, 
                                   double*, struct tree*, boolean, boolean, 
                                   boolean, double*);
boolean         generic_tree_addtraverse_1way(struct tree*, struct node*, 
                                   struct node*, traversetype, struct node*, 
                                   double*, struct tree*, boolean, 
                                   boolean, boolean*, double*);
#ifdef WIN32              /* if using screen attributes of a Windows system */
void 		phySaveConsoleAttributes(void);
void 		phySetConsoleAttributes(void);
void 		phyRestoreConsoleAttributes(void);
void 		phyFillScreenColor(void);
void 		phyClearScreen(void);
#endif

void            unrooted_tree_save_lr_nodes(struct tree*, 
                                             struct node*, struct node*);
void            unrooted_tree_restore_lr_nodes(struct tree*, 
                                                struct node*, struct node*);
void            generic_unrooted_locrearrange(struct tree*, struct node*, 
                                     boolean, double*, struct tree*, 
                                     struct tree*, boolean, double*);
boolean		unrooted_tree_locrearrange_recurs(struct tree*, struct node*, 
                                              double*, boolean, struct tree*, 
                                              struct tree*, boolean, double*);
void            generic_tree_save_traverses(struct tree*, 
                                             struct node*, struct node*);
void            generic_tree_restore_traverses(struct tree*, 
                                                  struct node*, struct node*);
void    	rooted_tryrearr(struct tree*, struct node*, boolean*);
void		rooted_repreorder(struct tree*, struct node*, boolean*);
void            rooted_locrearrange(struct tree*, struct node*, boolean, 
                                      double*, struct tree*, struct tree*, 
                                      boolean, double*);
void            generic_tree_save_lr_nodes(struct tree*, struct node*, 
                                             struct node*);
void            rooted_tree_restore_lr_nodes(struct tree*, 
                                                  struct node*, struct node*);
void*		pop(struct stack**);
struct stack* 	push(struct stack*,void*);
struct node*    generic_tree_get_fork(struct tree*, long);
void            generic_tree_release_fork(struct tree*, struct node*);
long		generic_tree_findemptyfork(struct tree*);
void            generic_tree_nuview(struct tree*, struct node*);
double          generic_tree_evaluate(struct tree*, struct node*, boolean);
void            generic_tree_insert_(struct tree*, struct node*, 
                                      struct node*, boolean);
void            generic_do_branchl_on_insert(struct tree*, 
                                                  struct node*, struct node*);
struct node*    generic_tree_get_forknode(struct tree*, long);
void            generic_tree_re_move(struct tree*, struct node*, 
                                       struct node**, boolean);
void            generic_re_move(struct tree*, struct node*, 
                                  struct node*, boolean);
void            allocdiscnontip(struct node*, long );
void            allocnode(struct node**, long);
void            allocdiscnode(struct node**, long);
void            gnudisctreenode(struct node**, struct node**, long, long);
void            generic_tree_restore_lr_nodes(struct tree*, struct node*, 
                                                                struct node*);
void            rooted_tree_save_lr_nodes(struct tree*, struct node*, 
                                            struct node*);
void            generic_tree_reinit_forknode(struct tree*, struct node*);
void            generic_initialvtrav(struct node*);
void            generic_treevaluate(struct tree*, boolean, boolean, boolean);
#endif /* OLDC */

#endif
/* end commenting out of whole header because it's been used before */

/* End. */
