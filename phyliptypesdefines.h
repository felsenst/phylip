/* Phylip typedefs and defines  (c) 2025 */
/* v4.0a by Joe Felsenstein, material from  phylip.h */

#define MAXNCH          20     /* extra-long names for Drawtree or Drawgram */

#define true    1                             /* messing with truth itself! */
#define false   0                              /* and now nothing is false! */

#define VERSION "4.0a"
#define true    1                             /* messing with truth itself! */
#define false   0                              /* and now nothing is false! */
#define NULL_EXPR                       ((void)(0))
#define phySaveConsoleAttributes()      NULL_EXPR
#define phySetConsoleAttributes()       NULL_EXPR
#define phyRestoreConsoleAttributes()   NULL_EXPR
#define phyFillScreenColor()            NULL_EXPR
#define phyClearScreen()                NULL_EXPR
#define DJGPP
#define DOS
#define MAC
#define MAC
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
#define SYSTEM_FIVE
#define SYSTEM_FIVE                      /* (which went extinct about 2002) */
#define UNIX
#define MACHINE_TYPE "BSD Unix C"
#define VMS                           /* running their VMS operating system */
#define MACHINE_TYPE "VAX/VMS C"                   /* also, now almost gone */
#define QUICKC                     /* which was multiplatform, now obsolete */
#define WATCOM
#define DOS
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#define DOS
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS                    /* DOS is always defined if on a DOS machine */
#define MSC                               /* MSC is defined for Microsoft C */
#define DOS
#define TURBOC
#define MACHINE_TYPE "ANSI C"
#define MALLOCRETURN void
#define MALLOCRETURN void   /* debug:  what ... ?  then why #ifdef ? */
#define signed                               /* signed doesn't exist in VMS */
#define IBMCRT true
#define ANSICRT false
#define IBMCRT false
#define ANSICRT false
#define IBMCRT true
#define ANSICRT false
#define IBMCRT false
#define ANSICRT true
#define MALLOCRETURN void
#define getch gettch
#define DELIMITER ':'                           /* In MacOS, use the colon? */
#define DELIMITER '\\'                             /* a backslash character */
#define DELIMITER '/'
#define FClose(file) if (file) fclose(file) ; file=NULL
#define Malloc(x) mymalloc((long)x)   /* mymalloc is our wrapper for malloc */
#define Signed    signed
#define Const     const
#define Volatile  volatile
#define Char      char                            /* Characters (not bytes) */
#define Static    static                   /* Private global funcs and vars */
#define Local     static                                /* Nested functions */
#define SETBITS 31
#define FNMLNGTH        200         /* length of array to store a file name */
#define nmlngth         10          /* number of characters in species name */
#define maxcategs       9                   /* maximum number of site types */
#define maxcategs2     11               /* maximum number of site types + 2 */
#define point           "."
#define pointe          '.'
#define down            2
#define MAXNUMTREES    10000000 /* greater than number of user trees can be */
#define MAXSHIMOTREES 100  /* SHT test. (Yes, he uses this as his nickname) */
#define smoothings      8   /* number of passes through smoothing algorithm */
#define iterations      8             /* number of iterates for each branch */
#define epsilon         0.0001             /* small number used in makenewv */
#define EPSILON         0.00001 /* small number used in Hermite rootfinding */
#define initialv        0.1      /* starting branch length unless otherwise */
#define over            60  /* maximum width all branches of tree on screen */
#define LIKE_EPSILON    1e-10   /* Estimate of round-off error in likelihood
#define SQRTPI 1.7724538509055160273       /* square root of Pi, for Normal */
#define SQRT2  1.4142135623730950488                   /* square root of 2. */
#define purset ((1 << (long)A) + (1 << (long)G))        /* the purine bases */
#define pyrset ((1 << (long)C) + (1 << (long)T))         /* the pyrimidines */
#define NLRSAVES 5    /* number of views that need to be saved during local
#define UNDEFINED -999999.99999
#define ebcdic EBCDIC                     /* IBM character set pre-ANSI/ISO */
#define TOO_MUCH_MEMORY 1000000000      /* debug: maybe should make bigger? */
#define MIN_DOUBLE 10e-100
#define MAX_32BITS        2147483647    /* max integer if 32-bit arithmetic */
#define MAX_32BITS_PLUS  +2147483647
#define maxuser        10000        /* maximum number of user-defined trees */
//#define node_init(n,b,l)        (((node*)(n))->init((node*)(n),(b),(l)))
//#define node_free(np)           (((node**)(np))->free((node*)(np)))
//#define node_copy(src,dst)      (((node*)(src))->copy((node*)(src),(node*)(dst)))
 * #define node_init(n,b,l)        (((node*)(n))->vtable->node_init_f((node*)(n),(b),(l)))
 * #define node_free(np)           (((node**)(np))->vtable->node_free_f((node*)(np)))
 * #define node_copy(src,dst)      (((node*)(src))->vtable->node_copy_f((node*)(src),(node*)(dst)))

typedef void *Anyptr;

#ifndef WIN32
typedef unsigned int boolean;
#endif

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

typedef enum { bottom, nonbottom, hslength, tip, iter, length,
                 hsnolength, treewt, unittrwt } initops ;

typedef struct node** pointarray; /* type is an array of pointers to nodes
                                  * and is the type of array nodep */

typedef struct tree tree;                            /* forward declaration */

typedef void (*initptr)(struct tree *, struct node **, long, long,
                         long *, long *, initops, pointarray,
                         Char *, Char *, FILE *);

long spp;                               /* global: number of species */
long chars;                 /* global: number of characters or sites */
long words, bits;    /* binary words, bit length for sets of species */
boolean ibmpc, ansi, tranvsp;     /* screens, transversion parsimony */
naym *nayme;                            /* array of names of species */
char progbuf[256];       /* string to display in the progress output */

typedef Char plotstring[MAXNCH];

/* The below pre-processor commands define the type used to store
   group arrays.  We can't use #elif for Metrowerks C, so we use
   cascaded if statements */
#if INT_MAX == MAX_32BITS                  /* If ints are 4 bytes, use them */
typedef int group_type;
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


typedef enum {  /* for local vs. not, how much further to go in addtraverse */
  nofurther,                                 /* have gone as far as we need */
  onestep,             /* go just one step further for local rearrangements */
  further                                   /* keep going as far as you can */
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

typedef enum {               /* names of discrete character states for Pars */
  zero = 0, one, two, three, four, five, six, seven
} discbases;

/* used by Protpars, Protdist, Proml and Promlk */
typedef enum {          /* the three-letter amino acid codes, extended */
  ala = 0, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, 
  pro, ser1, ser2, thr, trp, tyr, val, del, stop, asx, glx, ser, unk, quest
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
typedef unsigned char *discbaseptr;             /* discbaseptr used in Pars */
typedef double *phenotype3;                 /* for continuous char programs */

typedef double *vector;                     /* used in distance programs    */

typedef long nucarray[(long)O - (long)A + 1];
typedef long discnucarray[(long)seven - (long)zero + 1];

typedef double **transmatrix;     /* transition matrix for Restml, Restdist */
typedef transmatrix *transptr;                   /* transptr used in Restml */

typedef long sitearray[3];
typedef sitearray *seqptr;                       /* seqptr used in Protpars */

/* datastructure typedefs */
enum node_type { FORK_NODE = 0, TIP_NODE, FREE_NODE };
typedef enum node_type node_type;  /* debug:  maybe remove "enum"?  Needed at all? */

typedef enum nodetype {                                /* what kind of data */
  NODE_T_UNKNOWN,      /* debug:  maybe rename this type "nodedatatype"? */
  NODE_T_GENERIC,
  NODE_T_ML,
  NODE_T_DNA,
  NODE_T_PROT
} nodetype;

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


typedef void (*tree_re_move_t)(struct tree*, struct node*, struct node**, 
		                                             boolean);
typedef boolean (*tree_addtraverse_t)(struct tree*, struct node*, 
		           struct node*, traversetype, struct node*, double*, 
			   struct tree*, boolean, boolean, boolean, double*);
typedef boolean (*tree_addtraverse_1way_t)(struct tree*, struct node*, 
		   struct node*, traversetype, struct node**, double*, 
		   struct tree*, boolean, boolean, boolean*, double*);
typedef void (*tree_insert_t)(struct tree*, struct node*, struct node*, 
		                                                     boolean);
typedef boolean (*tree_try_insert_t)(struct tree*, struct node*, struct node*, 
                   struct node*, double*, struct tree*, boolean, 
		   boolean, boolean, double*);
typedef void (*tree_free_t)(struct tree*);
typedef void (*tree_globrearrange_t)(struct tree*, struct tree*, boolean, 
                                       boolean, double*);
typedef void (*tree_locrearrange_t)(struct tree*, struct node*, boolean, 
		                      double*, struct tree*,struct tree*, 
				      boolean, double*);
typedef void (*tree_smoothall_t)(struct tree*, struct node*);
typedef double (*tree_evaluate_t)(struct tree*, struct node*, boolean); /* debug: needed? */
typedef void (*tree_save_lr_nodes_t)(struct tree*, struct node*);
typedef void (*tree_restore_lr_nodes_t)(struct tree*, struct node*); 
typedef void (*tree_save_traverses_t)(struct tree*, struct node*);
typedef void (*tree_restore_traverses_t)(struct tree*, struct node*);
typedef void (*tree_release_fork_t)(struct tree*, long);
typedef struct node* (*tree_get_fork_t)(struct tree*,  long);
typedef struct node* (*tree_get_forknode_t)(struct tree*, long);
typedef void (*tree_release_forknode_t)(struct tree*, struct node*);
typedef void (*tree_reinit_forknode_t)(struct tree*, struct node*);
typedef void (*tree_nuview_t)(struct tree*, struct node*);
typedef void (*tree_makenewv_t)(struct tree*, struct node*);
typedef void (*tree_print_t)(struct tree*);

typedef boolean (*tree_good_t)(struct tree*);
typedef boolean (*node_good_t)(struct tree*, struct node*);

typedef struct tree_vtable tree_vtable;              /* forward declaration */

typedef enum {                           /* enum type which is type of tree */
  TREE_T_UNKNOWN,
  TREE_T_GENERIC,
  TREE_T_ROOTED,
  TREE_T_UNROOTED,
  TREE_T_ML
} treetype;

typedef void (*initfuncptr)(struct tree *, struct node **, long, long,
                         long *, long *, initops, pointarray,
                         Char *, Char *, FILE *);

/* some pointers to functions needed in the node and tree class hierarchies */

typedef struct initdata {
  tree_new_t tree_new;                                  /* makes a new tree */
  tree_init_t tree_init;                     /* initiates stuff in the tree */
  node_new_t node_new;                      /* makes a new node in the tree */
  node_init_t node_init;                     /* initiates stuff in the node */
} initdata;

/* End. */
