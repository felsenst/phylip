/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#define OVER              8
#define ADJACENT_PAIRS    1
#define CORR_IN_1_AND_2   2
#define ALL_IN_1_AND_2    3
#define NO_PAIRING        4
#define ALL_IN_FIRST      5
#define TREE1             8
#define TREE2             9

#define FULL_MATRIX       11
#define COMPUTER_READABLE_MATRIX  22
#define VERBOSE           33
#define SPARSE            44

/* Number of columns per block in a matrix output */
#define COLUMNS_PER_BLOCK 10

typedef struct cons_node {
  node node;
  group_type* nodeset;
} cons_node;

typedef struct pattern_elm {
  group_type *apattern;
  long *patternsize;
  double *length;
} pattern_elm;

#define NUM_BUCKETS 100

typedef struct namenode {
  struct namenode *next;
  plotstring naym;
  int hitCount;
} namenode;

typedef namenode **hashtype;

#ifndef OLDC
/* function prototypes */
void initconsnode(tree *, node **, long, long, long *, long *,
                  initops, pointarray, Char *, Char *, FILE *);
void   compress(long *);
void   sort(long);
void   eliminate(long *, long *);
void   printset(long);
void   bigsubset(group_type *, long);
void   recontraverse(tree*, node **, group_type *, long, long *);
void   reconstruct(tree*, long);
void   coordinates(node *, long *);
void   drawline(long i);

void   printree(void);
void   consensus(tree *, pattern_elm ***, long);
void   rehash(void);
void   enternodeset(tree * treep, node *r);
void   accumulate(tree * treep, node *);
void   dupname2(Char *, node *, node *);
void   dupname(node *);
void   missingname(node *);
void   initreenode(node *);
void   reroot(tree *, node *, long *);

void   store_pattern (pattern_elm ***, int);
boolean samename(naym, plotstring);
void   reordertips(tree * treep);
void   read_groups (pattern_elm ****, long , long, FILE *);
node*  cons_node_new(node_type type, long index);
void   adjustsupport(long );
void   censor(void);
boolean   compatible(long, long);
void   elimboth(long);
void   enternohash(group_type*, long*);
void   enterpartition (group_type*, long*);
void   reorient(tree * treep, node* n);
/* function prototypes */
#endif

extern long setsz;


// End.
