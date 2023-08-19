/* Version 4.0. (c) Copyright 1993-2023
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */
/*
  seq.h:  for dealing with molecular sequence data types
          included in dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnamove,
          dnapars, dnapenny, protdist, protpars, & restml
*/

#ifndef SEQ_H
#define SEQ_H

#include "phylip.h"

#define ebcdic          EBCDIC

/* All of this came over from cons.h    -plc*/
#define OVER              7
#define ADJACENT_PAIRS    1
#define CORR_IN_1_AND_2   2
#define ALL_IN_1_AND_2    3
#define NO_PAIRING        4
#define ALL_IN_FIRST      5
#define TREE1             8
#define TREE2             9

#define FULL_MATRIX       11
#define VERBOSE           22
#define SPARSE            33

/* Number of columns per block in a matrix output */
#define COLUMNS_PER_BLOCK 10

extern boolean transvp;
extern steptr weight, category, alias, location, ally;
extern sequence inputSequences;
/* debug:  extern FILE *infile, *outfile, *intree, *intree2, *outtree;  */
extern struct bl_node** lrsaves;

typedef void (*freex_t)(long, pointarray);       /* pointer to free fn type */

extern freex_t *freex_f;               /* forward: pointer to free function */

#ifndef OLDC
/* function prototypes.  Needed if not the old 
   original Kernighan & Ritchie compiler */
void inputdata(long);
void read_sequences(long nchars);
void output_sequences(long nchars);
void setuptree(pointarray, long, boolean);
void setuptree2(struct tree);
void alloctip(struct bl_node *);
void freetrans(transptr *, long, long );
void sitesort(long, steptr);
void sitecombine(long);
void sitescrunch(long);
void sitesort2(long, steptr);
void sitecombine2(long, steptr);
void sitescrunch2(long, long, long, steptr);
void drawline(long, double, struct bl_node *);
void treeout(struct node *, long, long *, struct node *);
void drawline2(long, double, struct tree *);
void standev(long, long, long, double, double *, long **, longer);
void standev2(long, long, long, long, double,
               double *, double **, steptr, longer);
void freex2(long, pointarray);
void inittrees(long, long);
void resetlrsaves(long, long);
/*function prototypes*/
#endif

#endif

/* end of #ifndef that conditions on this header file not already used */

/* End. */
