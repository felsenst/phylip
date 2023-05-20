/* Version 4.0. (c) Copyright 1993-2022
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


/*
  seq.h:  for dealing with molecular sequence data types
          included in dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnamove,
          dnapars, dnapenny, protdist, protpars, & restml
*/

#ifndef _SEQ_H_
#define _SEQ_H_

#include "phylip.h"
#include "bl.h"

/* move */
/* All the below moved here in the Great TreeRead Migration of '96 */

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


extern long endsite, outgrno, which;
extern boolean interleaved, printdata, outgropt, treeprint, dotdiff, transvp;
extern steptr weight, category, alias, location, ally;
extern sequence inputSequences;
extern struct bl_node** lrsaves;

#ifndef OLDC
/* function prototypes */
void inputdata(long);
void read_sequences(long nchars);
void output_sequences(long nchars);
void setuptree(pointarray, long, boolean);
void setuptree2(tree);
void alloctip(bl_node *);
void freetrans(transptr *, long, long );
void sitesort(long, steptr);
void sitecombine(long);
void sitescrunch(long);
void sitesort2(long, steptr);
void sitecombine2(long, steptr);
void sitescrunch2(long, long, long, steptr);
void drawline(long, double, struct bl_node *);
void treeout(struct node *, long, long *, struct node *);
void drawline2(long i, double scale, tree *curtree);
void standev(long, long, long, double, double *, long **, longer);
void standev2(long, long, long, long, double,
               double *, double **, steptr, longer);
void freex2(long, pointarray);
void inittrees(long, long);
void resetlrsaves(long param1, long param2);
/*function prototypes*/
#endif

#endif
/* above ifndef ... endif  if this header has not been #included before */

/* End. */
