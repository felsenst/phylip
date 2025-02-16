/* Version 4.0.  Copyright 2022.
   Written by Joe Felsenstein  */

/* specializing  ml_node  for DNA data type */

#ifndef MLDNA_H
#define MLDNA_H

/* end of ifdef block if have not previously defined the mldna.h stuff */

#ifndef SEQ_H
#include "seq.h"
#endif

#ifndef BL_H
#include "bl.h"
#endif

#ifndef ML_H
#include "ml.h"
#endif


typedef struct mldna_node{                          /* subclass of ml_node */
  struct ml_node ml_node;                     /* Base object, must be first */
  allocx_t allocx_f;
  freex_t freex_f;
  phenotype x;
} mldna_node;

typedef struct nuview_data {
  /* A big 'ol collection of pointers used in nuview */
  double *yy, *wwzz, *vvzz, *vzsumr, *vzsumy, *sum, *sumr, *sumy;
  sitelike *xx;
} nuview_data;

typedef struct basefreq {
  double        a, c, g, t, r, y;                       /* base frequencies */
  double        ar, cy, gr, ty;   /* base per purine or base per pyrimidine */
  double        xi, xv; /* rates: transition-like, transversion-like events */
  double        fracchange;      /* fraction of events that change the base */
  double        ttratio;                   /* transition/transversion ratio */
} basefreq;

#ifndef OLDC    /* Prototypes, if not original Kernighan & Ritchie compiler */
struct node* mldna_node_new(node_type, long, long);
void mldna_node_init(struct node*, node_type, long);
void mldna_node_copy(struct node*, struct node*);
void fix_x(struct mldna_node*, long, double, long);
void mldna_node_freex(struct node*);
void mldna_node_allocx(struct node*, long, long);
void makevalues2(long, pointarray, long, long, sequence, steptr);
void freex_notip(long, pointarray);
void freex(long, pointarray);
void print_basefreq(FILE *fp, basefreq *freq, boolean empirical);
void makebasefreq(basefreq *freq, double freqa, double freqc,
                   double freqg, double freqt, double ttratio);
void getbasefreqs(double, double, double, double, double *, double *,
                   double *, double *, double *, double *, double *,
                   double *, double *, double *, boolean, boolean);
void ttratio_warning(double ttratio);
void empiricalfreqs(double*, double*, double*, double*, steptr, pointarray);
#endif

#endif
