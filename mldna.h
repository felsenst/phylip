/* Version 4.0.  Copyright 2022.
   Written by Joe Felsenstein  */

/* specializing  ml_node  for DNA data type */

#ifndef ML_H
#include "ml.h"
#endif

#ifndef SEQ_H
#include "seq.h"
#endif

long rcategs;                    /* number of rate categories, default is 1 */

typedef struct mldna_node{                          /* subclass of ml_node */
  struct ml_node ml_node;                     /* Base object, must be first */
  double* underflows;
  phenotype x;
} mldna_node;

typedef struct nuview_data {
  /* A big 'ol collection of pointers used in nuview */
  double *yy, *wwzz, *vvzz, *vzsumr, *vzsumy, *sum, *sumr, *sumy;
  sitelike *xx;
} nuview_data;

typedef struct basefreq {
  double        a, c, g, t, r, y;       /* Base frequencies */
  double        ar, cy, gr, ty;         /* Base/class freq */
  double        xi, xv;
  double        fracchange;
  double        ttratio;                /* Transition/transversion ratio */
} basefreq;


#ifndef OLDC  /* Prototypes, if not original Kernighan & Ritchie compiler */
struct mldna_node* mldna_node_new(node_type, long, long);
void mldna_node_init(struct mldna_node *, node_type, long);
void mldna_node_copy(mldna_node*, mldna_node*);
void fix_x(mldna_node*, long, double, long);
typedef void (*freex_t)(struct bl_tree*);
void mldna_node_freex(mldna_node*);
void mldna_node_allocx(struct mldna_node*, long, long);
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

