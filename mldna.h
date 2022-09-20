/* Version 4.0.  Copyright 2022.
   Written by Joe Felsenstein  */

/* specializing  ml_node  for DNA data type */

long rcategs;                    /* number of rate categories, default is 1 */

typedef struct mldna_node{                          /* subclass of ml_node */
  struct ml_node ml_node;                     /* Base object, must be first */
  double* underflows;
  phenotype x;
} mldna_node;


#ifndef OLDC  /* Prototypes, if not original Kernighan & Ritchie compiler */
node*   mldna_node_new(node_type, long, long);
void    mldna_node_init(struct node *, node_type, long);
void    mldna_node_copy(node*, node*);
void    fix_x(mldna_node*, long, double, long);
void    mldna_node_freex(ml_node*);
void    mldna_node_allocx(struct node*, long, long);
void    makevalues2(long, pointarray, long, long, sequence, steptr);
void    freex_notip(long, pointarray);
void    freex(long, pointarray);
void    empiricalfreqs(double*, double*, double*, double*, steptr, pointarray);
#endif

