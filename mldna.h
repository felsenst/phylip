/* Copyright 2022 */

long rcategs;                    /* number of rate categories, default is 1 */

typedef struct ml_dna_node{                          /* subclass of ml_node */
  struct ml_node ml_node;                     /* Base object, must be first */
  phenotype x;
} ml_dna_node;

#ifndef OLDC  /* Prototypes, if not original Kernighan & Ritchie compiler */
node*   ml_dna_node_new(node_type, long, long);
void    ml_dna_node_init(struct node *, node_type, long);
void    ml_dna_node_allocx(struct node*, long, long);
void    fix_x(ml_dna_node*, long, double, long);
void    ml_dna_node_copy(node*, node*);
void    ml_dna_node_freex(ml_node*);
#endif

