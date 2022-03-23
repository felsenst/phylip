/* Version 4.0. (c) Copyright 1993-2022.
   Written mostly by Sean Lamont, Andrew Keeffe, Akiko Fuseki, and Joe Felsenstein */

/*
  cont.h: included in contml & contrast
*/


#ifndef _CONT_H_     /* debug:  why? */
#define _CONT_H_

typedef double* view;  /* debug: can we define this here? */

typedef struct cont_node_type
{
  struct node node;
  phenotype3 view;
  long totalleles;
} cont_node_type;


#ifndef OLDC      /* If not the old original Kernighan and Ritchie compiler */
/*function prototypes*/
node* cont_node_new(node_type, long, long);
void cont_node_init(struct node*, node_type, long);
void cont_node_copy(struct node* src, struct node* dst);
void alloctree(pointarray *, long);
void freetree(pointarray *, long);
void setuptree(tree *, long);
void allocview(tree *, long, long);
void freeview(tree *, long);
void standev2(long, long, long, long, double, double *, double **, longer);
#endif

#endif /* _CONT_H_ */


/* End. */
