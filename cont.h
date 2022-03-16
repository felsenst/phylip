/* Version 4.0. (c) Copyright 1993-2013.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe */


/*
  cont.h: included in contml & contrast
*/


#ifndef _CONT_H_
#define _CONT_H_

typedef double* view;  /* debig: can we define this here? */

typedef struct cont_node_type
{
  node node_var;
  phenotype3 view;
  long totalleles;
} cont_node_type;


#ifndef OLDC
/*function prototypes*/
node* cont_node_new(node_type, long);
void cont_node_init(cont_node_type*, node_type, long);
void cont_node_copy(node* src, node* dst);
void alloctree(pointarray *, long);
void freetree(pointarray *, long);
void setuptree(tree *, long);
void allocview(tree *, long, long);
void freeview(tree *, long);
void standev2(long, long, long, long, double, double *, double **, longer);
/*function prototypes*/
#endif

#endif /* _CONT_H_ */


// End.
