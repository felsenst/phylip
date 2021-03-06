
/* version 3.6. (c) Copyright 1993-2011 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
    cont.h: included in contml & contrast
*/

#include "phylip.h"

typedef struct cont_node_type 
{
  node node_var;
  phenotype3 view;
  long totalleles;

} cont_node_type;


#ifndef OLDC
/*function prototypes*/
node* cont_node_new(node_type, long);
void cont_node_init(cont_node_type* n, boolean tip, long index);
void cont_node_copy(node* src, node* dst);
void alloctree(pointarray *, long);
void freetree(pointarray *, long);
void setuptree(tree *, long);
void allocview(tree *, long, long);
void freeview(tree *, long);
void standev2(long, long, long, long, double, double *, double **, longer);
/*function prototypes*/
#endif
