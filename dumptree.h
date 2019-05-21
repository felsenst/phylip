/* Version 3.7. (c) Copyright 2012 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifndef DUMPTREE_H
#define DUMPTREE_H


#include "phylip.h"


void dumpNodeLoop(node *, char * indent);
void dumpNodeBack(node *, char * indent);
void dumpOneNode(node *, long index, char * indent);
void dumpSlist(Slist_ptr s, char * indent, char * name);
void dumppointarray(pointarray p, long numNodes, char * indent);
void dumptree(tree * t);

#endif /* DUMPTREE_H */


// End.
