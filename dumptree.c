/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#include "dumptree.h"
#include "phylip.h"

#include <stdio.h>


long maxindentch = 100;


void dumpNodeLoop(node * n, char * indent)
{
  /* copied from generic_fork_print */
  printf("%s",indent);
  boolean firstTime = true;
  boolean nulledOut = false;
  node * p = n;
  while((firstTime || (p != n)) && !nulledOut)
  {
    if(p == NULL)
    {
      nulledOut = true;
    }
    else
    {
      // might also instead use p->node_print_f(p)
      printf("%p [%p] :: ",p,p->back);
      p = p->next;
    }
    firstTime = false;
  }
  printf("\n");
}


void dumpNodeBack(node * n, char * indent)
{
  printf("%s%p\n",indent,n->back);
}


void dumpOneNode(node * n, long index, char * indent)
{
  char newIndent[maxindentch];
  sprintf(newIndent,"%s        ",indent);

  printf ("%s%3ld: %p %ld\n",indent,index,n,n->index);
  dumpNodeLoop(n,newIndent);
}


void dumpSlist(Slist_ptr s,char * indent, char * name)
{
  printf("%sSlist %s length %ld\n",indent,name,s->length);
}


void dumppointarray(pointarray p, long numNodes, char* indent)
{
  long pindex;
  char newIndent[maxindentch];

  printf("%sPOINTARRAY %p\n", indent,p);
  sprintf(newIndent,"%s\t",indent);

  for(pindex = 0; pindex < numNodes; pindex++)
  {
    dumpOneNode(p[pindex],pindex,newIndent);
  }
}


void dumptree(tree * t)
{
  printf ("Dumping tree %p\n", t);

  printf ("\ttree root %p\n", t->root);
  printf ("\tnum nodes %ld\n", t->nonodes);
  printf ("\tnum spp   %ld\n", t->spp);

  dumpSlist(t->free_forkrings, "\t", "forkrings");
  dumpSlist(t->free_forknodes, "\t", "forknodes");

  dumppointarray(t->nodep,t->nonodes, "\t");
  fflush(stdout);
}


// End.
