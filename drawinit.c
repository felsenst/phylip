/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#include "draw.h"


extern boolean haslengths;


void initdrawnode(tree * treep, node **p, long len,
                  long nodei, long *ntips, long *parens, initops whichinit,
                  pointarray treenode, Char *str, Char *ch,
                  FILE *intree, negValFix whatToDo)
{
  long i;
  boolean minusread;
  double valyew, divisor;

  (void)treenode;                       // RSGnote: Parameter never used.

  //printf("in initdrawnode name: %s ", str);

  switch (whichinit)
  {
    case bottom:
      //printf("creating bottom\n");
      *p = treep->get_forknode(treep, nodei);

      (*p)->index = nodei;
      (*p)->tip = false;
      for (i=0; i < MAXNCH; i++)
        (*p)->nayme[i] = '\0';
      treep->nodep[(*p)->index - 1] = (*p);
      break;

    case nonbottom:
      //printf("creating nonbottom\n");
      *p = treep->get_forknode(treep, nodei);

      (*p)->index = nodei;
      break;

    case tip:
      //printf("creating tip\n");
      (*ntips)++;
      *p = treep->get_forknode(treep, nodei);

      treep->nodep[(*ntips) - 1] = *p;
      setupnode(*p, *ntips);
      (*p)->tip        = true;
      (*p)->naymlength = len ;
      strncpy ((*p)->nayme, str, MAXNCH);
      break;

    case length:
      //printf("processing length\n");
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);

      if (!minusread)
        (*p)->oldlen = valyew / divisor;
      else
      {
        switch (whatToDo)
        {
          case negfix_AS_IS:
            (*p)->oldlen = valyew / divisor;
            break;
          case negfix_ZERO:
            (*p)->oldlen = 0.0;
            break;
          case negfix_FABS:
            (*p)->oldlen = fabs ( valyew / divisor );
            if ( (*p)->oldlen < epsilon )
            {
              (*p)->oldlen = epsilon;
            }
            break;
          default:
            // should not happen
            break;
        }
      }

      if ((*p)->back != NULL)
        (*p)->back->oldlen = (*p)->oldlen;

      break;

    case hsnolength:
      //printf("has no length\n");
      haslengths = false;
      break;

    case iter:
      // save internal node name if there is one
      (*p)->naymlength = len ;
      strncpy ((*p)->nayme, str, MAXNCH);
      break;

    default:        // cases hslength, treewt, unitrwt should not occur
      //printf("other whichinit: %i\n", whichinit);
      break;
  }
} // initdrawnode


void initdrawgramnode(tree * treep, node **p, long len,
                      long nodei, long *ntips, long *parens, initops whichinit,
                      pointarray treenode, Char *str, Char *ch,
                      FILE *intree)
{
  /* initialize a node for a tree for drawgram */

  initdrawnode(treep, p, len, nodei, ntips, parens, whichinit, treenode, str, ch, intree, negfix_ZERO);
} /* initdrawgramnode */


void initdrawtreenode(tree * treep, node **p, long len,
                      long nodei, long *ntips, long *parens, initops whichinit,
                      pointarray treenode, Char *str, Char *ch,
                      FILE *intree)
{
  /* initialize a node for a tree for drawtree */

  initdrawnode(treep, p, len, nodei, ntips, parens, whichinit, treenode, str, ch, intree, negfix_FABS);
} // initdrawtreenode


// End.
