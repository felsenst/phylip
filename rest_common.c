/* Version 4.0. (c) Copyright 2012-2013 by the University of Washington.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#include "phylip.h"
#include "rest_common.h"


extern steptr   weight, alias, aliasweight;
extern sequence inputSequences;

// Prototypes:
void    sitescrunch2(long, long, long, steptr);


void rest_makeweights(long sites, long * endsitePtr)
{
  /* make up weights vector to avoid duplicate computations */
  long i;

#ifdef BUG_968
  printf("BUG.968 at start, sites:%ld endsite:%ld\n", sites, *endsitePtr);
#endif

  for (i = 1; i <= sites; i++)
  {
    alias[i] = i;
    aliasweight[i] = weight[i];
  }

#ifdef BUG_968
  printf("BUG.968 alias : ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", alias[i]);
  }
  printf("\n");
  printf("BUG.968 aliasw: ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", aliasweight[i]);
  }
  printf("\n");
  printf("\n");
#endif

  rest_sitesort(sites);
#ifdef BUG_968
  printf("BUG.968 after sitesort:\n");
  printf("BUG.968 alias : ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", alias[i]);
  }
  printf("\n");
  printf("BUG.968 aliasw: ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", aliasweight[i]);
  }
  printf("\n");
  printf("\n");
#endif

  rest_sitecombine(sites);
#ifdef BUG_968
  printf("BUG.968 after sitecombine:\n");
  printf("BUG.968 alias : ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", alias[i]);
  }
  printf("\n");
  printf("BUG.968 aliasw: ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", aliasweight[i]);
  }
  printf("\n");
  printf("\n");
#endif

  sitescrunch2(sites + 1, 2, 3, aliasweight);
#ifdef BUG_968
  printf("BUG.968 after sitescrunch2:\n");
  printf("BUG.968 alias : ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", alias[i]);
  }
  printf("\n");
  printf("BUG.968 aliasw: ");
  for (i=0; i <= sites; i++)
  {
    printf("%3ld ", aliasweight[i]);
  }
  printf("\n");
  printf("\n");
#endif

  for (i = 1; i <= sites; i++)
  {
    weight[i] = aliasweight[i];
    if (weight[i] > 0)
      *endsitePtr = i;
  }
  weight[0] = 1;
}  /* makeweights */


void rest_sitecombine(long sites)
{
  /* combine sites that have identical patterns */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites)
  {
    j = i + 1;
    tied = true;
    while (j <= sites && tied)
    {
      k = 1;
      while (k <= spp && tied)
      {
        tied = (tied && inputSequences[k - 1][alias[i] - 1] == inputSequences[k - 1][alias[j] - 1]);
        k++;
      }
      if (tied && aliasweight[j] > 0)
      {
        aliasweight[i] += aliasweight[j];
        aliasweight[j] = 0;
        alias[j] = alias[i];
      }
      j++;
    }
    i = j - 1;
  }
}  /* rest_sitecombine */


void rest_sitesort(long sites)
{
  /* Shell sort keeping alias, aliasweight in same order */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0)
  {
    for (i = gap + 1; i <= sites; i++)
    {
      j = i - gap;
      flip = true;
      while (j > 0 && flip)
      {
        jj = alias[j];
        jg = alias[j + gap];
        flip = false;
        tied = true;
        k = 1;
        while (k <= spp && tied)
        {
          flip = (inputSequences[k - 1][jj - 1] > inputSequences[k - 1][jg - 1]);
          tied = (tied && inputSequences[k - 1][jj - 1] == inputSequences[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j];
        alias[j] = alias[j + gap];
        alias[j + gap] = itemp;
        itemp = aliasweight[j];
        aliasweight[j] = aliasweight[j + gap];
        aliasweight[j + gap] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* rest_sitesort */


// End.
