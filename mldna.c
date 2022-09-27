/* Copyright, 2022 */
/* functions for ML analysis on DNA/RNA sequence data */

#include "phylip.h"
#include "ml.h"
#include "mldna.h"


node * mldna_node_new(node_type type, long index, long nodesize) // RSGbugfix
{
  struct node* n;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  n = ml_node_new(type, index, nodesize);
  mldna_node_init(n, type, index);
  return n;
} /* mldna_node_new */


void mldna_node_init(node *node, node_type type, long index)
{
  /* initialize a node for an ml dna tree */

  mldna_node *n = (mldna_node *)node;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  ml_node_init(node, type, index);
  n->ml_node.allocx = (allocx_t)mldna_node_allocx;
  n->ml_node.bl_node.copy = mldna_node_copy;
  n->ml_node.bl_node.node_init = mldna_node_init;
  n->ml_node.freex = (freex_t)mldna_node_freex;
  n->x = NULL;

  if ( endsite != 0 && rcategs != 0 )
    n->ml_node.allocx((struct node*)n, endsite, rcategs);
} /* mldna_node_init */


void mldna_node_copy(node* srcn, node* destn)
{
  /* copy a node when DNA likelihoods are used */
  mldna_node * src  = (mldna_node *)srcn;
  mldna_node * dest = (mldna_node *)destn;
  long i, j;
  long oldendsite = dest->ml_node.endsite;

  ml_node_copy(srcn, destn);

  if ( oldendsite != 0 && oldendsite != src->ml_node.endsite )
  {
    dest->ml_node.freex((node*)dest);
    dest->ml_node.endsite = 0;
  }
  if ( oldendsite == 0 )
    ((ml_node*)dest)->allocx((node*)dest, ((ml_node*)src)->endsite, ((ml_node*)src)->categs);
  for (i = 0; i < ((ml_node*)src)->endsite; i++)
    for (j = 0; j < ((ml_node*)src)->categs; j++)
      memcpy(((mldna_node*)dest)->x[i][j], ((mldna_node*)src)->x[i][j], sizeof(sitelike));
} /* mldna_node_copy */


void fix_x(mldna_node* p, long site, double maxx, long rcategs)
{ /* used in  Dnaml, Dnamlk */
  long i, j;

  ((ml_node*)p)->underflows[site] += log(maxx);
  for ( i = 0 ; i < rcategs ; i++ )
  {
    for ( j = 0 ; j < ((long)T - (long)A + 1) ; j++)
      p->x[site][i][j] /= maxx;
  }
} /* fix_x */


void mldna_node_freex(ml_node* n)
{
  /* free a dna tree node */
  mldna_node *dn;
  long i;

  dn = (mldna_node *)n;
  for ( i = 0 ; i < n->endsite ; i++ )
  {
    free(dn->x[i]);
  }

  free(dn->x);
  dn->x = NULL;
  free(n->underflows);
  n->underflows = NULL;
} /* mldna_node_freex */


void mldna_node_allocx(node* n, long endsite, long rcategs)
{
  /* allocate space for sequences on a dna tree node */
  ml_node *mln = (ml_node *)n;
  mldna_node *dn = (mldna_node *)n;
  long i;

  dn->x = (phenotype)Malloc(endsite * sizeof(ratelike));
  for ( i = 0 ; i < endsite ; i++ )
    dn->x[i] = (ratelike)Malloc(rcategs * sizeof(sitelike));

  mln->categs = rcategs;
  mln->endsite = endsite;
  mln->underflows = Malloc(endsite * sizeof(double));
} /* mldna_node_allocx */


void makevalues2(long categs, pointarray nodep, long endsite, long spp, sequence y, steptr alias)
{
  /* set up fractional likelihoods at tips 
   * used by dnaml & dnamlk */
  long i, j, k, l;
  bases b;

  for (k = 0; k < endsite; k++)
  {
    j = alias[k];
    for (i = 0; i < spp; i++)
    {
      for (l = 0; l < categs; l++)
      {
        for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
          ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 0.0;

        switch (y[i][j - 1])
        {
          case 'A':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            break;

          case 'C':
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'G':
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'T':
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'U':
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'M':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'R':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'W':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'S':
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'Y':
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'K':
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'B':
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'D':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'H':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'V':
            ((mldna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((mldna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'N':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case 'X':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '?':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
          break;

          case 'O':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '-':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((mldna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;
        }
      }
    }
  }
}  /* makevalues2 */


void freex_notip(long nonodes, pointarray treenode)
{
  /* free interior fork nodes
   * used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++)
  {
    p = treenode[i];
    if ( p == NULL ) continue;
    do
    {
      for (j = 0; j < endsite; j++)
        free(((mldna_node*)p)->x[j]);
      free(((mldna_node*)p)->x);
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* freex_notip */


void freex(long nonodes, pointarray treenode)
{
  /* used in dnaml & dnamlk */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++)
  {
    for (j = 0; j < endsite; j++)
      free(((mldna_node*)treenode[i])->x[j]);
    free(((mldna_node*)treenode[i])->x);
  }

  for (i = spp; i < nonodes; i++)
  {
    if(treenode[i])
    {
      p = treenode[i];
      do {
        for (j = 0; j < endsite; j++)
          free(((mldna_node*)p)->x[j]);
        free(((mldna_node*)p)->x);
        p = p->next;
      } while (p != treenode[i]);
    }
  }
}  /* freex */


void empiricalfreqs(double *freqa, double *freqc, double *freqg, double *freqt, steptr weight, pointarray treenode)
{
  /* Get empirical base frequencies from the data */
  /* used in dnaml & dnamlk */
  /* this is kind of strange */
  long i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  *freqa = 0.25;
  *freqc = 0.25;
  *freqg = 0.25;
  *freqt = 0.25;

  for (k = 1; k <= 8; k++)
  {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;

    for (i = 0; i < spp; i++)
    {
      for (j = 0; j < endsite; j++)
      {
        w = weight[j];
        sum = (*freqa) * ((mldna_node*)treenode[i])->x[j][0][0];
        sum += (*freqc) * ((mldna_node*)treenode[i])->x[j][0][(long)C - (long)A];
        sum += (*freqg) * ((mldna_node*)treenode[i])->x[j][0][(long)G - (long)A];
        sum += (*freqt) * ((mldna_node*)treenode[i])->x[j][0][(long)T - (long)A];
        suma += w * (*freqa) * ((mldna_node*)treenode[i])->x[j][0][0] / sum;
        sumc += w * (*freqc) * ((mldna_node*)treenode[i])->x[j][0][(long)C - (long)A] / sum;
        sumg += w * (*freqg) * ((mldna_node*)treenode[i])->x[j][0][(long)G - (long)A] / sum;
        sumt += w * (*freqt) * ((mldna_node*)treenode[i])->x[j][0][(long)T - (long)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;
    *freqa = suma / sum;
    *freqc = sumc / sum;
    *freqg = sumg / sum;
    *freqt = sumt / sum;
  }
  if (*freqa <= 0.0)
    *freqa = 0.000001;
  if (*freqc <= 0.0)
    *freqc = 0.000001;
  if (*freqg <= 0.0)
    *freqg = 0.000001;
  if (*freqt <= 0.0)
    *freqt = 0.000001;
}  /* empiricalfreqs */


