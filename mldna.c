/* Copyright, 2022 */
/* functions for ML analysis on DNA/RNA sequence data */

void dna_node_copy(node* srcn, node* destn)
{
  /* copy a node when DNA likelihoods are used */
  dna_node * src  = (dna_node *)srcn;
  dna_node * dest = (dna_node *)destn;
  long i, j;
  long oldendsite = dest->ml_node.endsite;

  ml_node_copy(srcn, destn);

  if ( oldendsite != 0 && oldendsite != src->ml_node.endsite )
  {
    dest->ml_node.freex((node*)dest);
    dest->ml_node.endsite = 0;
  }
  if ( oldendsite == 0 )
    ((ml_node*)dest)->allocx(((node*)dest), ((ml_node*)src)->endsite, ((ml_node*)src)->categs);
  for (i = 0; i < ((ml_node*)src)->endsite; i++)
    for (j = 0; j < ((ml_node*)src)->categs; j++)
      memcpy(((dna_node*)dest)->x[i][j], ((dna_node*)src)->x[i][j], sizeof(sitelike));
}


void fix_x(dna_node* p, long site, double maxx, long rcategs)
{ /* used in  Dnaml, Dnamlk */
  long i, j;
  ((ml_node*)p)->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs ; i++ )
  {
    for ( j = 0 ; j < ((long)T - (long)A + 1) ; j++)
      p->x[site][i][j] /= maxx;
  }
} /* fix_x */


node * dna_node_new(node_type type, long index, long nodesize) // RSGbugfix
{
/* debug:  how does this work, is it in node_new hierarchy? */
  node* n = Malloc(sizeof(dna_node));  /* debug: where allocated? */

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

/* debug: instead do a call upwards here? */
  dna_node_init(n, type, index);
  return n;
} /* dna_node_new */


void dna_node_init(node *node, node_type type, long index)
{
  /* initialize a node for a dna tree */

  dna_node *n = (dna_node *)node;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  ml_node_init(node, type, index);
  n->ml_node.allocx = (allocx_t)dna_node_allocx;
  n->ml_node.node.copy = dna_node_copy;
  n->ml_node.node.init = dna_node_init;
  n->ml_node.freex = (freex_t)dna_node_freex;
  n->x = NULL;

  if ( endsite != 0 && rcategs != 0 )
    n->ml_node.allocx((struct node*)n, endsite, rcategs);
} /* dna_node_init */


void dna_node_freex(ml_node* n)
{
  /* free a dna tree node */
  dna_node *dn;
  long i;

  dn = (dna_node *)n;
  for ( i = 0 ; i < n->endsite ; i++ )
  {
    free(dn->x[i]);
  }

  free(dn->x);
  dn->x = NULL;
  free(n->underflows);
  n->underflows = NULL;
} /* dna_node_freex */


void dna_node_allocx(node* n, long endsite, long rcategs)
{
  /* allocate space for sequences on a dna tree node */
  ml_node *mln = (ml_node *)n;
  dna_node *dn = (dna_node *)n;
  long i;

  dn->x = (phenotype)Malloc(endsite * sizeof(ratelike));
  for ( i = 0 ; i < endsite ; i++ )
    dn->x[i] = (ratelike)Malloc(rcategs * sizeof(sitelike));

  mln->categs = rcategs;
  mln->endsite = endsite;
  mln->underflows = Malloc(endsite * sizeof(double));
} /* dna_node_allocx */


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
          ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 0.0;

        switch (y[i][j - 1])
        {
          case 'A':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            break;

          case 'C':
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'G':
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'T':
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'U':
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'M':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'R':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'W':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'S':
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'Y':
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'K':
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'B':
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'D':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'H':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'V':
            ((dna_node*)nodep[i])->x[k][l][0] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)C - (long)A] = 1.0;
            ((dna_node*)nodep[i])->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'N':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case 'X':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '?':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
          break;

          case 'O':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '-':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((dna_node*)nodep[i])->x[k][l][(long)b - (long)A] = 1.0;
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
        free(((dna_node*)p)->x[j]);
      free(((dna_node*)p)->x);
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
      free(((dna_node*)treenode[i])->x[j]);
    free(((dna_node*)treenode[i])->x);
  }

  for (i = spp; i < nonodes; i++)
  {
    if(treenode[i])
    {
      p = treenode[i];
      do {
        for (j = 0; j < endsite; j++)
          free(((dna_node*)p)->x[j]);
        free(((dna_node*)p)->x);
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
        sum = (*freqa) * ((dna_node*)treenode[i])->x[j][0][0];
        sum += (*freqc) * ((dna_node*)treenode[i])->x[j][0][(long)C - (long)A];
        sum += (*freqg) * ((dna_node*)treenode[i])->x[j][0][(long)G - (long)A];
        sum += (*freqt) * ((dna_node*)treenode[i])->x[j][0][(long)T - (long)A];
        suma += w * (*freqa) * ((dna_node*)treenode[i])->x[j][0][0] / sum;
        sumc += w * (*freqc) * ((dna_node*)treenode[i])->x[j][0][(long)C - (long)A] / sum;
        sumg += w * (*freqg) * ((dna_node*)treenode[i])->x[j][0][(long)G - (long)A] / sum;
        sumt += w * (*freqt) * ((dna_node*)treenode[i])->x[j][0][(long)T - (long)A] / sum;
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


