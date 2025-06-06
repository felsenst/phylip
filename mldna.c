/* Copyright, 2023 */
/* functions for ML analysis on DNA/RNA sequence data */

#include "mldna.h"

extern node_type type;
extern long endsite;
extern long rcategs;
extern FILE *outfile;
allocx_t allocx_f;
freex_t *freex_f;                      /* forward: pointer to free function */


struct node* mldna_node_new(node_type type, long index, long nodesize) // RSGbugfix
{
  struct node* n;

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  assert(index >= 0);

  nodesize = sizeof(mldna_node);
  n = (struct node*)ml_node_new(type, index, nodesize);
  return n;
} /* mldna_node_new */


void mldna_node_init(struct node* n, node_type type, long index)
{
  /* initialize a node for an ml dna tree */

  // RSGdebug: "index" should be > 0 if used for array access.  Can be 0 only
  // for initialization where it will be changed to > 0 before used for access.
  // Test here is for ">= 0", which allows both cases.
  struct mldna_node* mldn;
  assert(index >= 0);

  mldn = (struct mldna_node*)n;
  ml_node_init(n, type, index);
  n->node_copy = mldna_node_copy;
  mldn->allocx_f = (allocx_t)mldna_node_allocx;
  mldn->freex_f = (freex_t)mldna_node_freex;

  if ( endsite != 0 && rcategs != 0 )
    mldna_node_allocx(n, endsite, rcategs);
} /* mldna_node_init */


void mldna_node_copy(struct node* src, struct node* dest)
{
  /* copy a node when DNA likelihoods are used */
  long i, j, oldendsite;
  struct ml_node *srcmln, *destmln;
  struct mldna_node *srcmldn, *destmldn;

  ml_node_copy(src, dest);                               /* go up hierarchy */
  srcmln = (struct ml_node*)src;     /* src, dest considered as polymorphic */
  destmln = (struct ml_node*)dest;
  oldendsite = destmln->endsite;
  srcmldn = (struct mldna_node*)src;
  destmldn = (struct mldna_node*)dest;

  if ( oldendsite != 0 && oldendsite != srcmln->endsite )
  {
    mldna_node_freex(dest);
    destmln->endsite = 0;
  }
  if ( oldendsite == 0 )    /* debug: why would endsites ever be different? */
    mldna_node_allocx(dest, srcmln->endsite,  srcmln->categs);
  for (i = 0; i < srcmln->endsite; i++)
    for (j = 0; j < srcmln->categs; j++)
      memcpy(destmldn->x[i][j], srcmldn->x[i][j], sizeof(sitelike));
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


void mldna_node_freex(struct node* n)
{
  /* free a dna tree node */
  struct ml_node *mln;
  struct mldna_node *dn;
  long i;

  mln = (struct ml_node *)n;
  dn = (struct mldna_node *)n;

  for ( i = 0 ; i < mln->endsite ; i++ )
  {
    free(dn->x[i]);
  }

  free(dn->x);
  dn->x = NULL;
  free(mln->underflows);
  mln->underflows = NULL;   /* debug:  legal? */
} /* mldna_node_freex */


void mldna_node_allocx(struct node* n, long endsite, long rcategs)
{
  /* allocate space for sequences on a dna tree node */
  long i;
  struct ml_node *mln;
  struct mldna_node *mldn;

  mln = (struct ml_node*)n;           /* node considered as an ml_ node ... */
  mldn = (struct mldna_node*)mln;

  mldn->x = (phenotype)Malloc(endsite * sizeof(ratelike));
  for ( i = 0 ; i < endsite ; i++ )
    mldn->x[i] = (ratelike)Malloc(rcategs * sizeof(sitelike));

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
          ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 0.0;

        switch (y[i][j - 1])
        {
          case 'A':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            break;

          case 'C':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'G':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'T':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'U':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'M':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            break;

          case 'R':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'W':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'S':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'Y':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'K':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'B':
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'D':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'H':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)T - (long)A] = 1.0;
            break;

          case 'V':
            ((struct mldna_node*)(nodep[i]))->x[k][l][0] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)C - (long)A] = 1.0;
            ((struct mldna_node*)(nodep[i]))->x[k][l][(long)G - (long)A] = 1.0;
            break;

          case 'N':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case 'X':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '?':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 1.0;
          break;

          case 'O':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 1.0;
            break;

          case '-':
            for (b = A; (long)b <= (long)T; b = (bases)((long)b + 1))
              ((struct mldna_node*)(nodep[i]))->x[k][l][(long)b - (long)A] = 1.0;
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
  /* this is an EM algorithm valid for independently sampled
   * sequences. These are of course not actually independently sampled because
   * they are on a tree, but ... "abi gezint". */
  long i, j, k;
  double sum, suma, sumc, sumg, sumt, w;

  *freqa = 0.25;
  *freqc = 0.25;
  *freqg = 0.25;
  *freqt = 0.25;

  for (k = 1; k <= 8; k++)                   /* do the EM iteration 8 times */
  {
    suma = 0.0;
    sumc = 0.0;
    sumg = 0.0;
    sumt = 0.0;

    for (i = 0; i < spp; i++)                       /* for each species ... */
    {
      for (j = 0; j < endsite; j++)    /* for each aliased site pattern ... */
      {    /* count of inferred numbers of counts given current frequencies */
        w = weight[j];        /* ... taking into account the number aliased */
        sum = (*freqa) * ((mldna_node*)treenode[i])->x[j][0][0];
        sum += (*freqc) 
                 * ((mldna_node*)treenode[i])->x[j][0][(long)C - (long)A];
        sum += (*freqg) 
                 * ((mldna_node*)treenode[i])->x[j][0][(long)G - (long)A];
        sum += (*freqt) 
                 * ((mldna_node*)treenode[i])->x[j][0][(long)T - (long)A];
        suma += w * (*freqa) * ((mldna_node*)treenode[i])->x[j][0][0] / sum;
        sumc += w * (*freqc) 
               * ((mldna_node*)treenode[i])->x[j][0][(long)C - (long)A] / sum;
        sumg += w * (*freqg) 
               * ((mldna_node*)treenode[i])->x[j][0][(long)G - (long)A] / sum;
        sumt += w * (*freqt) 
               * ((mldna_node*)treenode[i])->x[j][0][(long)T - (long)A] / sum;
      }
    }
    sum = suma + sumc + sumg + sumt;  /* the denominators for the fractions */
    *freqa = suma / sum;        /* the inferred fractions in this iteration */
    *freqc = sumc / sum;
    *freqg = sumg / sum;
    *freqt = sumt / sum;
  }
  if (*freqa <= 0.0)            /* this is done to prevent underflows later */
    *freqa = 0.000001;
  if (*freqc <= 0.0)
    *freqc = 0.000001;
  if (*freqg <= 0.0)
    *freqg = 0.000001;
  if (*freqt <= 0.0)
    *freqt = 0.000001;
}  /* empiricalfreqs */


void makebasefreq(basefreq *freq, double freqa, double freqc,
                   double freqg, double freqt, double ttratio)
{
  /* Takes base frequencies and fills out a basefreq struct. If ttratio is
   * incompatible, a warning is printed and a reasonable value is returned in
   * freq. */
  /* (used by dnadist, dnaml, & dnamlk) */

  double aa, bb;

  assert(freq != NULL);

  freq->a = freqa;
  freq->c = freqc;
  freq->g = freqg;
  freq->t = freqt;

  freq->r = freqa + freqg;
  freq->y = freqc + freqt;

  freq->ar = freq->a / freq->r;
  freq->cy = freq->c / freq->y;
  freq->gr = freq->g / freq->r;
  freq->ty = freq->t / freq->y;

  aa = ttratio * freq->r * freq->y - freqa * freqg - freqc * freqt;
  bb = freqa * freq->gr + freqc * freq->ty;
  freq->xi = aa / (aa + bb);
  freq->xv = 1.0 - freq->xi;
  if (freq->xi < 0.0)
  {
    freq->xi = 0.0;
    freq->xv = 1.0;
    ttratio = (freq->a*freq->g+freq->c*freq->t)/(freq->r*freq->y);
    ttratio_warning(ttratio);
  }
  if (freqa <= 0.0)
    freqa = 0.000001;
  if (freqc <= 0.0)
    freqc = 0.000001;
  if (freqg <= 0.0)
    freqg = 0.000001;
  if (freqt <= 0.0)
    freqt = 0.000001;

  freq->fracchange = freq->xi * (2 * freqa * freq->gr + 2 * freqc * freq->ty) +
    freq->xv *
    (1.0 - freqa*freqa - freqc*freqc
     - freqg*freqg - freqt*freqt );
  freq->ttratio = ttratio;
} /* makebasefreq */


void print_basefreq(FILE *fp, basefreq *freq, boolean empirical)
{
  /* print out empirical base frequencies */

  putc('\n', fp);
  if (empirical)
    fprintf(outfile, "Empirical ");
  fprintf(fp, "Base Frequencies:\n\n");
  fprintf(fp, "   A    %10.5f\n", freq->a);
  fprintf(fp, "   C    %10.5f\n", freq->c);
  fprintf(fp, "   G    %10.5f\n", freq->g);
  fprintf(fp, "  T(U)  %10.5f\n", freq->t);
} /* print_basefreq */


void ttratio_warning(double ttratio)
{
  /* print warning that this ttratio is impossible */

  printf("\n WARNING: This transition/transversion ratio\n"
         "  is impossible with these base frequencies!\n"
         "  Using a transition/transversion ratio of %.6f\n\n", ttratio);
} /* ttratio_warning */


void getbasefreqs(double freqa, double freqc, double freqg, double freqt,
                   double *freqr, double *freqy, double *freqar,
                   double *freqcy, double *freqgr, double *freqty,
                   double *ttratio, double *xi, double *xv,
                   double *fracchange, boolean freqsfrom, boolean printdata)
{
  /* Inputs freq[acgt] and ttratio and calculates freq[ry], freq[acgt][ry],
   * and fracchange.  If ttratio is impossible, a warning is printed and a more
   * reasonable value is returned. If printdata is true, the base frequencies
   * are printed to outfile. If freqsfrom is also true, the output is
   * identified as "Empirical". */

  /* Deprecated in favor of direct use of struct freq, makebasefreq(),
   * and print_basefreq() */

  /* used by dnadist, dnaml, & dnamlk */
  basefreq freq;
  boolean isempirical = freqsfrom;

  freq.a = freqa;
  freq.c = freqc;
  freq.g = freqg;
  freq.t = freqt;

  if (printdata)
  {
    print_basefreq(outfile, &freq, isempirical);
  }
  makebasefreq(&freq, freq.a, freq.c, freq.g, freq.t, *ttratio);

  *freqr = freq.r;
  *freqy = freq.y;
  *freqar = freq.ar;
  *freqcy = freq.cy;
  *freqgr = freq.gr;
  *freqty = freq.ty;
  *xi     = freq.xi;
  *xv     = freq.xv;
  *fracchange = freq.fracchange;
  *ttratio    = freq.ttratio;

}  /* getbasefreqs */


/* End. */

