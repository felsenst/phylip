/* Version 4.0. (c) Copyright 1999-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cont.h"
#include "ml.h"

#define SAMPLES 1000


node * cont_node_new(node_type type, long index) // RSGbugfix
{ /* make a new node of type cont_node_type */
  node *n = Malloc(sizeof(cont_node_type));
  cont_node_init((cont_node_type *)n, type, index);
  return n;
} /* cont_node_new */


void cont_node_init(cont_node_type* n, node_type type, long index)
{  /* make and init variables for a new node of type cont_node_type */
  generic_node_init((node*)n, type, index);
  ((node*)n)->copy = cont_node_copy;
} /* cont_node_init */


void cont_node_copy(node* srcn, node* dstn)
{ /* copy a node of type cont_node */
  cont_node_type *src, *dst;

  src = (cont_node_type *)srcn;
  dst = (cont_node_type *)dstn;
  generic_node_copy(srcn, dstn);
  if ( dst->totalleles != 0 && dst->totalleles != src->totalleles )
  {
    free(dst->view);
    dst->view = NULL;
  }
  dst->totalleles = src->totalleles;
  if ( dst->view == NULL )
    dst->view = Malloc(src->totalleles * sizeof(double));
  memcpy(dst->view, src->view, dst->totalleles * sizeof(double));
} /* cont_node_copy */


void alloctree(pointarray *treenode, long nonodes)
{
  /* allocate array treenode dynamically */
  /* used in contml, contrast and threshml */
  long i;

  *treenode = (pointarray)Malloc(nonodes * sizeof(node *));
  for (i = 0; i < spp; i++)
  {
    (*treenode)[i] = functions.node_new(TIP_NODE, i + 1);
  }
  for (i = spp; i < nonodes; i++)
  {
    /* node allocation handled in treeread which calls initcontrastnode or initthreshmlnode */
    (*treenode)[i] = NULL;
  }
} /* alloctree */


void freetree(pointarray *treenode, long nonodes)
{
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
    free((*treenode)[i]);
  for (i = spp; i < nonodes; i++)
  {
    p = (*treenode)[i]->next;
    /* cannot guarantee always triplets, so travel around the ring */
    while (p != (*treenode)[i])
    {
      q = p;
      p = p->next;
      free(q);
    }
    free((*treenode)[i]);
  }
  free(*treenode);   /* debug:  should we also free the free_fork_nodes  list?  */
} /* freetree */


void setuptree(tree *a, long nonodes)
{
  (void)nonodes;                        // RSGnote: Parameter never used.

  /* initialize a tree */
  /* used in contrast and threshml */
  long i;

  for (i = 1; i <= spp; i++)
  {
    a->nodep[i - 1]->back = NULL;
    a->nodep[i - 1]->iter = true;
  }

  /* no setup for interior nodes -- handled in treeread
     which calls initcontnode */
  a->score = -99999.0;
  a->root = a->nodep[0];
  a->get_fork = generic_tree_get_fork;

}  /* setuptree */


void allocview(tree *a, long nonodes, long totalleles)
{
  /* allocate view array
     used in contml */
  long i, n;
  node *r, *s;
  Slist_node_ptr p;
  Slist_data_ptr q;

  for (i = 0; i < spp; i++)
    if (a->nodep[i] != NULL)
      {
        ((cont_node_type*)a->nodep[i])->view = (phenotype3)Malloc(totalleles
                                                            * sizeof(double));
        ((cont_node_type*)a->nodep[i])->totalleles = totalleles;
      }

  for (i = spp; i < nonodes; i++)
  {
    if (a->nodep[i] != NULL) {
      r = (node*)(a->nodep[i]);
      s = r;
      do {        /* go around circle */
        ((cont_node_type*)s)->view = (phenotype3)Malloc(totalleles 
                                                         * sizeof(double));
        ((cont_node_type*)s)->totalleles = totalleles;
        s = s->next;
      } while (s != r);
    }
  }
  p = (Slist_node_ptr)(a->free_fork_nodes->first);   /* go along free nodes
                                                        list as needed */
  q = (Slist_data_ptr)(a->free_fork_nodes->first->data);
  n = a->free_fork_nodes->length;
  for (i = 1; i <= n; i++) {
    ((cont_node_type *)q)->view = (phenotype3)Malloc(totalleles
                                                      * sizeof(double));
    ((cont_node_type *)q)->totalleles = totalleles;
    if (i < n) {
      p = p->next;
      q = p->data;
      }
    }
}  /* allocview */


void freeview(tree *a, long nonodes)
{
  /* deallocate view */
  /* used in contml and threshml */
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
  {
    free(((cont_node_type*)a->nodep[i])->view);
  }
  for (i = spp; i < nonodes; i++)
  {
    p = a->nodep[i];
    q = p;
    do {
      free(((cont_node_type*)q)->view);
      q = q->next;
    }
    while (q != p);
  }
}  /* freeview */


void standev2(long numtrees, long maxwhich, long a, long b, double maxlogl, double *l0gl, double **l0gf, longer seed)
{  /* do paired sites test (KHT or SH) on user-defined trees */
  /* used in contml */
  double **covar, *P, *f, *r;
  long i, j, k;
  double sumw, sum, sum2, sd;
  double temp;

  if (numtrees == 2)
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    logL    Diff logL    Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    i = 1;
    while (i <= numtrees)
    {
      fprintf(outfile, "%3ld%10.1f", i, l0gl[i - 1]);
      if (maxwhich == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (j = a; j <= b; j++)
        {
          sumw += 1;
          temp = l0gf[i - 1][j] - l0gf[maxwhich - 1][j];
          sum += temp;
          sum2 += temp * temp;
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - temp * temp));
        fprintf(outfile, "%10.1f%12.4f", (l0gl[i - 1])-maxlogl, sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      i++;
    }
    fprintf(outfile, "\n\n");
  }
  else            /* Shimodaira-Hasegawa test using normal approximation */
  {
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n", MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    }
    else
    {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees * sizeof(double *));
    sumw = b-a+1;
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees * sizeof(double));
    for (i = 0; i < numtrees; i++)          /* compute covariances of trees */
    {
      sum = l0gl[i]/sumw;
      for (j = 0; j <=i; j++)
      {
        sum2 = l0gl[j]/sumw;
        temp = 0.0;
        for (k = a; k <= b ; k++)
        {
          temp = temp + (l0gf[i][k]-sum)*(l0gf[j][k]-sum2);
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++)  /* in-place Cholesky decomposition of
                                       trees x trees covariance matrix */
    {
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++)
      {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees * sizeof(double)); /* resampled likelihoods */
    P = (double *)Malloc(numtrees * sizeof(double)); /* vector: P's of trees */
    r = (double *)Malloc(numtrees * sizeof(double)); /* put Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    for (i = 1; i <= SAMPLES; i++)           /* loop over resampled trees */
    {
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++)        /* compute vectors */
      {
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (maxlogl-l0gl[j] < sum-f[j])
          P[j] += 1.0 / SAMPLES;
    }
    fprintf(outfile, "Tree    logL    Diff logL    P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++)
    {
      fprintf(outfile, "%3ld%10.2f", i+1, l0gl[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        fprintf(outfile, " %9.2f  %10.2f", l0gl[i]-maxlogl, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
    fprintf(outfile, "\n");
    free(P);             /* free the variables we Malloc'ed */
    free(f);
    free(r);
    for (i = 0; i < numtrees; i++)
      free(covar[i]);
    free(covar);
  }
}  /* standev2 */


// End.
