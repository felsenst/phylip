/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "dist.h"

extern long nonodes;


node* dist_node_new(node_type type, long index)
{
  node *n = Malloc(sizeof(dist_node));
  dist_node_init(n, type, index);
  return n;
}


void dist_node_init(node* n, node_type type, long index)
{
  dist_node * dn = (dist_node *)n;

  generic_node_init(n, type, index);
  n->free = dist_node_free;
  n->copy = dist_node_copy;
  n->init = dist_node_init;
  dn->dist = 0;
  dn->d = (vector)Malloc(nonodes * sizeof(double));
  dn->w = (vector)Malloc(nonodes * sizeof(double));
}


void dist_node_free(node **np)
{
  dist_node *n = (dist_node *)*np;
  free(n->d);
  free(n->w);
  generic_node_free(np);
}


void dist_node_copy(node* srcn, node* dstn)
{
  dist_node *src = (dist_node *)srcn;
  dist_node *dst = (dist_node *)dstn;

  generic_node_copy(srcn, dstn);
  dst->dist = src->dist;
  memcpy(dst->d, src->d, nonodes * sizeof(double));
  memcpy(dst->w, src->w, nonodes * sizeof(double));
  dst->sametime = src->sametime;
  dst->t = src->t;
}


void alloctree(pointptr *treenode, long nonodes)
{
  /* allocate treenode dynamically */
  /* used in fitch, kitsch & neighbor */
  long i, j;
  node *p, *q;

  *treenode = (pointptr)Malloc(nonodes * sizeof(node *));
  for (i = 0; i < spp; i++)
    (*treenode)[i] = functions.node_new(TIP_NODE, i+1);
  for (i = spp; i < nonodes; i++)
  {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = functions.node_new(FORK_NODE, i+1);
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
} /* alloctree */


void freetree(pointptr *treenode, long nonodes)
{
  long i;
  node *p, *q;

  for (i = 0; i < spp; i++)
    (*treenode)[i]->free(&(*treenode)[i]);
  for (i = spp; i < nonodes; i++) {
    p = (*treenode)[i];
    q = p->next;
    while(q != p)
    {
      node * r = q;
      q = q->next;
      r->free(&r);
    }
    p->free(&p);
  }
  free(*treenode);
} /* freetree */


void allocd(long nonodes, pointptr treenode)
{
  /* used in fitch & kitsch */
  long i, j;
  dist_node *p;
  dist_node** dtreenode = (dist_node**)treenode;

  for (i = 0; i < spp; i++)
  {
    dtreenode[i]->d = (vector)Malloc(nonodes * sizeof(double));
  }
  for (i = spp; i < nonodes; i++) {
    p = dtreenode[i];
    for (j = 1; j <= 3; j++) {
      p->d = (vector)Malloc(nonodes * sizeof(double));
      p = (dist_node*)(p->node.next);
    }
  }
}


void freed(long nonodes, pointptr treenode)
{
  /* used in fitch */
  long i, j;
  dist_node *p;
  dist_node **dtreenode = (dist_node**)treenode;

  for (i = 0; i < spp; i++)
  {
    free(dtreenode[i]->d);
  }
  for (i = spp; i < nonodes; i++) {
    p = dtreenode[i];
    for (j = 1; j <= 3; j++) {
      free(((dist_node*)p)->d);
      p = (dist_node*)(p->node.next);
    }
  }
}


void allocw(long nonodes, pointptr treenode)
{
  /* used in fitch & kitsch */
  long i, j;
  dist_node *p;
  dist_node **dtreenode = (dist_node**)treenode;

  for (i = 0; i < spp; i++)
  {
    dtreenode[i]->w = (vector)Malloc(nonodes * sizeof(double));
  }
  for (i = spp; i < nonodes; i++) {
    p = dtreenode[i];
    for (j = 1; j <= 3; j++) {
      p->w = (vector)Malloc(nonodes * sizeof(double));
      p = (dist_node*)p->node.next;
    }
  }
}


void freew(long nonodes, pointptr treenode)
{
  /* free weights array
   * used in fitch */
  long i, j;
  dist_node *p;
  dist_node **dtreenode = (dist_node**)treenode;

  for (i = 0; i < spp; i++)
  {
    free(dtreenode[i]->w);
  }
  for (i = spp; i < nonodes; i++) {
    p = dtreenode[i];
    for (j = 1; j <= 3; j++) {
      free(p->w);
      p = (dist_node*)p->node.next;
    }
  }
} /* freew */


void setuptree(tree *a, long nonodes)
{
  /* initialize a tree
   * used in fitch, kitsch, & neighbor */
  long i=0;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    if (i > spp) {
      a->nodep[i-1] = generic_tree_get_fork(a, i-1); 
    }
    a->nodep[i - 1]->back = NULL;
    a->nodep[i - 1]->iter = true;
    ((dist_node*)a->nodep[i - 1])->t = 0.0;
    ((dist_node*)a->nodep[i - 1])->sametime = false;
    a->nodep[i - 1]->v = 0.0;
    if (i > spp) {
      p = a->nodep[i-1]->next;
      while (p != a->nodep[i-1]) {
        p->back = NULL;
        p->iter = true;
        ((dist_node*)p)->t = 0.0;
        ((dist_node*)p)->sametime = false;
        p = p->next;
      }
    }
  }
  /* Create garbage lists */ 
  a->free_forks = Slist_new();      /* debug:  Now unnecessary? */
  a->free_fork_nodes = Slist_new();
      
  /* Put all interior nodes on garbage lists by "releasing" them */
  for ( i = nonodes - 1 ; i >= spp ; i-- ) {
    a->release_fork(a, a->nodep[i]);
  } 
  a->score = -1.0;
  a->root = a->nodep[0];
}  /* setuptree */


void inputdata(boolean replicates, boolean printdata, boolean lower,
               boolean upper, vector *x, intvector *reps)
{
  /* read in distance matrix */
  /* used in fitch & neighbor */
  long i=0, j=0, k=0, columns=0;
  boolean skipit=false, skipother=false;

  if (replicates)
    columns = 4;
  else
    columns = 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  for (i = 0; i < spp; i++) {
    x[i][i] = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 0; j < spp; j++) {
      skipit = ((lower && j + 1 >= i + 1) || (upper && j + 1 <= i + 1));
      skipother = ((lower && i + 1 >= j + 1) || (upper && i + 1 <= j + 1));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][j]) != 1) {
          printf("The infile is of the wrong type\n");
          exxit(-1);
        }
        if (replicates)
        {
          if (eoln(infile))
            scan_eoln(infile);
          if(fscanf(infile, "%ld", &reps[i][j]) < 1)
          {
            printf("\nERROR reading input file.\n\n");
            exxit(-1);
          }
        }
        else
          reps[i][j] = 1;
      }
      if (!skipit && skipother) {
        x[j][i] = x[i][j];
        reps[j][i] = reps[i][j];
      }
      if ((i == j) && (fabs(x[i][j]) > 0.000000001)) {
        printf("\nERROR:  Diagonal element of row %ld of distance matrix ", i+1);
        printf("is not zero.\n");
        printf("        Is it a distance matrix?\n\n");
        exxit(-1);
      }
      if ((j < i) && (fabs(x[i][j]-x[j][i]) > 0.000000001)) {
        printf("ERROR:  Distance matrix is not symmetric:\n");
        printf("        (%ld,%ld) element and (%ld,%ld) element are unequal.\n", i+1, j+1, j+1, i+1);
        printf("        They are %10.6f and %10.6f, respectively.\n", x[i][j], x[j][i]);
        printf("        Is it a distance matrix?\n\n");
        exxit(-1);
      }
    }
  }
  scan_eoln(infile);
  checknames(spp);                      // Check NAYME array for duplicates.
  if (!printdata)
    return;
  for (i = 0; i < spp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[i][j], outfile);
    putc(' ', outfile);
    for (j = 1; j <= spp; j++) {
      fprintf(outfile, "%10.5f", x[i][j - 1]);
      if (replicates)
        fprintf(outfile, " (%3ld)", reps[i][j - 1]);
      if (j % columns == 0 && j < spp) {
        putc('\n', outfile);
        for (k = 1; k <= nmlngth + 1; k++)
          putc(' ', outfile);
      }
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


void coordinates(node *p, double lengthsum, long *tipy, double *tipmax, node *start)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  int nover;

  if (p->tip) {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > *tipmax)
      *tipmax = lengthsum;
    return;
  }
  q = p->next;
  do {
    if (q->back)
      coordinates(q->back, lengthsum + q->v, tipy, tipmax, start);
    q = q->next;
  } while ((p == start || p != q) && (p != start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p && q->next->back)  /* is this right ? */
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == start) {
    nover = count_sibs(p);
    if ( p->back != NULL ) nover++;
    if ( (nover % 2) == 0 ) {
      p->ycoord = (first->ycoord + last->ycoord) / 2;
    } else {
      nover = (nover / 2) + 1;
      q = p;
      while (nover) {
        q = q->next;
        nover--;
      }
      p->ycoord = q->back->ycoord;
    }
  }
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(long i, double scale, node *start, boolean rooted)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n=0, j=0;
  boolean extra=false, trif=false;
  node *r, *first =NULL, *last =NULL;
  boolean done=false;

  p = start;
  q = start;
  extra = false;
  trif = false;
  if (i == (long)p->ycoord && p == start) {  /* display the root */
    if (rooted) {
      if (p->index - spp >= 10)
        fprintf(outfile, "-");
      else
        fprintf(outfile, "--");
    }
    else {
      if (p->index - spp >= 10)
        fprintf(outfile, " ");
      else
        fprintf(outfile, "  ");
    }
    if (p->index - spp >= 10)
      fprintf(outfile, "%2ld", p->index - spp);
    else
      fprintf(outfile, "%ld", p->index - spp);
    extra = true;
    trif = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) { /* internal nodes */
      r = p->next;
      /* r->back here is going to the same node. */
      do {
        if (!r->back) {
          r = r->next;
          continue;
        }
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          break;
        }
        r = r->next;
      } while (!((p != start && r == p) || (p == start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (!rooted && (p == start))
        last = p->back;
    } /* end internal node case... */
    /* draw the line: */
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (!q->tip) {
      if ((n < 3) && (q->index - spp >= 10))
        n = 3;
      if ((n < 2) && (q->index - spp < 10))
        n = 2;
    }
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      if (trif) {
        n++;
        trif = false;
      }
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i
          && i != (long)p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
        trif = false;
      }
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(node *start, boolean treeprint, boolean rooted)
{
  /* prints out diagram of the tree */
  /* used in fitch & neighbor */
  long i;
  long tipy;
  double scale, tipmax;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(start, 0.0, &tipy, &tipmax, start);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, start, rooted);
  putc('\n', outfile);
}  /* printree */


void treeoutr(node *p, long *col, tree *curtree)
{
  /* write out file with representation of final tree.
   * Rooted case. Used in kitsch and neighbor. */
  long i, n, w;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    (*col) += n;
  } else {
    putc('(', outtree);
    (*col)++;
    treeoutr(p->next->back, col, curtree);
    putc(',', outtree);
    (*col)++;
    if ((*col) > 55) {
      putc('\n', outtree);
      (*col) = 0;
    }
    treeoutr(p->next->next->back, col, curtree);
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree->root)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    (*col) += w + 8;
  }
}  /* treeoutr */


void treeout(node *p, long *col, double m, boolean njoin, node *start)
{
  /* write out file with representation of final tree */
  /* used in fitch & neighbor */
  long i=0, n=0, w=0;
  Char c;
  double x=0.0;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    treeout(p->next->back, col, m, njoin, start);
    putc(',', outtree);
    (*col)++;
    if (*col > 55) {
      putc('\n', outtree);
      *col = 0;
    }
    treeout(p->next->next->back, col, m, njoin, start);
    if (p == start && njoin) {
      putc(',', outtree);
      treeout(p->back, col, m, njoin, start);
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)
    w = (long)(m * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(m * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == start)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int) w + 7, x);
    *col += w + 8;
  }
}  /* treeout */


// End.
