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


struct node* dist_node_new(node_type type, long index, long nodesize)
{
  /* make a new dist_node. */
  struct node* dn;

  dn =  ml_node_new(type, index, nodesize);        /* call upwards */
  return dn;
} /* dist_node_new */


void dist_node_init(struct node* n, node_type type, long index)
{
  /* initialize a new dist_node */
  struct dist_node* dn;

  ml_node_init(n, type, index);                  /* go up node hierarchy */
  n->free = dist_node_free;
  n->copy = dist_node_copy;
/* debug:  needed?    n->init = dist_node_init; */
  dn = (struct dist_node*)n;
  dn->dist = 0.0;                              /* used in  fitch, kitsch */
  dn->d = (vector)Malloc((nonodes+1) * sizeof(double));    /* extra room */
  dn->w = (vector)Malloc((nonodes+1) * sizeof(double));
} /* dist_node_init */


void dist_tree_init(struct tree *t, long nonodes, long spp)
{
  /* initialize a dist_tree.
   * used in fitch, kitsch, & neighbor
   * acts after phylip.c: generic_tree_init 
   * runs along nodep, circling each fork, and then runs
   * along  free_fork_nodes  list too, in each case calling dist_node_init */
 /* debug:  For that matter could we make a generic function that does running-along?
 * much of this function is supposed to be done in the generic_tree_init !! */
  long i;
  struct dist_node* p;
  struct node* pp;
  Slist_node_ptr q;

/* debug:  extra?  (yes)   ml_tree_init(mlt, nonodes+1, spp);         */
  for (i = 1; i <= nonodes; i++) {                       /* note  nonodes+1 */
    if (t->nodep[i - 1] != NULL) {           /* these will be NULL normally */
      t->nodep[i - 1]->back = NULL;  /* debug: why bother? */
      t->nodep[i - 1]->iter = true; /* debug: at what level initiate this? */
/* debug: needed? which type of node?      ((tree*)t)->nodep[i - 1]->t = 0.0; */
      ((dist_node*)(t->nodep[i - 1]))->sametime = false;
      t->nodep[i - 1]->v = 0.0;
      if (i > spp) {  /* go around fork circles initializing node variables */
        pp = t->nodep[i-1]->next;
        p = (dist_node*)pp;
        while (pp != t->nodep[i-1]) {  /* until get to where you entered */
          pp->back = NULL;   /* debug:  why bother? */
          pp->iter = true;   /* debug: where should  iter  be initialized?  ml_node_init? */
          p->t = 0.0;   /* debug: initialize where? */
          p->sametime = false;   /* debug: initialize where? */
          pp = pp->next;                /* go to next  node  in fork circle */
          p = (dist_node*)pp; /* make sure generic_node version is same one */
        }
      }
    }
  }
  q = t->free_fork_nodes->first;  /* go along list of forknode items too */
  while (q != NULL) {              
    p = q->data;                        /* p  is now the node that is there */
/* debug:     dist_node_init((dist_node*)p, FREE_NODE, 0);  wasn't this done from node_new? */
    q = q->next;                                    /* go to next list item */
  }
  t->score = -1.0;    /* debug: set this in  ml_tree_init? */
  t->root = t->nodep[0];   /* debug: set in which function? */
}  /* dist_tree_init */


void dist_tree_new(struct tree** treep, long nonodes, long spp, long treesize)
{ /* make a new pointer to dist_tree (which is itself a pointer to the tree
   * structure.  Calls to ml_tree_new, which calls up to generic_tree_new,
   * then after that it calls dist_tree_init */
  long tsz;

  tsz = (long)treesize;   /* debug:  need this cast or tsz variable? */
  ml_tree_new(treep, nonodes, spp, tsz); 
  dist_tree_init(*treep, nonodes, spp);
} /* dist_tree_new */


void dist_node_free(node **np)
{
  /* free a dist_node */
  dist_node *n = (dist_node *)*np;

  free(n->d);
  free(n->w);
  generic_node_free(np);
} /* dist_node_free */


void dist_node_copy(node* srcn, node* dstn)
{
  /* copy a dist_node */
   /* debug: copy functions too?  */
  dist_node *src = (dist_node *)srcn;
  dist_node *dst = (dist_node *)dstn;

  if ((dstn != NULL) && (srcn != NULL)) {
    generic_node_copy(srcn, dstn);
    dst->dist = src->dist;
    memcpy(dst->d, src->d, (nonodes+1) * sizeof(double));
    memcpy(dst->w, src->w, (nonodes+1) * sizeof(double));
    dst->sametime = src->sametime;
    dst->t = src->t;
  }
} /* dist_node_copy */


// void alloctree(tree *t, long nonodes)
// {
  /* allocate treenode dynamically
   * formerly used in fitch, kitsch & neighbor */
  /* debug: not used now in most programs*/
//  long i, j;
//  node *p = NULL, *q = NULL; /* debug to avoid warning message */
//
//  t->nodep = (node**)Malloc(nonodes * sizeof(node *));
//  for (i = 0; i < spp; i++) {
//    funcs.node_new(p, TIP_NODE, i+1, sizeof(dist_node));
//    t->nodep[i] = p;
//  }
//  for (i = spp; i < nonodes; i++)
//  {
//    q = NULL;
//    for (j = 1; j <= 3; j++) {
//      funcs.node_new(p, FORK_NODE, i+1, sizeof(dist_node));
//      t->nodep[i] = p;
//      p->next = q;
//      q = p;
//    }
//    p->next->next->next = p;
//    t->nodep[i] = p;
//  }
//} /* alloctree */


void allocd(long nonodes, pointptr treenode)
{
  /* formerly used in fitch & kitsch.  Not used now  */
  long i, j;
  dist_node *p;
  dist_node** dtreenode = (dist_node**)treenode;

  for (i = 0; i < spp; i++)
  {
    dtreenode[i]->d = (vector)Malloc((nonodes+1) * sizeof(double));
  }
  for (i = spp; i < nonodes; i++) {
    p = dtreenode[i];
    for (j = 1; j <= 3; j++) {
      p->d = (vector)Malloc((nonodes+1) * sizeof(double));
      p = (dist_node*)(((node*)p)->next);
    }
  }
} /* allocd */


void freed(long nonodes, pointptr treenode)
{
  /* formerly used in fitch, not used now */
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
      p = (dist_node*)(((node*)p)->next);
    }
  }
} /* freed */


void allocw(long nonodes, pointptr treenode)
{
  /* allocate weights array. used in fitch & kitsch */
/* debug:  formerly used in fitch, kitsch and neighbor, not used now */
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
      p = (dist_node*)(((node*)p)->next);
    }
  }
} /* allocw */


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
      p = (dist_node*)(((node*)p)->next);
    }
  }
} /* freew */


void inputdata(boolean replicates, boolean printdata, boolean lower,
               boolean upper, vector *x, intvector *reps)
{
  /* read in distance matrix
   * used in fitch, kitsch & neighbor */
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
        printf("        (%ld,%ld) element and (%ld,%ld) element are unequal.\n",
                 i+1, j+1, j+1, i+1);
        printf("        They are %10.6f and %10.6f, respectively.\n",
                 x[i][j], x[j][i]);
        printf("        Is it a distance matrix?\n\n");
        exxit(-1);
      }
    }
  }
  scan_eoln(infile);
  checknames(spp);       /* Check  nayme  array for duplicate species names */
  if (!printdata)
    return;
  for (i = 0; i < spp; i++) {         /* if printing out the input data too */
    for (j = 0; j < nmlngth; j++)           /* print the name of the species */
      putc(nayme[i][j], outfile);
    putc(' ', outfile);                                      /* then a blank */
    for (j = 1; j <= spp; j++) {
      fprintf(outfile, "%10.5f", x[i][j - 1]);         /* then its distances */
      if (replicates)        /* and, if relevant, number of replicates of it */
        fprintf(outfile, " (%3ld)", reps[i][j - 1]);
      if (j % columns == 0 && j < spp) { /* ... go to a new column as needed */
        putc('\n', outfile);
        for (k = 1; k <= nmlngth + 1; k++)    /* with blanksinstead of names */
          putc(' ', outfile);
      }
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* inputdata */


/* debug: can the following functions be promoted to  phylip.c ? */
void coordinates(node *p, double lengthsum, long *tipy, double *tipmax, node *start)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  int nover;

  if (p->tip) {                                         /* if at a tip ... */
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > *tipmax)
      *tipmax = lengthsum;
    return;
  }
  q = p->next;                    /* ... otherwise are at an interior node */
  do {
    if (q->back)
      coordinates(q->back, lengthsum + q->v, tipy, tipmax, start);
    q = q->next;
  } while ((p == start || p != q) && (p != start || p->next != q));
  q = p;
  do {                                                /* find leftmost furc */
    q = q->next; 
    first = q->back;
  } while (first == NULL);
  q = p;
  while ((q->next != p) && q->next->back)  /* debug: is this right ? */
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
      if (q->back != NULL)
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
  long n=0, j=0;
  boolean extra=false, trif=false, done=false;
  node *p, *q, *r, *s, *first =NULL, *last =NULL;

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
      } while (!(((p != start) && (r == p)) ||
                 ((p == start) && (r == p->next))));
      s = p;
      do {
        s = s->next;
      } while (s->back == NULL);
      first = s->back;
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
      if (((long)last->ycoord > i) && ((long)first->ycoord < i)
          && (i != (long)p->ycoord)) {
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
  /* write out file with representation of final tree recursively.
   * used in fitch & neighbor */
  /* debug:  Needs code for indenting. */
  /* debug:  Don't need the njoin argument, no longer use it */
  long i=0, n=0, w=0;
  Char c;
  double x=0.0;
  node *q;

  if (p->tip) {                               /* if at a tip, write out the name there */
    n = 0;
    for (i = 1; i <= nmlngth; i++) {                    /* its nonblank characters ... */
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {                         /* ... followed by enough blanks */
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  } else {                                                      /* if at a fork circle */
    q = p;                                           /* keep track of where we entered */
    putc('(', outtree);
    (*col)++;
    do {
      p = p->next;
      if (p != q) {
        if (p->back != NULL) {
          treeout(p->back, col, m, njoin, start);                 /* write the subtree */
          if (p->next != q) {              /* don't put in a comma if was last subtree */
            putc(',', outtree);
            (*col)++;
          }
          if (*col > 55) {
            putc('\n', outtree);
            *col = 0;
          }
        }
      }
    } while (p != q);
    if ((p->index == start->index) && (p->back != NULL)) {
      treeout(p->back, col, m, njoin, start);
    }
    putc(')', outtree);
    (*col)++;
  }
  x = p->v;
  if (x > 0.0)          /* computing width of number needed to print out branch length */
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
    fprintf(outtree, ":%*.5f", (int) w + 7, x);         /* print out the branch length */
    *col += w + 8;
  }
}  /* treeout */


/* End. */
