/* Version 4.0. Copyright 2021.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "seq.h"

#define SAMPLES 1000

extern long endsite, outgrno, which;
extern FILE *infile, *outfile, *outtree;

steptr weight, category, alias, location, ally;
boolean interleaved, printdata, outgropt, treeprint, dotdiff;


void inputdata(long chars)
{
  /* read in and if needed print out sequences */
  read_sequences(chars);
  if (printdata)
  {
    output_sequences(chars);
  }
}  /* inputdata */


void read_sequences(long nchars)
{
  /* input the names and sequences for each species
   * used by dnacomp, dnadist, dnainvar, dnaml, dnamlk, dnapars, & dnapenny */
  long i, j, basesread, basesnew=0;
  Char charstate;
  boolean allread, done;

  basesread = 0;
  allread = false;
  while (!(allread))
  {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp)
    {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = (interleaved) ? basesread : 0;
      done = false;
      while (!done && !eoff(infile))
      {
        if (interleaved)
          done = true;
        while (j < nchars && !(eoln(infile) || eoff(infile)))
        {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((strchr("ABCDGHKMNRSTUVWXY?O-", charstate)) == NULL)
          {
            printf("\nERROR:  Bad base: %c at site %5ld of species %3ld.\n",
                    charstate, j+1, i);
            if (charstate == '.')
            {
              printf(
               "        Periods (.) may not be used as gap characters.\n");
              printf("        The correct gap character is (-).\n");
            }
            exxit(-1);
          }
          j++;
          inputSequences[i - 1][j - 1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < nchars)
          scan_eoln(infile);
        else if (j == nchars)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;

      scan_eoln(infile);

      if ((interleaved && j != basesnew) || (!interleaved && j != nchars))
      {
        printf("\nERROR:  Sequences out of alignment at position %ld", j+1);
        printf(" of species %ld.\n\n", i);
        exxit(-1);
      }
      i++;
    }

    if (interleaved)
    {
      basesread = basesnew;
      allread = (basesread == nchars);
    }
    else
      allread = (i > spp);
  }
  checknames(spp);                      // Check NAYME array for duplicates.
} /* read_sequences */


void output_sequences(long nchars)
{
  /* print out sequences */
  long i, j, k, l;
  char charstate;

  headings(nchars, "Sequences", "---------");

  for (i = 1; i <= ((nchars - 1) / 60 + 1); i++)
  {
    for (j = 1; j <= spp; j++)
    {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j - 1][k], outfile);
      fprintf(outfile, "   ");
      l = i * 60;
      if (l > nchars)
        l = nchars;
      for (k = (i - 1) * 60 + 1; k <= l; k++)
      {
        if (dotdiff &&
            ((j > 1) &&
             inputSequences[j - 1][k - 1] == inputSequences[0][k - 1]))
          charstate = '.';
        else
          charstate = inputSequences[j - 1][k - 1];
        putc(charstate, outfile);
        if (k % 10 == 0 && k % 60 != 0)
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
} /* output_sequences */


void setuptree(pointarray treenode, long nonodes, boolean usertree)
{
  /* initialize treenodes */
  long i;
  struct node *p;

  for (i = 1; i <= nonodes; i++)
  {
    if (i <= spp || !usertree)
    {
      treenode[i-1]->back = NULL;
      treenode[i-1]->node_init(treenode[i-1], i<=spp, i);
/* debug:      treenode[i-1]->iter = true;    maybe initialize later? */
      treenode[i-1]->initialized = true;
    }
  }
  if (!usertree)
  {
    for (i = spp + 1; i <= nonodes; i++)
    {
      p = treenode[i-1]->next;
      while (p != treenode[i-1])
      {
        p->back = NULL;
        p->node_init(p, false, i);
/* debug         p->iter = true;   maybe initialize later? */
        p->initialized = false;
        p = p->next;
      }
    }
  }
} /* setuptree */


void freetrans(transptr *trans, long nonodes, long sitelength)
{
  /* free array trans */
  long i, j;
  for ( i = 0 ; i < nonodes ; i++ )
  {
    for ( j = 0 ; j < sitelength + 1; j++)
    {
      free ((*trans)[i][j]);
    }
    free ((*trans)[i]);
  }
  free(*trans);
} /* freetrans */


void sitesort(long chars, steptr weight)
{
  /* Shell sort keeping sites, weights in same order */
  /* used in dnainvar, dnapars, dnacomp & dnapenny */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = chars / 2;
  while (gap > 0)
  {
    for (i = gap + 1; i <= chars; i++)
    {
      j = i - gap;
      flip = true;
      while (j > 0 && flip)
      {
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = true;
        k = 1;
        while (k <= spp && tied)
        {
          flip = (inputSequences[k - 1][jj - 1] >
                   inputSequences[k - 1][jg - 1]);
          tied = (tied && inputSequences[k - 1][jj - 1]
                           == inputSequences[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = weight[j - 1];
        weight[j - 1] = weight[j + gap - 1];
        weight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort */


void sitecombine(long chars)
{
  /* combine sites that have identical patterns */
  /* used in dnapars, dnapenny, & dnacomp */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < chars)
  {
    j = i + 1;
    tied = true;
    while (j <= chars && tied)
    {
      k = 1;
      while (k <= spp && tied)
      {
        tied = (tied && inputSequences[k - 1][alias[i - 1] - 1] == inputSequences[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (tied)
      {
        weight[i - 1] += weight[j - 1];
        weight[j - 1] = 0;
        ally[alias[j - 1] - 1] = alias[i - 1];
      }
      j++;
    }
    i = j - 1;
  }
}  /* sitecombine */


void sitescrunch(long chars)
{
  /* move so one representative of each pattern of sites comes first */
  /* used in dnapars & dnacomp */
  long i, j, itemp;
  boolean done, found;

  done = false;
  i = 1;
  j = 2;
  while (!done)
  {
    if (ally[alias[i - 1] - 1] != alias[i - 1])
    {
      if (j <= i)
        j = i + 1;
      if (j <= chars)
      {
        do {
          found = (ally[alias[j - 1] - 1] == alias[j - 1]);
          j++;
        } while (!(found || j > chars));
        if (found)
        {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = weight[i - 1];
          weight[i - 1] = weight[j - 1];
          weight[j - 1] = itemp;
        }
        else
          done = true;
      }
      else
        done = true;
    }
    i++;
    done = (done || i >= chars);
  }
}  /* sitescrunch */


void sitesort2(long sites, steptr aliasweight)
{
  /* Shell sort keeping sites, weights in same order
     alias[i-1] is the number of the site that comes ith in lexicographical
     order by rate category and pattern (sites with category 1 and pattern
     AAA ... AAA are first, etc.)
     aliasweight[i-1] will end up as the weight of all sites that have
     the same pattern as the site that is ith in lexicographical order */
  /* used in dnaml, dnamlk, proml, promlk, and codml */
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
        jj = alias[j - 1];
        jg = alias[j + gap - 1];
        tied = (category[jj - 1] == category[jg - 1]);
        flip = (category[jj - 1] > category[jg - 1]);
        k = 1;
        while (k <= spp && tied)
        {
          flip = (inputSequences[k - 1][jj - 1] > inputSequences[k - 1][jg - 1]);
          tied = (tied && inputSequences[k - 1][jj - 1] == inputSequences[k - 1][jg - 1]);
          k++;
        }
        if (!flip)
          break;
        itemp = alias[j - 1];
        alias[j - 1] = alias[j + gap - 1];
        alias[j + gap - 1] = itemp;
        itemp = aliasweight[j - 1];
        aliasweight[j - 1] = aliasweight[j + gap - 1];
        aliasweight[j + gap - 1] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* sitesort2 */


void sitecombine2(long sites, steptr aliasweight)
{
  /* combine sites that have identical rate categories and patterns
   * by making one of them which "represents" them have the total weight
   * and the rest have weight 0.  They are tied if they have the same
   * rate category number and the same site pattern
   * ally[i-1] ends up as the number of the site (in the original site
   * numbering that represents the site which is  ith  in its original order.
   * used in dnaml, dnamlk, proml, promlk, and codml */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites)
  {
    j = i + 1;
    tied = true;
    while (j <= sites && tied)
    {
      tied = (category[alias[i - 1] - 1] == category[alias[j - 1] - 1]);
      k = 1;
      while (k <= spp && tied)
      {
        tied = (tied && inputSequences[k - 1][alias[i - 1] - 1]
                         == inputSequences[k - 1][alias[j - 1] - 1]);
        k++;
      }
      if (!tied)
        break;
      aliasweight[i - 1] += aliasweight[j - 1];
      aliasweight[j - 1] = 0;
      ally[alias[j - 1] - 1] = alias[i - 1];
      j++;
    }
    i = j;
  }
}  /* sitecombine2 */


void sitescrunch2(long sites, long i, long j, steptr aliasweight)
{
  /* move so positively weighted sites come first
   * used by dnaml, dnamlk, proml, promlk, codml, dnainvar, and restml */
  long itemp;
  boolean done, found;

  done = false;
  while (!done)
  {
    if (aliasweight[i - 1] > 0)
      i++;
    else
    {
      if (j <= i)
        j = i + 1;
      if (j <= sites)
      {
        do {
          found = (aliasweight[j - 1] > 0);
          j++;
        } while (!(found || j > sites));
        if (found)
        {
          j--;
          itemp = alias[i - 1];
          alias[i - 1] = alias[j - 1];
          alias[j - 1] = itemp;
          itemp = aliasweight[i - 1];
          aliasweight[i - 1] = aliasweight[j - 1];
          aliasweight[j - 1] = itemp;
        }
        else
          done = true;
      }
      else
        done = true;
    }
    done = (done || i >= sites);
  }
}  /* sitescrunch2 */


void drawline(long i, double scale, struct bl_node *rt)
{
  /* draws one row of the tree diagram by moving up tree */
  struct node *root, *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done, noplus;

  root = (struct node*)rt;
  p = root;
  q = root;
  assert(p->index > 0);                 // RSGdebug

  extra = false;
  noplus = false;
  if (i == (long)p->ycoord && p == root)
  {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
    noplus = true;
  }
  else
    fprintf(outfile, "  ");
  do {
    if (!p->tip)
    {
      r = p->next;
      done = false;
      do {
        if (r->back != NULL) {
          if ((i >= r->back->ymin) && (i <= r->back->ymax))
          {
            q = r->back;
            done = true;
          }
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if (((long)q->ycoord == i) && !done)
    {
      if (noplus)
      {
        putc('-', outfile);
        noplus = false;
      }
      else
        putc('+', outfile);
      if (!q->tip)
      {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        assert(q->index > 0);           // RSGdebug
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
        noplus = true;
      }
      else
      {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    }
    else if (!p->tip)
    {
      if (((long)last->ycoord > i) && ((long)first->ycoord < i) 
          && (i != (long)p->ycoord))
      {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      }
      else
      {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
      noplus = false;
    }
    else
    {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
      noplus = false;
    }
    if (p != q)
      p = q;
  } while (!done);
  if (((long)p->ycoord == i) && p->tip)
  {
    assert(p->index > 0);               // RSGdebug
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void treeout(struct node *p, long nextree, long *col, struct node *root)
{
  /* write out file with representation of final tree
   * used in dnacomp, dnamove, dnapars, & dnapenny */
  struct node *q;
  long i, n;
  Char c;

  assert(p->index > 0);                 // RSGdebug

  if (p->tip)
  {
    n = 0;
    for (i = 1; i <= nmlngth; i++)
    {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++)
    {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    *col += n;
  }
  else
  {
    putc('(', outtree);
    (*col)++;
    q = p->next;
    while (q != p)
    {
      treeout(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 60)
      {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout */


void drawline2(long i, double scale, struct node *p, struct tree* curtree)
{
	/* debug:  the newer version, older one follows */
  /* draws one row of the tree diagram by moving up tree
   * the argument  i  is the vertical number (y) of the row we draw,
   * numbered from top (1) to bottom
   * used in Dnaml, Proml, & Restml */

  struct node *r, *q;
  long n, j;
  boolean itoleft, iequal, iinsubtree, iatitsroot;
  boolean printedbar, done;

  itoleft = i < (long)p->ycoord;         /* Is  i  to left, right or at ... */
  iequal = i == (long)p->ycoord;               /* ... the coordinate of  p  */                
  q = curtree->root;
  if (q->tip)
    q = curtree->root->back;
  if (iequal && p->tip) {                           /* if now at a tip, ... */
    for (j = 0; j < nmlngth; j++)                 /* ... write the name ... */
      putc(nayme[p->index-1][j], outfile);
    return;             /* exit: we're all done after printing species name */
  }
  if (iequal) {                           /* if at an interior node instead */
    if (p->index - spp >= 100)           /* print out a number for the node */
      fprintf(outfile, "%3ld", p->index - spp);
    else {  
      if (p->index - spp >= 10)
        fprintf(outfile, "-%2ld", p->index - spp);
      else
        fprintf(outfile, "--%ld", p->index - spp);
    }
  }
  else {
      fprintf(outfile, "  ");              /* start by indenting two spaces */
  }
  if ((p->back != 0) && (p == q))     /* if at root and nonempty descendant */
     r = p;
  else                                /* otherwise move to first descendant */
     r = p->next;
  done = false;
  printedbar = false;         /* not (yet) printed a vertical bar character */
  do {  /* now check for each of  p's  descendants if  i  is in subtree ... */
    n = (long)(scale * ((long)r->back->xcoord - (long)p->xcoord) + 0.5);
    iinsubtree = (i >= r->back->ymin) && (i <= r->back->ymax);
    if (iinsubtree) {
      iatitsroot = (i == (long)r->back->ycoord);
      if (iatitsroot) {
        if (itoleft)                    /* print any turn-corner characters */
          putc(',', outfile);
        else {
          if (!iequal) {                   /* i.e., "itoright", so to speak */
            putc('\'', outfile);           /* "quoting" a single apostrophe */
          }
        }
        for (j = 1; j <= n - 3; j++)    /* ...  print dashes out to subtree */
          putc('-', outfile);
      }
      else {                           /* if in subtree but not at its root */
        for (j = 1; j <= n - 3; j++)    /* ...  print spaces out to subtree */
          putc(' ', outfile);
      } 
    }
    if (itoleft && (i > (long)r->back->ycoord)) {
        putc('|', outfile);           /* if branch to left crosses this row */
	printedbar = true;
    } else {
      if ((!iequal) && (!itoleft) && (i < (long)r->back->ycoord)) {
        putc('|', outfile);          /* if branch to right crosses this row */
        printedbar = true;
      }
      if ((!printedbar) && (!iequal))
        putc(' ', outfile);
    } 
    if (iinsubtree) {
      if (r->back != 0) {                     /* if branch is not empty ... */
        drawline2(i, scale, r->back, curtree);          /* ... start out it */
      }
    }
    r = r->next;                         /* move to next descendant, if any */
    if (!done) {
      if (r->back == 0) {              /* making sure not at bottom of tree */
        done = true;
      } else {
        if (r == p)        /* done if finished with all descendant branches */
          done = true;
	}
      }
  } while (!done);
}  /* drawline2 */


/* the previous version */
void drawline3(long i, double scale, struct node *p, struct tree* curtree)
{
  /* draws one row of the tree diagram by moving up tree
   * the argument  i  is the vertical number (y) of the row we draw,
   * numbered from top (1) to bottom
   * used in Dnaml, Proml, & Restml */
  struct node *pprev,  *q,  *r,  *rnext;
  long n, j;
  boolean extra, done, done2;

  q = p;                                               /* ... and so is  q  */
  extra = false;
  if (i == (long)p->ycoord)         /* if  i  is rootmost node's coordinate */
  {                                     /* write out the number of the node */
    if (p->index - spp >= 100)   /* can be changed to go beyond 999 species */
      fprintf(outfile, " %3ld", p->index - spp);
    else {
      if (p->index - spp >= 10)
        fprintf(outfile, "  %2ld", p->index - spp);
      else
        fprintf(outfile, "   %ld", p->index - spp);
    }
    extra = true;
  }
  else
    fprintf(outfile, "   ");                /* start by indenting two spaces */
  do {                                   /* working our way up the tree ... */
    if (!p->tip)
    {
      if (p->back != 0)      /* start with first nonempty descendant branch */
        r = p;
      else
        r = p->next;
      done2 = false;
      do {                                         /* go around fork circle */
        if (r->back != 0) {
          if ((i >= r->back->ymin) && (i <= r->back->ymax))
          {                            /* if this row intersects that clade */
            drawline3(i, scale, q, curtree);
            q = r->back;    /* ... then move to next node out that branch */
            done2 = true;  /* ... and note that are done circling that fork */
          }
	  rnext = r->next;                           /* next in fork circle */
          if ((i > r->back->ymax) && (i < rnext->back->ymin)) {
            done2 = true;     /* if in the gap between consecutive subtrees */
	  }
        }
        if (!done2)
          r = rnext;                         /* ... otherwise keep circling */
	done2 = done2 || (r == p) || (r->back == 0);  /* till where started */
      } while (!done2);                /* finished going around fork circle */
      pprev = p;                                /* pointer to where  p  was */
      p = q;                    /* ... and set  p  to next step up the tree */
    }
    done = (p->tip) || (pprev == q); /* done if at a tip, or not moved node */
    n = (long)(scale * (q->xcoord - pprev->xcoord) + 0.5); /* it's how far? */
    if ((n < 3) && !q->tip)    /* if interior branch, at least 3 chars long */
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i)                     /* if on row of next node */  
    {
      if (i < (long)pprev->ycoord)      /* print any turn-corner characters */
        putc(',', outfile);
      else {
        if (i > (long)(pprev)->ycoord)
          putc('\'', outfile);
      }
      for (j = 1; j <= n - 3; j++)         /* print line of "-" out to node */
        putc('-', outfile);
      assert(q->index > 0);           // RSGdebug
      if (!q->tip) {                               /*   if not at a tip ... */
        if (q->index - spp >= 100)       /* print out a number for the node */
          fprintf(outfile, "%3ld", q->index - spp);
	else {  
          if (q->index - spp >= 10)
            fprintf(outfile, "-%2ld", q->index - spp);
          else
            fprintf(outfile, "--%ld", q->index - spp);
        }
      }
    }
    if ((i > (long)r->back->ycoord) && ((long)pprev->ycoord > i)) {
      putc('|', outfile);                   /* if branch crosses this row */
    }
    if ((i < (long)r->back->ycoord) && ((long)pprev->ycoord < i)) {
      putc('|', outfile);                   /* if branch crosses this row */
    }
    if (!q->tip) {
      for (j = 1; j < n; j++)
        putc(' ', outfile);
    }
    extra = true;
  } while (!done);
  if (((long)q->ycoord == i) && q->tip)             /* if now at a tip, ... */
  {
    for (j = 0; j < nmlngth; j++)                 /* ... write the name ... */
      putc(nayme[q->index-1][j], outfile);
  }
  putc('\n', outfile);                               /* ... and end the row */
}  /* drawline2 */


void standev(long chars, long numtrees, long minwhich, double minsteps,
              double *nsteps, long **fsteps, longer seed)
{
  /* do paired sites test (KHT or SH test) on user-defined trees
   * used in dnapars & protpars */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

  (void)chars;                          // RSGnote: Parameter never used.

  if (numtrees == 2)
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees)
    {
      fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1]);
      if (minwhich == which)
        fprintf(outfile, "  <------ best\n");
      else
      {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < endsite; i++)
        {
          if (weight[i] > 0)
          {
            wt = weight[i] ;
            sumw += wt;
            temp = (fsteps[which - 1][i] - fsteps[minwhich - 1][i]);
            sum += wt * temp;
            sum2 += wt * temp * temp;
          }
        }
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw));
        fprintf(outfile, "%10.1f%12.4f", nsteps[which - 1] - minsteps, sd);
        if ((sum > 0.0) && (sum > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  }
  else
  {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile,
                   "Shimodaira-Hasegawa test on first %d of %ld trees\n\n",
                   MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    }
    else
    {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees * sizeof(double *));
    sumw = 0.0;
    for (i = 0; i < endsite; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees * sizeof(double));
    for (i = 0; i < numtrees; i++)          /* compute covariances of trees */
    {
      sum = nsteps[i]/(sumw);
      for (j = 0; j <=i; j++)
      {
        sum2 = nsteps[j]/sumw;
        temp = 0.0;
        for (k = 0; k < endsite; k++)
        {
          if (weight[k] > 0)
          {
            wt = weight[k];
            temp = temp + wt*(fsteps[i][k]-sum)
              *(fsteps[j][k]-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    /* in-place Cholesky decomposition of trees x trees covariance matrix */
    for (i = 0; i < numtrees; i++)
    {
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
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
    f = (double *)Malloc(numtrees * sizeof(double)); /* resampled sums */
    P = (double *)Malloc(numtrees * sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees * sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0];               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i])
        sum2 = nsteps[i];
    for (i = 1; i <= SAMPLES; i++)            /* loop over resampled trees */
    {
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++)          /* compute vectors */
      {
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get min of vector */
        if (f[j] < sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (nsteps[j]-sum2 <= f[j] - sum)
          P[j] += 1.0 / SAMPLES;
    }
    fprintf(outfile, "Tree    Steps   Diff Steps   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++)
    {
      fprintf(outfile, "%3ld%10.1f", i+1, nsteps[i]);
      if ((minwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        fprintf(outfile, " %9.1f  %10.3f", nsteps[i]-sum2, P[i]);
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
}  /* standev */


void standev2(long numtrees, long maxwhich, long a, long b, double maxlogl,
               double *l0gl, double **l0gf, steptr aliasweight, longer seed)
{
  /* do paired sites test (KHT or SH) for user-defined trees
   * used in dnaml, dnamlk, proml, promlk, and restml */
  double **covar, *P, *f, *r;
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;

  if (numtrees == 2)
  {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    logL    Diff logL    Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees)
    {
      fprintf(outfile, "%3ld %9.1f", which, l0gl[which - 1]);
      if (maxwhich == which)
        fprintf(outfile, "  <------ best\n");
      else
      {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = a; i <= b; i++)
        {
          if (aliasweight[i] > 0)
          {
            wt = aliasweight[i];
            sumw += wt;
            temp = l0gf[which - 1][i] - l0gf[maxwhich - 1][i];
            sum += temp * wt;
            sum2 += wt * temp * temp;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum / sumw ));
        fprintf(outfile, "%10.1f %11.4f", (l0gl[which - 1])-maxlogl, sd);
        if ((sum < 0.0) && ((-sum) > 1.95996 * sd))
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  }
  else
  {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES)
    {
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n",
                        MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    }
    else
    {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees * sizeof(double *));
    sumw = 0.0;
    for (i = a; i <= b; i++)
      sumw += aliasweight[i];
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
          if (aliasweight[k] > 0)
          {
            wt = aliasweight[k];
            temp = temp + wt * (l0gf[i][k]-sum) * (l0gf[j][k]-sum2);
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++)  /* in-place Cholesky decomposition
                                       of trees x trees covariance matrix */
    {
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i] <= sum)
        temp = 0.0;
      else
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
    f = (double *)Malloc(numtrees * sizeof(double)); /*resampled likelihoods*/
    P = (double *)Malloc(numtrees * sizeof(double)); /* vector of P's */
    r = (double *)Malloc(numtrees * sizeof(double)); /* Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    for (i = 1; i <= SAMPLES; i++)            /* loop over resampled trees */
    {
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++)          /* compute vectors */
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
        if (maxlogl-l0gl[j] <= sum-f[j])
          P[j] += 1.0 / SAMPLES;
    }
    fprintf(outfile, "Tree    logL    Diff logL    P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++)
    {
      fprintf(outfile, "%3ld%10.1f", i+1, l0gl[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else
      {
        fprintf(outfile, " %9.1f  %10.3f", l0gl[i]-maxlogl, P[i]);
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


/* End. */
