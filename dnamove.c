/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "moves.h"
#include "seq.h"
#include "dnaparsimony.h"

#define overr           4
#define which           1
#define maxsz           999   /* size of pointer array for the undo trees */
                              /* this can be large without eating memory */

tree* treesets[2];

node **treeone, **treetwo;

typedef enum {
  horiz, vert, up, overt, upcorner, midcorner, downcorner, aa, cc, gg, tt, question
} chartype;

typedef enum {
  rearr, flipp, reroott, none
} rearrtype;

typedef struct gbase2 {
  baseptr base2;
  struct gbase2 *next;
} gbase2;

typedef enum {
  arb, use, spec
} howtree;

typedef enum {beforenode, atnode} movet;

movet fromtype;

typedef node **pointptr;
tree* curtree;

#ifndef OLDC
/* function prototypes */
void getoptions(void);
void inputoptions(void);
void allocrest(void);
void doinput(void);
void configure(void);
void prefix(chartype);
void postfix(chartype);

void makechar(chartype);
void dnamove_reroot(node *);
void firstrav(node *, long);
void dnamove_hyptrav(node *, long *, long, boolean *);

void grwrite(chartype, long, long *);
void dnamove_drawline(long);
void dnamove_printree(void);
void arbitree(void);
void yourtree(void);
void initdnamovenode(node **, node **, node *, long, long, long *, long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void buildtree(void);
void setorder(void);
void mincomp(void);

void rearrange(void);
void dnamove_nextinc(void);
void dnamove_nextchar(void);
void dnamove_prevchar(void);
void dnamove_show(void);
void tryadd(node *, node **, node **, double *);
void addpreorder(node *, node *, node *, double *);
void try(void);
void undo(void);
void treewrite(boolean);

void clade(void);
void flip(long);
void changeoutgroup(void);
void redisplay(void);
void treeconstruct(void);
void prepare_node(node *p);
void numdesctrav(node *p);
void makeweights(void);
void add_at(node *below, node *newtip, node *newfork);
void add_before(node *atnode, node *newtip);
void add_child(node *parent, node *newchild);
void consolidatetree(long index);
void fliptrav(node *p, boolean recurse);
/* function prototypes */
#endif

char basechar[32]="ACMGRSVTWYHKDBNO???????????????";
char infilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
long chars, screenlines, col, treelines, leftedge, topedge, vmargin, hscroll, vscroll, scrollinc, screenwidth, whichtree, othertree, nonodes = 0;
double farthest;
boolean weights, thresh, waswritten;
boolean usertree, goteof, firsttree, haslengths;   /* treeread variables */

#if 0                                   // RSGbugfix: Global variable never used (slot in TREE object used instead).
pointarray nodep;
#endif

double threshold;                       /* treeread variables */
double *threshwt;
boolean reversed[(long)question - (long)horiz + 1];
boolean graphic[(long)question - (long)horiz + 1];
unsigned char chh[(long)question - (long)horiz + 1];
howtree how;
char *progname;

/* Local variables for treeconstruct, propagated global for C version: */
long dispchar, atwhat, what, fromwhere, towhere, oldoutgrno, compatible;
double like, bestyet, gotlike;
boolean display, newtree, changed, subtree, written, oldwritten, restoring, wasleft, oldleft, earlytree;
steptr necsteps;
boolean *in_tree;
long sett[31];
steptr numsteps;
node *nuroot;
rearrtype lastop;
Char  ch;
boolean *names;


void numdesctrav(node *p)
{
  node *q;
  long childcount = 0;

  if (p->tip)
  {
    return;
  }

  q = p->next;

  do {
    numdesctrav(q->back);
    childcount++;
    q = q->next;
  } while (q != p);

} /* numdesctrav */


void getoptions(void)
{
  /* interactively set options */
  Char ch;
  boolean done, gotopt;
  long loopcount;

  how = arb;
  usertree = false;
  goteof = false;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  weights = false;
  interleaved = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nInteractive DNA parsimony, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  O                             Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  W                            Sites weighted?  %s\n", (weights ? "Yes" : "No"));
    printf("  T                   Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  I               Input sequences interleaved?  %s\n", (interleaved ? "Yes" : "No, sequential"));
    printf("  U   Initial tree (arbitrary, user, specify)?  %s\n", (how == arb) ? "Arbitrary"                : (how == use) ? "User tree from tree file" : "Tree you specify");
    printf("  0        Graphics type (IBM PC, ANSI, none)?  %s\n", ibmpc ? "IBM PC" : ansi  ? "ANSI"     : "(none)");
    printf("  S                  Width of terminal screen?");
    printf("%4ld\n", screenwidth);
    printf("  L                 Number of lines on screen?%4ld\n", screenlines);
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    gotopt = (strchr("SOTIU0WL", ch) != NULL) ? true : false;
    if (gotopt)
    {
      switch (ch)
      {
        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
          break;

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
          break;

        case 'I':
          interleaved = !interleaved;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'U':
          if (how == arb)
          {
            how = use;
            usertree = 1;}
          else if (how == use)
          {
            how = spec;
            usertree = 0;}
          else
            how = arb;
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case 'S':
          screenwidth= readlong("Width of terminal screen (in characters)?\n");
          break;

        case 'L':
          initnumlines(&screenlines);
          break;
      }
    }
    if (!(gotopt || done))
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  } while (!done);
  if (scrollinc < screenwidth / 2.0)
    hscroll = scrollinc;
  else
    hscroll = screenwidth / 2;
  if (scrollinc < screenlines / 2.0)
    vscroll = scrollinc;
  else
    vscroll = screenlines / 2;
}  /* getoptions */


void inputoptions(void)
{
  /* input the information on the options */
  long i;

  for (i = 0; i < chars; i++)
    weight[i] = 1;
  if (weights)
  {
    inputweights(chars, weight, &weights);
    printweights(stdout, 0, chars, weight, "Sites");
  }
  if (!thresh)
    threshold = spp;
  for (i = 0; i < chars; i++)
    threshwt[i] = threshold * weight[i];
}  /* inputoptions */


void allocrest(void)
{
  long i;

  nayme = (naym *)Malloc(spp * sizeof(naym));
  in_tree = (boolean *)Malloc(nonodes * sizeof(boolean));
  weight = (steptr)Malloc(chars * sizeof(long));
  numsteps = (steptr)Malloc(chars * sizeof(long));
  necsteps = (steptr)Malloc(chars * sizeof(long));
  threshwt = (double *)Malloc(chars * sizeof(double));
  alias = (long *)Malloc(chars * sizeof(long));     /* from dnapars */
  ally = (long *)Malloc(chars * sizeof(long));      /* from dnapars */
  inputSequences = (Char **)Malloc(spp * sizeof(Char *));        /* from dnapars */
  for (i = 0; i < spp; i++)                       /* from dnapars */
    inputSequences[i] = (Char *)Malloc(chars * sizeof(Char));    /* from dnapars */
  location = (long *)Malloc(chars * sizeof(long));  /* from dnapars */
}  /* allocrest */


void makeweights(void)
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++)
  {
    alias[i - 1] = i;
    ally[i - 1] = i;
  }
  endsite = 0;
  for (i = 1; i <= chars; i++)
  {
    if (ally[i - 1] == i)
      endsite++;
  }

  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  if (!thresh)
    threshold = spp;
}  /* makeweights */


void doinput(void)
{
  /* reads the input data */
  inputnumbers(&spp, &chars, &nonodes, 1);
  printf("%2ld species, %3ld  sites\n", spp, chars);
  getoptions();
  printf("\nReading input file ...\n\n");
  if (weights)
    openfile(&weightfile, WEIGHTFILE, "weights file", "r", progname, weightfilename);
  allocrest();
  inputoptions();
  endsite = chars;                                 // needed for allocation ??
  curtree = funcs->tree_new(nonodes, spp);      // FIXME: usertree or spp?
  inputdata(chars);
  makeweights();
}  /* doinput */


void configure(void)
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)question; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)question; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc)
  {
    chh[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    chh[(long)vert] = 186;
    graphic[(long)vert] = true;
    chh[(long)up] = 186;
    graphic[(long)up] = true;
    chh[(long)overt] = 205;
    graphic[(long)overt] = true;
    chh[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    chh[(long)midcorner] = 204;
    graphic[(long)midcorner] = true;
    chh[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 176;
    chh[(long)cc] = 178;
    chh[(long)gg] = 177;
    chh[(long)tt] = 219;
    chh[(long)question] = '\001';
    return;
  }
  if (ansi)
  {
    chh[(long)horiz] = ' ';
    reversed[(long)horiz] = true;
    chh[(long)vert] = chh[(long)horiz];
    reversed[(long)vert] = true;
    chh[(long)up] = 'x';
    graphic[(long)up] = true;
    chh[(long)overt] = 'q';
    graphic[(long)overt] = true;
    chh[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    chh[(long)midcorner] = 't';
    graphic[(long)midcorner] = true;
    chh[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 'a';
    reversed[(long)aa] = true;
    chh[(long)cc] = 'c';
    reversed[(long)cc] = true;
    chh[(long)gg] = 'g';
    reversed[(long)gg] = true;
    chh[(long)tt] = 't';
    reversed[(long)tt] = true;
    chh[(long)question] = '?';
    reversed[(long)question] = true;
    return;
  }
  chh[(long)horiz] = '=';
  chh[(long)vert] = ' ';
  chh[(long)up] = '!';
  chh[(long)upcorner] = '`';
  chh[(long)midcorner] = '+';
  chh[(long)downcorner] = ',';
  chh[(long)overt] = '-';
  chh[(long)aa] = 'a';
  chh[(long)cc] = 'c';
  chh[(long)gg] = 'g';
  chh[(long)tt] = 't';
  chh[(long)question] = '.';
}  /* configure */


void prefix(chartype a)
{
  /* give prefix appropriate for this character */
  if (reversed[(long)a])
    prereverse(ansi);
  if (graphic[(long)a])
    pregraph2(ansi);
}  /* prefix */


void postfix(chartype a)
{
  /* give postfix appropriate for this character */
  if (reversed[(long)a])
    postreverse(ansi);
  if (graphic[(long)a])
    postgraph2(ansi);
}  /* postfix */


void makechar(chartype a)
{
  /* print out a character with appropriate prefix or postfix */
  prefix(a);
  putchar(chh[(long)a]);
  postfix(a);
}  /* makechar */


void add_at(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  node *leftdesc, *rtdesc;

  if (below != curtree->nodep[below->index - 1])
    below = curtree->nodep[below->index - 1];

  if (newfork == NULL)
  {
    nonodes++;
    newfork = curtree->get_forkring(curtree);
  }
  if (below->back != NULL)
  {
    below->back->back = newfork;
  }
  newfork->back = below->back;
  leftdesc = newtip;
  rtdesc = below;
  rtdesc->back = newfork->next->next;
  newfork->next->next->back = rtdesc;
  newfork->next->back = leftdesc;
  leftdesc->back = newfork->next;
  if (curtree->root == below)
    curtree->root = newfork;
  curtree->root->back = NULL;
}  /* add_at */


void add_before(node *atnode, node *newtip)
{
  /* Inserts the node newtip together with its ancestral fork into the tree next to the node atnode. */
  node *q;

  if (atnode != curtree->nodep[atnode->index - 1])
    atnode = curtree->nodep[atnode->index - 1];
  q = curtree->nodep[newtip->index-1]->back;
  if (q != NULL)
  {
    q = curtree->nodep[q->index-1];
    if (newtip == q->next->next->back)
    {
      q->next->back = newtip;
      newtip->back = q->next;
      q->next->next->back = NULL;
    }
  }
  if (newtip->back != NULL)
  {
    add_at(atnode, newtip, curtree->nodep[newtip->back->index-1]);
  }
  else
  {
    add_at(atnode, newtip, NULL);
  }
}  /* add_before */


void add_child(node *parent, node *newchild)
{
  /* adds the node newchild into the tree as the last child of parent */

  int i;
  node *newnode, *q;

  if (parent != curtree->nodep[parent->index - 1])
    parent = curtree->nodep[parent->index - 1];
  newnode = curtree->get_forknode(curtree, parent->index);
  newnode->tip = false;
  newnode->deleted=false;
  newnode->deadend=false;
  for (i=0;i<MAXNCH;i++)
    newnode->nayme[i] = '\0';
  q = parent;
  do {
    q = q->next;
  } while (q->next != parent);
  newnode->next = parent;
  q->next = newnode;
  newnode->back = newchild;
  newchild->back = newnode;
  if (newchild->haslength)
  {
    newnode->length = newchild->length;
    newnode->haslength = true;
  }
  else
    newnode->haslength = false;
}  /* add_child */


void dnamove_reroot(node *outgroup)
{
  /* Reorient tree so that outgroup is by itself on the left of the root */
  node *p, *q, *r;
  long nodecount = 0;
  double templen;

  if(outgroup->back->index == curtree->root->index)
    return;

  q = curtree->root->next;
  do {                    /* when this loop exits, p points to the internal */
    p = q;                /* node to the right of root */
    nodecount++;
    q = p->next;
  } while (q != curtree->root);
  r = p;

  /* reorient nodep array

     The nodep array must point to the ring member of each ring
     that is closest to the root.  The while loop changes the ring member
     pointed to by treenode[] for those nodes that will have their
     orientation changed by the reroot operation.
  */
  p = outgroup->back;
  while (p->index != curtree->root->index)
  {
    q = curtree->nodep[p->index - 1]->back;
    curtree->nodep[p->index - 1] = p;
    p = q;
  }
  if (nodecount > 2)
    curtree->nodep[p->index - 1] = p;

  /* If nodecount > 2, the current node ring to which root is pointing
     will remain in place and root will point somewhere else. */
  /* detach root from old location */
  if (nodecount > 2)
  {
    r->next = curtree->root->next;
    curtree->root->next = NULL;
    nonodes++;
    curtree->root = curtree->get_forkring(curtree);

    if (haslengths)
    {
      /* root->haslength remains false, or else treeout() will generate
         a bogus extra length */
      curtree->root->next->haslength = true;
      curtree->root->next->next->haslength = true;
    }
  }
  else  /* if (nodecount > 2) else */
  {
    q = curtree->root->next;
    q->back->back = r->back;
    r->back->back = q->back;

    if (haslengths)
    {
      r->back->length = r->back->length + q->back->length;
      q->back->length = r->back->length;
    }
  } /* if (nodecount > 2) endif */

  /* tie root into new location */
  curtree->root->next->back = outgroup;
  curtree->root->next->next->back = outgroup->back;
  outgroup->back->back = curtree->root->next->next;
  outgroup->back = curtree->root->next;

  /* place root equidistant between left child (outgroup) and
     right child by deviding outgroup's length */
  if (haslengths)
  {
    templen = outgroup->length / 2.0;
    outgroup->length = templen;
    outgroup->back->length = templen;
    curtree->root->next->next->length = templen;
    curtree->root->next->next->back->length = templen;
  }
} /* dnamove_reroot */


void grwrite(chartype c, long num, long *pos)
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++)
  {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(chh[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */


void dnamove_drawline(long i)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j, pos;
  boolean extra, done;
  Char st;
  chartype c, d;

  pos = 1;
  p = nuroot;
  q = nuroot;
  extra = false;
  if (i == p->ycoord && (p == curtree->root || subtree))
  {
    extra = true;
    c = overt;
    if (display)
    {
      switch (p->state)
      {
        case 'A':
          c = aa;
          break;

        case 'C':
          c = cc;
          break;

        case 'G':
          c = gg;
          break;

        case 'T':
          c = tt;
          break;

        case '?':
          c = question;
          break;
      }
    }
    if ((subtree))
      stwrite("Subtree:", 8, &pos, leftedge, screenwidth);
    if (p->index >= 100)
      nnwrite(p->index, 3, &pos, leftedge, screenwidth);
    else if (p->index >= 10)
    {
      grwrite(c, 1, &pos);
      nnwrite(p->index, 2, &pos, leftedge, screenwidth);
    }
    else
    {
      grwrite(c, 2, &pos);
      nnwrite(p->index, 1, &pos, leftedge, screenwidth);
    }
  }
  else
  {
    if ((subtree))
      stwrite("          ", 10, &pos, leftedge, screenwidth);
    else
      stwrite("  ", 2, &pos, leftedge, screenwidth);
  }
  do {
    if (!p->tip)
    {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax)
        {
          q = r->back;
          done = true;
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
    n = p->xcoord - q->xcoord;
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done)
    {
      c = overt;

      if (q == first)
        d = downcorner;
      else if (q == last)
        d = upcorner;
      else if ((long)q->ycoord == (long)p->ycoord)
        d = c;
      else
        d = midcorner;

      if (display)
      {
        switch (q->state)
        {
          case 'A':
            c = aa;
            break;

          case 'C':
            c = cc;
            break;

          case 'G':
            c = gg;
            break;

          case 'T':
            c = tt;
            break;

          case '?':
            c = question;
            break;
        }
        d = c;
      }
      if (n > 1)
      {
        grwrite(d, 1, &pos);
        grwrite(c, n - 3, &pos);
      }
      if (q->index >= 100)
        nnwrite(q->index, 3, &pos, leftedge, screenwidth);
      else if (q->index >= 10)
      {
        grwrite(c, 1, &pos);
        nnwrite(q->index, 2, &pos, leftedge, screenwidth);
      }
      else
      {
        grwrite(c, 2, &pos);
        nnwrite(q->index, 1, &pos, leftedge, screenwidth);
      }
      extra = true;
    }
    else if (!q->tip)
    {
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord)
      {
        c = up;
        if (i < p->ycoord)
          st = p->next->back->state;
        else
          st = p->next->next->back->state;
        if (display)
        {
          switch (st)
          {
            case 'A':
              c = aa;
              break;

            case 'C':
              c = cc;
              break;

            case 'G':
              c = gg;
              break;

            case 'T':
              c = tt;
              break;

            case '?':
              c = question;
              break;
          }
        }
        grwrite(c, 1, &pos);
        chwrite(' ', n - 1, &pos, leftedge, screenwidth);
      }
      else
        chwrite(' ', n, &pos, leftedge, screenwidth);
    }
    else
      chwrite(' ', n, &pos, leftedge, screenwidth);
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip)
  {
    n = 0;
    for (j = 1; j <= nmlngth; j++)
    {
      if (nayme[p->index - 1][j - 1] != '\0')
        n = j;
    }
    chwrite(':', 1, &pos, leftedge, screenwidth);
    for (j = 0; j < n; j++)
      chwrite(nayme[p->index - 1][j], 1, &pos, leftedge, screenwidth);
  }
  putchar('\n');
}  /* dnamove_drawline */


void dnamove_printree(void)
{
  /* prints out diagram of the tree */
  long tipy;
  long i, dow;

  if (!subtree)
    nuroot = curtree->root;
  if (changed || newtree)
    curtree->evaluate(curtree, curtree->root, 0);
  if (display)
  {
    outfile = stdout;
    dna_hypstates(curtree, chars, basechar);
  }
#ifdef WIN32
  if(ibmpc || ansi)
    phyClearScreen();
  else
    printf("\n");
#else
  printf((ansi || ibmpc) ? "\033[2J\033[H" : "\n");
#endif
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree)
    dow--;

  printf("  (unrooted)");
  if (display)
  {
    printf(" ");
    makechar(aa);
    printf(":A ");
    makechar(cc);
    printf(":C ");
    makechar(gg);
    printf(":G ");
    makechar(tt);
    printf(":T ");
    makechar(question);
    printf(":?");
  }
  else
    printf("                    ");
  if (!earlytree)
  {
    printf("%10.1f Steps", -like);
  }
  if (display)
    printf(" SITE%4ld", dispchar);
  else
    printf("         ");
  if (!earlytree)
  {
    printf("  %3ld sites compatible\n", compatible);
  }

  printf("                            ");
  if (changed && !earlytree)
  {
    if (-like < bestyet)
    {
      printf("     BEST YET!");
      bestyet = -like;
    }
    else if (fabs(-like - bestyet) < 0.000001)
      printf("     (as good as best)");
    else
    {
      if (-like < gotlike)
        printf("     better");
      else if (-like > gotlike)
        printf("     worse!");
    }
  }
  printf("\n");

  farthest = 0;
  coordinates(curtree, nuroot, 1.5, &tipy,  &farthest);
  vmargin = 4;
  treelines = tipy - dow;
  if (topedge != 1)
  {
    printf("** %ld lines above screen **\n", topedge - 1);
    vmargin++;
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin))
    vmargin++;
  for (i = 1; i <= treelines; i++)
  {
    if (i >= topedge && i < topedge + screenlines - vmargin)
      dnamove_drawline(i);
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin))
  {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
  gotlike = -like;
  changed = false;
}  /* dnamove_printree */


void arbitree(void)
{
  long i;
  curtree->root = curtree->nodep[0];
  hookup(curtree->nodep[0], curtree->nodep[1]);
  for (i = 3; i <= spp; i++)
  {
    curtree->insert_(curtree, curtree->nodep[i - 1], curtree->nodep[i - 2], false, false);
  }
}  /* arbitree */


void yourtree(void)
{
  long i, j;
  boolean ok;

  curtree->root = curtree->nodep[0];
  hookup(curtree->nodep[0], curtree->nodep[1]);
  i = 2;
  do {
    i++;
    dnamove_printree();
    printf("Add species%3ld: ", i);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[i - 1][j]);
    do {
      printf("\n at or before which node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((j >= 1 && j < i) || (j > spp && j < spp + i - 1)));
      if (!ok)
        printf("Impossible number. Please try again:\n");
    } while (!ok);

    if (j >= i)   /* has user chosen a non-tip? if so, offer choice */
    {
      do {
        printf(" Insert at node (A) or before node (B)? ");
        phyFillScreenColor();
        if(scanf("%c%*[^\n]", &ch)) {}  // Read char and scan to EOL.
        (void)getchar();
        if (ch == '\n')
          ch = ' ';
        ch = isupper(ch) ? ch : toupper(ch);
      } while (ch != 'A' && ch != 'B');
    }
    else ch = 'B';   /* if user has chosen a tip, set Before */

    if (j != 0)
    {
      if (ch == 'A')
      {
        if (!curtree->nodep[j - 1]->tip)
        {
          add_child(curtree->nodep[j - 1], curtree->nodep[i - 1]);
        }
      }
      else
      {
        printf("dnamove_add(below %ld, newtip %ld, newfork %ld)\n", j-1, i-1, spp+i-2);
        curtree->insert_(curtree, curtree->nodep[i-1], curtree->nodep[j-1], false, false);
      } /* endif (before or at node) */
    }
  } while (i != spp);
}  /* yourtree */


void buildtree(void)                    // RSGbugfix
{
  long i, nextnode;
  node *p;
  long j;

  treesets[0] = (tree*)dnapars_tree_new(nonodes, spp);
  treesets[1] = (tree*)dnapars_tree_new(nonodes, spp);
  changed = false;
  newtree = false;
  switch (how)
  {
    case arb:
      arbitree();
      break;

    case use:
      /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
      openfile(&intree, intreename, "input tree file", "rb", progname, intreename);
      names = (boolean *)Malloc(spp * sizeof(boolean));
      firsttree = true;

#if 0                                   // RSGbugfix: Global variable never used.
      nodep = NULL;
#endif

      nextnode = 0;
      haslengths = 0;
      treeread(curtree, intree, &curtree->root, curtree->nodep, &goteof, &firsttree, &nextnode, &haslengths, initparsnode, true, nonodes);

      for (i = spp; i < (nextnode); i++)
      {
        p = curtree->nodep[i];
        for (j = 1; j <= 3; j++)
        {
          ((dnapars_node*)p)->base = (baseptr)Malloc(chars * sizeof(long));
          p = p->next;
        }
      } /* debug: see comment at initdnamovenode() */

      free(names);
      FClose(intree);
      break;

    case spec:
      yourtree();
      break;
  }
  if (!outgropt)
  {
    if (curtree->root->tip)
      // other code seems to believe a tip isn't a root, yet
      // it appears to be true. it is possible that fixing
      // another bug will make the if half of this case moot
    {
      outgrno = curtree->root->index;
    }
    else
    {
      outgrno = curtree->root->next->back->index;
    }
  }
  if (outgropt)
    dnamove_reroot(curtree->nodep[outgrno - 1]);

}  /* buildtree */


void setorder(void)
{
  /* sets in order of number of members */
  sett[0] = 1L << ((long)A);
  sett[1] = 1L << ((long)C);
  sett[2] = 1L << ((long)G);
  sett[3] = 1L << ((long)T);
  sett[4] = 1L << ((long)O);
  sett[5] = (1L << ((long)A)) | (1L << ((long)C));
  sett[6] = (1L << ((long)A)) | (1L << ((long)G));
  sett[7] = (1L << ((long)A)) | (1L << ((long)T));
  sett[8] = (1L << ((long)A)) | (1L << ((long)O));
  sett[9] = (1L << ((long)C)) | (1L << ((long)G));
  sett[10] = (1L << ((long)C)) | (1L << ((long)T));
  sett[11] = (1L << ((long)C)) | (1L << ((long)O));
  sett[12] = (1L << ((long)G)) | (1L << ((long)T));
  sett[13] = (1L << ((long)G)) | (1L << ((long)O));
  sett[14] = (1L << ((long)T)) | (1L << ((long)O));
  sett[15] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G));
  sett[16] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)T));
  sett[17] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)O));
  sett[18] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)T));
  sett[19] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[20] = (1L << ((long)A)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[21] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)T));
  sett[22] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[23] = (1L << ((long)C)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[24] = (1L << ((long)G)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[25] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)T));
  sett[26] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)O));
  sett[27] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[28] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[29] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[30] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)T)) | (1L << ((long)O));
}  /* setorder */


void mincomp(void)
{
  /* computes for each site the minimum number of steps necessary to accomodate those species already in the analysis */
  long i, j, k;
  boolean done;

  for (i = 0; i < chars; i++)
  {
    done = false;
    j = 0;
    while (!done)
    {
      j++;
      done = true;
      k = 1;
      do {
        if (k < nonodes)
        {
          // RSGdebug: Potential infinite loop.  When crashing on memory-protection, j has value 9005.
          done = (done && (((dnapars_node*)curtree->nodep[k - 1])->base[i] & sett[j - 1]) != 0);
        }
        k++;
      } while (k <= spp && done);
    }
    if (j == 31)
      necsteps[i] = 4;
    if (j <= 30)
      necsteps[i] = 3;
    if (j <= 25)
      necsteps[i] = 2;
    if (j <= 15)
      necsteps[i] = 1;
    if (j <= 5)
      necsteps[i] = 0;
    necsteps[i] *= weight[i];
  }
}  /* mincomp */


void consolidatetree(long index)
{
  node *start, *r, *q;
  int i;

  i = 0;

  start = curtree->nodep[index - 1];
  q = start->next;

  while (q != start)
  {
    r = q;
    q = q->next;
    r->next = NULL; // FIXME: is garbage collection needed?
  }
  q->next = NULL; // FIXME: is garbage collection needed?

  i = index;

  while (i <= nonodes)
  {
    r = curtree->nodep[i - 1];
    if (!(r->tip))
      r->index--;
    if (!(r->tip))
    {
      q = r->next;
      do {
        q->index--;
        q = q->next;
      } while (r != q && q != NULL);
    }
    curtree->nodep[i - 1] = curtree->nodep[i];
    i++;
  }

  nonodes--;
} /* consolidatetree */


void rearrange(void)
{
  long i, j, maxinput;
  boolean ok1, ok2;
  node *p, *q;
  char ch;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i <= (spp * 2 - 1) && i != curtree->root->index);
  if (ok1)
    ok1 = !curtree->nodep[i - 1]->deleted;
  if (ok1)
  {
    printf("Add at or before which node? ");
    inpnum(&j, &ok2);
    ok2 = (ok2 && j >= 1 && j <= (spp * 2 - 1));
    if (ok2)
    {
      if (j != curtree->root->index)
        ok2 = !curtree->nodep[curtree->nodep[j - 1]->back->index - 1]->deleted;
    }
    if (ok2)
    {
      // xx This edit says "j must not be i's parent".  Is this necessary anymore?
#if 0
      ok2 = (nodep[j - 1] != nodep[nodep[i - 1]->back->index - 1]);
#endif

      p = curtree->nodep[j - 1];
      /* make sure that j is not a descendent of i */
      while (p != curtree->root)
      {
        ok2 = (ok2 && p != curtree->nodep[i - 1]);
        p = curtree->nodep[p->back->index - 1];
      }
      if (ok1 && ok2)
      {
        maxinput = 1;
        do {
          printf("Insert at node (A) or before node (B)? ");
          phyFillScreenColor();
          if(scanf("%c%*[^\n]", &ch)) {} // Read char and scan to EOL.
          (void)getchar();
          if (ch == '\n')
            ch = ' ';
          ch = isupper(ch) ? ch : toupper(ch);
          maxinput++;
          if (maxinput == 100)
          {
            printf("\nERROR:  Too many tries at choosing option.\n");
            exxit(-1);
          }
        } while (ch != 'A' && ch != 'B');
        if (ch == 'A')
        {
          if (!(curtree->nodep[j - 1]->deleted) && !curtree->nodep[j - 1]->tip)
          {
            changed = true;
            curtree->copy(curtree, treesets[othertree]);
            curtree->re_move(curtree, curtree->nodep[i - 1], &q, false);
            add_child(curtree->nodep[j - 1], curtree->nodep[i - 1]);
            if (fromtype == beforenode)
              consolidatetree(q->index);
          }
          else
            ok2 = false;
        }
        else
        {
          if (j != curtree->root->index) /* can't insert at root */
          {
            changed = true;
            curtree->copy(curtree, treesets[othertree]);
            curtree->re_move(curtree, curtree->nodep[i-1], &q, true);
            if (q != NULL)
            {
              curtree->nodep[q->index-1]->next->back = curtree->nodep[i-1];
              curtree->nodep[i-1]->back = curtree->nodep[q->index-1]->next;
            }
            add_before(curtree->nodep[j - 1], curtree->nodep[i - 1]);
          }
          else
            ok2 = false;
        } /* endif (before or at node) */
      } /* endif (ok to do move) */
    } /* endif (destination node ok) */
  } /* endif (from node ok) */
  dnamove_printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else
  {
    written = false;
  }
}  /* rearrange */


void dnamove_nextinc(void)
{
  /* show next incompatible site */
  long disp0;
  boolean done;

  display = true;
  disp0 = dispchar;
  done = false;
  do {
    dispchar++;
    if (dispchar > chars)
    {
      dispchar = 1;
      done = (disp0 == 0);
    }
  } while (!(necsteps[dispchar - 1] != numsteps[dispchar - 1] ||
             dispchar == disp0 || done));
  dnamove_printree();
}  /* dnamove_nextinc */


void dnamove_nextchar(void)
{
  /* show next site */
  display = true;
  dispchar++;
  if (dispchar > chars)
    dispchar = 1;
  dnamove_printree();
}  /* dnamove_nextchar */


void dnamove_prevchar(void)
{
  /* show previous site */
  display = true;
  dispchar--;
  if (dispchar < 1)
    dispchar = chars;
  dnamove_printree();
}  /* dnamove_prevchar */


void dnamove_show(void)
{
  long i;
  boolean ok;

  do {
    printf("SHOW: (Character number or 0 to see none)? ");
    inpnum(&i, &ok);
    ok = (ok && (i == 0 || (i >= 1 && i <= chars)));
    if (ok && i != 0)
    {
      display = true;
      dispchar = i;
    }
    if (ok && i == 0)
      display = false;
  } while (!ok);
  dnamove_printree();
}  /* dnamove_show */


void tryadd(node *p, node **item, node **nufork, double *place)
{
  /* Temporarily adds one fork and one tip to the tree.
     Records scores in ARRAY place */
  curtree->insert_(curtree, *item, p, false, false);
  curtree->evaluate(curtree, curtree->root, 0);
  place[p->index - 1] = -like;
  curtree->re_move(curtree, *item, nufork, false);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_, double *place)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *item, *nufork, *q;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item, &nufork, place);
  if (!p->tip)
  {
    q = p->next;
    do {
      addpreorder(q->back, item, nufork, place);
      q = q->next;
    } while (q != p);
  }
}  /* addpreorder */


void try(void)
{
  /* Remove node, try it in all possible places */
  double *place;
  long i, j, oldcompat, saveparent;
  double current;
  node *q, *dummy, *rute;
  boolean tied, better, ok, madenode;

  madenode = false;
  printf("Try other positions for which node? ");
  inpnum(&i, &ok);
  if (!(ok && i >= 1 && i <= nonodes && i != curtree->root->index))
  {
    printf("Not a possible choice! ");
    return;
  }
  curtree->copy(curtree, treesets[othertree]);
  printf("WAIT ...\n");
  place = (double *)Malloc(nonodes * sizeof(double));
  for (j = 0; j < nonodes; j++)
    place[j] = -1.0;
  curtree->evaluate(curtree, curtree->root, 0);
  current = -like;

  oldcompat = compatible;
  what = i;
  /* q = ring base of i's parent */
  q = curtree->nodep[curtree->nodep[i - 1]->back->index - 1];
  saveparent = q->index;
  /* if i is a left child, fromwhere = index of right sibling (binary) */
  /* if i is a right child, fromwhere = index of left sibling (binary) */
  if (q->next->back->index == i)
    fromwhere = q->next->next->back->index;
  else
    fromwhere = q->next->back->index;
  rute = curtree->root;

  /* if root is i's parent ... */
  if (q->next->next->next == q)
  {
    if (curtree->root == curtree->nodep[curtree->nodep[i - 1]->back->index - 1])
    {
      /* if i is left child then rute becomes right child, and vice-versa */
      if (curtree->nodep[curtree->nodep[i - 1]->back->index - 1]->next->back == curtree->nodep[i - 1])
        rute = curtree->nodep[curtree->nodep[i - 1]->back->index - 1]->next->next->back;
      else
        rute = curtree->nodep[curtree->nodep[i - 1]->back->index - 1]->next->back;
    }
  }

  /* Remove i and perhaps its parent node from the tree.  If i is part of a
     multifurcation, *dummy will come back null.  If so, make a new internal
     node to be i's parent as it is inserted in various places around the
     tree.
  */
  curtree->re_move(curtree, curtree->nodep[i-1], &dummy, false);
  if (dummy == NULL)
  {
    madenode = true;
    nonodes++;
    curtree->get_forkring(curtree);
  }
  oldleft = wasleft;
  curtree->root = rute;
  addpreorder(curtree->root, curtree->nodep[i - 1], dummy, place);
  wasleft = oldleft;
  restoring = true;
  if (madenode)
  {
    add_child(curtree->nodep[saveparent - 1], curtree->nodep[i - 1]);
    nonodes--;
  }
  else
    curtree->insert_(curtree, curtree->nodep[what - 1], curtree->nodep[fromwhere- 1], true, false);
  like = -current;
  compatible = oldcompat;
  restoring = false;
  better = false;
  printf("       BETTER: ");
  for (j = 1; j <= nonodes; j++)
  {
    if (place[j - 1] < current && place[j - 1] >= 0.0)
    {
      printf("%3ld:%6.2f", j, place[j - 1]);
      better = true;
    }
  }
  if (!better)
    printf(" NONE");
  printf("\n       TIED:    ");
  tied = false;
  for (j = 1; j <= nonodes; j++)
  {
    if (fabs(place[j - 1] - current) < 1.0e-6 && j != fromwhere)
    {
      if (j < 10)
        printf("%2ld", j);
      else
        printf("%3ld", j);
      tied = true;
    }
  }
  if (tied)
    printf(":%6.2f\n", current);
  else
    printf("NONE\n");
  changed = true;
  free(place);
}  /* try */


void undo(void)
{
  boolean btemp;

  /* don't undo to an uninitialized tree */
  if (treesets[othertree] == NULL)
  {
    dnamove_printree();
    printf("Nothing to undo.\n");
    return;
  }

  if (othertree == 0)
    othertree = 1;
  else
    othertree = 0;

  changed = true;
  dnamove_printree();
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* undo */


void treewrite(boolean done)
{
  /* write out tree to a file */
  Char ch;

  treeoptions(waswritten, &ch, &outtree, outtreename, progname);
  if (!done)
    dnamove_printree();
  if (waswritten && ch != 'A' && ch != 'R')
    return;
  col = 0;
  treeout(curtree->root, 1, &col, curtree->root);
  printf("\nTree written to file \"%s\".\n\n", outtreename);
  waswritten = true;
  written = true;
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
}
/* treewrite */


void clade(void)
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && ((unsigned)i) <= ((unsigned)nonodes));
  if (ok)
  {
    subtree = (i > 0);
    if (subtree)
      nuroot = curtree->nodep[i - 1];
    else
      nuroot = curtree->root;
  }
  dnamove_printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


void fliptrav(node *p, boolean recurse)
{
  node *q, *temp, *r =NULL, *rprev =NULL, *l, *lprev;
  boolean lprevflag;
  int nodecount, loopcount, i;

  if (p->tip)
    return;

  q = p->next;
  l = q;
  lprev = p;
  nodecount = 0;

  do {
    nodecount++;
    if (q->next->next == p)
    {
      rprev = q;
      r = q->next;
    }
    q = q->next;
  } while (p != q);

  if (nodecount == 1)
    return;
  loopcount = nodecount / 2;

  for (i=0; i<loopcount; i++)
  {
    lprev->next = r;
    rprev->next = l;
    temp = r->next;
    r->next = l->next;
    l->next = temp;
    if (i < (loopcount - 1))
    {
      lprevflag = false;
      q = p->next;
      do {
        if (q == lprev->next && !lprevflag)
        {
          lprev = q;
          l = q->next;
          lprevflag = true;
        }
        if (q->next == rprev)
        {
          rprev = q;
          r = q->next;
        }
        q = q->next;
      } while (p != q);
    }
  }
  if (recurse)
  {
    q = p->next;
    do {
      fliptrav(q->back, true);
      q = q->next;
    } while (p != q);
  }
}  /* fliptrav */


void flip(long atnode)
{
  /* flip at a node left-right */
  long i;
  boolean ok;

  if (atnode == 0)
  {
    printf("Flip branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
  }
  else
  {
    i = atnode;
    ok = true;
  }
  if (ok)
  {
    curtree->copy(curtree, treesets[othertree]);
    fliptrav(curtree->nodep[i - 1], true);
  }
  if (atnode == 0)
    dnamove_printree();
  if (ok)
  {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) ||
      (i > spp && i <= nonodes))
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */


void changeoutgroup(void)
{
  long i;
  boolean ok;

  oldoutgrno = outgrno;
  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && i >= 1 && i <= nonodes &&
          i != curtree->root->index);
    if (ok)
      outgrno = i;
  } while (!ok);

  curtree->copy(curtree, treesets[othertree]);
  dnamove_reroot(curtree->nodep[outgrno - 1]);
  changed = true;
  lastop = reroott;
  dnamove_printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */


void redisplay(void)
{
  boolean done = false;
  waswritten = false;
  do {
    printf("NEXT? (Options: R # + - S . T U W O F H J K L C ? X Q) ");
    printf("(? for Help) ");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    uppercase(&ch);
    if (strchr("HJKLCFORSTUXQ+#-.W?", ch) != NULL)
    {
      switch (ch)
      {
        case 'R':
          rearrange();
          break;

        case '#':
          dnamove_nextinc();
          break;

        case '+':
          dnamove_nextchar();
          break;

        case '-':
          dnamove_prevchar();
          break;

        case 'S':
          dnamove_show();
          break;

        case '.':
          dnamove_printree();
          break;

        case 'T':
          try();
          break;

        case 'U':
          undo();
          break;

        case 'W':
          treewrite(done);
          break;

        case 'O':
          changeoutgroup();
          break;

        case 'F':
          flip(0);
          break;

        case 'H':
          window(left, &leftedge, &topedge, hscroll, vscroll, treelines, screenlines, screenwidth, farthest, subtree);
          dnamove_printree();
          break;

        case 'J':
          window(downn, &leftedge, &topedge, hscroll, vscroll, treelines, screenlines, screenwidth, farthest, subtree);
          dnamove_printree();
          break;

        case 'K':
          window(upp, &leftedge, &topedge, hscroll, vscroll, treelines, screenlines, screenwidth, farthest, subtree);
          dnamove_printree();
          break;

        case 'L':
          window(right, &leftedge, &topedge, hscroll, vscroll, treelines, screenlines, screenwidth, farthest, subtree);
          dnamove_printree();
          break;

        case 'C':
          clade();
          break;

        case '?':
          help("site");
          dnamove_printree();
          break;

        case 'X':
          done = true;
          break;

        case 'Q':
          done = true;
          break;
      }
    }
  } while (!done);
  if (written)
    return;
  do {
    printf("Do you want to write out the tree to a file? (Y or N) ");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == 'Y' || ch == 'y')
      treewrite(done);
  } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
}  /* redisplay */


void treeconstruct(void)
{
  /* constructs a binary tree from the pointers in treenode. */
  int i;

  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  earlytree = true;
  waswritten = false;
  buildtree();

  /* get an accurate value for nonodes by finding out where the nodes really stop */
  for (i=0;i<nonodes;i++)
  {
    if (curtree->nodep[i]==NULL)
      break;
  }
  nonodes = i;

  printf("\nComputing steps needed for compatibility in sites ...\n\n");
  setorder();
  mincomp();
  newtree = true;
  earlytree = false;
  dnamove_printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  redisplay();
}  /* treeconstruct */


int main(int argc, Char *argv[])
{ /* Interactive DNA parsimony */
  /* reads in spp, chars, and the data. Then calls treeconstruct to
     construct the tree and query the user */
  initdata funcs;
#ifdef MAC
  argc = 1;                /* macsetup("Dnamove", "");        */
  argv[0] = "Dnamove";
#endif
  memset(&funcs, 0, sizeof(funcs));
  funcs.node_new = dnapars_node_new;
  funcs.tree_new = dnapars_tree_new;

  phylipinit(argc, argv, &funcs, false);
  progname = argv[0];
  strcpy(infilename, INFILE);
  strcpy(intreename, INTREE);
  strcpy(outtreename, OUTTREE);

  openfile(&infile, infilename, "input file", "r", argv[0], infilename);
  openfile(&outtree, outtreename, "output file", "w", argv[0], outtreename);

  whichtree = 0;
  othertree = 1;
  treesets[whichtree] = NULL;
  treesets[othertree] = NULL;
  screenlines = 24;
  scrollinc = 20;
  screenwidth = 80;
  printdata = 0;
  topedge = 1;
  leftedge = 1;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinput();
  configure();
  treeconstruct();
  FClose(infile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
  printf("\nDone.\n\n");
  phyRestoreConsoleAttributes();
  return 0;
}  /* Interactive DNA parsimony */


// End.
