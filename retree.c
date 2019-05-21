/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joseph Felsenstein and Andrew Keeffe.  Permission is granted to
   copy and use this program provided no fee is charged for it and provided
   that this copyright notice is not removed. */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "moves.h"
//#include "draw.h"

#define maxsp           50000   /* maximum number of species               */
#define maxsz           99999  /* size of pointer array.  >= 2*maxsp - 1  */
                               /* these can be large without eating memory */

#define overr           4
#define which           1


typedef enum {valid, remoov, quit} reslttype;

typedef enum {
  horiz, vert, up, updel, ch_over, upcorner, midcorner, downcorner, aa, cc,
  gg, tt, deleted
} chartype;

typedef struct treeset_t {
  tree * tree_p;
  boolean waswritten, hasmult, haslengths, nolengths, initialized;
} treeset_t;

treeset_t treesets[2];
treeset_t simplifiedtree;

typedef enum { arb, use, spec } howtree;

typedef enum {beforenode, atnode} movet;

movet fromtype;

#ifndef OLDC
/* function prototypes */
void   initretreenode(tree *, node **, long, long, long *, long *, initops, pointarray, Char *, Char *, FILE *);
void   maketriad(tree *, node **, long);
void   maketip(tree *, node **, long);
void   copytree(void);
void   getoptions(void);
void   configure(void);
void   prefix(chartype);

void   postfix(chartype);
void   ltrav(node *, boolean *);
boolean ifhaslengths(void);
void   add_at(tree *, node *, node *, node *);
void   add_before(tree *, node *, node *);
void   add_child(tree *, node *, node *);
void   re_move(tree *, node **, node **);
void   reroot(tree *, node *);
void   ltrav_(tree *, node *, double, double, double *, long *, long *);

void   precoord(node *, boolean *, double *, long *);
void   coordinates(node *, double, long *, long *, double *);
void   flatcoordinates(node *, long *);
void   grwrite(chartype, long, long *);
void   drawline(long, node *, boolean *);
void   printree(void);
void   togglelengths(void);
void   arbitree(tree *);
void   yourtree(void);
void   buildtree(void);
void   unbuildtree(void);
void   retree_help(void);
void   consolidatetree(long);
void   rearrange(void);
boolean any_deleted(node *);
void   fliptrav(node *, boolean);
void   flip(long);
void   transpose(long);
void   ifdeltrav(node *, boolean *);
double oltrav(node *);
void   outlength(void);
void   midpoint(void);
void   deltrav(node *, boolean );
void   reg_del(node *, boolean);
boolean isdeleted(long);
void   deletebranch(void);
void   restorebranch(void);
void   del_or_restore(void);
void   undo(void);
void   treetrav(node *);
void   simcopynode(node *, node *);
node  *simcopytrav(tree * src, node * srcRoot, tree * dest);
void   simcopytree(void);
void   writebranchlength(double);
void   treeout(node *, boolean, double, long, boolean);
void   maketemptriad(node **, long);
void   roottreeout(boolean *);
void   notrootedtorooted(void);
void   rootedtonotrooted(void);
void   treewrite(boolean *);
void   retree_window(adjwindow);
void   getlength(double *, reslttype *, boolean *);
void   changelength(void);
void   changename(void);
void   cladecollapse(void);
void   clade(void);
void   changeoutgroup(void);
void   redisplay(void);
void   treeconstruct(void);
void   retreeread(char * intreename, int intreenum, char * usefont, int usebranchlengths);
void   retree(char * intreename, char * intree, char * outtreename, char * outtreeopt, char * outtreefmt, char * usefont,
              char * plotfilename, char * plotfileopt, int usebranchlengths, char * option, int activenode, int refnode,
              double branchlength, int rearrbefore, char * newname, int rootnode, char * outtreeformat, int readtree,
              int writetree, int update, int undo, int doplot, char * finalplotkind);
void   fill_del(node*p);
/* function prototypes */
#endif

long outgrno, screenwidth, vscreenwidth, screenlines, col, treenumber, leftedge, topedge, treelines,
  hscroll, vscroll, scrollinc, whichtree, othertree, numtrees, treesread, nonodes = 0;

double     trweight;
boolean    anywritten, waswritten, onfirsttree, hasmult, haslengths, nolengths, nexus, xmltree;

tree * treeone, * treetwo, * curtree;

boolean    reversed[14];
boolean    graphic[14];
unsigned char       cch[14];
howtree    how;

char       intreename[FNMLNGTH], outtreename[FNMLNGTH];

boolean    subtree, written, readnext, writeindented, indentnotchanged;
long       indentlevel;
node      *nuroot;
Char      ch;

boolean delarray[maxsz];


void initretreenode(tree *treep, node **p, long len,
                    long nodei, long *ntips, long *parens, initops whichinit,
                    pointarray treenode, Char *str,
                    Char *ch, FILE *intree)
{
  /* initializes a node */
  long i;
  boolean minusread;
  double valyew, divisor;

  (void)len;                            // RSGnote: Parameter never used.
  (void)treenode;                       // RSGnote: Parameter never used.

  switch (whichinit)
  {
    case bottom:
      *p = treep->get_forknode(treep, nodei);
      (*p)->index = nodei;
      (*p)->tip = false;
      (*p)->deleted=false;
      (*p)->deadend=false;
      (*p)->onebranch=false;
      (*p)->onebranchhaslength=false;
      for (i=0;i<MAXNCH;i++)
        (*p)->nayme[i] = '\0';
      treep->nodep[(*p)->index - 1] = (*p);
      break;

    case nonbottom:
      *p = treep->get_forknode(treep, nodei);
      (*p)->index = nodei;
      break;

    case hslength:
      if ((*p)->back)
      {
        (*p)->back->back = *p;
        (*p)->haslength = (*p)->back->haslength;
        if ((*p)->haslength)
          (*p)->length = (*p)->back->length;
      }
      break;

    case tip:
      (*ntips)++;
      *p = treep->get_forknode(treep, nodei);
      treep->nodep[(*ntips) - 1] = *p;
      (*p)->index = *ntips;
      (*p)->tip = true;
      (*p)->hasname = true;
      strncpy ((*p)->nayme, str, MAXNCH);
      break;

    case length:
      (*p)->haslength = true;
      processlength(&valyew, &divisor, ch, &minusread, intree, parens);
      if (!minusread)
        (*p)->length = valyew / divisor;
      else
        (*p)->length = 0.0;

      if ( (*p)->back != NULL)
      {
        (*p)->back->haslength = (*p)->haslength;
        (*p)->back->length = (*p)->length;
      }
      break;

    case hsnolength:
      (*p)->haslength = false;
      if ( (*p)->back != NULL )
      {
        (*p)->back->haslength  = false;
      }
      break;

    default:        /*cases iter, treewt, unttrwt         */
      break;        /*should not occur                */
  }
} /* initretreenode */


void maketriad(tree *treep, node **p, long index)
{
  /* Initiate an internal node with stubs for two children */
  long i, j;
  node *q;
  q = NULL;
  for (i = 1; i <= 3; i++)
  {
    *p = treep->get_forknode(treep, index);
    (*p)->index = index;
    (*p)->hasname = false;
    (*p)->haslength = false;
    (*p)->deleted=false;
    (*p)->deadend=false;
    (*p)->onebranch=false;
    (*p)->onebranchhaslength=false;
    for (j=0;j<MAXNCH;j++)
      (*p)->nayme[j] = '\0';
    (*p)->next = q;
    q = *p;
  }
  (*p)->next->next->next = *p;
  q = (*p)->next;
  while (*p != q)
  {
    (*p)->back = NULL;
    (*p)->tip = false;
    *p = (*p)->next;
  }
  treep->nodep[index - 1] = *p;
}  /* maketriad */


void maketip(tree *treep, node **p, long index)
{
  /*  Initiate a tip node */
  *p = treep->get_forknode(treep, index);
  (*p)->index = index;
  (*p)->tip = true;
  (*p)->hasname = false;
  (*p)->haslength = false;
  treep->nodep[index - 1] = *p;
}  /* maketip */


void copytree(void)
{
  /* Make a complete copy of the current tree for undo purposes */
  if (whichtree == 1)
    othertree = 0;
  else
    othertree = 1;

  generic_tree_copy(treesets[whichtree].tree_p, treesets[othertree].tree_p);

  treesets[othertree].waswritten = waswritten;
  treesets[othertree].hasmult = hasmult;
  treesets[othertree].haslengths = haslengths;
  treesets[othertree].nolengths = nolengths;
  treesets[othertree].initialized = true;

} /* copytree */


void getoptions(void)
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done, gotopt;

  how = use;
  outgrno = 1;
  indentnotchanged = true;
  loopcount = 0;
  onfirsttree = true;
  do {
    cleerhome();
    printf("\nTree Rearrangement, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
    printf("  U          Initial tree (arbitrary, user, specify)?");
    if (how == arb)
      printf("  Arbitrary\n");
    else if (how == use)
      printf("  User tree from tree file\n");
    else
      printf("  Tree you specify\n");
    printf("  F   Format to write out trees (PHYLIP, Nexus, XML)?");
    if (nexus)
      printf("  Nexus\n");
    else
    {
      if (xmltree)
        printf("  PhyloXML\n");
      else
        printf("  PHYLIP\n");
    }
    if (indentnotchanged)   /* default indented if XML; once changed, stays */
      writeindented = xmltree;
    printf("  I      Indent when writing out trees (for clarity)?");
    if (writeindented)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  0               Graphics type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi )
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  W       Width of terminal screen, of plotting area?");
    printf("%4ld, %2ld\n", screenwidth, vscreenwidth);
    printf("  L                        Number of lines on screen?");
    printf("%4ld\n", screenlines);
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    ch = (isupper(ch)) ? ch : toupper(ch);
    done = (ch == 'Y');
    gotopt = (ch == 'U' || ch == 'F' || ch == 'I'
              || ch == '0' || ch == 'L' || ch == 'W');
    if (gotopt)
    {
      switch (ch)
      {
        case 'U':
          if (how == arb)
            how = use;
          else if (how == use)
            how = spec;
          else
            how = arb;
          break;

        case 'F':
          if (nexus)
          {
            nexus = false;
            xmltree = true;
          }
          else if (xmltree)
            xmltree = false;
          else
            nexus = true;
          break;

        case 'I':
          writeindented = !writeindented;
          indentnotchanged = false;
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case 'L':
          initnumlines(&screenlines);
          break;

        case 'W':
          screenwidth= readlong("Width of terminal screen (in characters)?\n");
          vscreenwidth=readlong("Width of plotting area (in characters)?\n");
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


void configure(void)
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)deleted; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)deleted; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  cch[(long)deleted] = '.';
  cch[(long)updel] = ':';
  if (ibmpc)
  {
    cch[(long)horiz] = '>';
    cch[(long)vert] = 186;
    graphic[(long)vert] = true;
    cch[(long)up] = 186;
    graphic[(long)up] = true;
    cch[(long)ch_over] = 205;
    graphic[(long)ch_over] = true;
    cch[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    cch[(long)midcorner] = 204;
    graphic[(long)midcorner] = true;
    cch[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    return;
  }
  if (ansi)
  {
    cch[(long)horiz] = '>';
    cch[(long)vert] = cch[(long)horiz];
    reversed[(long)vert] = true;
    cch[(long)up] = 'x';
    graphic[(long)up] = true;
    cch[(long)ch_over] = 'q';
    graphic[(long)ch_over] = true;
    cch[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    cch[(long)midcorner] = 't';
    graphic[(long)midcorner] = true;
    cch[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  cch[(long)horiz] = '>';
  cch[(long)vert] = ' ';
  cch[(long)up] = '!';
  cch[(long)upcorner] = '`';
  cch[(long)midcorner] = '+';
  cch[(long)downcorner] = ',';
  cch[(long)ch_over] = '-';
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


void ltrav(node *p, boolean *localhl)
{
  /* Traversal function for ifhaslengths() */
  node *q;

  if (p->tip)
  {
    (*localhl) = ((*localhl) && p->haslength);
    return;
  }
  q = p->next;
  do {
    (*localhl) = ((*localhl) && q->haslength);
    if ((*localhl))
      ltrav(q->back, localhl);
    q = q->next;
  } while (p != q);
}  /* ltrav */


boolean ifhaslengths(void)
{
  /* return true if every branch in tree has a length */
  boolean localhl;
  localhl = true;
  ltrav(curtree->root, &localhl);
  return localhl;
}  /* ifhaslengths */


void add_at(tree *treep, node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  node *leftdesc, *rtdesc;
  double length;

  if (below != treep->nodep[below->index - 1])
    below = treep->nodep[below->index - 1];

  if (newfork == NULL)
  {
    nonodes++;
    maketriad (treep, &newfork, nonodes);
    if (haslengths)
    {
      newfork->haslength = true;
      newfork->next->haslength = true;
      newfork->next->next->haslength = true;
    }
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
  if (treep->root == below)
    treep->root = newfork;
  treep->root->back = NULL;
  if (!haslengths)
    return;
  if (newfork->back != NULL)
  {
    length = newfork->back->length / 2.0;
    newfork->length = length;
    newfork->back->length = length;
    below->length = length;
    below->back->length = length;
  }
  else
  {
    length = newtip->length / 2.0;
    newtip->length = length;
    newtip->back->length = length;
    below->length = length;
    below->back->length = length;
    below->haslength = true;
  }
  newtip->back->length = newtip->length;
}  /* add_at */


void add_before(tree *treep, node *atnode, node *newtip)
{
  /* Inserts the node newtip together with its ancestral fork into the tree next to the node atnode. */
  /*xx ?? debug what to do if no ancestral node -- have to create one */
  /*xx this case is handled by add_at.  However, add_at does not account for
    when there is more than one sibling for the relocated newtip */
  node *q;

  if (atnode != treep->nodep[atnode->index - 1])
    atnode = treep->nodep[atnode->index - 1];
  q = treep->nodep[newtip->index-1]->back;
  if (q != NULL)
  {
    q = treep->nodep[q->index-1];
    if (newtip == q->next->next->back)
    {
      q->next->back = newtip;
      newtip->back = q->next;
      q->next->next->back = NULL;
    }
  }
  if (newtip->back != NULL)
  {
    add_at(treep, atnode, newtip, treep->nodep[newtip->back->index-1]);
  }
  else
  {
    add_at(treep, atnode, newtip, NULL);
  }
}  /* add_before */


void add_child(tree *treep, node *parent, node *newchild)
{
  /* adds the node newchild into the tree as the last child of parent */

  int i;
  node *newnode, *q;

  if (parent != treep->nodep[parent->index - 1])
    parent = treep->nodep[parent->index - 1];

  newnode = treep->get_forknode(treep, parent->index);
  newnode->tip = false;
  newnode->deleted=false;
  newnode->deadend=false;
  newnode->onebranch=false;
  newnode->onebranchhaslength=false;

  for (i=0;i<MAXNCH;i++)
    newnode->nayme[i] = '\0';
  newnode->index = parent->index;
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


void re_move(tree *treep, node **item, node **forknode)
{
  /* Removes node item from the tree.  If item has one sibling,
     removes its ancestor, forknode, from the tree as well and attach
     item's sib to forknode's ancestor.  In this case, it returns a pointer
     to the removed forknode which is still attached to item.
  */
  node *p =NULL, *q;
  int nodecount;

  if ((*item)->back == NULL)
  {
    *forknode = NULL;
    return;
  }
  *forknode = treep->nodep[(*item)->back->index - 1];
  nodecount = 0;
  if ((*forknode)->next->back == *item)
    p = *forknode;
  q = (*forknode)->next;
  do {
    nodecount++;
    if (q->next->back == *item)
      p = q;
    q = q->next;
  } while (*forknode != q);

  if (nodecount > 2)
  {
    fromtype = atnode;
    p->next = (*item)->back->next;
    (*item)->back = NULL;
    /*xx*/ *forknode = NULL;
  }
  else
  {
    /* traditional (binary tree) remove code */
    if (*item == (*forknode)->next->back)
    {
      if (treep->root == *forknode)
        treep->root = (*forknode)->next->next->back;
    }
    else
    {
      if (treep->root == *forknode)
        treep->root = (*forknode)->next->back;
    }
    fromtype = beforenode;
    /* stitch nodes together, leaving out item */
    p = (*item)->back->next->back;
    q = (*item)->back->next->next->back;
    if (p != NULL)
      p->back = q;
    if (q != NULL)
      q->back = p;
    if (haslengths)
    {
      if (p != NULL && q != NULL)
      {
        p->length += q->length;
        q->length = p->length;
      }
      else
        (*item)->length = (*forknode)->next->length + (*forknode)->next->next->length;
    }
    (*forknode)->back = NULL;
    p = (*forknode)->next;
    while (p != *forknode)
    {
      p->back = NULL;
      p = p->next;
    }
    (*item)->back = NULL;
  } /* endif nodecount > 2 else */
}  /* re_move */


void reroot(tree *treep, node *outgroup)
{
  /* Reorient tree so that outgroup is by itself on the left of the root */
  node *p, *q, *r;
  long nodecount = 0;
  double templen;

  q = treep->root->next;
  do {                    /* when this loop exits, p points to the internal */
    p = q;                /* node to the right of root */
    nodecount++;
    q = p->next;
  } while (q != treep->root);
  r = p;

  /* There is no point in proceeding if
       1. outgroup is a child of root, and
       2. the tree bifurcates at the root.
  */
  if((outgroup->back->index == treep->root->index) && !(nodecount > 2))
    return;

  /* reorient nodep array

     The nodep array must point to the ring member of each ring
     that is closest to the root.  The while loop changes the ring member
     pointed to by nodep[] for those nodes that will have their
     orientation changed by the reroot operation.
  */
  p = outgroup->back;
  while (p->index != treep->root->index)
  {
    q = treep->nodep[p->index - 1]->back;
    treep->nodep[p->index - 1] = p;
    p = q;
  }
  if (nodecount > 2)
    treep->nodep[p->index - 1] = p;

  /* If nodecount > 2, the current node ring to which root is pointing
     will remain in place and root will point somewhere else. */
  /* detach root from old location */
  if (nodecount > 2)
  {
    r->next = treep->root->next;
    treep->root->next = NULL;
    nonodes++;
    maketriad(treep, &treep->root, nonodes);

    if (haslengths)
    {
      /* root->haslength remains false, or else treeout() will generate a bogus extra length */
      treep->root->next->haslength = true;
      treep->root->next->next->haslength = true;
    }
  }
  else
  { /* if (nodecount > 2) else */
    q = treep->root->next;
    q->back->back = r->back;
    r->back->back = q->back;

    if (haslengths)
    {
      r->back->length = r->back->length + q->back->length;
      q->back->length = r->back->length;
    }
  } /* if (nodecount > 2) endif */

  /* tie root into new location */
  treep->root->next->back = outgroup;
  treep->root->next->next->back = outgroup->back;
  outgroup->back->back = treep->root->next->next;
  outgroup->back = treep->root->next;

  /* place root equidistant between left child (outgroup) and
     right child by dividing outgroup's length */
  if (haslengths)
  {
    templen = outgroup->length / 2.0;
    outgroup->length = templen;
    outgroup->back->length = templen;
    treep->root->next->next->length = templen;
    treep->root->next->next->back->length = templen;
  }
} /* reroot */


void ltrav_(tree *treep, node *p, double lengthsum, double lmin, double *tipmax, long *across, long *maxchar)
{
  node *q;
  long rchar, nl;
  double sublength;

  if (p->tip)
  {
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    if (lmin == 0.0)
      return;
    rchar = (long)(lengthsum / (*tipmax) * (*across) + 0.5);

    nl = strlen(treep->nodep[p->index - 1]->nayme);
    if (rchar + nl > (*maxchar))
      (*across) = (*maxchar) - (long)(nl * (*tipmax) / lengthsum + 0.5);
    return;
  }
  q = p->next;
  do {
    if (q->length >= lmin)
      sublength = q->length;
    else
      sublength = lmin;
    ltrav_(treep, q->back, lengthsum + sublength, lmin, tipmax, across, maxchar);
    q = q->next;
  } while (p != q);
}  /* ltrav */


void precoord(node *nuroot, boolean *subtree, double *tipmax, long *across)
{
  /* set tipmax and across so that tree is scaled to screenwidth */

  double oldtipmax, minimum;
  long i, maxchar;

  (*tipmax) = 0.0;
  if ((*subtree))
    maxchar = screenwidth - 13;
  else
    maxchar = screenwidth - 5;
  (*across) = maxchar;
  ltrav_(curtree, nuroot, 0.0, 0.0, tipmax, across, &maxchar);
  i = 0;
  do {
    oldtipmax = (*tipmax);
    minimum = 3.0 / (*across) * (*tipmax);
    ltrav_(curtree, nuroot, 0.0, minimum, tipmax, across, &maxchar);
    i++;
  } while (fabs((*tipmax) - oldtipmax) > 0.01 * oldtipmax && i <= 40);
}  /* precoord */


void coordinates(node *p, double lengthsum, long *across, long *tipy,
                 double *tipmax)
{
  /* establishes coordinates of nodes for display with lengths */
  node *q, *first, *last;

  if (p->tip)
  {
    p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, lengthsum + q->length, across, tipy, tipmax);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
  if (p == curtree->root)
  {
    if (curtree->root->next->next->next == curtree->root)
      p->ycoord = (first->ycoord + last->ycoord) / 2;
    else
      p->ycoord = p->next->next->back->ycoord;
  }
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void flatcoordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes for display without lengths */
  node *q, *first, *last;

  if (p->tip)
  {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy) += down;
    return;
  }
  q = p->next;
  do {
    flatcoordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (last->ymax - first->ymin) * 3 / 2;
  if (p == curtree->root)
  {
    if (curtree->root->next->next->next == curtree->root)
      p->ycoord = (first->ycoord + last->ycoord) / 2;
    else
      p->ycoord = p->next->next->back->ycoord;
  }
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* flatcoordinates */


void grwrite(chartype c, long num, long *pos)
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++)
  {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(cch[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */


void drawline(long i, node *nuroot, boolean *subtree)
{
  /* draws one row of the tree diagram by moving up tree */
  long pos;
  node *p, *q, *r, *s, *first =NULL, *last =NULL;
  long n, j;
  long up_nondel, down_nondel;
  boolean extra, done;
  chartype c, d;
  pos = 1;
  p = nuroot;
  q = nuroot;

  extra = false;
  if (i == (long)p->ycoord && (p == curtree->root || (*subtree)))
  {
    c = ch_over;
    if ((*subtree))
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
    extra = true;
  }
  else
  {
    if ((*subtree))
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
    if (haslengths && !nolengths)
      n = (long)(q->xcoord - p->xcoord);
    else
      n = (long)(p->xcoord - q->xcoord);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra)
    {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done)
    {
      c = ch_over;
      if (!haslengths && !q->haslength)
        c = horiz;
      if (q->deleted)
        c = deleted;
      if (q == first)
        d = downcorner;
      else if (q == last)
        d = upcorner;
      else if ((long)q->ycoord == (long)p->ycoord)
        d = c;
      else
        d = midcorner;
      if (n > 1 || q->tip)
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
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          i != (long)p->ycoord)
      {
        c = up;
        if(p->deleted)
          c = updel;
        if (!p->tip)
        {
          up_nondel = 0;
          down_nondel = 0;
          r = p->next;
          do {
            s = r->back;
            if ((long)s->ycoord < (long)p->ycoord && !s->deleted)
              up_nondel = (long)s->ycoord;
            if (s->ycoord > p->ycoord && !s->deleted &&
                (down_nondel == 0))
              down_nondel = (long)s->ycoord;
            if (i < (long)p->ycoord && s->deleted && i > (long)s->ycoord)
              c = updel;
            if (i > (long)p->ycoord && s->deleted && i < (long)s->ycoord)
              c = updel;
            r = r->next;
          } while (r != p);

          if ((up_nondel != 0) && i < (long)p->ycoord && i > up_nondel)
            c = up;
          if ((down_nondel != 0) && i > (long)p->ycoord && i < down_nondel)
            c = up;
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
  if ((long)p->ycoord == i && p->tip)
  {
    if (p->hasname)
    {
      n = 0;
      for (j = 1; j <= MAXNCH; j++)
      {
        if (curtree->nodep[p->index - 1]->nayme[j - 1] != '\0')
          n = j;
      }
      chwrite(':', 1, &pos, leftedge, screenwidth);
      for (j = 0; j < n; j++)
        chwrite(curtree->nodep[p->index - 1]->nayme[j], 1, &pos, leftedge, screenwidth);
    }
  }
  putchar('\n');
}  /* drawline */


void printree(void)
{
  /* prints out diagram of the tree */
  long across;
  long tipy;
  double tipmax;
  long i, dow, vmargin;

  haslengths = ifhaslengths();
  if (!subtree)
    nuroot = curtree->root;
  cleerhome();
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree)
  {
    dow--;
  }
  if (haslengths && !nolengths)
  {
    precoord(nuroot, &subtree, &tipmax, &across);
    /* protect coordinates() from div/0 errors if user decides to
       examine a tip as a subtree */
    if (tipmax == 0)
      tipmax = 0.01;
    coordinates(nuroot, 0.0, &across, &tipy, &tipmax);
  }
  else
    flatcoordinates(nuroot, &tipy);
  vmargin = 2;
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
      drawline(i, nuroot, &subtree);
  }
  if (leftedge > 1)
    printf("** %ld characters to left of screen ", leftedge);
  if ((treelines - topedge + 1) > (screenlines - vmargin))
  {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
}  /* printree */


void togglelengths(void)
{
  nolengths = !nolengths;
  printree();
}  /* togglengths */


void arbitree(tree *treep)
{
  long i, maxinput;
  node *newtip, *newfork;

  maxinput = 1;
  do {
    spp = readlong("How many species?\n");
    maxinput++;
    if (maxinput == 100)
    {
      printf("ERROR:  Too many tries at choosing species.\n");
      exxit(-1);
    }
  } while (spp <= 0);
  nonodes = spp * 2 - 1;
  maketip(treep, &treep->root, 1);
  maketip(treep, &newtip, 2);
  maketriad(treep, &newfork, spp + 1);
  add_at(treep, treep->root, newtip, newfork);
  for (i = 3; i <= spp; i++)
  {
    maketip(treep, &newtip, i);
    maketriad(treep, &newfork, spp + i - 1);
    add_at(treep, treep->nodep[spp + i - 3], newtip, newfork);
  }
}  /* arbitree */


void yourtree(void)
{
  long uniquearray[maxsz];
  long uniqueindex = 0;
  long i, j, k, k_max, maxinput;
  boolean ok, done;
  node *newtip, *newfork;
  Char ch;

  uniquearray[0] = 0;
  spp = 2;
  nonodes = spp * 2 - 1;
  maketip(curtree, &curtree->root, 1);
  maketip(curtree, &newtip, 2);
  maketriad(curtree, &newfork, spp + 3);
  add_at(curtree, curtree->root, newtip, newfork);
  i = 2;
  maxinput = 1;
  k_max = 5;
  do {
    i++;
    printree();
    printf("Enter 0 to stop building tree.\n");
    printf("Add species%3ld", i);
    do {
      printf("\n at or before which node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && (j < i || (j > spp + 2 && j < spp + i + 1)));
      if (!ok)
        printf("Impossible number. Please try again:\n");
      maxinput++;
      if (maxinput == 100)
      {
        printf("ERROR:  Too many tries at choosing number.\n");
        exxit(-1);
      }
    } while (!ok);
    maxinput = 1;
    if (j >= i)                         /* has user chosen a non-tip? if so, offer choice */
    {
      do {
        printf(" Insert at node (A) or before node (B)? ");
        phyFillScreenColor();
        if(scanf("%c%*[^\n]", &ch)) {}  // Read char and scan to EOL.
        (void)getchar();
        if (ch == '\n')
          ch = ' ';
        ch = isupper(ch) ? ch : toupper(ch);
        maxinput++;
        if (maxinput == 100)
        {
          printf("ERROR:  Too many tries at choosing option.\n");
          exxit(-1);
        }
      } while (ch != 'A' && ch != 'B');
    }
    else ch = 'B';   /* if user has chosen a tip, set Before */

    if (j != 0)
    {
      if (ch == 'A')
      {
        if (!curtree->nodep[j - 1]->tip)
        {
          maketip(curtree, &newtip, i);
          add_child(curtree, curtree->nodep[j - 1], curtree->nodep[i - 1]);
        }
      }
      else
      {
        maketip(curtree, &newtip, i);
        maketriad(curtree, &newfork, spp + i + 1);
        curtree->nodep[i-1]->back = newfork;
        newfork->back = curtree->nodep[i-1];
        add_before(curtree, curtree->nodep[j - 1], curtree->nodep[i - 1]);
      } /* endif (before or at node) */
    }

    done = (j == 0);
    if (!done)
    {
      if (ch == 'B')
        k = spp * 2 + 3;
      else
        k = spp * 2 + 2;

      k_max = k;
      do {
        if (curtree->nodep[k - 2] != NULL)
        {
          curtree->nodep[k - 1] = curtree->nodep[k - 2];
          curtree->nodep[k - 1]->index = k;
          curtree->nodep[k - 1]->next->index = k;
          curtree->nodep[k - 1]->next->next->index = k;
        }
        k--;
      } while (k != spp + 3);
      if (j > spp + 1)
        j++;
      spp++;
      nonodes = spp * 2 - 1;
    }
  } while (!done);

  for (i = spp + 1; i <= k_max; i++)
  {
    if ((curtree->nodep[i - 1] != curtree->nodep[i]) && (curtree->nodep[i - 1] != NULL))
    {
      uniquearray[uniqueindex++] = i;
      uniquearray[uniqueindex] = 0;
    }
  }

  for ( i = 0; uniquearray[i] != 0; i++)
  {
    curtree->nodep[spp + i] = curtree->nodep[uniquearray[i] - 1];
    curtree->nodep[spp + i]->index = spp + i + 1;
    curtree->nodep[spp + i]->next->index = spp + i + 1;
    curtree->nodep[spp + i]->next->next->index = spp + i + 1;
  }
  for (i = spp + uniqueindex; i <= k_max; i++)
    curtree->nodep[i] = NULL;

  nonodes = spp * 2 - 1;
}  /* yourtree */


void buildtree(void)
{
  /* variables needed to be passed to treeread() */
  long    nextnode   = 0;
  pointarray dummy_treenode = NULL;  /* Ignore what happens to this */
  boolean goteof     = false;
  boolean haslengths = false;
  boolean firsttree;
  node *p, *q;
  long nodecount = 0;


  /* These assignments moved from treeconstruct -- they seem to happen only here. */
  /*xx treeone & treetwo assignments should probably happen in treeconstruct.  Memory leak if user reads multiple trees. */
  /*
    treeone = (node **)Malloc(maxsz * sizeof(node *));
    treetwo = (node **)Malloc(maxsz * sizeof(node *));
  */
  // BUG.967 -- these are default values below, this should be
  // done more carefully
  nonodes = 100;
  spp = 50;

  // BUG.967 these tree_new replace the Malloc's above
  treeone = functions.tree_new(nonodes, spp);
  treetwo = functions.tree_new(nonodes, spp);
  treesets[whichtree].tree_p = treeone;
  treesets[othertree].tree_p = treetwo;
  curtree = treeone;

  simplifiedtree.tree_p = functions.tree_new(nonodes, spp);
  subtree     = false;
  topedge     = 1;
  leftedge    = 1;
  switch (how)
  {
    case arb:
      treesets[whichtree].tree_p = treeone;
      treesets[othertree].tree_p = treetwo;
      curtree = treeone;
      arbitree(curtree);
      break;

    case use:
      printf("\nReading tree file ...\n\n");
      fflush(stdout);
      if (!readnext)
      {
        /* This is the first time through here, act accordingly */
        firsttree = true;
        /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
        if (!javarun)
        {
            openfile(&intree, INTREE, "input tree file", "rb", "retree", intreename);
        }
        numtrees = countsemic(intree);
        treesread = 0;
      }
      else
      {
        /* This isn't the first time through here ... */
        firsttree = false;
      }
      // BUG.967 these may need to be re-incorporated in a more modern way
      // allocate_nodep(&nodep, intree, &spp);
      // treesets[whichtree].nodep = nodep;
      printf("numtrees: %li\n", numtrees);
      fflush(stdout);

      if (firsttree)
        nayme = (naym *)Malloc(spp * sizeof(naym));
      printf("calling treeread\n");
      fflush(stdout);
      treeread(treesets[whichtree].tree_p, intree, &(treesets[whichtree].tree_p->root), dummy_treenode, &goteof, &firsttree, &nextnode, &haslengths, initretreenode, true, -1);
      nonodes = nextnode;
      treesread++;
      treesets[othertree].tree_p = treetwo;
      break;

    case spec:
      // BUG.967 nodep = treeone; // can this be removed?
      treesets[othertree].tree_p = treetwo;
      yourtree();
      break;
  }

  q = treesets[whichtree].tree_p->root->next;
  do {
    p = q;
    nodecount++;
    q = p->next;
  } while (q != treesets[whichtree].tree_p->root);

  outgrno = treesets[whichtree].tree_p->root->next->back->index;
  if(!(nodecount > 2))
  {
    reroot(curtree, curtree->nodep[outgrno - 1]);
  }
}  /* buildtree */


void unbuildtree(void)
{
  // BUG.967 -- need to test and possibly replace
  // treesets[0].tree_p->free(treesets[0].tree_p);
  // treesets[0].tree_p = NULL;
  // treesets[1].tree_p->free(treesets[1].tree_p);
  // treesets[1].tree_p = NULL;
  // simplifiedtree.tree_p->free(treesets[1].tree_p);
  // simplifiedtree.tree_p = NULL;
}  // unbuildtree


void retree_help(void)
{
  /* display help information */
  char tmp[100];
  printf("\n\n . Redisplay the same tree again\n");
  if (haslengths)
  {
    printf(" = Redisplay the same tree with");
    if (!nolengths)
      printf("out/with");
    else
      printf("/without");
    printf(" lengths\n");
  }
  printf(" U Undo the most recent change in the tree\n");
  printf(" W Write tree to a file\n");
  printf(" + Read next tree from file (may blow up if none is there)\n");
  printf("\n");
  printf(" R Rearrange a tree by moving a node or group\n");
  printf(" O select an Outgroup for the tree\n");
  if (haslengths)
    printf(" M Midpoint root the tree\n");
  printf(" T Transpose immediate branches at a node\n");
  printf(" F Flip (rotate) subtree at a node\n");
  printf(" D Delete or restore nodes\n");
  printf(" B Change or specify the length of a branch\n");
  printf(" N Change or specify the name(s) of tip(s)\n");
  printf(" H Move viewing window to the left\n");
  printf(" J Move viewing window downward\n");
  printf(" K Move viewing window upward\n");
  printf(" L Move viewing window to the right\n");
  printf(" C collapse/uncollapse a clade so it looks as if it were a tip\n");
  printf(" Y show only one clade (subtree) (might be useful if tree is ");
  printf("too big)\n");
  printf(" ? Help (this screen)\n");
  printf(" Q (Quit) Exit from program\n");
  printf(" X Exit from program\n\n");
  printf(" TO CONTINUE, PRESS ON THE Return OR Enter KEY");
  getstryng(tmp);
  printree();
}  /* retree_help */


void consolidatetree(long index)
{
  node *start, *r, *q;
  int i;

  start = curtree->nodep[index - 1];
  q = start->next;

  while (q != start)
  {
    r = q;
    q = q->next;
  }

  i = index;
  while (curtree->nodep[i-1] != NULL)
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
  boolean ok;
  node *p, *q;
  char ch;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok);
  if ( ok == false )
  {
    /* fall through */
  }
  else if ( i < 1 || i > spp*2 - 1 )
  {
    /* i is not in range */
    ok = false;
  }
  else if (i == curtree->root->index )
  {
    /* i is root */
    ok = false;
  }
  else if ( curtree->nodep[i-1]->deleted )
  {
    /* i has been deleted */
    ok = false;
  }
  else
  {
    printf("Add at or before which node? ");
    inpnum(&j, &ok);
    if ( ok == false )
    {
      /* fall through */
    }
    else if ( j < 1 || j > spp*2 - 1 )
    {
      /* j is not in range */
      ok = false;
    }
    else if ( curtree->nodep[j-1]->deleted )
    {
      /* j has been deleted */
      ok = false;
    }
    else if (j != curtree->root->index && curtree->nodep[curtree->nodep[j-1]->back->index - 1]->deleted )
    {
      /* parent of j has been deleted */
      ok = false;
    }
    else if ( curtree->nodep[j-1] == curtree->nodep[curtree->nodep[i-1]->back->index -1] )
    {
      /* i is j's parent */
      ok = false;
    }
    else
    {
      /* make sure that j is not a descendant of i */
      for ( p = curtree->nodep[j-1]; p != curtree->root; p = curtree->nodep[p->back->index - 1] )
      {
        if ( p == curtree->nodep[i-1] )
        {
          ok = false;
          break;
        }
      }
      if ( ok )
      {
        maxinput = 1;
        do {
          printf("Insert at node (A) or before node (B)? ");
          phyFillScreenColor();
          if(scanf("%c%*[^\n]", &ch)) {} // Read char and scan to EOL.
          (void)getchar();
          if (ch == '\n')
            ch = ' ';
          ch = toupper(ch);
          maxinput++;
          if (maxinput == 100)
          {
            printf("ERROR:  Input failed too many times.\n");
            exxit(-1);
          }
        } while (ch != 'A' && ch != 'B');

        if (ch == 'A')
        {
          if ( curtree->nodep[j - 1]->deleted || curtree->nodep[j - 1]->tip )
          {
            /* If j is a tip or has been deleted */
            ok = false;
          }
          else if ( curtree->nodep[j-1] == curtree->nodep[curtree->nodep[i-1]->back->index -1] )
          {
            /* If j is i's parent */
            ok = false;
          }
          else
          {
            copytree();
            re_move(curtree, &curtree->nodep[i - 1], &q);
            add_child(curtree, curtree->nodep[j - 1], curtree->nodep[i - 1]);
            if (fromtype == beforenode)
              consolidatetree(q->index);
          }
        }
        else
        { /* ch == 'B' */
          if (j == curtree->root->index) /* can't insert at root */
          {
            ok = false;
          }
          else
          {
            copytree();
            printf("Insert before node %ld.\n", j);
            re_move(curtree, &(curtree->nodep[i - 1]), &q);
            if (q != NULL)
            {
              curtree->nodep[q->index-1]->next->back = curtree->nodep[i-1];
              curtree->nodep[i-1]->back = curtree->nodep[q->index-1]->next;
            }
            add_before(curtree, curtree->nodep[j - 1], curtree->nodep[i - 1]);
          }
        } /* endif (before or at node) */
      } /* endif (ok to do move) */
    } /* endif (destination node ok) */
  } /* endif (from node ok) */

  printree();

  if ( !ok )
    printf("Not a possible rearrangement.  Try again: \n");
  else
  {
    written = false;
  }
}  /* rearrange */


boolean any_deleted(node *p)
{
  /* return true if there are any deleted branches from branch on down */
  boolean localdl;
  localdl = false;
  ifdeltrav(p, &localdl);
  return localdl;
}  /* any_deleted */


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
    if (ok)
      ok = !any_deleted(curtree->nodep[i - 1]);
  }
  else
  {
    i = atnode;
    ok = true;
  }
  if (ok)
  {
    copytree();
    fliptrav(curtree->nodep[i - 1], true);
  }
  if (atnode == 0)
    printree();
  if (ok)
  {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) || (i > spp && i <= nonodes && any_deleted(curtree->nodep[i - 1])))
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */


void transpose(long atnode)
{
  /* flip at a node left-right */
  long i;
  boolean ok;

  if (atnode == 0)
  {
    printf("Transpose branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
    if (ok)
      ok = !curtree->nodep[i - 1]->deleted;
  }
  else
  {
    i = atnode;
    ok = true;
  }
  if (ok)
  {
    copytree();
    fliptrav(curtree->nodep[i - 1], false);
  }
  if (atnode == 0)
    printree();
  if (ok)
  {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) || (i > spp && i <= nonodes && curtree->nodep[i - 1]->deleted))
    printf("Can't transpose there. ");
  else
    printf("No such node. ");
}  /* transpose */


void ifdeltrav(node *p, boolean *localdl)
{
  node *q;

  if (*localdl)
    return;

  if (p->tip)
  {
    (*localdl) = ((*localdl) || p->deleted);
    return;
  }
  q = p->next;
  do {
    (*localdl) = ((*localdl) || q->deleted);
    ifdeltrav(q->back, localdl);
    q = q->next;
  } while (p != q);
}  /* ifdeltrav */


double oltrav(node *p)
{
  node *q;
  double maxlen, templen;
  if (p->deleted)
    return 0.0;
  if (p->tip)
  {
    p->beyond = 0.0;
    return 0.0;
  }
  else
  {
    q = p->next;
    maxlen = 0;
    do {
      templen = q->back->deleted ? 0.0 : q->length + oltrav(q->back);
      maxlen = (maxlen > templen) ? maxlen : templen;
      q->beyond = templen;
      q = q->next;
    } while (p != q);
    p->beyond = maxlen;
    return (maxlen);
  }
}  /* oltrav */


void outlength(void)
{
  /* compute the farthest combined length out from each node */
  oltrav(curtree->root);
}  /* outlength */


void midpoint(void)
{
  /* midpoint root the tree */
  double balance, greatlen, lesslen, grlen, maxlen;
  node *maxnode, *grnode, *lsnode =NULL;
  boolean ok = true;
  boolean changed = false;
  node *p, *q;
  long nodecount = 0;
  boolean multi = false;

  copytree();
  p = curtree->root;
  outlength();
  q = p->next;
  greatlen = 0;
  grnode = q->back;
  lesslen = 0;

  q = curtree->root->next;
  do {
    p = q;
    nodecount++;
    q = p->next;
  } while (q != curtree->root);
  if (nodecount > 2)
    multi = true;

  /* Find the two greatest lengths reaching from root to tips.
     Also find the lengths and node pointers of the first nodes in the
     direction of those two greatest lengths.  */
  p = curtree->root;
  q = curtree->root->next;
  do {
    if (greatlen <= q->beyond)
    {
      lesslen = greatlen;
      lsnode = grnode;
      greatlen = q->beyond;
      grnode = q->back;
    }
    if ((greatlen > q->beyond) && (q->beyond >= lesslen))
    {
      lesslen = q->beyond;
      lsnode = q->back;
    }
    q = q->next;
  } while (p != q);

  /* If we don't have two non-deleted nodes to balance between then we can't midpoint root the tree */
  if (grnode->deleted || lsnode->deleted || grnode == lsnode)
    ok = false;
  balance = greatlen - (greatlen + lesslen) / 2.0;
  grlen = grnode->length;

  while ((balance - grlen > 1e-10) && ok)
  {
    /* First, find the most distant immediate child of grnode and reroot to it. */
    p = grnode;
    q = p->next;
    maxlen = 0;
    maxnode = q->back;
    do {
      if (maxlen <= q->beyond)
      {
        maxlen = q->beyond;
        maxnode = q->back;
      }
      q = q->next;
    } while (p != q);
    reroot(curtree, maxnode);
    changed = true;

    /* Reassess the situation, using the same "find the two greatest
       lengths" code as occurs before the while loop.  If another reroot
       is necessary, this while loop will repeat. */
    p = curtree->root;
    outlength();
    q = p->next;
    greatlen = 0;
    grnode = q->back;
    lesslen = 0;
    do {
      if (greatlen <= q->beyond)
      {
        lesslen = greatlen;
        lsnode = grnode;
        greatlen = q->beyond;
        grnode = q->back;
      }
      if ((greatlen > q->beyond) && (q->beyond >= lesslen))
      {
        lesslen = q->beyond;
        lsnode = q->back;
      }
      q = q->next;
    } while (p != q);
    if (grnode->deleted || lsnode->deleted || grnode == lsnode)
      ok = false;
    balance = greatlen - (greatlen + lesslen) / 2.0;
    grlen = grnode->length;
  }; /* end of while ((balance > grlen) && ok) */

  if (ok)
  {
    /*xx the following ignores deleted nodes */
    /*   this may be ok because deleted nodes are omitted from length calculations */
    if (multi)
    {
      reroot(curtree, grnode); /*xx need length corrections */

      p = curtree->root;
      outlength();
      q = p->next;
      greatlen = 0;
      grnode = q->back;
      lesslen = 0;

      do {
        if (greatlen <= q->beyond)
        {
          lesslen = greatlen;
          lsnode = grnode;
          greatlen = q->beyond;
          grnode = q->back;
        }
        if ((greatlen > q->beyond) && (q->beyond >= lesslen))
        {
          lesslen = q->beyond;
          lsnode = q->back;
        }
        q = q->next;
      } while (p != q);
      balance = greatlen - (greatlen + lesslen) / 2.0;
    }
    grnode->length -= balance;
    if (((grnode->length) < 0.0) && (grnode->length > -1.0e-10))
      grnode->length = 0.0;
    grnode->back->length = grnode->length;
    if (((lsnode->length) < 0.0) && (lsnode->length > -1.0e-10))
      lsnode->length = 0.0;
    lsnode->length += balance;
    lsnode->back->length = lsnode->length;
  }
  printree();
  if (ok)
  {
    if (any_deleted(curtree->root))
      printf("Deleted nodes were not used in midpoint calculations.\n");
  }
  else
  {
    printf("Can't perform midpoint because of deleted branches.\n");
    if (changed)
    {
      undo();
      printf("Tree restored to original state.  Undo information lost.\n");
    }
  }
} /* midpoint */


void deltrav(node *p, boolean value)
{
  /* register p and p's children as deleted or extant, depending on value */
  node *q;

  p->deleted = value;

  if (p->tip)
    return;

  q = p->next;
  do {
    deltrav(q->back, value);
    q = q->next;
  } while (p != q);
}  /* deltrav */


void fill_del(node*p)
{
  int alldell;
  node *q = p;

  if ( p->next == NULL) return;

  q=p->next;
  while ( q != p)
  {
    fill_del(q->back);
    q=q->next;
  }

  alldell = 1;

  q=p->next;
  while ( q != p)
  {
    if ( !q->back->deleted )
    {
      alldell = 0;
    }
    q=q->next;
  }

  p->deleted = alldell;
}


void reg_del(node *delp, boolean value)
{
  /* register delp and all of delp's children as deleted */
  deltrav(delp, value);
  fill_del(curtree->root);
}  /* reg_del */


boolean isdeleted(long nodenum)
{
  /* true if nodenum is a node number in a deleted branch */
  return(curtree->nodep[nodenum - 1]->deleted);
} /* isdeleted */


void deletebranch(void)
{
  /* delete a node */
  long i;
  boolean ok1;

  printf("Delete everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i <= nonodes && i != curtree->root->index && !isdeleted(i));
  if (ok1)
  {
    copytree();
    reg_del(curtree->nodep[i - 1], true);
  }
  printree();
  if (!ok1)
    printf("Not a possible deletion.  Try again.\n");
  else
  {
    written = false;
  }
}  /* deletebranch */


void restorebranch(void)
{
  /* restore deleted branches */
  long i;
  boolean ok1;

  printf("Restore everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i < spp * 2 && isdeleted(i) && ( i == curtree->root->index || !curtree->nodep[curtree->nodep[i - 1]->back->index - 1]->deleted));

  if (ok1)
  {
    reg_del(curtree->nodep[i - 1], false);
  }
  printree();
  if (!ok1)
    printf("Not a possible restoration.  Try again: \n");
  else
  {
    written = false;
  }
} /* restorebranch */


void del_or_restore(void)
{
  /* delete or restore a branch */
  long maxinput;
  Char ch;

  if (any_deleted(curtree->root))
  {
    maxinput = 1;
    do {
      printf("Enter D to delete a branch\n");
      printf("OR enter R to restore a branch: ");
      phyFillScreenColor();
      if(scanf("%c%*[^\n]", &ch)) {}    // Read char and scan to EOL.
      (void)getchar();
      if (ch == '\n')
        ch = ' ';
      ch = (isupper(ch)) ? ch : toupper(ch);
      maxinput++;
      if (maxinput == 100)
      {
        printf("ERROR:  Too many tries at choosing option.\n");
        exxit(-1);
      }
    } while (ch != 'D' && ch != 'R');
    if (ch == 'R')
      restorebranch();
    else
      deletebranch();
  }
  else
    deletebranch();
} /* del_or_restore */


void undo(void)
{
  /* don't undo to an uninitialized tree */
  if (!treesets[othertree].initialized)
  {
    printree();
    printf("Nothing to undo.\n");
    return;
  }

  treesets[whichtree].tree_p = curtree;

  treesets[whichtree].waswritten = waswritten;
  treesets[whichtree].hasmult = hasmult;
  treesets[whichtree].haslengths = haslengths;
  treesets[whichtree].nolengths = nolengths;
  treesets[whichtree].initialized = true;

  whichtree = othertree;

  curtree = treesets[whichtree].tree_p;

  waswritten = treesets[whichtree].waswritten;
  hasmult = treesets[whichtree].hasmult;
  haslengths = treesets[whichtree].haslengths;
  nolengths = treesets[whichtree].nolengths;

  if (othertree == 0)
    othertree = 1;
  else
    othertree = 0;

  printree();
}  /* undo */


/*
   These attributes of nodes in the tree are modified by treetrav()
   in preparation for writing a tree to disk.

   boolean deadend      This node is not deleted but all of its children
     are, so this node will be treated as such when
     the tree is written or displayed.

   boolean onebranch    This node has only one valid child, so that this
     node will not be written and its child will be
     written as a child of its grandparent with the
     appropriate summing of lengths.

   nodep *onebranchnode
     Used if onebranch is true.  Onebranchnode points
     to the one valid child.  This child may be one or
     more generations down from the current node.

   double onebranchlength
     Used if onebranch is true.  Onebranchlength is
     the length from the current node to the valid
     child.
*/


void treetrav(node *p)
{
  long branchcount = 0;
  node *q, *onebranchp =NULL;

  /* Count the non-deleted branches hanging off of this node into branchcount.
     If there is only one such branch, onebranchp points to that branch. */

  if (p->tip)
    return;

  q = p->next;
  do {
    if (!q->back->deleted)
    {
      if (!q->back->tip)
        treetrav(q->back);

      if (!q->back->deadend && !q->back->deleted)
      {
        branchcount++;
        onebranchp = q->back;
      }
    }
    q = q->next;
  } while (p != q);

  if (branchcount == 0)
    p->deadend = true;
  else
    p->deadend = false;
  p->onebranch = false;
  if (branchcount == 1 && onebranchp->tip)
  {
    p->onebranch = true;
    p->onebranchnode = onebranchp;
    p->onebranchhaslength = (p->haslength || (p == curtree->root)) && onebranchp->haslength;
    if (p->onebranchhaslength)
      p->onebranchlength = onebranchp->length + p->length;
  }
  if (branchcount == 1 && !onebranchp->tip)
  {
    p->onebranch = true;
    if (onebranchp->onebranch)
    {
      p->onebranchnode = onebranchp->onebranchnode;
      p->onebranchhaslength = (p->haslength || (p == curtree->root)) && onebranchp->onebranchhaslength;
      if (p->onebranchhaslength)
        p->onebranchlength = onebranchp->onebranchlength + p->length;
    }
    else
    {
      p->onebranchnode = onebranchp;
      p->onebranchhaslength = p->haslength && onebranchp->haslength;
      if (p->onebranchhaslength)
        p->onebranchlength = onebranchp->length + p->length;
    }
  }
} /* treetrav */


void simcopynode(node *fromnode, node *tonode)
{
  /* Copy the contents of a node from fromnode to tonode. */
  int i;

  tonode->index   = fromnode->index;
  tonode->deleted = fromnode->deleted;
  tonode->tip     = fromnode->tip;
  tonode->hasname = fromnode->hasname;
  if (fromnode->hasname)
    for (i=0;i<MAXNCH;i++)
      tonode->nayme[i] = fromnode->nayme[i];
  tonode->haslength = fromnode->haslength;
  if (fromnode->haslength)
    tonode->length = fromnode->length;

} /* simcopynode */


// BUG.967 -- we are mixing two trees here -- careful
node *simcopytrav(tree * srctree, node *srcNode, tree * desttree)
{
  /* Traverse the tree from p on down, copying nodes to the other tree */
  node *q, *newnode, *newnextnode, *temp;
  long lastnodeidx = 0;

  newnode = desttree->get_forknode(desttree, srcNode->index); // BUG.967 worried about index
  simcopynode(srcNode, newnode);

  if (srctree->nodep[srcNode->index - 1] == srcNode)
    desttree->nodep[srcNode->index - 1] = newnode;

  /* if this is a tip, return now */
  if (srcNode->tip)
  {
    return newnode;
  }
  if (srcNode->onebranch && srcNode->onebranchnode->tip)
  {
    simcopynode(srcNode->onebranchnode, newnode);
    if (srcNode->onebranchhaslength)
      newnode->length = srcNode->onebranchlength;
    return newnode;
  }
  else if (srcNode->onebranch && !srcNode->onebranchnode->tip)
  {
    /* recurse down p->onebranchnode */
    srcNode->onebranchnode->length = srcNode->onebranchlength;
    srcNode->onebranchnode->haslength = srcNode->onebranchnode->haslength;
    return simcopytrav(srctree, srcNode->onebranchnode, desttree);
  }
  else
  {
    /* Multiple non-deleted branch case:  go round the node recursing
       down the branches. Don't go down deleted branches or dead ends. */
    q = srcNode->next;
    while (q != srcNode)
    {
      if (!q->back->deleted && !q->back->deadend)
        lastnodeidx = q->back->index;
      q = q->next;
    }

    q = srcNode->next;
    newnextnode = desttree->get_forknode(desttree, q->index); // BUG.967 worried about index
    simcopynode(q, newnextnode);
    newnode->next = newnextnode;
    do {
      /* If branch is deleted or is a dead end, do not recurse
         down the branch. */
      if (!q->back->deleted && !q->back->deadend)
      {
        newnextnode->back = simcopytrav(srctree, q->back, desttree);
        newnextnode->back->back = newnextnode;
        q = q->next;
        if (newnextnode->back->index == lastnodeidx)
        {
          newnextnode->next = newnode;
          break;
        }
        if (q == srcNode)
        {
          newnextnode->next = newnode;
        }
        else
        {
          temp = newnextnode;
          newnextnode = desttree->get_forknode(desttree, q->index); // BUG.967 worried about index
          simcopynode(q, newnextnode);
          temp->next = newnextnode;
        }
      }
      else
      { /*xx this else and q=q->next are experimental
                 (seems to be working) */
        q = q->next;
      }

    } while (q != srcNode);
  }
  return newnode;
}  /* simcopytrav */


void simcopytree(void)
{
  /* Make a simplified copy of the current tree for rooting/unrooting
     on output.  Deleted notes are removed and lengths are consolidated. */

  simplifiedtree.tree_p->root = simcopytrav(curtree, curtree->root, simplifiedtree.tree_p);

  simplifiedtree.waswritten = waswritten;
  simplifiedtree.hasmult = hasmult;
  simplifiedtree.haslengths = haslengths;
  simplifiedtree.nolengths = nolengths;
  simplifiedtree.initialized = true;
} /* simcopytree */


void writebranchlength(double x)
{
  long w;

  /* write branch length onto output file, keeping track of what
     column of line you are in, and writing to correct precision */

  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if ((long)(100000*x) == 100000*(long)x)
  {
    if (!xmltree)
      putc(':', outtree);
    fprintf(outtree, "%*.1f", (int)(w + 2), x);
    col += w + 3;
  }
  else
  {
    if ((long)(100000*x) == 10000*(long)(10*x))
    {
      if (!xmltree)
        putc(':', outtree);
      fprintf(outtree, "%*.1f", (int)(w + 3), x);
      col += w + 4;
    }
    else
    {
      if ((long)(100000*x) == 1000*(long)(100*x))
      {
        if (!xmltree)
          putc(':', outtree);
        fprintf(outtree, "%*.2f", (int)(w + 4), x);
        col += w + 5;
      }
      else
      {
        if ((long)(100000*x) == 100*(long)(1000*x))
        {
          if (!xmltree)
            putc(':', outtree);
          fprintf(outtree, "%*.3f", (int)(w + 5), x);
          col += w + 6;
        }
        else
        {
          if ((long)(100000*x) == 10*(long)(10000*x))
          {
            if (!xmltree)
              putc(':', outtree);
            fprintf(outtree, "%*.4f", (int)(w + 6), x);
            col += w + 7;
          }
          else
          {
            if (!xmltree)
              putc(':', outtree);
            fprintf(outtree, "%*.5f", (int)(w + 7), x);
            col += w + 8;
          }
        }
      }
    }
  }
} /* writebranchlength */


void treeout(node *p, boolean writeparens, double addlength, long indent, boolean rooted)
{
  /* write out file with representation of final tree */
  long i, n, lastnodeidx = 0;
  Char c;
  double x;
  boolean comma;
  node *q;

  /* If this is a tip or there are no non-deleted branches from this node, render this node as a tip (write its name).  */
  if (p == curtree->root)
  {
    indent = 0;
    if (xmltree)
    {
      putc('\n', outtree);
      indent += 2;
      if (writeindented)
      {
        for (i = 1; i <= indent; i++)
          putc(' ', outtree);
      }
      fprintf(outtree, "<phylogeny");
      if (rooted)
        fprintf(outtree, " rooted=\"true\"");
      else
        fprintf(outtree, " rooted=\"false\"");
      fprintf(outtree, ">");  /* assumes no length at root! */
      indent += 2;
    }
    else
    {
      putc('(', outtree);
      if (writeindented)
      {
        putc('\n', outtree);
        indent += 2;
      }
    }
  }
  if (p->tip)
  {
    if (p->hasname)
    {
      n = 0;
      for (i = 1; i <= MAXNCH; i++)
      {
        if ((curtree->nodep[p->index - 1]->nayme[i - 1] != '\0') && (curtree->nodep[p->index - 1]->nayme[i - 1] != ' '))
          n = i;
      }
      if (xmltree)
      {
        putc('\n', outtree);
        if (writeindented)
        {
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
        }
        fprintf(outtree, "<clade");
        if (p->haslength)
        {
          fprintf(outtree, " branch_length=\"");
          x = p->length;
          writebranchlength(x);
          fprintf(outtree, "\"");
        }
        putc('>', outtree);
        fprintf(outtree, "<name>");
      }
      if ((!xmltree) && writeindented)
      {
        for (i = 1; i <= indent; i++)
          putc(' ', outtree);
      }
      for (i = 0; i < n; i++)
      {
        c = curtree->nodep[p->index - 1]->nayme[i];
        if (c == ' ')
          c = '_';
        putc(c, outtree);
      }
      col += n;
      if (xmltree)
        fprintf(outtree, "</name></clade>");
    }
  }
  else if (p->onebranch && p->onebranchnode->tip)
  {
    if (p->onebranchnode->hasname)
    {
      n = 0;
      for (i = 1; i <= MAXNCH; i++)
      {
        if ((curtree->nodep[p->index - 1]->nayme[i - 1] != '\0') && (curtree->nodep[p->index - 1]->nayme[i - 1] != ' '))
          n = i;
      }
      if (writeindented)
        indent += 2;
      if (xmltree)
      {
        putc('\n', outtree);
        if (writeindented)
        {
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
        }
        fprintf(outtree, "<clade");
        if ((p->haslength && writeparens) || p->onebranch)
        {
          if (!(p->onebranch && !p->onebranchhaslength))
          {
            fprintf(outtree, " branch_length=");
            if (p->onebranch)
              x = p->onebranchlength;
            else
              x = p->length;
            x += addlength;
            writebranchlength(x);
            fprintf(outtree, "\"");
          }
          putc('>', outtree);
          fprintf(outtree, "<name>");
        }
      }
      for (i = 0; i < n; i++)
      {
        c = p->onebranchnode->nayme[i];
        if (c == '_')
          c = ' ';
        putc(c, outtree);
      }
      col += n;
      if (xmltree)
        fprintf(outtree, "</name></clade>");
    }
  }
  else if (p->onebranch && !p->onebranchnode->tip)
  {
    treeout(p->onebranchnode, true, 0.0, indent, true);
  }
  else
  {
    /* Multiple non-deleted branch case:  go round the node
       recursing down the branches.  */
    if (xmltree)
    {
      putc('\n', outtree);
      if (writeindented)
      {
        for (i = 1; i <= indent; i++)
          putc(' ', outtree);
      }
      if (p == curtree->root)
      {
        fprintf(outtree, "<clade>");
        if (writeindented)
          indent += 2;
      }
    }
    if (p != curtree->root)
    {
      if (xmltree)
      {
        fprintf(outtree, "<clade");
        if ((p->haslength && writeparens) || p->onebranch)
        {
          if (!(p->onebranch && !p->onebranchhaslength))
          {
            fprintf(outtree, " branch_length=\"");
            if (p->onebranch)
              x = p->onebranchlength;
            else
              x = p->length;
            x += addlength;
            writebranchlength(x);
          }
          fprintf(outtree, "\">");
        }
        else fprintf(outtree, ">");
        if (writeindented)
          indent += 2;
      }
      else
      {
        if (writeindented)
        {
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
        }
        putc('(', outtree);
        if (writeindented)
        {
          indent += 2;
          putc('\n', outtree);
        }
      }
    }
    (col)++;
    q = p->next;
    while (q != p)
    {
      if (!q->back->deleted && !q->back->deadend)
        lastnodeidx = q->back->index;
      q = q->next;
    }
    q = p->next;
    while (q != p)
    {
      comma = true;
      /* If branch is deleted or is a dead end, do not recurse
         down the branch and do not write a comma afterwards.
      */
      if (!q->back->deleted && !q->back->deadend)
        treeout(q->back, true, 0.0, indent, true);
      else
        comma = false;
      if (q->back->index == lastnodeidx)
        comma = false;
      q = q->next;
      if (q == p)
        break;
      if ((q->next == p) && (q->back->deleted || q->back->deadend))
        break;
      if (comma && !xmltree)
      {
        putc(',', outtree);
        if (writeindented)
          fprintf(outtree, "\n");
      }
      (col)++;
      if ((!xmltree) && (col > 65) && !writeindented)
      {
        putc('\n', outtree);
        col = 0;
      }
    }
    /* The right paren ')' closes off this level of recursion. */
    if (p != curtree->root)
    {
      if (xmltree)
      {
        fprintf(outtree, "\n");
        if (writeindented)
        {
          indent -= 2;
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
        }
        fprintf(outtree, "</clade>");
      }
      else
      {
        if (writeindented)
        {
          putc('\n', outtree);
          indent -= 2;
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
        }
        putc(')', outtree);
      }
    }
    (col)++;
  }
  if (!xmltree)
    if ((p->haslength && writeparens) || p->onebranch)
    {
      if (!(p->onebranch && !p->onebranchhaslength))
      {
        if (p->onebranch)
          x = p->onebranchlength;
        else
          x = p->length;
        x += addlength;
        writebranchlength(x);
      }
    }
  if (p == curtree->root)
  {
    if (xmltree)
    {
      if (writeindented)
        fprintf(outtree, "\n   </clade>\n  </phylogeny>\n");
      else
        fprintf(outtree, "\n</clade>\n</phylogeny>\n");
      if (writeindented)
        putc('\n', outtree);
    }
    else
    {
      if (writeindented)
        putc('\n', outtree);
      putc(')', outtree);
    }
  }
}  /* treeout */


void maketemptriad(node **p, long index)
{
  /* Initiate an internal node with stubs for two children */
  long i, j;
  node *q;
  q = NULL;
  for (i = 1; i <= 3; i++)
  {
    *p = curtree->get_forknode(curtree, index);

    (*p)->index = index;
    (*p)->hasname = false;
    (*p)->haslength = false;
    (*p)->deleted=false;
    (*p)->deadend=false;
    (*p)->onebranch=false;
    (*p)->onebranchhaslength=false;
    for (j=0;j<MAXNCH;j++)
      (*p)->nayme[j] = '\0';
    (*p)->next = q;
    q = *p;
  }
  (*p)->next->next->next = *p;
  q = (*p)->next;
  while (*p != q)
  {
    (*p)->back = NULL;
    (*p)->tip = false;
    *p = (*p)->next;
  }
}  /* maketemptriad */


void roottreeout(boolean *userwantsrooted)
{
  /* write out file with representation of final tree */
  long trnum, trnumwide;
  boolean treeisrooted = false;

  treetrav(curtree->root);
  simcopytree(); /* Prepare a copy of the going tree without deleted branches */
  treesets[whichtree].tree_p->root = curtree->root;        /* Store the current root */

  if (nexus)
  {
    trnum = treenumber;
    trnumwide = 1;
    while (trnum >= 10)
    {
      trnum /= 10;
      trnumwide++;
    }
    fprintf(outtree, "TREE PHYLIP_%*ld = ", (int)trnumwide, treenumber);
    if (!(*userwantsrooted))
      fprintf(outtree, "[&U] ");
    else
      fprintf(outtree, "[&R] ");
    col += 15;
  }
  curtree->root = simplifiedtree.tree_p->root;                /* Point root at simplified tree */
  curtree->root->haslength = false;              /* Root should not have a length */
  if (curtree->root->tip)
    treeisrooted = true;
  else
  {
    if (curtree->root->next->next->next == curtree->root)
      treeisrooted = true;
    else
      treeisrooted = false;
  }
  if (*userwantsrooted && !treeisrooted)
    notrootedtorooted();
  if (!(*userwantsrooted) && treeisrooted)
    rootedtonotrooted();
  if ((*userwantsrooted && treeisrooted) ||
      (!(*userwantsrooted) && !treeisrooted))
  {
    treeout(curtree->root, true, 0.0, 0, true);
  }
  curtree->root = treesets[whichtree].tree_p->root;     /* Point root at original (real) tree */
  if (!xmltree)
  {
    if (hasmult)
      fprintf(outtree, "[%6.4f];\n", trweight);
    else
      fprintf(outtree, ";\n");
  }
}  /* roottreeout */


void notrootedtorooted(void)
{
  node *newbase, *temp;

  /* root halfway along leftmost branch of unrooted tree */
  /* create a new triad for the new base */
  maketemptriad(&newbase, nonodes+1);

  /* Take left branch and make it the left branch of newbase */
  newbase->next->back = curtree->root->next->back;
  newbase->next->next->back = curtree->root;
  /* If needed, divide length between left and right branches */
  if (newbase->next->back->haslength)
  {
    newbase->next->back->length /= 2.0;
    newbase->next->next->back->length =
      newbase->next->back->length;
    newbase->next->next->back->haslength = true;
  }
  /* remove leftmost ring node from old base ring */
  temp = curtree->root->next->next;
  curtree->root->next = temp;
  /* point root at new base and write the tree */
  curtree->root = newbase;
  treeout(curtree->root, true, 0.0, 0, true);
  /* (since tree mods are to simplified tree and will not be used
     for general purpose tree editing, much initialization can be
     skipped.) */
} /* notrootedtorooted */


void rootedtonotrooted(void)
{
  node *q, *r, *temp, *newbase;
  boolean sumhaslength = false;
  double sumlength = 0;

  /* Use the leftmost non-tip immediate descendant of the root,
     root at that, write a multifurcation with that as the base.
     If both descendants are tips, write tree as is. */
  curtree->root = simplifiedtree.tree_p->root;
  /* first, search for leftmost non-tip immediate descendent of root */
  q = curtree->root->next->back;
  r = curtree->root->next->next->back;
  if (q->tip && r->tip)
  {
    treeout(curtree->root, true, 0.0, 0, false);
  }
  else if (!(q->tip))
  {
    /* allocate new base pointer */

    newbase = curtree->get_forknode(curtree, 0); // BUG.967 there a better index than 0 ??

    newbase->next = q->next;
    q->next = newbase;
    q->back = r;
    r->back = q;
    if (q->haslength && r->haslength)
    {
      sumlength = q->length + r->length;
      sumhaslength = true;
    }
    if (sumhaslength)
    {
      q->length = sumlength;
      q->back->length = sumlength;
    }
    else
    {
      q->haslength = false;
      r->haslength = false;
    }
    curtree->root = newbase;
    treeout(curtree->root, true, 0.0, 0, false);
  }
  else if (q-tip && !(r->tip))
  {
    temp = r;
    do {
      temp = temp->next;
    } while (temp->next != r);

    newbase = curtree->get_forknode(curtree, 0); // BUG.967 there a better index than 0 ??

    newbase->next = temp->next;
    temp->next = newbase;
    q->back = r;
    r->back = q;
    if (q->haslength && r->haslength)
    {
      sumlength = q->length + r->length;
      sumhaslength = true;
    }
    if (sumhaslength)
    {
      q->length = sumlength;
      q->back->length = sumlength;
    }
    else
    {
      q->haslength = false;
      r->haslength = false;
    }
    curtree->root = newbase;
    treeout(curtree->root, true, 0.0, 0, false);
  }
} /* rootedtonotrooted */


void treewrite(boolean *done)
{
  /* write out tree to a file */
  long maxinput;
  boolean rooted;

  if ( curtree->root->deleted )
  {
    printf("Cannot write tree because every branch in the tree is deleted\n");
    return;
  }
  if (onfirsttree)
    openfile(&outtree, OUTTREE, "output tree file", "w", "retree", outtreename);
  if (nexus && onfirsttree)
  {
    fprintf(outtree, "#NEXUS\n");
    fprintf(outtree, "BEGIN TREES\n");
    fprintf(outtree, "TRANSLATE;\n"); /* MacClade needs this */
  }
  if (xmltree && onfirsttree)
  {
    fprintf(outtree, "<phyloxml xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\">\n");
  }
  onfirsttree = false;
  maxinput = 1;
  do {
    printf("Enter R if the tree is to be rooted\n");
    printf("OR enter U if the tree is to be unrooted: ");
    phyFillScreenColor();
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
    ch = (isupper(ch)) ? ch : toupper(ch);
    maxinput++;
    if (maxinput == 100)
    {
      printf("ERROR:  Too many tries at choosing option.\n");
      exxit(-1);
    }
  } while (ch != 'R' && ch != 'U');
  col = 0;
  rooted = (ch == 'R');
  roottreeout(&rooted);
  treenumber++;
  printf("\nTree written to file \"%s\".\n\n", outtreename);
  waswritten = true;
  anywritten = true;
  written = true;
  if (!(*done))
    printree();
}  /* treewrite */


void retree_window(adjwindow action)
{
  /* move viewing window of tree */
  switch (action)
  {
    case left:
      if (leftedge != 1)
        leftedge -= hscroll;
      break;

    case downn:
      /* The 'topedge + 3' is needed to allow downward scrolling
         when part of the tree is above the screen and only 1 or 2 lines
         are below it. */
      if (treelines - topedge + 3 >= screenlines)
        topedge += vscroll;
      break;

    case upp:
      if (topedge != 1)
        topedge -= vscroll;
      break;

    case right:
      if (leftedge < vscreenwidth+2)
      {
        if (hscroll > leftedge - vscreenwidth + 1)
          leftedge = vscreenwidth;
        else
          leftedge += hscroll;
      }
      break;
  }
  printree();
}  /* retree_window */


void getlength(double *length, reslttype *reslt, boolean *hslngth)
{
  long maxinput;
  double valyew;
  char tmp[100];

  valyew = 0.0;
  maxinput = 1;
  do {
    printf("\nEnter the new branch length\n");
    printf("OR enter U to leave the length unchanged\n");
    if (*hslngth)
      printf("OR enter R to remove the length from this branch: \n");
    getstryng(tmp);

    if (tmp[0] == 'u' || tmp[0] == 'U')
    {
      *reslt = quit;
      break;
    }
    else if (tmp[0] == 'r' || tmp[0] == 'R')
    {
      (*reslt) = remoov;
      break;}
    else if (sscanf(tmp, "%lf", &valyew) == 1)
    {
      (*reslt) = valid;
      break;}
    maxinput++;
    if (maxinput == 100)
    {
      printf("ERROR:  Too many tries at choosing option.\n");
      exxit(-1);
    }
  } while (1);
  (*length) = valyew;
}  /* getlength */


void changelength(void)
{
  /* change or specify the length of a tip */
  boolean hslngth;
  boolean ok;
  long i, w, maxinput;
  double length, x;
  Char ch;
  reslttype reslt;
  node *p;

  maxinput = 1;
  do {
    printf("Specify length of which branch (0 = all branches)? ");
    inpnum(&i, &ok);
    ok = (ok && i <= nonodes);
    if (ok && (i != 0))
      ok = (ok && !curtree->nodep[i - 1]->deleted);
    if (i == 0)
      ok = (curtree->nodep[i - 1] != curtree->root);
    maxinput++;
    if (maxinput == 100)
    {
      printf("ERROR:  Too many tries at choosing option.\n");
      exxit(-1);
    }
  } while (!ok);
  if (i != 0)
  {
    p = curtree->nodep[i - 1];
    putchar('\n');
    if (p->haslength)
    {
      x = p->length;
      if (x > 0.0)
        w = (long)(0.43429448222 * log(x));
      else if (x == 0.0)
        w = 0;
      else
        w = (long)(0.43429448222 * log(-x)) + 1;
      if (w < 0)
        w = 0;
      printf("The current length of this branch is %*.5f\n", (int)(w + 7), x);
    }
    else
      printf("This branch does not have a length\n");
    hslngth = p->haslength;
    getlength(&length, &reslt, &hslngth);
    switch (reslt)
    {
      case valid:
        copytree();

        p->length = length;
        p->haslength = true;
        if (p->back != NULL)
        {
          p->back->length = length;
          p->back->haslength = true;
        }
        break;

      case remoov:
        copytree();

        p->haslength = false;
        if (p->back != NULL)
          p->back->haslength = false;
        break;

      case quit:
        /* blank case */
        break;
    }
  }
  else
  {
    printf("\n (this operation cannot be undone)\n");
    maxinput = 1;
    do {
      printf("\n   enter U to leave the lengths unchanged\n");
      printf("OR enter R to remove the lengths from all branches: \n");
      phyFillScreenColor();
      if(scanf("%c%*[^\n]", &ch)) {}    // Read char and scan to EOL.
      (void)getchar();
      if (ch == '\n')
        ch = ' ';
      maxinput++;
      if (maxinput == 100)
      {
        printf("ERROR:  Too many tries at choosing option.\n");
        exxit(-1);
      }
    } while (ch != 'U' && ch != 'u' && ch != 'R' && ch != 'r');
    if (ch == 'R' || ch == 'r')
    {
      copytree();
      for (i = 0; i < spp; i++)
        curtree->nodep[i]->haslength = false;
      for (i = spp; i < nonodes; i++)
      {
        if (curtree->nodep[i] != NULL)
        {
          curtree->nodep[i]->haslength = false;
          curtree->nodep[i]->next->haslength = false;
          curtree->nodep[i]->next->next->haslength = false;
        }
      }
    }
  }
  printree();
}  /* changelength */


void changename(void)
{
  /* change or specify the name of a tip */
  boolean ok;
  long i, n, tipno;
  char tipname[100];

  for(;;)
  {
    for(;;)
    {
      printf("Specify name of which tip? (enter its number or 0 to quit): ");
      inpnum(&i, &ok);
      if (i > 0 && (i <= spp) && ok)
        if (!curtree->nodep[i - 1]->deleted)
        {
          tipno = i;
          break;
        }
      if (i == 0)
      {
        tipno = 0;
        break;
      }
    }
    if (tipno == 0)
      break;
    if (curtree->nodep[tipno - 1]->hasname)
    {
      n = 0;

      /* this is valid because names are padded out to MAXNCH with nulls */
      for (i = 1; i <= MAXNCH; i++)
      {
        if (curtree->nodep[tipno - 1]->nayme[i - 1] != '\0')
          n = i;
      }
      printf("The current name of tip %ld is \"", tipno);
      for (i = 0; i < n; i++)
        putchar(curtree->nodep[tipno - 1]->nayme[i]);
      printf("\"\n");
    }
    copytree();
    for (i = 0; i < MAXNCH; i++)
      curtree->nodep[tipno - 1]->nayme[i] = ' ';
    printf("Enter new tip name: ");
    i = 1;
    getstryng(tipname);
    strncpy(curtree->nodep[tipno-1]->nayme, tipname, MAXNCH);
    curtree->nodep[tipno - 1]->hasname = true;
    printree();
  }
  printree();
}  /* changename */


void cladecollapse(void)
{
  /* put code here to
     (1) ask for a node number of an interior node
     (2) check a boolean attribute of the node "collapsed"
     (3) if true, uncollapse it (just toggle that boolean off)
     (4) if false, change the boolean and ...
     (5)  ... also ask what the name of the node is to be.  Default to
          any string stored in the naym field for that node, and if
          there is none ask the user to choose between:
          (i) the name of the leftmost descendant, and
          (ii) the names of the leftmost and rightmost ones, hyphenated, or
          (iii) something the user enters.
     Make sure that the tree displays only out as far on any branch as the
     first "collapsed"
     Make sure that in writing out the tree the user is asked, if any
     "collapsed" boolean is true, whether they want to write out the
     collapsed tree or the full tree.
  */
  printf("\nFor now this routine is not yet written so nothing has been done.\n");
}  /* cladecollapse */


void clade(void)
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && i <= nonodes);
  if (ok)
  {
    subtree = (i > 0);
    if (subtree)
      nuroot = curtree->nodep[i - 1];
    else
      nuroot = curtree->root;
  }
  printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


void changeoutgroup(void)
{
  long i, maxinput;
  boolean ok;

  maxinput = 1;
  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && i >= 1 && i <= nonodes && i != curtree->root->index);
    if (ok)
      ok = (ok && !curtree->nodep[i - 1]->deleted);
    if (ok)
      ok = !curtree->nodep[curtree->nodep[i - 1]->back->index - 1]->deleted;
    if (ok)
      outgrno = i;
    maxinput++;
    if (maxinput == 100)
    {
      printf("ERROR:  Too many tries at choosing option.\n");
      exxit(-1);
    }
  } while (!ok);
  copytree();
  reroot(curtree, curtree->nodep[outgrno - 1]);
  printree();
  written = false;
}  /* changeoutgroup */


void redisplay(void)
{
  long maxinput;
  boolean done;
  char    ch;

  done = false;
  maxinput = 1;
  do {
    printf("\nNEXT? (Options: R . ");
    if (haslengths)
      printf("= ");
    printf("U W O ");
    if (haslengths)
      printf("M ");
    printf("T F D B N H J K L C Y + ? X Q) (? for Help) ");
    if(scanf("%c%*[^\n]", &ch)) {}      // Read char and scan to EOL.
    (void)getchar();
#ifdef WIN32
    phyFillScreenColor();
#endif
    if (ch == '\n')
      ch = ' ';
    ch = isupper(ch) ? ch : toupper(ch);
    if (ch == 'C' || ch == 'F' || ch == 'O' || ch == 'R' ||
        ch == 'U' || ch == 'X' || ch == 'Q' || ch == '.' ||
        ch == 'W' || ch == 'B' || ch == 'N' || ch == '?' ||
        ch == 'H' || ch == 'J' || ch == 'K' || ch == 'L' ||
        ch == '+' || ch == 'T' || ch == 'Y' || ch == 'D' ||
        (haslengths && ch == 'M') ||
        (haslengths && ch == '='))
    {
      switch (ch)
      {
        case 'R':
          rearrange();
          break;

        case '.':
          printree();
          break;

        case '=':
          togglelengths();
          break;

        case 'U':
          undo();
          break;

        case 'W':
          treewrite(&done);
          break;

        case 'O':
          changeoutgroup();
          break;

        case 'M':
          midpoint();
          break;

        case 'T':
          transpose(0);
          break;

        case 'F':
          flip(0);
          break;

        case 'C':
          cladecollapse();
          break;

        case 'Y':
          clade();
          break;

        case 'D':
          del_or_restore();
          break;

        case 'B':
          changelength();
          break;

        case 'N':
          changename();
          break;

        case 'H':
          retree_window(left);
          break;

        case 'J':
          retree_window(downn);
          break;

        case 'K':
          retree_window(upp);
          break;

        case 'L':
          retree_window(right);
          break;

        case '?':
          retree_help();
          break;

        case '+':
          if (treesread <numtrees)
          {
            readnext = true;
            done = true;
          }
          else
          {
            printf("No more trees to read in input file.\n");
          }
          break;

        case 'X':
          done = true;
          break;

        case 'Q':
          done = true;
          break;
      }
    }
    else
    {
      printf("Not a possible option!\n");
      maxinput++;
      if (maxinput == 100)
      {
        printf("ERROR:  Too many tries at choosing option.\n");
        exxit(-1);
      }
    }
  } while (!done);
  if (!written)
  {
    maxinput = 1;
    do {
      printf("Do you want to write out the tree to a file? (Y or N) ");
      phyFillScreenColor();
      if(scanf("%c%*[^\n]", &ch)) {}    // Read char and scan to EOL.
      (void)getchar();
      if (ch == '\n')
        ch = ' ';
      if (ch == 'Y' || ch == 'y')
        treewrite(&done);
      maxinput++;
      if (maxinput == 100)
      {
        printf("ERROR:  Too many tries at choosing option.\n");
        exxit(-1);
      }
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
}  /* redisplay */


void treeconstruct(void)
{
  /* constructs a tree from the pointers in nodep. */
  waswritten = false;
  readnext = false;

  do {
    whichtree = 0;
    othertree = 1;
    treesets[whichtree].initialized = false;
    treesets[othertree].initialized = false;
    /*xx    nodep = treesets[whichtree].nodep;   debug */
    subtree = false;
    topedge = 1;
    leftedge = 1;
    buildtree();
    curtree->root = treesets[whichtree].tree_p->root;
    readnext = false;
    written = false;
    waswritten = false;
    printree();
    redisplay();
    unbuildtree();
  } while (readnext);
}  /* treeconstruct */


void retreeread(
  char * intreename,
  int    intreenum,
  char * usefont,
  int    usebranchlengths)
{
  (void)intreenum;                      // RSGnote: unused.
  (void)usefont;                        // RSGnote: unused.

  printf("Hello from RetreeRead!\n");   // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Retree";
  phylipinit(argc, argv, NULL, true);

  if (usebranchlengths != 0)
  {
    //weights = true;
  }
  else
  {
    //weights = false;
  }

  intree = fopen(intreename, "r");
  how = use;
  readnext = false;
  whichtree = 0;
  othertree = 1;
  treesets[whichtree].initialized = false;
  treesets[othertree].initialized = false;
  /*xx    nodep = treesets[whichtree].nodep;   debug */
  subtree = false;
  topedge = 1;
  leftedge = 1;
  printf("calling buildtree\n");
  fflush(stdout);
  buildtree();
  curtree->root = treesets[whichtree].tree_p->root;
  readnext = false;
  written = false;
  waswritten = false;

  // save the initial tree for later modification
  outtree = fopen("JavaTree.txt", "w");
  nexus = false;
  xmltree = false;
  boolean roottest = true;
  printf("calling roottreeout\n");
  roottreeout(&roottest);  // hardwired for testing JRM

  // display initial tree
  workingplot = fopen("JavaPreview.ps", "wb");
  //previewing = false;
  //previewer = lw;  // hardwired to ps
  //initplotter(spp, fontname);
  //drawit(fontname, &xoffset, &yoffset, numlines, curtree->root);
  //finishplotter();

  fclose(intree);
  fclose(outtree);
  fclose(workingplot);

  printf("Done.\n\n");

}


void retree(
  char * intreename,
  char * intree,
  char * outtreename,
  char * outtreeopt,
  char * outtreefmt,
  char * usefont,
  char * plotfilename,
  char * plotfileopt,
  int    usebranchlengths,
  char * option,
  int    activenode,
  int    refnode,
  double branchlength,
  int    rearrbefore,
  char * newname,
  int    rootnode,
  char * outtreeformat,
  int    readtree,
  int    writetree,
  int    update,
  int    undo,
  int    doplot,
  char * finalplotkind)
{
  (void)intreename;                     // RSGnote: unused so far.
  (void)intree;                         // RSGnote: unused so far.
  (void)outtreename;                    // RSGnote: unused so far.
  (void)outtreeopt;                     // RSGnote: unused so far.
  (void)outtreefmt;                     // RSGnote: unused so far.
  (void)usefont;                        // RSGnote: unused so far.
  (void)plotfilename;                   // RSGnote: unused so far.
  (void)plotfileopt;                    // RSGnote: unused so far.
  (void)usebranchlengths;               // RSGnote: unused so far.
  (void)option;                         // RSGnote: unused so far.
  (void)activenode;                     // RSGnote: unused so far.
  (void)refnode;                        // RSGnote: unused so far.
  (void)branchlength;                   // RSGnote: unused so far.
  (void)rearrbefore;                    // RSGnote: unused so far.
  (void)newname;                        // RSGnote: unused so far.
  (void)rootnode;                       // RSGnote: unused so far.
  (void)outtreeformat;                  // RSGnote: unused so far.
  (void)readtree;                       // RSGnote: unused so far.
  (void)writetree;                      // RSGnote: unused so far.
  (void)update;                         // RSGnote: unused so far.
  (void)undo;                           // RSGnote: unused so far.
  (void)doplot;                         // RSGnote: unused so far.
  (void)finalplotkind;                  // RSGnote: unused so far.

  printf("Hello from Retree!\n"); // JRMdebug

  int argc;
  Char *argv[1];
  argc = 1;
  argv[0] = "Retree";
  phylipinit(argc, argv, NULL, true);

  /*
    char * intreename,
    char * intree,
    char * outtreename,
    char * outtreeopt,
    char * outtreefmt,
    char * usefont,
    char * plotfilename,
    char * plotfileopt,
    int    usebranchlengths,
    char * option,
    int    activenode,
    int    refnode,
    double branchlength,
    int    rearrbefore,
    char * newname,
    int    rootnode,
    char * outtreeformat,
    int    readtree,
    int    writetree,
    int    update,
    int    undo,
    int    doplot,
    char * finalplotkind);
  */
}


int main(int argc, Char *argv[])
{
  /* Reads in spp.  Then calls treeconstruct() to construct the tree and query the user. */
  int i;

#ifdef MAC
  argc = 1;                /* macsetup("Retree", "");        */
  argv[0] = "Retree";
#endif
  phylipinit(argc, argv, NULL, false);
#if 0
  treesets[0].nodep = treeone;
  treesets[1].nodep = treetwo;
#endif
  nexus     = false;
  nolengths = false;
  scrollinc = 20;
  screenlines = 24;
  screenwidth = 80;
  vscreenwidth = 80;
  ibmpc = IBMCRT;
  ansi  = ANSICRT;
  for(i = 0; i < maxsz; i++)
    delarray[i] = false;
  treenumber = 1;
  getoptions();
  configure();
  treeconstruct();
  if (anywritten)
  {
    if (xmltree)
      fprintf(outtree, "</phyloxml>\n");
    else if (nexus)
      fprintf(outtree, "END;\n");
  }
  FClose(intree);
  FClose(outtree);
  printf("\nDone.\n\n");
#ifdef MAC
  fixmacfile(outtreename);
#endif
  phyRestoreConsoleAttributes();
  return 0;
}  /* Retree */


// End.
