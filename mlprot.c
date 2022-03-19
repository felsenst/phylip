/* Copyright, 2022 */
/* Routines for protein ML */

#ifndef OLDC  /* prototypes, if not original Kernighan & Ritchie compiler */
node*   ml_prot_node_new(node_type, long, long);
void    ml_prot_node_init(struct node *, node_type, long);
void    ml_prot_node_allocx(ml_node*, long, long);
void    ml_prot_node_copy(node*, node*);
void    ml_prot_node_freex(ml_node*);
void    ml_prot_freex_notip(long, pointarray);
void    fix_protx(ml_ml_prot_node*, long, double, long);
#endif


typedef struct ml_prot_node {
  struct ml_node ml_node;                     /* Base object, must be first */
  pphenotype x;
} ml_prot_node;

void ml_prot_node_copy(node* srcn, node* destn)
{
  /* Copy a node when data is protein and method is likelihood */
  ml_prot_node *src = (ml_prot_node *)srcn;
  ml_prot_node *dest = (ml_prot_node *)destn;
  long i, j;
  long oldendsite = dest->ml_node.endsite;

  ml_node_copy(srcn, destn);
  if ( oldendsite != 0 && oldendsite != src->ml_node.endsite )
  {
    ((ml_node*)destn)->freex((node*)dest);
    oldendsite = 0;
  }
  if ( oldendsite == 0 )
    ((ml_node*)dest)->allocx((node*)dest, ((ml_node*)src)->endsite,
                              ((ml_node*)src)->categs);

  for (i = 0; i < src->ml_node.endsite; i++)
    for (j = 0; j < src->ml_node.categs; j++)
      memcpy(dest->x[i][j], src->x[i][j], sizeof(psitelike));
} /* ml_prot_node_copy */


void fix_protx(ml_prot_node* p, long site, double maxx, long rcategs)
{ /* used in Proml, Promlk */
  long i, m;

  ((ml_node*)p)->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs  ; i++ )
    for (m = 0; m <= 19; m++)
      p->x[site][i][m] /= maxx;
} /* fix_protx */


node * ml_prot_node_new(node_type type, long index, long nodesize) // RSGbugfix
{
  /* create a node for a protein tree */

  node *n = Malloc(sizeof(struct ml_prot_node));
/* debug: instead do a call upwards here? */
  ml_prot_node_init(n, type, index);
  return n;
} /* ml_prot_node_new */


void ml_prot_node_init(node *n, node_type type, long index)
{
  /* initialize a node for a protein tree */

  ml_prot_node *pn = (ml_prot_node *)n;

  ml_node_init(n, type, index);
  pn->ml_node.allocx = (allocx_t)ml_prot_node_allocx;
  pn->ml_node.node.copy = ml_prot_node_copy;
  pn->ml_node.node.init = ml_prot_node_init;
  pn->ml_node.freex = (freex_t)ml_prot_node_freex;
  pn->x = NULL;
  if ( endsite != 0 && rcategs != 0 )
    pn->ml_node.allocx((struct node*)(&(pn->ml_node.node)), endsite, rcategs);
} /* ml_prot_node_init */


void ml_prot_node_freex(ml_node* n)
{
  /* free a protein tree node */
  ml_prot_node *pn;
  long i;

  pn = (ml_prot_node *)n;
  for ( i = 0 ; i < n->endsite ; i++ )
  {
    free(pn->x[i]);
  }

  free(pn->x);
  pn->x = NULL;
  free(n->underflows);
  n->underflows = NULL;
} /* ml_prot_node_freex */


void ml_prot_node_allocx(ml_node* nn, long endsite, long rcategs)
{
  /* allocate space for sequences on a protein tree node */
  ml_prot_node *n = (ml_prot_node *)nn;
  long i;

  n->ml_node.categs = rcategs;
  n->ml_node.endsite = endsite;

  n->x = (pphenotype)Malloc(endsite * sizeof(pratelike));
  for ( i = 0 ; i < endsite ; i++ )
    n->x[i] = (pratelike)Malloc(rcategs * sizeof(psitelike));
  n->ml_node.underflows= Malloc(endsite * sizeof(double));
} /* ml_prot_node_allocx */


void ml_prot_freex_notip(long nonodes, pointarray treenode)
{
  /* free interior fork nodes
   * used in proml */
  long i, j;
  node *p;

  for (i = spp; i < nonodes; i++)
  {
    p = treenode[i];
    if ( p == NULL ) continue;
    do {
      for (j = 0; j < endsite; j++)
      {
        free(((ml_prot_node*)p)->x[j]);
        ((ml_prot_node*)p)->x[j] = NULL;
      }
      free(((ml_prot_node*)p)->x);
      ((ml_prot_node*)p)->x = NULL;
      p = p->next;
    } while (p != treenode[i]);
  }
}  /* ml_prot_freex_notip */


