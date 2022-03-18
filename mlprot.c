/* Copyright, 2022 */
/* Routines for protein ML */

void prot_node_copy(node* srcn, node* destn)
{
  /* Copy a node when data is protein and method is likelihood */
  prot_node *src = (prot_node *)srcn;
  prot_node *dest = (prot_node *)destn;
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
} /* prot_node_copy */


void fix_protx(prot_node* p, long site, double maxx, long rcategs)
{ /* used in Proml, Promlk */
  long i, m;

  ((ml_node*)p)->underflows[site] += log(maxx);

  for ( i = 0 ; i < rcategs  ; i++ )
    for (m = 0; m <= 19; m++)
      p->x[site][i][m] /= maxx;
} /* fix_protx */


node * prot_node_new(node_type type, long index, long nodesize) // RSGbugfix
{
  /* create a node for a protein tree */

  node *n = Malloc(sizeof(struct prot_node));
/* debug: instead do a call upwards here? */
  prot_node_init(n, type, index);
  return n;
} /* prot_node_new */


void prot_node_init(node *n, node_type type, long index)
{
  /* initialize a node for a protein tree */

  prot_node *pn = (prot_node *)n;

  ml_node_init(n, type, index);
  pn->ml_node.allocx = (allocx_t)prot_node_allocx;
  pn->ml_node.node.copy = prot_node_copy;
  pn->ml_node.node.init = prot_node_init;
  pn->ml_node.freex = (freex_t)prot_node_freex;
  pn->x = NULL;
  if ( endsite != 0 && rcategs != 0 )
    pn->ml_node.allocx((struct node*)(&(pn->ml_node.node)), endsite, rcategs);
} /* prot_node_init */


