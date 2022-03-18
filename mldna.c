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



