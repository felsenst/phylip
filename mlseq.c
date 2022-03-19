/* Copyright 2022
 * Functions for ML on sequences, applicable to DNA, protein, codon, and restriction sites data
 * Written by Mike Palczewski and Joe Felsenstein */


void allocx(long nonodes, long endsite, long param, ml_node** treenode)
{
  /* allocate sequences */
  /* param =  sitelength in restml */
  /* param =  rcategs in dnaml/proml */
  long i;
  ml_node *p;
  ml_node *q;

  for (i = 0; i < spp; i++)
    treenode[i]->allocx((node*)treenode[i], endsite, param);
  for (i = spp; i < nonodes; i++)
  {
    p = treenode[i];
    q = p;
    do
    {
      q->allocx((node*)q, endsite, param);
      q = (ml_node*)q->node.next;
    } while ( q != p);
  }
}  /* allocx */


