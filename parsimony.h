/* Version 4.0. (c) Copyright 1993-2013 by the University of Washington.
   Written by Joe Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Michal Palczewski.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


typedef struct pars_node{
  node node;
  steptr numsteps;               /* bookkeeps steps                        */

} pars_node;

typedef long (*pars_supplement_t)(tree*, long);
typedef boolean (*branchcollapsible_t)(tree*, node*);

typedef struct pars_tree {
  tree tree;
  pars_supplement_t supplement;
  branchcollapsible_t branchcollapsible;

} pars_tree;


#ifndef OLDC
/* prototypes */
void    pars_node_init(node* p, node_type type, long index);
void    pars_node_copy(node* src, node* dst);
void    pars_node_print(node* n);
void    pars_node_reinit(node* p);
tree*   pars_tree_new(long nonodes, long spp);
void    pars_tree_init(tree* t, long nonodes, long spp);

#if 0                                   // RSGbugfix: Never used.
void    pars_tree_re_move(tree *, node **, node **, boolean);
#endif

void    updatenumdesc(node *, node *, long);
void    preorder(tree*, node *, node *, node *, node *, node *, long);
void    branchlength(node *, node *, double *, pointarray);
void    minpostorder(node *, pointarray);
void    collapsebestrees(tree *, bestelm *, long *, long, boolean, long *);
void    reroot(node *, node *);
void    reroot2(node *, node *);
void    collapsetree(tree * t, node * n);
void    savetree(tree *, long *);
void    addnsave(tree *, node *, node *, boolean, long *);
void    inittreetrav(node *, long);
void    printbranchlengths(node *);
void    branchlentrav(node *, node *, long, long, double *, pointarray);
void    initbranchlen(node *p);
boolean allcommonbases(node *, node *, boolean *);
void    flipindexes(long, pointarray);
long    sibsvisited(node *, long *);
boolean parentinmulti(node *, node*);
long    smallest(node *, long *);
boolean alltips(node *, node *);
void    newindex(long, node *);
void    load_tree(tree*, long, bestelm*);
void    pars_node_free(node **);
void    pars_globrearrange(tree*, boolean, boolean);
boolean treecollapsible(tree*, node*);
void    collapsebranch(tree*, node*);
void    writesteps(tree*, long, boolean, steptr);
void    addbestever(long, long *, long, boolean, long *, bestelm *, double);
void    addtiedtree(long *, long *, long, boolean, long *, bestelm *, double);
void    add_to_besttrees(tree*, long, bestelm*);
boolean	pars_addtraverse(tree*, node*, node*, boolean, node*,
                                   double*, bestelm*, boolean, boolean);
void    treeout3(node *, long, long *, node *);
void    disc_treelength(node *, long, pointarray);
void    grandrearr(tree*, boolean, boolean);
void    drawline3(long, double, node *);
boolean pars_tree_try_insert_(tree*, node *, node *, node *, double *, tree*,
                               tree*, boolean, boolean, boolean *);
void    coordinates(tree*, node *, double, long *, double *);
void    printree(tree*);
node*   root_tree(tree*, node*);
void    reroot_tree(tree*, node*); // RSGbugfix: Name change.
void    initparsnode(tree *, node **, long, long, long *, long *, initops,
                      pointarray, Char *, Char *, FILE *);
double  pars_tree_evaluate(tree*, node*, boolean);
bestelm*  allocbestree(void) ;
bestelm** allocbestrees(void);
#endif


// End.
